use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs;
use std::io::Write;
use std::time::Instant;
use anyhow::{bail, Result};


use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use statrs::function::gamma::ln_gamma;

use crate::cr_bam_to_mtx::{write_h5, write_matrix};
use super::h5::{filtered_matrix, read_raw, resolve_h5, totals};
use super::types::round_count;
use super::types::{
    Cli, Metrics, RawMatrix, MIN_RECOVERED_CELLS, ORDMAG_BOOTSTRAPS, ORDMAG_QUANTILE,
};

use super::helpers::{
    bootstrap_counts, estimate_recovered, resolve_tagseq_barcodes, within_ordmag,
};



fn ordmag(
    totals: &[u64],
    expected: Option<usize>,
    rng: &mut ChaCha8Rng,
    max_expected: usize,
) -> (Vec<usize>, usize, u64) {
    let nonzero: Vec<u64> = totals.iter().copied().filter(|&x| x > 0).collect();
    if nonzero.is_empty() {
        return (Vec::new(), 0, 0);
    }
    let recovered = expected
        .unwrap_or_else(|| estimate_recovered(&nonzero, max_expected.min(nonzero.len()), rng))
        .max(MIN_RECOVERED_CELLS);
    let baseline =
        ((recovered as f64 * (1.0 - ORDMAG_QUANTILE)).round() as usize).min(nonzero.len() - 1);
    let mut called_counts = Vec::with_capacity(ORDMAG_BOOTSTRAPS);
    for _ in 0..ORDMAG_BOOTSTRAPS {
        called_counts.push(within_ordmag(&bootstrap_counts(&nonzero, rng), baseline));
    }
    let target =
        (called_counts.iter().sum::<usize>() as f64 / called_counts.len() as f64).round() as usize;
    let mut ranked: Vec<usize> = (0..totals.len()).filter(|&i| totals[i] > 0).collect();
    ranked.sort_by_key(|&i| (std::cmp::Reverse(totals[i]), i));
    let target = target.min(ranked.len());
    if target == 0 {
        return (Vec::new(), recovered, 0);
    }
    let cutoff = totals[ranked[target - 1]];
    let called = ranked
        .into_iter()
        .take_while(|&i| totals[i] >= cutoff)
        .collect();
    (called, recovered, cutoff)
}

fn linear_slope(x: &[f64], y: &[f64]) -> f64 {
    let xm = x.iter().sum::<f64>() / x.len() as f64;
    let ym = y.iter().sum::<f64>() / y.len() as f64;
    let numerator: f64 = x.iter().zip(y).map(|(a, b)| (a - xm) * (b - ym)).sum();
    let denominator: f64 = x.iter().map(|a| (a - xm).powi(2)).sum();
    numerator / denominator
}

fn sgt_proportions(freq: &[u64]) -> Option<(Vec<f64>, f64)> {
    let mut nr = BTreeMap::<u64, u64>::new();
    for &r in freq {
        *nr.entry(r).or_default() += 1;
    }
    if nr.len() < 10 {
        return None;
    }
    let r: Vec<f64> = nr.keys().map(|&x| x as f64).collect();
    let n: Vec<f64> = nr.values().map(|&x| x as f64).collect();
    let total: f64 = r.iter().zip(&n).map(|(a, b)| a * b).sum();
    let mut z = Vec::with_capacity(r.len());
    for i in 0..r.len() {
        let left = if i == 0 { 0.0 } else { r[i - 1] };
        let right = if i + 1 == r.len() {
            2.0 * r[i] - left
        } else {
            r[i + 1]
        };
        z.push(n[i] / (0.5 * (right - left)));
    }
    let mut slope = linear_slope(
        &r.iter().map(|x| x.ln()).collect::<Vec<_>>(),
        &z.iter().map(|x| x.ln()).collect::<Vec<_>>(),
    );
    if slope >= -1.0 {
        slope = -1.0;
    }
    let smooth: Vec<f64> = r
        .iter()
        .map(|&x| x * (1.0 + 1.0 / x).powf(1.0 + slope))
        .collect();
    let p0 = n[0] / total;
    let denominator: f64 = smooth.iter().zip(&n).map(|(a, b)| a * b).sum();
    let adjusted: HashMap<u64, f64> = nr
        .keys()
        .copied()
        .zip(smooth.iter().map(|x| (1.0 - p0) * x / denominator))
        .collect();
    Some((freq.iter().map(|r| adjusted[r]).collect(), p0))
}

struct AliasSampler {
    prob: Vec<f64>,
    alias: Vec<usize>,
}

impl AliasSampler {
    fn new(weights: &[f64]) -> Self {
        let n = weights.len();
        let mut scaled: Vec<f64> = weights.iter().map(|x| x * n as f64).collect();
        let (mut small, mut large) = (Vec::new(), Vec::new());
        for (i, &x) in scaled.iter().enumerate() {
            if x < 1.0 {
                small.push(i)
            } else {
                large.push(i)
            }
        }
        let mut prob = vec![1.0; n];
        let mut alias = (0..n).collect::<Vec<_>>();
        while let (Some(s), Some(l)) = (small.pop(), large.pop()) {
            prob[s] = scaled[s];
            alias[s] = l;
            scaled[l] += scaled[s] - 1.0;
            if scaled[l] < 1.0 {
                small.push(l)
            } else {
                large.push(l)
            }
        }
        Self { prob, alias }
    }
    fn draw(&self, rng: &mut ChaCha8Rng) -> usize {
        let i = rng.gen_range(0..self.prob.len());
        if rng.gen::<f64>() < self.prob[i] {
            i
        } else {
            self.alias[i]
        }
    }
}

pub fn bh_adjust(p: &[f64]) -> Vec<f64> {
    let mut order: Vec<usize> = (0..p.len()).collect();
    order.sort_by(|&a, &b| p[a].total_cmp(&p[b]));
    let mut adjusted = vec![1.0; p.len()];
    let mut running: f64 = 1.0;
    for (rank0, &index) in order.iter().enumerate().rev() {
        running = running.min(p[index] * p.len() as f64 / (rank0 + 1) as f64);
        adjusted[index] = running.min(1.0);
    }
    adjusted
}

fn emptydrops(
    raw: &RawMatrix,
    totals: &[u64],
    ordmag: &[usize],
    cli: &Cli,
    rng: &mut ChaCha8Rng,
) -> (Vec<usize>, Vec<(usize, f64, f64)>, usize, u64, u64) {
    let mut ranked: Vec<usize> = (0..totals.len()).collect();
    ranked.sort_by_key(|&i| (std::cmp::Reverse(totals[i]), i));
    let empty: Vec<usize> = ranked
        .iter()
        .skip(cli.ambient_lower)
        .take(cli.ambient_upper.saturating_sub(cli.ambient_lower))
        .copied()
        .collect();
    let ambient: Vec<usize> = empty.iter().copied().filter(|&i| totals[i] > 0).collect();
    let max_ambient = empty.iter().map(|&i| totals[i]).max().unwrap_or(0);
    let threshold = cli.emptydrops_minimum_umis.max(max_ambient + 1);
    if ambient.is_empty() {
        return (Vec::new(), Vec::new(), 0, max_ambient, threshold);
    }
    let ordmag_set: BTreeSet<usize> = ordmag.iter().copied().collect();
    let candidates: Vec<usize> = (0..totals.len())
        .filter(|i| !ordmag_set.contains(i) && totals[*i] >= threshold)
        .collect();
    if candidates.is_empty() {
        return (
            Vec::new(),
            Vec::new(),
            ambient.len(),
            max_ambient,
            threshold,
        );
    }

    // Summed in f64 so fractional matrices keep every partial molecule;
    // the ambient profile is rounded per feature (SGT needs integer
    // frequencies), which is exact for integer matrices.
    let mut ambient_sums = vec![0f64; raw.features.len()];
    let mut global_counts = vec![0f64; raw.features.len()];
    for bc in 0..raw.barcodes.len() {
        for pos in raw.indptr[bc] as usize..raw.indptr[bc + 1] as usize {
            global_counts[raw.indices[pos] as usize] += raw.data[pos];
        }
    }
    for &bc in &ambient {
        for pos in raw.indptr[bc] as usize..raw.indptr[bc + 1] as usize {
            ambient_sums[raw.indices[pos] as usize] += raw.data[pos];
        }
    }
    let ambient_counts: Vec<u64> = ambient_sums.into_iter().map(round_count).collect();
    let eval_features: Vec<usize> = (0..raw.features.len())
        .filter(|&i| global_counts[i] > 0.0)
        .collect();
    let positive: Vec<u64> = eval_features
        .iter()
        .map(|&i| ambient_counts[i])
        .filter(|&x| x > 0)
        .collect();
    let Some((smoothed, p0)) = sgt_proportions(&positive) else {
        return (
            Vec::new(),
            Vec::new(),
            ambient.len(),
            max_ambient,
            threshold,
        );
    };
    let zeros = eval_features
        .iter()
        .filter(|&&i| ambient_counts[i] == 0)
        .count();
    let mut probabilities = vec![0.0; raw.features.len()];
    let mut cursor = 0;
    for &feature in &eval_features {
        probabilities[feature] = if ambient_counts[feature] == 0 {
            p0 / zeros.max(1) as f64
        } else {
            let p = smoothed[cursor];
            cursor += 1;
            p
        };
    }
    let eval_prob: Vec<f64> = eval_features.iter().map(|&i| probabilities[i]).collect();
    let sampler = AliasSampler::new(&eval_prob);
    let feature_to_eval: HashMap<usize, usize> = eval_features
        .iter()
        .enumerate()
        .map(|(i, &f)| (f, i))
        .collect();
    let logp: Vec<f64> = eval_prob.iter().map(|p| p.ln()).collect();
    let mut observed = Vec::with_capacity(candidates.len());
    for &bc in &candidates {
        let n = totals[bc];
        let mut ll = ln_gamma(n as f64 + 1.0);
        for pos in raw.indptr[bc] as usize..raw.indptr[bc + 1] as usize {
            let count = raw.data[pos];
            let e = feature_to_eval[&(raw.indices[pos] as usize)];
            ll -= ln_gamma(count + 1.0);
            ll += count * logp[e];
        }
        observed.push(ll);
    }
    let mut by_n = HashMap::<u64, Vec<usize>>::new();
    for (i, &bc) in candidates.iter().enumerate() {
        by_n.entry(totals[bc]).or_default().push(i);
    }
    let max_n = *by_n.keys().max().unwrap() as usize;
    let mut extreme = vec![0usize; candidates.len()];
    let mut seen = vec![0u32; eval_features.len()];
    let mut touched = Vec::new();
    for _ in 0..cli.simulations {
        let mut ll = 0.0;
        for draw_number in 1..=max_n {
            let feature = sampler.draw(rng);
            if seen[feature] == 0 {
                touched.push(feature);
            }
            seen[feature] += 1;
            ll += (draw_number as f64).ln() - (seen[feature] as f64).ln() + logp[feature];
            if let Some(indexes) = by_n.get(&(draw_number as u64)) {
                for &index in indexes {
                    if ll <= observed[index] {
                        extreme[index] += 1;
                    }
                }
            }
        }
        for &feature in &touched {
            seen[feature] = 0;
        }
        touched.clear();
    }
    let pvalues: Vec<f64> = extreme
        .iter()
        .map(|&x| (x + 1) as f64 / (cli.simulations + 1) as f64)
        .collect();
    let adjusted = bh_adjust(&pvalues);
    let calls: Vec<usize> = candidates
        .iter()
        .zip(&adjusted)
        .filter_map(|(&bc, &q)| (q <= cli.fdr).then_some(bc))
        .collect();
    let diagnostics = candidates
        .into_iter()
        .zip(pvalues)
        .zip(adjusted)
        .map(|((bc, p), q)| (bc, p, q))
        .collect();
    (calls, diagnostics, ambient.len(), max_ambient, threshold)
}

pub fn run(cli: &Cli) -> Result<()> {
    if cli.ambient_lower >= cli.ambient_upper {
        bail!("--ambient-lower must be less than --ambient-upper");
    }
    if cli.simulations == 0 {
        bail!("--simulations must be greater than zero");
    }
    if !(0.0..=1.0).contains(&cli.fdr) {
        bail!("--fdr must lie between 0 and 1");
    }
    let started = Instant::now();
    let input = resolve_h5(&cli.raw_matrix);
    let raw = read_raw(&input)?;
    let totals = totals(&raw);
    let nonzero = totals.iter().filter(|&&x| x > 0).count();
    let mut rng = ChaCha8Rng::seed_from_u64(cli.seed);
    let tagseq_calls: Option<Vec<usize>> = cli
        .tagseq_barcodes
        .as_ref()
        .map(|path| resolve_tagseq_barcodes(path, &raw.barcodes))
        .transpose()?;
    let (ordmag_calls, recovered, cutoff) = if let Some(calls) = &tagseq_calls {
        let cutoff = calls.iter().map(|&bc| totals[bc]).min().unwrap_or(0);
        (calls.clone(), calls.len(), cutoff)
    } else if let Some(force) = cli.force_cells {
        let mut ranked: Vec<usize> = (0..totals.len()).filter(|&i| totals[i] > 0).collect();
        ranked.sort_by_key(|&i| (std::cmp::Reverse(totals[i]), i));
        let n = force.min(ranked.len());
        let cutoff = n.checked_sub(1).map(|i| totals[ranked[i]]).unwrap_or(0);
        (ranked.into_iter().take(n).collect(), force, cutoff)
    } else {
        ordmag(&totals, cli.expect_cells, &mut rng, cli.ambient_lower)
    };
    let (empty_calls, diagnostics, ambient_n, max_ambient, ed_threshold) =
        if tagseq_calls.is_some() || cli.force_cells.is_some() {
            (Vec::new(), Vec::new(), 0, 0, cli.emptydrops_minimum_umis)
        } else {
            emptydrops(&raw, &totals, &ordmag_calls, cli, &mut rng)
        };
    let mut selected = ordmag_calls.clone();
    selected.extend(&empty_calls);
    selected.sort_unstable();
    selected.dedup();
    let filtered = filtered_matrix(&raw, &selected);
    fs::create_dir_all(&cli.output)?;
    let mex = cli.output.join("filtered_feature_bc_matrix");
    let h5 = cli.output.join("filtered_feature_bc_matrix.h5");
    write_matrix(&filtered, &mex, true)?;
    write_h5(&filtered, &h5, &raw.genome)?;
    let ordmag_set: BTreeSet<usize> = ordmag_calls.iter().copied().collect();
    let empty_set: BTreeSet<usize> = empty_calls.iter().copied().collect();
    let diagnostic_map: HashMap<usize, (f64, f64)> = diagnostics
        .into_iter()
        .map(|(bc, p, q)| (bc, (p, q)))
        .collect();
    let mut calls = fs::File::create(cli.output.join("cell_calls.csv"))?;
    writeln!(calls, "barcode,total_umis,call,pvalue,adjusted_pvalue")?;
    for &bc in &selected {
        let call = if tagseq_calls.is_some() {
            "tagseq"
        } else if cli.force_cells.is_some() {
            "forced"
        } else if ordmag_set.contains(&bc) {
            "ordmag"
        } else if empty_set.contains(&bc) {
            "emptydrops"
        } else {
            "forced"
        };
        let (p, q) = diagnostic_map
            .get(&bc)
            .copied()
            .unwrap_or((f64::NAN, f64::NAN));
        writeln!(
            calls,
            "{},{},{},{},{}",
            raw.barcodes[bc], totals[bc], call, p, q
        )?;
    }
    let mut filtered_barcodes = fs::File::create(cli.output.join("filtered_barcodes.csv"))?;
    writeln!(filtered_barcodes, "genome,barcode")?;
    for &bc in &selected {
        writeln!(filtered_barcodes, "{},{}", raw.genome, raw.barcodes[bc])?;
    }
    let metrics = Metrics {
        algorithm: if tagseq_calls.is_some() {
            "explicit tag-seq sample list"
        } else {
            "OrdMag + multinomial EmptyDrops"
        },
        raw_matrix: input.display().to_string(),
        filtered_h5: h5.display().to_string(),
        raw_barcodes: raw.barcodes.len(),
        nonzero_barcodes: nonzero,
        recovered_cells: recovered,
        ordmag_cells: ordmag_calls.len(),
        ordmag_cutoff_umis: cutoff,
        ambient_rank_lower: cli.ambient_lower,
        ambient_rank_upper: cli.ambient_upper,
        ambient_barcodes_nonzero: ambient_n,
        max_ambient_umis: max_ambient,
        emptydrops_minimum_umis: ed_threshold,
        emptydrops_candidates: diagnostic_map.len(),
        emptydrops_cells: empty_calls.len(),
        emptydrops_simulations: cli.simulations,
        emptydrops_fdr: cli.fdr,
        filtered_cells: selected.len(),
        filtered_molecules: filtered.molecules,
        elapsed_secs: started.elapsed().as_secs_f64(),
    };
    fs::write(
        cli.output.join("cell_calling_metrics.json"),
        serde_json::to_string_pretty(&metrics)?,
    )?;
    eprintln!(
        "[sbxpress cells] OrdMag={} EmptyDrops={} filtered={} molecules={} -> {}",
        metrics.ordmag_cells,
        metrics.emptydrops_cells,
        metrics.filtered_cells,
        metrics.filtered_molecules,
        cli.output.display()
    );
    Ok(())
}
