//! 10x-style GEX cell calling: OrdMag followed by multinomial EmptyDrops.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, Context, Result};
use clap::Parser;
use hdf5::types::VarLenAscii;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use serde::Serialize;
use statrs::function::gamma::ln_gamma;

use crate::cr_bam_to_mtx::{write_h5, write_matrix, Feature, Matrix};

const ORDMAG_BOOTSTRAPS: usize = 100;
const ORDMAG_QUANTILE: f64 = 0.99;
const MIN_RECOVERED_CELLS: usize = 50;

#[derive(Debug, Parser)]
#[command(about = "Call cells with 10x OrdMag + EmptyDrops and write filtered matrices")]
pub struct Cli {
    /// raw_feature_bc_matrix.h5, or an outs directory containing it.
    #[arg(value_name = "RAW_MATRIX")]
    pub raw_matrix: PathBuf,

    /// 10x-style outs directory receiving filtered outputs.
    #[arg(short, long, default_value = "outs")]
    pub output: PathBuf,

    /// Expected recovered cells. If omitted, estimate with the 10x OrdMag grid search.
    #[arg(long, conflicts_with = "force_cells")]
    pub expect_cells: Option<usize>,

    /// Select exactly the top N nonzero barcodes and skip EmptyDrops.
    #[arg(long, conflicts_with = "expect_cells")]
    pub force_cells: Option<usize>,

    /// Skip OrdMag and EmptyDrops entirely and call exactly these barcodes
    /// as cells: one per line, matching the extraction whitelist format
    /// (a gem-group "-1" suffix, if present, is ignored). OrdMag's
    /// order-of-magnitude knee and EmptyDrops' ambient-profile test are
    /// both shaped around droplet single-cell data -- a few real cells
    /// against a large background of empty droplets with a shared ambient
    /// RNA noise profile. A pooled tag-seq run has a different noise
    /// mechanism (index hopping into a small closed barcode whitelist) and
    /// the real barcodes are already known from the experiment design, so
    /// there's nothing to infer statistically.
    #[arg(long, conflicts_with_all = ["expect_cells", "force_cells"])]
    pub tagseq_barcodes: Option<PathBuf>,

    /// First barcode rank used for the ambient background (NextGEM default).
    #[arg(long, default_value_t = 45_000)]
    pub ambient_lower: usize,

    /// Exclusive last barcode rank used for the ambient background.
    #[arg(long, default_value_t = 90_000)]
    pub ambient_upper: usize,

    /// Minimum EmptyDrops candidate UMIs; raised above the ambient maximum.
    #[arg(long, default_value_t = 500)]
    pub emptydrops_minimum_umis: u64,

    /// Monte Carlo multinomial simulations per tested UMI total.
    #[arg(long, default_value_t = 100_000)]
    pub simulations: usize,

    /// Benjamini-Hochberg FDR for 5' / NextGEM EmptyDrops calls.
    #[arg(long, default_value_t = 0.001)]
    pub fdr: f64,

    /// Reproducible bootstrap and Monte Carlo seed.
    #[arg(long, default_value_t = 0)]
    pub seed: u64,
}

/// 10x-style CSC feature-barcode matrix (features x barcodes, `indptr`
/// column/barcode-indexed, `indices` = feature row per nonzero entry).
/// `pub(crate)`: also read directly by `protospacer_calling.rs` rather than
/// duplicating this H5-reading logic a third time.
pub(crate) struct RawMatrix {
    pub(crate) features: Vec<Feature>,
    pub(crate) barcodes: Vec<String>,
    pub(crate) genome: String,
    pub(crate) data: Vec<u32>,
    pub(crate) indices: Vec<u32>,
    pub(crate) indptr: Vec<u32>,
}

#[derive(Serialize)]
struct Metrics {
    algorithm: &'static str,
    raw_matrix: String,
    filtered_h5: String,
    raw_barcodes: usize,
    nonzero_barcodes: usize,
    recovered_cells: usize,
    ordmag_cells: usize,
    ordmag_cutoff_umis: u64,
    ambient_rank_lower: usize,
    ambient_rank_upper: usize,
    ambient_barcodes_nonzero: usize,
    max_ambient_umis: u64,
    emptydrops_minimum_umis: u64,
    emptydrops_candidates: usize,
    emptydrops_cells: usize,
    emptydrops_simulations: usize,
    emptydrops_fdr: f64,
    filtered_cells: usize,
    filtered_molecules: u64,
    elapsed_secs: f64,
}

/// Resolve a --tagseq-barcodes file into raw-matrix column indices. The
/// gem-group suffix ("-1") that the matrix carries but a plain whitelist
/// file doesn't is ignored on both sides of the match.
fn resolve_tagseq_barcodes(path: &Path, barcodes: &[String]) -> Result<Vec<usize>> {
    let by_prefix: HashMap<&str, usize> = barcodes
        .iter()
        .enumerate()
        .map(|(i, bc)| (bc.split('-').next().unwrap_or(bc.as_str()), i))
        .collect();
    let text = fs::read_to_string(path)
        .with_context(|| format!("reading --tagseq-barcodes {}", path.display()))?;
    let mut indices = Vec::new();
    for (lineno, line) in text.lines().enumerate() {
        let bc = line.trim();
        if bc.is_empty() {
            continue;
        }
        let bare = bc.split('-').next().unwrap_or(bc);
        let &index = by_prefix.get(bare).with_context(|| {
            format!(
                "{}:{}: barcode {bare} is not present in the raw matrix",
                path.display(),
                lineno + 1
            )
        })?;
        indices.push(index);
    }
    if indices.is_empty() {
        bail!("--tagseq-barcodes {} contained no barcodes", path.display());
    }
    indices.sort_unstable();
    indices.dedup();
    Ok(indices)
}

fn resolve_h5(path: &Path) -> PathBuf {
    if path.is_dir() {
        path.join("raw_feature_bc_matrix.h5")
    } else {
        path.to_owned()
    }
}

fn ascii_strings(dataset: hdf5::Dataset) -> Result<Vec<String>> {
    Ok(dataset
        .read_raw::<VarLenAscii>()?
        .into_iter()
        .map(|value| value.as_str().to_owned())
        .collect())
}

pub(crate) fn read_raw(path: &Path) -> Result<RawMatrix> {
    let file = hdf5::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let matrix = file.group("matrix")?;
    let features_group = matrix.group("features")?;
    let ids = ascii_strings(features_group.dataset("id")?)?;
    let names = ascii_strings(features_group.dataset("name")?)?;
    let genomes = ascii_strings(features_group.dataset("genome")?)?;
    let barcodes = ascii_strings(matrix.dataset("barcodes")?)?;
    let data = matrix.dataset("data")?.read_raw::<u32>()?;
    let indices = matrix.dataset("indices")?.read_raw::<u32>()?;
    let indptr = matrix.dataset("indptr")?.read_raw::<u32>()?;
    if ids.len() != names.len() || data.len() != indices.len() || indptr.len() != barcodes.len() + 1
    {
        bail!("invalid 10x CSC matrix dimensions");
    }
    Ok(RawMatrix {
        features: ids
            .into_iter()
            .zip(names)
            .map(|(id, name)| Feature { id, name })
            .collect(),
        barcodes,
        genome: genomes.first().cloned().unwrap_or_else(|| "unknown".into()),
        data,
        indices,
        indptr,
    })
}

fn totals(raw: &RawMatrix) -> Vec<u64> {
    (0..raw.barcodes.len())
        .map(|bc| {
            let start = raw.indptr[bc] as usize;
            let end = raw.indptr[bc + 1] as usize;
            raw.data[start..end].iter().map(|&x| x as u64).sum()
        })
        .collect()
}

fn within_ordmag(sorted_ascending: &[u64], baseline_idx: usize) -> usize {
    let idx = baseline_idx.min(sorted_ascending.len() - 1);
    let baseline = sorted_ascending[sorted_ascending.len() - idx - 1];
    let cutoff = ((baseline as f64 * 0.1).round() as u64).max(1);
    let first = sorted_ascending.partition_point(|&x| x < cutoff);
    sorted_ascending.len() - first
}

fn bootstrap_counts(source: &[u64], rng: &mut ChaCha8Rng) -> Vec<u64> {
    let mut sample = Vec::with_capacity(source.len());
    for _ in 0..source.len() {
        sample.push(source[rng.gen_range(0..source.len())]);
    }
    sample.sort_unstable();
    sample
}

fn recovered_grid(maximum: usize) -> Vec<usize> {
    let max_log = (maximum.max(2) as f64).log2();
    let mut values = BTreeSet::new();
    for i in 0..2000 {
        let exponent = 1.0 + (max_log - 1.0) * i as f64 / 1999.0;
        values.insert(2f64.powf(exponent).round() as usize);
    }
    values.into_iter().collect()
}

fn estimate_recovered(nonzero: &[u64], max_expected: usize, rng: &mut ChaCha8Rng) -> usize {
    let grid = recovered_grid(max_expected);
    let mut estimates = Vec::with_capacity(ORDMAG_BOOTSTRAPS);
    for _ in 0..ORDMAG_BOOTSTRAPS {
        let sample = bootstrap_counts(nonzero, rng);
        let mut best = (f64::INFINITY, MIN_RECOVERED_CELLS);
        for &expected in &grid {
            let baseline = ((expected as f64 * (1.0 - ORDMAG_QUANTILE)).round() as usize)
                .min(sample.len() - 1);
            let observed = within_ordmag(&sample, baseline);
            let loss = (observed as f64 - expected as f64).powi(2) / expected as f64;
            if loss < best.0 {
                best = (loss, expected);
            }
        }
        estimates.push(best.1);
    }
    ((estimates.iter().sum::<usize>() as f64 / estimates.len() as f64).round() as usize)
        .max(MIN_RECOVERED_CELLS)
}

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

/// Simple Good-Turing proportions, following the Gale-Sampson implementation
/// used by 10x reference. Returns None when fewer than ten frequency classes exist.
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

fn bh_adjust(p: &[f64]) -> Vec<f64> {
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

    let mut ambient_counts = vec![0u64; raw.features.len()];
    let mut global_counts = vec![0u64; raw.features.len()];
    for bc in 0..raw.barcodes.len() {
        for pos in raw.indptr[bc] as usize..raw.indptr[bc + 1] as usize {
            global_counts[raw.indices[pos] as usize] += raw.data[pos] as u64;
        }
    }
    for &bc in &ambient {
        for pos in raw.indptr[bc] as usize..raw.indptr[bc + 1] as usize {
            ambient_counts[raw.indices[pos] as usize] += raw.data[pos] as u64;
        }
    }
    let eval_features: Vec<usize> = (0..raw.features.len())
        .filter(|&i| global_counts[i] > 0)
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
            let count = raw.data[pos] as f64;
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

fn filtered_matrix(raw: &RawMatrix, selected: &[usize]) -> Matrix {
    let selected_set: BTreeSet<usize> = selected.iter().copied().collect();
    let ordered: Vec<usize> = (0..raw.barcodes.len())
        .filter(|i| selected_set.contains(i))
        .collect();
    let mut counts = BTreeMap::new();
    let mut molecules = 0u64;
    for (new_bc, &old_bc) in ordered.iter().enumerate() {
        for pos in raw.indptr[old_bc] as usize..raw.indptr[old_bc + 1] as usize {
            let count = raw.data[pos] as u64;
            counts.insert((new_bc, raw.indices[pos] as usize), count);
            molecules += count;
        }
    }
    Matrix {
        features: raw.features.clone(),
        barcodes: ordered.iter().map(|&i| raw.barcodes[i].clone()).collect(),
        counts,
        records: 0,
        molecules,
        barcodes_skipped: 0,
        probe_barcodes_skipped: 0,
    }
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

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn bh_is_monotone_in_rank() {
        let q = bh_adjust(&[0.01, 0.04, 0.03, 0.002]);
        assert!((q[3] - 0.008).abs() < 1e-12);
        assert!(q[0] <= q[2] && q[2] <= q[1]);
    }
    #[test]
    fn ordmag_threshold_uses_top_one_percent() {
        let mut x = vec![1; 1000];
        x.extend([1000; 10]);
        x.sort_unstable();
        assert_eq!(within_ordmag(&x, 1), 10);
    }

    #[test]
    fn force_cells_writes_filtered_outputs() {
        let temp = tempfile::tempdir().unwrap();
        let raw_path = temp.path().join("raw.h5");
        let matrix = Matrix {
            features: vec![Feature {
                id: "g1".into(),
                name: "G1".into(),
            }],
            barcodes: vec!["bc1-1".into(), "bc2-1".into()],
            counts: [((0, 0), 10), ((1, 0), 2)].into_iter().collect(),
            records: 0,
            molecules: 12,
            barcodes_skipped: 0,
            probe_barcodes_skipped: 0,
        };
        write_h5(&matrix, &raw_path, "GRCh38").unwrap();
        let output = temp.path().join("outs");
        run(&Cli {
            raw_matrix: raw_path,
            output: output.clone(),
            expect_cells: None,
            force_cells: Some(1),
            tagseq_barcodes: None,
            ambient_lower: 1,
            ambient_upper: 2,
            emptydrops_minimum_umis: 1,
            simulations: 10,
            fdr: 0.001,
            seed: 0,
        })
        .unwrap();
        assert!(output.join("filtered_feature_bc_matrix.h5").is_file());
        assert_eq!(
            fs::read_to_string(output.join("filtered_barcodes.csv")).unwrap(),
            "genome,barcode\nGRCh38,bc1-1\n"
        );
    }

    #[test]
    fn tagseq_barcodes_selects_exact_list_ignoring_umi_rank() {
        let temp = tempfile::tempdir().unwrap();
        let raw_path = temp.path().join("raw.h5");
        let matrix = Matrix {
            features: vec![Feature {
                id: "g1".into(),
                name: "G1".into(),
            }],
            barcodes: vec!["bc1-1".into(), "bc2-1".into()],
            // bc1 has far more UMIs than bc2 -- OrdMag/force_cells would
            // pick bc1, proving a --tagseq-barcodes selection of bc2 alone
            // comes from the list, not the counts.
            counts: [((0, 0), 10), ((1, 0), 2)].into_iter().collect(),
            records: 0,
            molecules: 12,
            barcodes_skipped: 0,
            probe_barcodes_skipped: 0,
        };
        write_h5(&matrix, &raw_path, "GRCh38").unwrap();
        let barcodes_path = temp.path().join("samples.txt");
        // Bare barcode, no gem-group suffix, matching the whitelist format.
        fs::write(&barcodes_path, "bc2\n").unwrap();
        let output = temp.path().join("outs");
        run(&Cli {
            raw_matrix: raw_path,
            output: output.clone(),
            expect_cells: None,
            force_cells: None,
            tagseq_barcodes: Some(barcodes_path),
            ambient_lower: 1,
            ambient_upper: 2,
            emptydrops_minimum_umis: 1,
            simulations: 10,
            fdr: 0.001,
            seed: 0,
        })
        .unwrap();
        assert_eq!(
            fs::read_to_string(output.join("filtered_barcodes.csv")).unwrap(),
            "genome,barcode\nGRCh38,bc2-1\n"
        );
        let calls = fs::read_to_string(output.join("cell_calls.csv")).unwrap();
        assert!(calls.contains("bc2-1,2,tagseq"));
        let metrics_json = fs::read_to_string(output.join("cell_calling_metrics.json")).unwrap();
        assert!(metrics_json.contains("explicit tag-seq sample list"));
    }

    #[test]
    fn tagseq_barcodes_rejects_unknown_barcode() {
        let temp = tempfile::tempdir().unwrap();
        let raw_path = temp.path().join("raw.h5");
        let matrix = Matrix {
            features: vec![Feature {
                id: "g1".into(),
                name: "G1".into(),
            }],
            barcodes: vec!["bc1-1".into()],
            counts: [((0, 0), 10)].into_iter().collect(),
            records: 0,
            molecules: 10,
            barcodes_skipped: 0,
            probe_barcodes_skipped: 0,
        };
        write_h5(&matrix, &raw_path, "GRCh38").unwrap();
        let barcodes_path = temp.path().join("samples.txt");
        fs::write(&barcodes_path, "not_in_the_matrix\n").unwrap();
        let output = temp.path().join("outs");
        let err = run(&Cli {
            raw_matrix: raw_path,
            output,
            expect_cells: None,
            force_cells: None,
            tagseq_barcodes: Some(barcodes_path),
            ambient_lower: 1,
            ambient_upper: 2,
            emptydrops_minimum_umis: 1,
            simulations: 10,
            fdr: 0.001,
            seed: 0,
        })
        .unwrap_err();
        assert!(err.to_string().contains("not_in_the_matrix"));
    }
}
