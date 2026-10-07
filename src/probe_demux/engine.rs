use std::collections::{BTreeMap, HashMap};
use std::fs::{self, File};
use std::io::{BufWriter, Write as _};
use std::path::Path;
use std::time::Instant;

use anyhow::{Context, Result};
use rust_htslib::bam::{self, Read};
use tracing::info;

use crate::probe_sample_sheet::ProbeSampleSheet;
use super::types::{Assignment, Cli, DemuxCfg, Metrics, OverlapRow, SampleCount, Status};

/// Negative binomial NB(r, p) 99.9th-percentile quantile, via the standard
/// PMF recurrence `P(k) = P(k-1) * p * (k-1+r) / k`, `P(0) = (1-p)^r` — exact
/// for real-valued `r` without a gamma-function dependency. Returns the
/// smallest `k` whose CDF clears `quantile`. `r <= 0.0` (no occupied GEMs to
/// fit against) disables the filter by returning `u64::MAX`.
pub fn negbinom_quantile(r: f64, p: f64, quantile: f64) -> u64 {
    if r <= 0.0 {
        return u64::MAX;
    }
    let mut pmf = (1.0 - p).powf(r); // P(0)
    let mut cdf = pmf;
    let mut k: u64 = 0;
    while cdf < quantile && k < 100_000 {
        k += 1;
        pmf *= p * (k as f64 - 1.0 + r) / k as f64;
        cdf += pmf;
    }
    k
}

/// Assign every cell barcode (GEM) in `tally` (keyed `(cell_barcode,
/// probe_barcode_id) -> molecule count`) to every sample from `sheet` whose
/// declared probe-barcode set it genuinely carries.
pub fn assign(
    sheet: &ProbeSampleSheet,
    tally: &BTreeMap<(String, String), u64>,
    cfg: &DemuxCfg,
) -> Vec<Assignment> {
    let mut per_cb: BTreeMap<&str, BTreeMap<&str, u64>> = BTreeMap::new();
    for ((cb, pb), &n) in tally {
        *per_cb
            .entry(cb.as_str())
            .or_default()
            .entry(pb.as_str())
            .or_insert(0) += n;
    }

    struct CbInfo<'a> {
        cb: &'a str,
        total: u64,
        probe_barcodes: Vec<String>,
        /// Samples whose declared set clears `min_molecules`, desc by sum.
        occupants: Vec<(&'a str, u64)>,
    }

    let mut infos: Vec<CbInfo> = Vec::with_capacity(per_cb.len());
    for (cb, pb_counts) in per_cb {
        let total: u64 = pb_counts.values().sum();
        let probe_barcodes: Vec<String> = pb_counts.keys().map(|s| s.to_string()).collect();

        let mut sample_sums: Vec<(&str, u64)> = sheet
            .samples
            .iter()
            .map(|s| {
                let sum: u64 = s
                    .probe_barcode_ids
                    .iter()
                    .filter_map(|p| pb_counts.get(p.as_str()))
                    .sum();
                (s.id.as_str(), sum)
            })
            .collect();
        sample_sums.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(b.0)));
        let occupants: Vec<(&str, u64)> = sample_sums
            .into_iter()
            .filter(|&(_, n)| n >= cfg.min_molecules)
            .collect();

        infos.push(CbInfo {
            cb,
            total,
            probe_barcodes,
            occupants,
        });
    }

    let occupied: Vec<u64> = infos
        .iter()
        .filter(|i| !i.occupants.is_empty())
        .map(|i| i.occupants.len() as u64)
        .collect();
    let r_hat = if occupied.is_empty() {
        0.0
    } else {
        occupied.iter().sum::<u64>() as f64 / occupied.len() as f64
    };
    let occupancy_cutoff = negbinom_quantile(r_hat, 0.5, 0.999);

    let mut out = Vec::with_capacity(infos.len());
    for info in infos {
        let occupancy = info.occupants.len();
        if occupancy == 0 {
            out.push(Assignment {
                cell_barcode: info.cb.to_string(),
                sample_id: None,
                probe_barcodes: info.probe_barcodes,
                molecules: info.total,
                purity: 0.0,
                occupancy,
                status: Status::Unassigned,
            });
        } else if occupancy as u64 > occupancy_cutoff {
            let best_sum = info.occupants[0].1;
            let purity = if info.total > 0 {
                best_sum as f64 / info.total as f64
            } else {
                0.0
            };
            out.push(Assignment {
                cell_barcode: info.cb.to_string(),
                sample_id: None,
                probe_barcodes: info.probe_barcodes,
                molecules: info.total,
                purity,
                occupancy,
                status: Status::HighOccupancyGem,
            });
        } else {
            for (sample_id, sum) in &info.occupants {
                out.push(Assignment {
                    cell_barcode: info.cb.to_string(),
                    sample_id: Some(sample_id.to_string()),
                    probe_barcodes: info.probe_barcodes.clone(),
                    molecules: *sum,
                    purity: if info.total > 0 {
                        *sum as f64 / info.total as f64
                    } else {
                        0.0
                    },
                    occupancy,
                    status: Status::Assigned,
                });
            }
        }
    }
    out
}

pub fn gem_overlap(
    tally: &BTreeMap<(String, String), u64>,
    declared: &std::collections::HashSet<&str>,
) -> Vec<OverlapRow> {
    let mut pbs_per_cb: BTreeMap<&str, std::collections::BTreeSet<&str>> = BTreeMap::new();
    for (cb, pb) in tally.keys() {
        if !declared.contains(pb.as_str()) {
            continue;
        }
        pbs_per_cb
            .entry(cb.as_str())
            .or_default()
            .insert(pb.as_str());
    }
    let mut gems_per_pb: HashMap<&str, u64> = HashMap::new();
    let mut common: HashMap<(&str, &str), u64> = HashMap::new();
    for pbs in pbs_per_cb.values() {
        for &pb in pbs {
            *gems_per_pb.entry(pb).or_insert(0) += 1;
        }
        let list: Vec<&str> = pbs.iter().copied().collect();
        for i in 0..list.len() {
            for j in (i + 1)..list.len() {
                let key = if list[i] < list[j] {
                    (list[i], list[j])
                } else {
                    (list[j], list[i])
                };
                *common.entry(key).or_insert(0) += 1;
            }
        }
    }
    let mut rows: Vec<OverlapRow> = common
        .into_iter()
        .map(|((a, b), c)| {
            let ga = gems_per_pb[a];
            let gb = gems_per_pb[b];
            OverlapRow {
                barcode1: a.to_string(),
                barcode2: b.to_string(),
                barcode1_gems: ga as usize,
                barcode2_gems: gb as usize,
                common_gems: c as usize,
                overlap: c as f64 / ga.min(gb) as f64,
            }
        })
        .collect();
    rows.sort_by(|x, y| y.overlap.partial_cmp(&x.overlap).unwrap());
    rows
}

use crate::bam_helpers::{aux_i64, aux_string};

/// Tally molecule counts per `(CB, PB)` from an annotated+corrected BAM,
/// gated on the `xf` representative-molecule bit (avoids PCR-duplicate
/// inflation — matches `cr_bam_to_mtx.rs`'s own `xf & 8` gate). Records
/// missing `CB` or with an empty/missing `PB` (no confident probe-barcode
/// call) are skipped.
pub fn read_bam_tally(path: &Path) -> Result<BTreeMap<(String, String), u64>> {
    let mut reader =
        bam::Reader::from_path(path).with_context(|| format!("opening BAM {}", path.display()))?;
    let mut tally = BTreeMap::new();
    for result in reader.records() {
        let record = result.context("reading BAM record")?;
        if aux_i64(&record, b"xf")?.unwrap_or(0) & 8 == 0 {
            continue;
        }
        let Some(cb) = aux_string(&record, b"CB")? else {
            continue;
        };
        let Some(pb) = aux_string(&record, b"PB")? else {
            continue;
        };
        if pb.is_empty() {
            continue;
        }
        *tally.entry((cb, pb)).or_insert(0) += 1;
    }
    Ok(tally)
}

fn write_sample_assignment(assignments: &[Assignment], path: &Path) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(
        w,
        "cell_barcode,sample_id,probe_barcodes,molecules,purity,occupancy,status"
    )?;
    for a in assignments {
        writeln!(
            w,
            "{},{},{},{},{:.4},{},{}",
            a.cell_barcode,
            a.sample_id.as_deref().unwrap_or(""),
            a.probe_barcodes.join("|"),
            a.molecules,
            a.purity,
            a.occupancy,
            a.status.as_str(),
        )?;
    }
    Ok(())
}

/// `matrix --cell-barcodes` (see `cr_bam_to_mtx.rs::read_cell_barcodes`)
/// expects *bare* ACGT barcodes and appends its own `-{gem_group}` suffix —
/// it rejects a barcode that already carries one (not pure ACGT). Our own
/// `CB` tag convention always carries `-1` already (see `lib.rs`'s
/// `minimap2_header` emission), so strip it back off here to keep
/// `<sample>.barcodes.txt` directly usable, matching the plan's own
/// "feeds `matrix --cell-barcodes` directly" requirement. Strips any
/// trailing `-<digits>`, not just `-1`, in case a future multi-gem-group run
/// carries a different suffix.
pub fn strip_gem_suffix(cb: &str) -> &str {
    match cb.rfind('-') {
        Some(i) if !cb[i + 1..].is_empty() && cb[i + 1..].bytes().all(|b| b.is_ascii_digit()) => {
            &cb[..i]
        }
        _ => cb,
    }
}

pub fn write_per_sample_barcodes(
    sheet: &ProbeSampleSheet,
    assignments: &[Assignment],
    out_dir: &Path,
) -> Result<()> {
    let mut by_sample: HashMap<&str, Vec<&str>> = HashMap::new();
    for a in assignments {
        if let (Status::Assigned, Some(sid)) = (a.status, a.sample_id.as_deref()) {
            by_sample
                .entry(sid)
                .or_default()
                .push(a.cell_barcode.as_str());
        }
    }
    for s in &sheet.samples {
        let path = out_dir.join(format!("{}.barcodes.txt", s.id));
        let mut w = BufWriter::new(File::create(&path)?);
        if let Some(cbs) = by_sample.get(s.id.as_str()) {
            for cb in cbs {
                writeln!(w, "{}", strip_gem_suffix(cb))?;
            }
        }
    }
    Ok(())
}

/// `<sample>.probe_barcodes.txt`: this sample's declared probe-barcode ids,
/// one per line, straight from the sample sheet (not derived from
/// assignments — a sample's declared set is static). Feeds
/// `cr-bam-to-mtx --probe-barcodes` so that a GEM shared with another sample
/// (see the module doc comment) only contributes *this* sample's own
/// molecules to *this* sample's matrix, not the shared GEM's whole total.
pub fn write_per_sample_probe_barcodes(sheet: &ProbeSampleSheet, out_dir: &Path) -> Result<()> {
    for s in &sheet.samples {
        let path = out_dir.join(format!("{}.probe_barcodes.txt", s.id));
        let mut w = BufWriter::new(File::create(&path)?);
        for pb in &s.probe_barcode_ids {
            writeln!(w, "{pb}")?;
        }
    }
    Ok(())
}

fn write_gem_overlap(rows: &[OverlapRow], path: &Path) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(
        w,
        "barcode1,barcode2,barcode1_gems,barcode2_gems,common_gems,overlap"
    )?;
    for r in rows {
        writeln!(
            w,
            "{},{},{},{},{},{:.6}",
            r.barcode1, r.barcode2, r.barcode1_gems, r.barcode2_gems, r.common_gems, r.overlap
        )?;
    }
    Ok(())
}

pub fn build_metrics(
    sheet: &ProbeSampleSheet,
    assignments: &[Assignment],
    elapsed_secs: f64,
) -> Metrics {
    let mut per_sample_counts: HashMap<&str, usize> = HashMap::new();
    let mut distinct_cbs: std::collections::HashSet<&str> = std::collections::HashSet::new();
    let mut occupancy_by_cb: HashMap<&str, usize> = HashMap::new();
    let mut assigned = 0;
    let mut unassigned = 0;
    // One row per GEM for this status (unlike Assigned, never duplicated
    // across samples), so the row count is already the distinct-GEM count.
    let mut high_occupancy_gems = 0;
    for a in assignments {
        distinct_cbs.insert(a.cell_barcode.as_str());
        occupancy_by_cb.insert(a.cell_barcode.as_str(), a.occupancy);
        match a.status {
            Status::Assigned => {
                assigned += 1;
                if let Some(sid) = &a.sample_id {
                    *per_sample_counts.entry(sid.as_str()).or_insert(0) += 1;
                }
            }
            Status::Unassigned => unassigned += 1,
            Status::HighOccupancyGem => high_occupancy_gems += 1,
        }
    }
    let occupied: Vec<u64> = occupancy_by_cb
        .values()
        .filter(|&&occ| occ > 0)
        .map(|&occ| occ as u64)
        .collect();
    let occupancy_r = if occupied.is_empty() {
        0.0
    } else {
        occupied.iter().sum::<u64>() as f64 / occupied.len() as f64
    };
    let occupancy_cutoff = negbinom_quantile(occupancy_r, 0.5, 0.999);

    Metrics {
        elapsed_secs,
        cell_barcodes: distinct_cbs.len(),
        assigned,
        unassigned,
        high_occupancy_gems,
        occupancy_r,
        occupancy_cutoff,
        per_sample: sheet
            .samples
            .iter()
            .map(|s| SampleCount {
                sample: s.id.clone(),
                cell_barcodes: per_sample_counts.get(s.id.as_str()).copied().unwrap_or(0),
            })
            .collect(),
    }
}

pub fn run(cli: &Cli) -> Result<()> {
    let started = Instant::now();
    let sheet = ProbeSampleSheet::load(&cli.probe_sample_sheet)?;
    info!(samples = sheet.samples.len(), path = %cli.probe_sample_sheet.display(), "loaded probe sample sheet");

    let tally = read_bam_tally(&cli.bam)?;
    let cfg = DemuxCfg {
        min_molecules: cli.min_molecules,
    };
    let assignments = assign(&sheet, &tally, &cfg);
    let declared: std::collections::HashSet<&str> = sheet
        .samples
        .iter()
        .flat_map(|s| s.probe_barcode_ids.iter().map(String::as_str))
        .collect();
    let overlap = gem_overlap(&tally, &declared);

    fs::create_dir_all(&cli.output)
        .with_context(|| format!("creating {}", cli.output.display()))?;
    write_sample_assignment(&assignments, &cli.output.join("sample_assignment.csv"))?;
    write_per_sample_barcodes(&sheet, &assignments, &cli.output)?;
    write_per_sample_probe_barcodes(&sheet, &cli.output)?;
    write_gem_overlap(&overlap, &cli.output.join("gem_barcode_overlap.csv"))?;

    let metrics = build_metrics(&sheet, &assignments, started.elapsed().as_secs_f64());
    let metrics_path = cli.output.join("probe_demux_metrics.json");
    fs::write(&metrics_path, serde_json::to_string_pretty(&metrics)?)
        .with_context(|| format!("writing {}", metrics_path.display()))?;

    info!(
        cell_barcodes = metrics.cell_barcodes,
        assigned = metrics.assigned,
        unassigned = metrics.unassigned,
        high_occupancy_gems = metrics.high_occupancy_gems,
        occupancy_r = format!("{:.3}", metrics.occupancy_r),
        occupancy_cutoff = metrics.occupancy_cutoff,
        elapsed = format!("{:.1}s", metrics.elapsed_secs),
        out = %cli.output.display(),
        "RUN SUMMARY"
    );
    Ok(())
}
