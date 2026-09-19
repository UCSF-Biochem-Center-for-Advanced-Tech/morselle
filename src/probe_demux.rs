//! Assign cell barcodes to samples by their declared probe-barcode set
//! (FLEX_SUPPORT_PLAN.md Phase 2).
//!
//! Published-algorithm parity notes (see FLEX_SUPPORT_PLAN.md "Where this
//! plan corrects the initial design"): assignment sums molecule counts across
//! a sample's *whole declared set* (which can be >1 probe barcode — real
//! `config.csv` example: `Pos,D-A01|D-A02`) before comparing samples, not a
//! naive top-single-barcode vote.
//!
//! Composite cell barcode (verified against 10x's own "Cell Ranger's Flex
//! Algorithm" support doc, "Multiplex workflow" + "High occupancy GEM
//! filtering"): in real Flex multiplex pooling, the true cell barcode is the
//! *(GEM barcode, Probe Barcode)* pair, not the GEM barcode alone — multiple
//! fixed, permeabilized samples are pooled before loading, so it is normal
//! and expected for one physical GEM droplet to legitimately hold cells from
//! more than one sample. `assign` therefore emits one `Assigned` row *per
//! sample* whose declared set clears `min_molecules` for a GEM (not a single
//! winner-take-all call), and only discards a GEM outright (`HighOccupancyGem`)
//! when its number of clearing samples exceeds a negative-binomial-fit
//! 99.9th-percentile cutoff on the observed per-GEM occupancy distribution —
//! mirroring 10x's real "High occupancy GEM filtering" step, which estimates
//! NB(r, p=0.5) from the data (r from the observed mean, p fixed at 0.5) and
//! discards only that rare chip-loading-artifact tail (their own worked
//! examples: >7 probe barcodes/GEM for a 16-plex run, >19 for 96-plex).
//! Earlier revision of this module treated *any* second sample clearing the
//! floor as an all-or-nothing reject — real CAT06 data showed that fires on
//! 60%+ of GEMs with enough signal to assign, which cannot be real
//! multiplet/doublet biology at that rate and was discarding the majority of
//! otherwise-real, correctly multiplexed cells.

use std::collections::{BTreeMap, HashMap};
use std::fs::{self, File};
use std::io::{BufWriter, Write as _};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result};
use clap::Parser;
use rust_htslib::bam::{self, record::Aux, Read};
use serde::Serialize;
use tracing::info;

use crate::probe_sample_sheet::ProbeSampleSheet;

/// Assign cell barcodes to samples by probe barcode (10x Flex).
#[derive(Debug, Parser)]
#[command(
    name = "sbxpress-probe-demux",
    version,
    about = "Assign cell barcodes to samples by probe barcode (10x Flex)"
)]
pub struct Cli {
    /// GEX-library BAM after `annotate`+`molecules` (must carry `CB`, `PB`,
    /// and the `xf` representative-molecule bit).
    #[arg(value_name = "BAM")]
    pub bam: PathBuf,

    /// Probe-barcode sample sheet: a real cellranger `multi` `config.csv`
    /// (only `[samples]` is read) or a flat `sample_id,probe_barcode_ids` CSV.
    #[arg(long)]
    pub probe_sample_sheet: PathBuf,

    /// Output directory.
    #[arg(short, long, default_value = "outs")]
    pub output: PathBuf,

    /// Minimum molecules a sample's declared probe-barcode set must capture
    /// for a GEM to be Assigned to that sample. A GEM can independently clear
    /// this floor for more than one sample — see [`assign`]'s module-level
    /// doc comment on composite `(GEM, ProbeBarcode)` cell barcodes.
    #[arg(long, default_value_t = 10)]
    pub min_molecules: u64,
}

/// Tuning for [`assign`].
#[derive(Debug, Clone, Copy)]
pub struct DemuxCfg {
    pub min_molecules: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    /// This GEM independently cleared `min_molecules` for `sample_id`'s
    /// declared set. Ordinary, expected outcome — including when another
    /// sample also has its own `Assigned` row for the same `cell_barcode`
    /// (real Flex multiplex pooling: see the module doc comment).
    Assigned,
    /// Below `min_molecules` total, or no sample's declared set captured
    /// enough of it to call — low-confidence/noisy.
    Unassigned,
    /// This GEM's number of independently-clearing samples exceeds the
    /// negative-binomial 99.9th-percentile cutoff fit to the whole run's
    /// observed per-GEM occupancy distribution — the rare chip-loading-
    /// artifact tail 10x's own "High occupancy GEM filtering" step targets,
    /// not ordinary multi-sample pooling. All of this GEM's molecules are
    /// excluded from every sample. See [`gem_overlap`] for the
    /// whole-experiment version of the underlying overlap diagnostic.
    HighOccupancyGem,
}

impl Status {
    pub fn as_str(&self) -> &'static str {
        match self {
            Status::Assigned => "Assigned",
            Status::Unassigned => "Unassigned",
            Status::HighOccupancyGem => "HighOccupancyGem",
        }
    }
}

/// One cell barcode's sample-assignment result. A single `cell_barcode` can
/// produce more than one `Assigned` row (one per sample whose declared set
/// independently clears the floor) — see the module doc comment.
#[derive(Debug, Clone)]
pub struct Assignment {
    pub cell_barcode: String,
    /// `Some` only when `status == Assigned`.
    pub sample_id: Option<String>,
    /// Distinct probe-barcode ids observed for this CB, sorted.
    pub probe_barcodes: Vec<String>,
    /// This row's own molecule count: `sample_id`'s declared-set sum when
    /// `Assigned`, the whole GEM's total otherwise.
    pub molecules: u64,
    /// This row's sample's molecules / the whole GEM's total molecules
    /// (dominant-sample fraction when not `Assigned`).
    pub purity: f64,
    /// Number of samples whose declared set independently cleared
    /// `min_molecules` for this GEM (10x's own "Probe barcodes per GEM").
    pub occupancy: usize,
    pub status: Status,
}

/// Negative binomial NB(r, p) 99.9th-percentile quantile, via the standard
/// PMF recurrence `P(k) = P(k-1) * p * (k-1+r) / k`, `P(0) = (1-p)^r` — exact
/// for real-valued `r` without a gamma-function dependency. Returns the
/// smallest `k` whose CDF clears `quantile`. `r <= 0.0` (no occupied GEMs to
/// fit against) disables the filter by returning `u64::MAX`.
fn negbinom_quantile(r: f64, p: f64, quantile: f64) -> u64 {
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
///
/// Two passes, mirroring 10x's real "Assign reads to samples" then "Filter
/// high occupancy GEMs" steps (see the module doc comment):
///
/// 1. For each GEM: sum molecule counts across each sample's *whole* declared
///    probe-barcode set, and collect the set of samples whose sum clears
///    `cfg.min_molecules` ("occupancy"). Fit NB(r, p=0.5) to the observed
///    occupancy distribution across all occupied GEMs (r = the observed
///    mean, matching NB's mean = r at p=0.5) and take its 99.9th-percentile
///    quantile as the high-occupancy-GEM cutoff.
/// 2. For each GEM: zero occupancy is `Unassigned`; occupancy above the
///    cutoff is `HighOccupancyGem` (one row, no sample, all molecules
///    excluded); otherwise emit one `Assigned` row *per* clearing sample,
///    each carrying that sample's own molecule sum and purity.
pub fn assign(sheet: &ProbeSampleSheet, tally: &BTreeMap<(String, String), u64>, cfg: &DemuxCfg) -> Vec<Assignment> {
    let mut per_cb: BTreeMap<&str, BTreeMap<&str, u64>> = BTreeMap::new();
    for ((cb, pb), &n) in tally {
        *per_cb.entry(cb.as_str()).or_default().entry(pb.as_str()).or_insert(0) += n;
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
        let occupants: Vec<(&str, u64)> =
            sample_sums.into_iter().filter(|&(_, n)| n >= cfg.min_molecules).collect();

        infos.push(CbInfo { cb, total, probe_barcodes, occupants });
    }

    let occupied: Vec<u64> = infos.iter().filter(|i| !i.occupants.is_empty()).map(|i| i.occupants.len() as u64).collect();
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
            let purity = if info.total > 0 { best_sum as f64 / info.total as f64 } else { 0.0 };
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
                    purity: if info.total > 0 { *sum as f64 / info.total as f64 } else { 0.0 },
                    occupancy,
                    status: Status::Assigned,
                });
            }
        }
    }
    out
}

/// Whole-experiment probe-barcode-pair co-occurrence diagnostic, mirroring
/// cellranger's real `frp_gem_barcode_overlap.csv`: for every pair of probe
/// barcodes with at least one shared cell barcode, `overlap =
/// common_gems / min(barcode1_gems, barcode2_gems)`.
///
/// Restricted to `declared` (a sample sheet's real probe-barcode set):
/// without that filter, a single-read background/error probe-barcode call
/// (common — `tally` is unfiltered raw calls, not whitelist-corrected ids)
/// produces a degenerate `barcode1_gems == 1` row that hits `overlap == 1.0`
/// merely by coinciding with any one real barcode's gem, drowning out the
/// real cross-sample signal this diagnostic exists to surface. Verified
/// against real CAT06 data: on an 800k-read sample, the unfiltered version's
/// top-10-by-overlap rows were *all* such 1-gem noise artifacts.
#[derive(Debug, Clone, Serialize)]
pub struct OverlapRow {
    pub barcode1: String,
    pub barcode2: String,
    pub barcode1_gems: usize,
    pub barcode2_gems: usize,
    pub common_gems: usize,
    pub overlap: f64,
}

pub fn gem_overlap(tally: &BTreeMap<(String, String), u64>, declared: &std::collections::HashSet<&str>) -> Vec<OverlapRow> {
    let mut pbs_per_cb: BTreeMap<&str, std::collections::BTreeSet<&str>> = BTreeMap::new();
    for (cb, pb) in tally.keys() {
        if !declared.contains(pb.as_str()) {
            continue;
        }
        pbs_per_cb.entry(cb.as_str()).or_default().insert(pb.as_str());
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
                let key = if list[i] < list[j] { (list[i], list[j]) } else { (list[j], list[i]) };
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

fn aux_string(record: &bam::Record, tag: &[u8; 2]) -> Result<Option<String>> {
    match record.aux(tag) {
        Ok(Aux::String(value)) => Ok(Some(value.to_owned())),
        Ok(_) => anyhow::bail!("tag {}{} is not a string", tag[0] as char, tag[1] as char),
        Err(_) => Ok(None),
    }
}

fn aux_i64(record: &bam::Record, tag: &[u8; 2]) -> Result<Option<i64>> {
    let value = match record.aux(tag) {
        Ok(Aux::I8(v)) => v as i64,
        Ok(Aux::U8(v)) => v as i64,
        Ok(Aux::I16(v)) => v as i64,
        Ok(Aux::U16(v)) => v as i64,
        Ok(Aux::I32(v)) => v as i64,
        Ok(Aux::U32(v)) => v as i64,
        Ok(_) => anyhow::bail!("tag {}{} is not an integer", tag[0] as char, tag[1] as char),
        Err(_) => return Ok(None),
    };
    Ok(Some(value))
}

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
        let Some(cb) = aux_string(&record, b"CB")? else { continue };
        let Some(pb) = aux_string(&record, b"PB")? else { continue };
        if pb.is_empty() {
            continue;
        }
        *tally.entry((cb, pb)).or_insert(0) += 1;
    }
    Ok(tally)
}

/// Run-level metrics, written to `probe_demux_metrics.json`.
#[derive(Debug, Serialize)]
pub struct Metrics {
    pub elapsed_secs: f64,
    /// Distinct GEM barcodes seen (not row count — a GEM legitimately
    /// present in more than one sample produces more than one row).
    pub cell_barcodes: usize,
    /// `Assigned` rows (one GEM present in 2 samples counts as 2).
    pub assigned: usize,
    pub unassigned: usize,
    /// Distinct GEMs discarded by the negative-binomial high-occupancy-GEM
    /// tail filter — see the module doc comment. Expect this to be small; a
    /// large value here (or `occupancy_cutoff` far below the sheet's sample
    /// count) means the run doesn't fit 10x's own "usually >90% of cells
    /// pass" expectation and is worth a closer look.
    pub high_occupancy_gems: usize,
    /// NB(r, p=0.5) fit's `r` (the observed mean occupancy across occupied
    /// GEMs) and the resulting 99.9th-percentile cutoff actually applied.
    pub occupancy_r: f64,
    pub occupancy_cutoff: u64,
    pub per_sample: Vec<SampleCount>,
}

#[derive(Debug, Serialize)]
pub struct SampleCount {
    pub sample: String,
    pub cell_barcodes: usize,
}

fn write_sample_assignment(assignments: &[Assignment], path: &Path) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w, "cell_barcode,sample_id,probe_barcodes,molecules,purity,occupancy,status")?;
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
fn strip_gem_suffix(cb: &str) -> &str {
    match cb.rfind('-') {
        Some(i) if !cb[i + 1..].is_empty() && cb[i + 1..].bytes().all(|b| b.is_ascii_digit()) => {
            &cb[..i]
        }
        _ => cb,
    }
}

fn write_per_sample_barcodes(sheet: &ProbeSampleSheet, assignments: &[Assignment], out_dir: &Path) -> Result<()> {
    let mut by_sample: HashMap<&str, Vec<&str>> = HashMap::new();
    for a in assignments {
        if let (Status::Assigned, Some(sid)) = (a.status, a.sample_id.as_deref()) {
            by_sample.entry(sid).or_default().push(a.cell_barcode.as_str());
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
fn write_per_sample_probe_barcodes(sheet: &ProbeSampleSheet, out_dir: &Path) -> Result<()> {
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
    writeln!(w, "barcode1,barcode2,barcode1_gems,barcode2_gems,common_gems,overlap")?;
    for r in rows {
        writeln!(
            w,
            "{},{},{},{},{},{:.6}",
            r.barcode1, r.barcode2, r.barcode1_gems, r.barcode2_gems, r.common_gems, r.overlap
        )?;
    }
    Ok(())
}

fn build_metrics(sheet: &ProbeSampleSheet, assignments: &[Assignment], elapsed_secs: f64) -> Metrics {
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
    let cfg = DemuxCfg { min_molecules: cli.min_molecules };
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::probe_sample_sheet::ProbeSample;

    fn sheet() -> ProbeSampleSheet {
        ProbeSampleSheet {
            samples: vec![
                ProbeSample { id: "Pos".into(), probe_barcode_ids: vec!["D-A01".into(), "D-A02".into()] },
                ProbeSample { id: "RZ211".into(), probe_barcode_ids: vec!["D-A03".into(), "D-A04".into()] },
                ProbeSample { id: "RZ213NR".into(), probe_barcode_ids: vec!["D-A07".into()] },
            ],
        }
    }

    fn tally(entries: &[(&str, &str, u64)]) -> BTreeMap<(String, String), u64> {
        entries
            .iter()
            .map(|&(cb, pb, n)| ((cb.to_string(), pb.to_string()), n))
            .collect()
    }

    #[test]
    fn sums_across_a_samples_whole_declared_set() {
        // Pos's declared set is D-A01|D-A02; a CB split across both should
        // sum to a single Pos call, not lose to a naive top-single-barcode vote.
        let t = tally(&[("CB1", "D-A01", 30), ("CB1", "D-A02", 25)]);
        let out = assign(&sheet(), &t, &DemuxCfg { min_molecules: 10 });
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].status, Status::Assigned);
        assert_eq!(out[0].sample_id.as_deref(), Some("Pos"));
        assert_eq!(out[0].molecules, 55);
        assert!((out[0].purity - 1.0).abs() < 1e-9);
    }

    #[test]
    fn below_floor_is_unassigned() {
        let t = tally(&[("CB1", "D-A01", 3)]);
        let out = assign(&sheet(), &t, &DemuxCfg { min_molecules: 10 });
        assert_eq!(out[0].status, Status::Unassigned);
        assert_eq!(out[0].sample_id, None);
    }

    #[test]
    fn two_samples_each_clearing_the_floor_both_get_assigned() {
        // Real Flex multiplex pooling: a GEM legitimately holding cells from
        // two different pooled samples is the expected case, not a doublet
        // to reject outright -- see the module doc comment. With only 3
        // possible samples in this toy sheet, occupancy=2 is nowhere near
        // any realistic high-occupancy-GEM cutoff.
        let t = tally(&[("CB1", "D-A01", 50), ("CB1", "D-A03", 40)]);
        let out = assign(&sheet(), &t, &DemuxCfg { min_molecules: 10 });
        assert_eq!(out.len(), 2);
        let pos = out.iter().find(|a| a.sample_id.as_deref() == Some("Pos")).unwrap();
        assert_eq!(pos.status, Status::Assigned);
        assert_eq!(pos.molecules, 50);
        assert_eq!(pos.occupancy, 2);
        assert!((pos.purity - 50.0 / 90.0).abs() < 1e-9);
        let rz = out.iter().find(|a| a.sample_id.as_deref() == Some("RZ211")).unwrap();
        assert_eq!(rz.status, Status::Assigned);
        assert_eq!(rz.molecules, 40);
        assert!((rz.purity - 40.0 / 90.0).abs() < 1e-9);
    }

    #[test]
    fn a_little_bleed_from_a_second_sample_does_not_trigger_high_occupancy() {
        // second sample present but below the floor -> ordinary Assigned, not
        // a doublet call.
        let t = tally(&[("CB1", "D-A01", 50), ("CB1", "D-A03", 2)]);
        let out = assign(&sheet(), &t, &DemuxCfg { min_molecules: 10 });
        assert_eq!(out[0].status, Status::Assigned);
        assert_eq!(out[0].sample_id.as_deref(), Some("Pos"));
        assert!((out[0].purity - 50.0 / 52.0).abs() < 1e-9);
    }

    #[test]
    fn molecules_on_an_undeclared_probe_barcode_count_against_purity_but_not_toward_any_sample() {
        let t = tally(&[("CB1", "D-A01", 20), ("CB1", "A-F04", 15)]); // A-F04 not declared by anyone
        let out = assign(&sheet(), &t, &DemuxCfg { min_molecules: 10 });
        assert_eq!(out[0].status, Status::Assigned);
        // molecules is this row's own sample sum (20), not the GEM total
        // (35) -- A-F04's 15 molecules aren't Pos's to claim.
        assert_eq!(out[0].molecules, 20);
        assert!((out[0].purity - 20.0 / 35.0).abs() < 1e-9);
    }

    #[test]
    fn negbinom_quantile_matches_hand_computed_values() {
        // r=2.0, p=0.5: cdf crosses 0.999 at k=12 (hand-verified via the
        // same P(k)=P(k-1)*p*(k-1+r)/k recurrence, P(0)=(1-p)^r=0.25).
        assert_eq!(negbinom_quantile(2.0, 0.5, 0.999), 12);
        // r<=0 (no occupied GEMs to fit) disables the filter.
        assert_eq!(negbinom_quantile(0.0, 0.5, 0.999), u64::MAX);
    }

    #[test]
    fn high_occupancy_gem_tail_is_filtered_but_ordinary_multi_sample_gems_are_not() {
        // 40 samples, one probe barcode each. 39 "normal" GEMs each hold
        // only one sample's worth of signal (occupancy=1); one "artifact"
        // GEM holds every sample's probe barcode at once (occupancy=40) --
        // the kind of GEM 10x's real high-occupancy-GEM filter targets
        // (chip clog / wetting failure), not ordinary pooling.
        let samples: Vec<ProbeSample> = (1..=40)
            .map(|i| ProbeSample { id: format!("S{i}"), probe_barcode_ids: vec![format!("P{i}")] })
            .collect();
        let big_sheet = ProbeSampleSheet { samples };

        let mut entries: Vec<(String, String, u64)> = Vec::new();
        for i in 1..=39 {
            entries.push((format!("NORMAL{i}"), format!("P{i}"), 20));
        }
        for i in 1..=40 {
            entries.push(("ARTIFACT".to_string(), format!("P{i}"), 20));
        }
        let t: BTreeMap<(String, String), u64> =
            entries.into_iter().map(|(cb, pb, n)| ((cb, pb), n)).collect();

        let out = assign(&big_sheet, &t, &DemuxCfg { min_molecules: 10 });

        let normal_rows: Vec<&Assignment> = out.iter().filter(|a| a.cell_barcode.starts_with("NORMAL")).collect();
        assert_eq!(normal_rows.len(), 39);
        assert!(normal_rows.iter().all(|a| a.status == Status::Assigned));

        let artifact_rows: Vec<&Assignment> = out.iter().filter(|a| a.cell_barcode == "ARTIFACT").collect();
        assert_eq!(artifact_rows.len(), 1, "the artifact GEM must be filtered as one row, not 40 Assigned rows");
        assert_eq!(artifact_rows[0].status, Status::HighOccupancyGem);
        assert_eq!(artifact_rows[0].sample_id, None);
        assert_eq!(artifact_rows[0].occupancy, 40);
    }

    #[test]
    fn gem_overlap_reports_only_pairs_with_shared_cells() {
        let t = tally(&[
            ("CB1", "D-A01", 10),
            ("CB1", "D-A02", 10),
            ("CB2", "D-A01", 10),
            ("CB3", "D-A03", 10),
        ]);
        let declared: std::collections::HashSet<&str> = ["D-A01", "D-A02", "D-A03"].into();
        let rows = gem_overlap(&t, &declared);
        assert_eq!(rows.len(), 1);
        let r = &rows[0];
        assert_eq!((r.barcode1.as_str(), r.barcode2.as_str()), ("D-A01", "D-A02"));
        assert_eq!(r.barcode1_gems, 2);
        assert_eq!(r.barcode2_gems, 1);
        assert_eq!(r.common_gems, 1);
        assert!((r.overlap - 1.0).abs() < 1e-9); // common / min(2,1) = 1/1
    }

    #[test]
    fn gem_overlap_ignores_undeclared_noise_barcodes() {
        // A-F04 is a single-read background call, not a real declared probe
        // barcode -- without filtering it'd produce a spurious overlap=1.0
        // row against any real barcode it happens to coincide with.
        let t = tally(&[("CB1", "D-A01", 10), ("CB1", "A-F04", 1)]);
        let declared: std::collections::HashSet<&str> = ["D-A01", "D-A02"].into();
        let rows = gem_overlap(&t, &declared);
        assert!(rows.is_empty());
    }

    #[test]
    fn bam_tally_gates_on_xf_bit_8_and_skips_missing_or_empty_pb() {
        let temp = tempfile::tempdir().unwrap();
        let bam_path = temp.path().join("in.bam");
        let mut header = bam::Header::new();
        header.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "guide_contig")
                .push_tag(b"LN", 100),
        );
        let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
        let cigar = bam::record::CigarString(vec![bam::record::Cigar::Match(10)]);

        let mut push = |cb: Option<&str>, pb: Option<&str>, xf: i32| {
            let mut r = bam::Record::new();
            r.set(b"read", Some(&cigar), b"ACGTACGTAC", &[30; 10]);
            r.set_tid(0);
            r.set_pos(0);
            if let Some(cb) = cb {
                r.push_aux(b"CB", Aux::String(cb)).unwrap();
            }
            if let Some(pb) = pb {
                r.push_aux(b"PB", Aux::String(pb)).unwrap();
            }
            r.push_aux(b"xf", Aux::I32(xf)).unwrap();
            writer.write(&r).unwrap();
        };
        push(Some("CB1-1"), Some("D-A01"), 25); // representative (bit 8 set, 25 = 0b11001)
        push(Some("CB1-1"), Some("D-A01"), 17); // NOT representative (17 = 0b10001, no bit 8)
        push(Some("CB1-1"), Some(""), 25); // representative but empty PB -> skipped
        push(Some("CB2-1"), None, 25); // representative but no PB tag -> skipped
        drop(writer);

        let tally = read_bam_tally(&bam_path).unwrap();
        assert_eq!(tally.len(), 1);
        assert_eq!(tally[&("CB1-1".to_string(), "D-A01".to_string())], 1);
    }

    #[test]
    fn strip_gem_suffix_removes_trailing_digits_only() {
        assert_eq!(strip_gem_suffix("ACGTACGTACGTACGT-1"), "ACGTACGTACGTACGT");
        assert_eq!(strip_gem_suffix("ACGTACGTACGTACGT-12"), "ACGTACGTACGTACGT");
        assert_eq!(strip_gem_suffix("ACGTACGTACGTACGT"), "ACGTACGTACGTACGT"); // no suffix
        assert_eq!(strip_gem_suffix("D-A01"), "D-A01"); // trailing non-digit -> untouched
    }

    #[test]
    fn per_sample_barcode_files_are_bare_acgt_matrix_cell_barcodes_ready() {
        // The real bug this guards: matrix --cell-barcodes rejects a barcode
        // that already carries a gem-group suffix (it appends its own), so
        // writing "CB1-1" instead of "CB1" would break the Phase 2 -> 3 handoff.
        let t = tally(&[
            ("AAACCCAAGAAACCCA-1", "D-A01", 50),
            ("TTTGGGTTCTTTGGGT-1", "D-A02", 50),
        ]);
        let assignments = assign(&sheet(), &t, &DemuxCfg { min_molecules: 10 });
        let dir = tempfile::tempdir().unwrap();
        write_per_sample_barcodes(&sheet(), &assignments, dir.path()).unwrap();

        let pos_barcodes = std::fs::read_to_string(dir.path().join("Pos.barcodes.txt")).unwrap();
        let lines: Vec<&str> = pos_barcodes.lines().collect();
        assert_eq!(lines, vec!["AAACCCAAGAAACCCA", "TTTGGGTTCTTTGGGT"]);
        for line in &lines {
            assert!(
                line.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')),
                "barcode {line} is not bare ACGT"
            );
        }
    }
}
