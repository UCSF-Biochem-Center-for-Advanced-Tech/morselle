//! Assign cell barcodes to samples by their dominant declared probe-barcode
//! set (FLEX_SUPPORT_PLAN.md Phase 2).
//!
//! Published-algorithm parity notes (see FLEX_SUPPORT_PLAN.md "Where this
//! plan corrects the initial design"): assignment sums molecule counts across
//! a sample's *whole declared set* (which can be >1 probe barcode — real
//! `config.csv` example: `Pos,D-A01|D-A02`) before comparing samples, not a
//! naive top-single-barcode vote.

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
    /// for a cell barcode to be Assigned (also the floor a *second* sample
    /// must clear for HighOccupancy — see [`assign`]).
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
    Assigned,
    /// Below `min_molecules` total, or no sample's declared set captured
    /// enough of it to call — low-confidence/noisy, not a doublet signature.
    Unassigned,
    /// A *second* sample's declared set also independently clears
    /// `min_molecules` for this cell barcode — a real cross-sample doublet
    /// signature (two samples' probe barcodes both genuinely present in the
    /// same droplet), distinct from ordinary background noise. See
    /// [`gem_overlap`] for the whole-experiment version of this diagnostic.
    HighOccupancy,
}

impl Status {
    pub fn as_str(&self) -> &'static str {
        match self {
            Status::Assigned => "Assigned",
            Status::Unassigned => "Unassigned",
            Status::HighOccupancy => "HighOccupancy",
        }
    }
}

/// One cell barcode's sample-assignment result.
#[derive(Debug, Clone)]
pub struct Assignment {
    pub cell_barcode: String,
    /// `Some` only when `status == Assigned`.
    pub sample_id: Option<String>,
    /// Distinct probe-barcode ids observed for this CB, sorted.
    pub probe_barcodes: Vec<String>,
    /// Total molecules across all probe barcodes observed for this CB.
    pub molecules: u64,
    /// Dominant sample's molecules / total molecules.
    pub purity: f64,
    pub status: Status,
}

/// Assign every cell barcode in `tally` (keyed `(cell_barcode, probe_barcode_id)
/// -> molecule count`) to a sample from `sheet`.
///
/// For each CB: sum molecule counts across each sample's *whole* declared
/// probe-barcode set, rank samples by that sum, and call the top sample
/// Assigned when its sum clears `cfg.min_molecules`. HighOccupancy fires when
/// a *second* sample's sum also independently clears the same floor — using
/// one symmetric threshold for "a sample is really present" keeps the rule
/// simple and self-contained per CB (the separate [`gem_overlap`] diagnostic
/// covers the whole-experiment probe-barcode-pair view, e.g. for comparison
/// against a real `frp_gem_barcode_overlap.csv`).
pub fn assign(sheet: &ProbeSampleSheet, tally: &BTreeMap<(String, String), u64>, cfg: &DemuxCfg) -> Vec<Assignment> {
    let mut per_cb: BTreeMap<&str, BTreeMap<&str, u64>> = BTreeMap::new();
    for ((cb, pb), &n) in tally {
        *per_cb.entry(cb.as_str()).or_default().entry(pb.as_str()).or_insert(0) += n;
    }

    let mut out = Vec::with_capacity(per_cb.len());
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

        let (best_id, best_sum) = sample_sums.first().copied().unwrap_or(("", 0));
        let second_sum = sample_sums.get(1).map_or(0, |&(_, n)| n);
        let purity = if total > 0 { best_sum as f64 / total as f64 } else { 0.0 };

        let (status, sample_id) = if best_sum >= cfg.min_molecules && second_sum >= cfg.min_molecules {
            (Status::HighOccupancy, None)
        } else if best_sum >= cfg.min_molecules {
            (Status::Assigned, Some(best_id.to_string()))
        } else {
            (Status::Unassigned, None)
        };

        out.push(Assignment {
            cell_barcode: cb.to_string(),
            sample_id,
            probe_barcodes,
            molecules: total,
            purity,
            status,
        });
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
    pub cell_barcodes: usize,
    pub assigned: usize,
    pub unassigned: usize,
    pub high_occupancy: usize,
    pub per_sample: Vec<SampleCount>,
}

#[derive(Debug, Serialize)]
pub struct SampleCount {
    pub sample: String,
    pub cell_barcodes: usize,
}

fn write_sample_assignment(assignments: &[Assignment], path: &Path) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w, "cell_barcode,sample_id,probe_barcodes,molecules,purity,status")?;
    for a in assignments {
        writeln!(
            w,
            "{},{},{},{},{:.4},{}",
            a.cell_barcode,
            a.sample_id.as_deref().unwrap_or(""),
            a.probe_barcodes.join("|"),
            a.molecules,
            a.purity,
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
    let mut assigned = 0;
    let mut unassigned = 0;
    let mut high_occupancy = 0;
    for a in assignments {
        match a.status {
            Status::Assigned => {
                assigned += 1;
                if let Some(sid) = &a.sample_id {
                    *per_sample_counts.entry(sid.as_str()).or_insert(0) += 1;
                }
            }
            Status::Unassigned => unassigned += 1,
            Status::HighOccupancy => high_occupancy += 1,
        }
    }
    Metrics {
        elapsed_secs,
        cell_barcodes: assignments.len(),
        assigned,
        unassigned,
        high_occupancy,
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
    write_gem_overlap(&overlap, &cli.output.join("gem_barcode_overlap.csv"))?;

    let metrics = build_metrics(&sheet, &assignments, started.elapsed().as_secs_f64());
    let metrics_path = cli.output.join("probe_demux_metrics.json");
    fs::write(&metrics_path, serde_json::to_string_pretty(&metrics)?)
        .with_context(|| format!("writing {}", metrics_path.display()))?;

    info!(
        cell_barcodes = metrics.cell_barcodes,
        assigned = metrics.assigned,
        unassigned = metrics.unassigned,
        high_occupancy = metrics.high_occupancy,
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
    fn two_samples_each_clearing_the_floor_is_high_occupancy() {
        let t = tally(&[("CB1", "D-A01", 50), ("CB1", "D-A03", 40)]);
        let out = assign(&sheet(), &t, &DemuxCfg { min_molecules: 10 });
        assert_eq!(out[0].status, Status::HighOccupancy);
        assert_eq!(out[0].sample_id, None);
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
        assert_eq!(out[0].molecules, 35);
        assert!((out[0].purity - 20.0 / 35.0).abs() < 1e-9);
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
