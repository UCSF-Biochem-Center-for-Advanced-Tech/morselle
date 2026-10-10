use anyhow::{bail, Context, Result};
use std::path::{Path, PathBuf};

use crate::anchors::revcomp;
use crate::barcode::{encode16, Whitelist};
use crate::barcode_metrics::{NO_BARCODE, UNRESOLVED_BARCODE};
use crate::complexity::ComplexityCfg;
use crate::extract::{extract, ExtractCfg};
use crate::gzin::open_fastx;
use crate::homopolymer;
use crate::overlap::{overlap_in, InsertClass, Merged, OverlapCfg, Scratch};
use crate::yc::YcTag;

use super::types::{
    canonical_id, AsyncFailedWriter, AsyncWriter, PairMetrics, CB_LEN, UMI_LEN, WRITER_BATCH_BYTES,
};
use super::writer::{write_failed_record, write_record};

#[allow(clippy::too_many_arguments)]
pub fn convert_pair(
    r1_path: &Path,
    r2_path: &Path,
    output_paths: &[PathBuf],
    whitelist: &Whitelist,
    gem_group: u32,
    overlap_cfg: Option<&OverlapCfg>,
    extract_cfg: &ExtractCfg,
    emit_yc: bool,
    r2_polya_min_run: Option<usize>,
    failed_paths: &[(PathBuf, PathBuf)],
) -> Result<PairMetrics> {
    convert_pair_with(
        r1_path,
        r2_path,
        output_paths,
        whitelist,
        gem_group,
        overlap_cfg,
        extract_cfg,
        emit_yc,
        r2_polya_min_run,
        failed_paths,
        false,
    )
}

/// [`convert_pair`] that also fills `PairMetrics::per_barcode` when
/// `barcode_metrics` is set (`--barcode-metrics`).
#[allow(clippy::too_many_arguments)]
pub fn convert_pair_with(
    r1_path: &Path,
    r2_path: &Path,
    output_paths: &[PathBuf],
    whitelist: &Whitelist,
    gem_group: u32,
    overlap_cfg: Option<&OverlapCfg>,
    extract_cfg: &ExtractCfg,
    emit_yc: bool,
    r2_polya_min_run: Option<usize>,
    failed_paths: &[(PathBuf, PathBuf)],
    barcode_metrics: bool,
) -> Result<PairMetrics> {
    convert_pair_full(
        r1_path,
        r2_path,
        output_paths,
        whitelist,
        gem_group,
        overlap_cfg,
        extract_cfg,
        emit_yc,
        r2_polya_min_run,
        failed_paths,
        barcode_metrics,
        &ComplexityCfg::default(),
    )
}

/// [`convert_pair_with`] with explicit low-complexity settings
/// (`--low-complexity-*`, `--drop-low-complexity`).
#[allow(clippy::too_many_arguments)]
pub fn convert_pair_full(
    r1_path: &Path,
    r2_path: &Path,
    output_paths: &[PathBuf],
    whitelist: &Whitelist,
    gem_group: u32,
    overlap_cfg: Option<&OverlapCfg>,
    extract_cfg: &ExtractCfg,
    emit_yc: bool,
    r2_polya_min_run: Option<usize>,
    failed_paths: &[(PathBuf, PathBuf)],
    barcode_metrics: bool,
    complexity: &ComplexityCfg,
) -> Result<PairMetrics> {
    let writers = output_paths
        .iter()
        .map(|p| AsyncWriter::create(p))
        .collect::<Result<Vec<_>>>()?;
    let failed_writers = failed_paths
        .iter()
        .map(|(a, b)| AsyncFailedWriter::create(a, b))
        .collect::<Result<Vec<_>>>()?;
    let mut buffers = (0..writers.len())
        .map(|_| Vec::with_capacity(WRITER_BATCH_BYTES))
        .collect::<Vec<_>>();
    let mut failed_buffers = (0..failed_writers.len())
        .map(|_| {
            (
                Vec::with_capacity(WRITER_BATCH_BYTES),
                Vec::with_capacity(WRITER_BATCH_BYTES),
            )
        })
        .collect::<Vec<_>>();
    let mut metrics = PairMetrics {
        per_barcode: crate::barcode_metrics::BarcodeTable::new(barcode_metrics),
        ..PairMetrics::default()
    };
    let mut r1 =
        open_fastx(r1_path).with_context(|| format!("cannot open R1 {}", r1_path.display()))?;
    let mut r2 =
        open_fastx(r2_path).with_context(|| format!("cannot open R2 {}", r2_path.display()))?;
    let mut pair_number = 0u64;
    let mut scratch = Scratch::default();
    let gem_suffix = format!("-{gem_group}");
    loop {
        let a = r1.next();
        let b = r2.next();
        match (a, b) {
            (None, None) => break,
            (None, Some(_)) | (Some(_), None) => bail!(
                "R1/R2 record-count mismatch between {} and {} after {} pairs",
                r1_path.display(),
                r2_path.display(),
                pair_number
            ),
            (Some(a), Some(b)) => {
                let a = a.with_context(|| format!("invalid FASTQ in {}", r1_path.display()))?;
                let b = b.with_context(|| format!("invalid FASTQ in {}", r2_path.display()))?;
                pair_number += 1;
                metrics.records_read += 1;
                let mut barcode_row = NO_BARCODE;
                if canonical_id(a.id()) != canonical_id(b.id()) {
                    bail!(
                        "R1/R2 name mismatch at pair {}: '{}' versus '{}'",
                        pair_number,
                        String::from_utf8_lossy(a.id()),
                        String::from_utf8_lossy(b.id())
                    );
                }
                let r1_seq = a.seq();
                let r2_seq = b.seq();
                let Some(r1_qual) = a.qual() else {
                    bail!("R1 record {} has no qualities", pair_number)
                };
                let Some(r2_qual) = b.qual() else {
                    bail!("R2 record {} has no qualities", pair_number)
                };

                /// Record a failure, emit the untouched pair to the side channel,
                /// and move on. Every `continue` below goes through this, so
                /// `records_read == records_emitted + records_failed` holds.
                macro_rules! fail_pair {
                    ($reason:expr) => {{
                        let reason: &str = $reason;
                        metrics.records_failed += 1;
                        *metrics
                            .failed_by_reason
                            .entry(reason.to_string())
                            .or_insert(0) += 1;
                        let row = metrics.per_barcode.row(&barcode_row, ".");
                        row[0] += 1;
                        let failure_index = match reason {
                            "merged_too_short" => 4,
                            "no_cdna" => 5,
                            "short_r1" => 6,
                            "r2_polya_short" => 7,
                            "cb_invalid" => 8,
                            "cb_unmatched" => 9,
                            "low_complexity" => 11,
                            _ => unreachable!("new failure reason needs a barcode column"),
                        };
                        row[failure_index] += 1;
                        if !failed_writers.is_empty() {
                            let shard =
                                (metrics.records_failed as usize - 1) % failed_writers.len();
                            write_failed_record(
                                &mut failed_buffers[shard].0,
                                reason,
                                a.id(),
                                &r1_seq,
                                r1_qual,
                            )?;
                            write_failed_record(
                                &mut failed_buffers[shard].1,
                                reason,
                                b.id(),
                                &r2_seq,
                                r2_qual,
                            )?;
                            if failed_buffers[shard].0.len() >= WRITER_BATCH_BYTES {
                                let batch = (
                                    std::mem::take(&mut failed_buffers[shard].0),
                                    std::mem::take(&mut failed_buffers[shard].1),
                                );
                                failed_writers[shard]
                                    .sender
                                    .send(batch)
                                    .map_err(|_| anyhow::anyhow!("failed writer stopped"))?;
                            }
                        }
                        continue;
                    }};
                }

                // Step 1: try to reconstruct the whole fragment from the RAW
                // mates. This must happen before anything is trimmed — the RT
                // primer's dT block is the polyA complement, so R1 and rc(R2)
                // describe the same span only while both are intact.
                let merged: Option<Merged> = overlap_cfg.and_then(|cfg| {
                    overlap_in(&mut scratch, &r1_seq, r1_qual, &r2_seq, r2_qual, cfg)
                });

                // Buffers for the merged path, kept alive across the branch.
                let mut merged_seq: Vec<u8> = Vec::new();
                let mut merged_qual: Vec<u8> = Vec::new();
                let mut merged_tag: Option<YcTag> = None;
                // Non-merged path only: end of the emitted R2 window, trimmed
                // back from r2_seq.len() when --trim-r2-polya finds read-through.
                let mut r2_trim_end = r2_seq.len();

                if let Some(m) = &merged {
                    // (a) sins / mins. The merged molecule is
                    // `CB | UMI | dT | rc(cDNA)` — the same canonical layout the
                    // SBX arm produces — so it goes through the *same* extractor
                    // rather than bespoke logic. Running identical code on both
                    // arms is what keeps the comparison honest.
                    match m.class {
                        InsertClass::Short => metrics.overlap_short += 1,
                        InsertClass::Medium => metrics.overlap_medium += 1,
                        InsertClass::Long => {}
                    }
                    metrics.merged_bases += m.seq.len() as u64;
                    metrics.discordant_bases += m.discordant as u64;
                    metrics.yc_n_columns_coerced += m.yc_stats.n_columns_coerced;
                    if emit_yc && m.yc.is_none() {
                        metrics.yc_unencodable += 1;
                    }
                    if m.seq.len() < CB_LEN + UMI_LEN {
                        metrics.short_r1 += 1;
                        metrics.barcode_rejected += 1;
                        fail_pair!("merged_too_short");
                    }

                    let e = extract(
                        &m.seq,
                        Some((0, CB_LEN + UMI_LEN)),
                        Some(CB_LEN + UMI_LEN),
                        extract_cfg,
                        false, // no Flex support in this binary; see the extract_cfg literal above
                    );
                    let (cdna_start, cdna_end) = match e.r2 {
                        Some(range) => range,
                        None => {
                            metrics.merged_without_cdna += 1;
                            barcode_row = UNRESOLVED_BARCODE;
                            fail_pair!("no_cdna");
                        }
                    };
                    if e.tso_trimmed {
                        metrics.tso_trimmed += 1;
                    }
                    *metrics
                        .insert_length_histogram
                        .entry(cdna_end - cdna_start)
                        .or_insert(0) += 1;

                    // Emit the cDNA reverse-complemented into R2 orientation, so
                    // every downstream stage sees the convention it already
                    // expects. The tag is trimmed to the emitted window and
                    // flipped with it — reversing the string would be wrong.
                    merged_seq = revcomp(&m.seq[cdna_start..cdna_end]);
                    merged_qual = m.qual[cdna_start..cdna_end].iter().rev().copied().collect();
                    if emit_yc {
                        merged_tag = m.yc.clone().map(|mut tag| {
                            tag.trim(cdna_start as u64, (m.seq.len() - cdna_end) as u64);
                            tag.reverse_complement();
                            tag
                        });
                        debug_assert!(merged_tag
                            .as_ref()
                            .is_none_or(|t| t.sequence_length() as usize == merged_seq.len()));
                    }
                } else {
                    // (b) lins, or overlap disabled: the standard Illumina path.
                    if overlap_cfg.is_some() {
                        metrics.overlap_long += 1;
                    }
                    if r1_seq.len() < CB_LEN + UMI_LEN || r1_qual.len() < CB_LEN + UMI_LEN {
                        metrics.short_r1 += 1;
                        metrics.barcode_rejected += 1;
                        fail_pair!("short_r1");
                    }
                    // R2 alone can still read through a short fragment into its
                    // poly-A tail (and beyond, into the reverse-complemented
                    // UMI/CB/primer) even when there is no R1 to merge against.
                    // --overlap catches this by reconstructing the fragment from
                    // both mates; this is the same check for reads that never
                    // get that treatment.
                    if let Some(min_run) = r2_polya_min_run {
                        if let Some(run) = homopolymer::first_run_at_least(&r2_seq, b'A', min_run) {
                            metrics.r2_polya_bases_trimmed += (r2_seq.len() - run.start) as u64;
                            r2_trim_end = run.start;
                            if r2_trim_end < extract_cfg.min_cdna_len {
                                metrics.r2_polya_discarded += 1;
                                barcode_row = UNRESOLVED_BARCODE;
                                fail_pair!("r2_polya_short");
                            }
                            metrics.r2_polya_trimmed += 1;
                        }
                    }
                }

                // The barcode comes from the merged molecule when we have one —
                // on a short fragment R2 read back through it, so CB and UMI
                // carry summed duplex quality.
                let (raw_cb, raw_umi, cb_qual, umi_qual, out_seq, out_qual) = match &merged {
                    Some(m) => (
                        &m.seq[..CB_LEN],
                        &m.seq[CB_LEN..CB_LEN + UMI_LEN],
                        &m.qual[..CB_LEN],
                        &m.qual[CB_LEN..CB_LEN + UMI_LEN],
                        &merged_seq[..],
                        &merged_qual[..],
                    ),
                    None => (
                        &r1_seq[..CB_LEN],
                        &r1_seq[CB_LEN..CB_LEN + UMI_LEN],
                        &r1_qual[..CB_LEN],
                        &r1_qual[CB_LEN..CB_LEN + UMI_LEN],
                        &r2_seq[..r2_trim_end],
                        &r2_qual[..r2_trim_end],
                    ),
                };

                let cb_invalid = encode16(raw_cb).is_none();
                let cb_exact = !cb_invalid && whitelist.contains(raw_cb);
                let corrected = if cb_invalid {
                    None
                } else if cb_exact {
                    metrics.barcode_exact += 1;
                    whitelist.correct_within1(raw_cb)
                } else {
                    let value = whitelist.correct_within1(raw_cb);
                    if value.is_some() {
                        metrics.barcode_hamming1_corrected += 1;
                    }
                    value
                };
                let corrected: [u8; 16] = match corrected {
                    Some(value) => value,
                    None => {
                        metrics.barcode_rejected += 1;
                        barcode_row = UNRESOLVED_BARCODE;
                        fail_pair!(if cb_invalid {
                            "cb_invalid"
                        } else {
                            "cb_unmatched"
                        });
                    }
                };
                let mut corrected_barcode = [0u8; CB_LEN + 11];
                corrected_barcode[..CB_LEN].copy_from_slice(&corrected);
                corrected_barcode[CB_LEN..CB_LEN + gem_suffix.len()]
                    .copy_from_slice(gem_suffix.as_bytes());
                let barcode = std::str::from_utf8(&corrected_barcode[..CB_LEN + gem_suffix.len()])?;
                // Checked after the barcode so the count is per cell; see
                // crate::complexity for why (DZEC05 dead read 2).
                let dominant = complexity.check(out_seq);
                let low_complexity = dominant.is_some();
                if let Some(base) = dominant {
                    metrics.low_complexity += 1;
                    metrics.low_complexity_by_base.0[base] += 1;
                    metrics.per_barcode.row(barcode, ".")[10] += 1;
                    if complexity.drop_low_complexity {
                        barcode_row = barcode;
                        fail_pair!("low_complexity");
                    }
                }
                let row = metrics.per_barcode.row(barcode, ".");
                row[0] += 1;
                row[1] += u64::from(cb_exact);
                row[2] += u64::from(!cb_exact);
                row[3] += 1;

                let shard = metrics.records_emitted as usize % writers.len();
                write_record(
                    &mut buffers[shard],
                    b.id(),
                    raw_cb,
                    raw_umi,
                    cb_qual,
                    umi_qual,
                    &corrected,
                    gem_group,
                    out_seq,
                    out_qual,
                    merged_tag.as_ref(),
                    low_complexity,
                )?;
                metrics.records_emitted += 1;
                if buffers[shard].len() >= WRITER_BATCH_BYTES {
                    writers[shard]
                        .sender
                        .send(std::mem::take(&mut buffers[shard]))
                        .map_err(|_| anyhow::anyhow!("output writer stopped"))?;
                }
            }
        }
    }
    for (i, buf) in buffers.iter_mut().enumerate() {
        if !buf.is_empty() {
            writers[i]
                .sender
                .send(std::mem::take(buf))
                .map_err(|_| anyhow::anyhow!("output writer stopped"))?;
        }
    }
    for (i, (a, b)) in failed_buffers.iter_mut().enumerate() {
        if !a.is_empty() {
            failed_writers[i]
                .sender
                .send((std::mem::take(a), std::mem::take(b)))
                .map_err(|_| anyhow::anyhow!("failed writer stopped"))?;
        }
    }
    for w in writers {
        w.finish()?;
    }
    for w in failed_writers {
        w.finish()?;
    }
    Ok(metrics)
}
