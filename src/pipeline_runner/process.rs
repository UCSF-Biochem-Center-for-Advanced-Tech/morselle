use anyhow::Result;

use crate::demux::{index_anchor_starts, DemuxVia, Outcome};
use crate::extract::extract;
use crate::orient::{orient, OrientCfg};
use crate::pipeline::{CbStatus, Orientation, ReadClass};

use super::barcode_metrics::{bump, observe};
use super::formatting::{
    append_probe_tags, append_processing_tags, flush_failed_buf, minimap2_header, push_fastq,
    push_fastq_r2_rc, FLUSH_BYTES, MIN_POLYT_FLANK,
};
use super::relock::relock_r1;
use super::types::{ProcessingTags, Stats};
use super::worker::Shared;

/// Process one read: orient → extract → CB classify/rescue → gate → demux → emit
/// into `bufs` (flushing the crossing batch to `chunk`), accumulating into `stats`.
/// `orient_cfg` carries the enforced per-file orientation. Shared by both the
/// one-file-per-thread path ([`process_file`]) and the pipelined path
/// ([`worker_loop`]).
#[allow(clippy::too_many_arguments)]
pub(crate) fn process_one(
    id: &[u8],
    seq: &[u8],
    qual: &[u8],
    orient_cfg: &OrientCfg,
    sh: &Shared,
    bufs: &mut [(Vec<u8>, Vec<u8>, u64)],
    failed_buf: &mut Vec<u8>,
    chunk: usize,
    stats: &mut Stats,
) -> Result<()> {
    stats.reads += 1;

    {
        let o = orient(seq, qual, orient_cfg);
        if o.orientation != Orientation::Unknown {
            stats.oriented += 1;
        }

        // Early, i5-ONLY library-type hint -- routes ExtractCfg::flex_r2_boundary's
        // hard-gate default away from libraries the sample sheet's own i5 index
        // confidently says aren't Flex GEX (e.g. CRISPR/Feature Barcoding sharing
        // this invocation). Deliberately separate from `early_outcome` below: i5
        // alone needs only `r1s` (`index_anchor_starts`'s i5 arm), so unlike the
        // full i5+i7 dual lookup -- which needs `e.adapter_pos`, i.e. extract()'s
        // own output -- this can run *before* extract() decides which boundary
        // strategy to use. A read this doesn't confidently resolve (no sheet, or
        // i5 itself unreadable) defaults to `false` -- `--flex` assumes every read
        // is Flex GEX unless the sheet affirmatively says otherwise.
        let likely_non_gex = sh
            .demuxer
            .and_then(|d| {
                let r1s = o.cb_umi.map(|(s, _)| s)?;
                let (i5s, _) = index_anchor_starts(r1s, None, sh.index_len);
                let i5_region = i5s.and_then(|s| o.seq.get(s..s + sh.index_len));
                d.assign(i5_region, None).sample
            })
            .map(|idx| !sh.gex_by_sample[idx])
            .unwrap_or(false);
        let e = extract(
            &o.seq,
            o.cb_umi,
            o.polyt_start,
            &sh.extract_cfg,
            likely_non_gex,
        );

        // Demux (early): i5/i7 assignment the moment R1 is located, before
        // any CB-specific work and independent of whether R2/cDNA is ever
        // found. Both index anchors depend only on `e.r1`/`e.adapter_pos`
        // (both already set by `extract()` above), never on CB classification
        // or gating -- so this is the earliest point the assignment can run,
        // and it's cached (not recomputed) for the actual routing decision
        // further down, so output is unaffected by moving it here.
        let early_detail = sh.demuxer.and_then(|d| {
            let r1s = e.r1.map(|(s, _)| s)?;
            let (i5s, i7s) = index_anchor_starts(r1s, e.adapter_pos, sh.index_len);
            Some(d.assign_windowed_detail(&o.seq, i5s, i7s, sh.index_window))
        });
        if let Some((outcome, ev)) = &early_detail {
            stats
                .index_evidence
                .observe(ev, outcome.via == DemuxVia::Dual);
        }
        let early_outcome: Option<Outcome> = early_detail.map(|(o, _)| o);
        if sh.demuxer.is_some() {
            match e.r1 {
                None => stats.index_unreachable += 1,
                Some(_) => {
                    if let Some(o) = early_outcome {
                        match o.via {
                            DemuxVia::Dual => stats.via_dual_early += 1,
                            DemuxVia::I5Only => stats.via_i5_only_early += 1,
                            DemuxVia::I7Only => stats.via_i7_only_early += 1,
                            DemuxVia::Conflict => stats.via_conflict_early += 1,
                            DemuxVia::None => {}
                        }
                    }
                    match early_outcome.and_then(|o| o.sample) {
                        Some(i) => {
                            if i < stats.per_sample_early.len() {
                                stats.per_sample_early[i] += 1;
                                if e.adapter_trimmed
                                    && i < stats.per_sample_early_boundary_found.len()
                                {
                                    stats.per_sample_early_boundary_found[i] += 1;
                                }
                            }
                        }
                        None => stats.undetermined_early += 1,
                    }
                }
            }
        }

        // Classify the cell barcode (CB = first cb_len bp of R1) if a whitelist.
        let cb_status = match (sh.whitelist, e.r1) {
            (Some(w), Some((r1s, _))) => Some(w.classify(&o.seq[r1s..r1s + sh.cb_len])),
            _ => None,
        };
        match cb_status {
            Some(CbStatus::OnWhitelist) => stats.cb_on_whitelist += 1,
            Some(CbStatus::Within1) => stats.cb_within1 += 1,
            Some(CbStatus::Fail) => stats.cb_fail += 1,
            _ => {}
        }

        // Location confidence: the 16 bp CB is definitively bracketed when both
        // flanks are present — the R1-primer tail (5') is lock-confirmed and the
        // dT-primer run (3') is long enough. Then a whitelist miss is a true base
        // error, not a mis-location.
        let located_confident =
            o.r1_tail_confirmed && e.polyt_run_len >= MIN_POLYT_FLANK && e.r1.is_some();
        if cb_status.is_some() && located_confident {
            stats.cb_located_confident += 1;
        }
        if cb_status == Some(CbStatus::Fail) && located_confident {
            stats.cb_fail_confident += 1;
        }

        // CB rescue (synthesizes a corrected R1). M5 indel realign (--realign);
        // then CB reframe (--cb-reframe); then lenient Hamming-2 (--lenient).
        // Reframe and lenient are gated on confident location.
        // A unique Hamming-1 match is a correction, not merely a classification:
        // emit the canonical whitelist sequence in R1 so downstream CB tags and
        // fixed barcode templates refer to the same cell identifier.
        let mut corrected_r1: Option<(Vec<u8>, Vec<u8>)> = match (sh.whitelist, e.r1, cb_status) {
            (Some(w), Some((r1s, r1e)), Some(CbStatus::Within1)) => {
                w.correct_within1(&o.seq[r1s..r1s + sh.cb_len]).map(|cb| {
                    let mut s = cb.to_vec();
                    s.extend_from_slice(&o.seq[r1s + sh.cb_len..r1e]);
                    (s, o.qual[r1s..r1e].to_vec())
                })
            }
            _ => None,
        };

        // XB: how the emitted barcode passed or was rescued. Values 0/1 are
        // whitelist exact/within-1; 2/3/4 are indel/reframe/Hamming-2 rescue;
        // 5 means no whitelist and 6 means an unrescued miss emitted with gating off;
        // 7 is the R1-tail re-lock rescue.
        let mut barcode_route = match cb_status {
            Some(CbStatus::OnWhitelist) => 0,
            Some(CbStatus::Within1) => 1,
            Some(CbStatus::Fail) => 6,
            _ => 5,
        };
        if cb_status == Some(CbStatus::Fail) {
            if let (Some(w), Some((r1s, r1e))) = (sh.whitelist, e.r1) {
                if sh.realign {
                    if let Some(pts) = o.polyt_start {
                        let delta = (pts as i64 - r1s as i64) - (sh.cb_len + sh.umi_len) as i64;
                        if (delta == 1 || delta == -1) && pts >= sh.umi_len && pts <= o.seq.len() {
                            if let Some(cb) = w.realign_cb(&o.seq[r1s..], delta as i32) {
                                let mut s = cb.to_vec();
                                s.extend_from_slice(&o.seq[pts - sh.umi_len..pts]);
                                let mut q =
                                    o.qual[r1s..(r1s + sh.cb_len).min(o.qual.len())].to_vec();
                                q.resize(sh.cb_len, b'I');
                                q.extend_from_slice(&o.qual[pts - sh.umi_len..pts]);
                                corrected_r1 = Some((s, q));
                                barcode_route = 2;
                                stats.realigned += 1;
                            }
                        }
                    }
                }
                // CB reframe rescue: a homopolymer indel shifted the 16 bp CB grab;
                // try nearby frames (+1 first) for a unique within-1 whitelist hit,
                // re-registering the UMI off the dT anchor (so R1 stays 28 bp).
                // Gated on confident location to keep precision.
                if corrected_r1.is_none() && sh.cb_reframe && located_confident {
                    if let Some(pts) = o.polyt_start {
                        if pts >= sh.umi_len {
                            if let Some((cb, off)) = w.reframe_cb(&o.seq, r1s, sh.reframe_window) {
                                let cb_s = (r1s as i64 + off).max(0) as usize;
                                let mut s = cb.to_vec();
                                s.extend_from_slice(&o.seq[pts - sh.umi_len..pts]);
                                let mut q =
                                    o.qual[cb_s..(cb_s + sh.cb_len).min(o.qual.len())].to_vec();
                                q.resize(sh.cb_len, b'I');
                                q.extend_from_slice(&o.qual[pts - sh.umi_len..pts]);
                                corrected_r1 = Some((s, q));
                                barcode_route = 3;
                                stats.reframed += 1;
                            }
                        }
                    }
                }
                // R1-tail re-lock: an indel in the R1-primer tail can leave the
                // mismatch-only lock a base or two off; re-find the tail with
                // indels allowed and read the CB + UMI after it.
                if corrected_r1.is_none() && sh.r1_tail_relock {
                    if let Some(r1) = relock_r1(w, &o.seq, &o.qual, r1s, sh.cb_len, sh.umi_len) {
                        corrected_r1 = Some(r1);
                        barcode_route = 7;
                        stats.relocked += 1;
                    }
                }
                // Lenient within-2 substitution rescue: only when the CB is
                // confidently located (so the 2 mismatches are real base errors),
                // keeping the original (correctly-positioned) UMI.
                if corrected_r1.is_none() && sh.lenient && located_confident {
                    if let Some(cb2) = w.correct_within2(&o.seq[r1s..r1s + sh.cb_len]) {
                        let mut s = cb2.to_vec();
                        s.extend_from_slice(&o.seq[r1s + sh.cb_len..r1e]);
                        corrected_r1 = Some((s, o.qual[r1s..r1e].to_vec()));
                        barcode_route = 4;
                        stats.corrected2 += 1;
                    }
                }
            }
        }
        let cb_ok = corrected_r1.is_some() || !matches!(cb_status, Some(CbStatus::Fail));

        let class = classify_read(o.orientation, e.r1.is_some(), e.r2.is_some(), cb_ok);
        match class {
            ReadClass::Full => stats.class_full += 1,
            ReadClass::Partial => stats.class_partial += 1,
            ReadClass::Unparseable => stats.class_unparseable += 1,
        }
        let sample = early_outcome
            .and_then(|outcome| outcome.sample)
            .and_then(|i| sh.demuxer.and_then(|d| d.sample_ids().get(i)))
            .map(String::as_str)
            .unwrap_or(".");
        let barcode_row = observe(
            sh.barcodes,
            &o.seq,
            e.r1.map(|(s, _)| s),
            corrected_r1.as_ref().map(|(s, _)| s.as_slice()),
            sh.cb_len,
            cb_status,
            sh.cb_gate,
            class,
            sample,
        );

        // Emit the ORIGINAL read (pre-orientation, untrimmed) to the side
        // channel on every path that does not produce a pair, so no input read
        // is silently dropped.
        let mut push_failed = |reason: &str, stats: &mut Stats| -> Result<()> {
            if sh.failed.is_none() {
                return Ok(());
            }
            stats.failed_emitted += 1;
            failed_buf.push(b'@');
            failed_buf.extend_from_slice(id);
            failed_buf.extend_from_slice(format!("\tXF:Z:{reason}").as_bytes());
            failed_buf.push(b'\n');
            failed_buf.extend_from_slice(seq);
            failed_buf.extend_from_slice(b"\n+\n");
            failed_buf.extend_from_slice(qual);
            failed_buf.push(b'\n');
            if failed_buf.len() >= FLUSH_BYTES {
                flush_failed_buf(failed_buf, sh)?;
            }
            Ok(())
        };

        let (Some((r1s, r1e)), Some((r2s, r2e))) = (e.r1, e.r2) else {
            // need both R1 and R2 to emit a pair
            let reason = if o.orientation == Orientation::Unknown {
                "unoriented"
            } else if e.r1.is_none() {
                "no_barcode"
            } else {
                "no_cdna"
            };
            push_failed(reason, stats)?;
            return Ok(());
        };

        // Gate: drop CB-fail pairs that weren't rescued (unless gating is off).
        if sh.cb_gate && !cb_ok {
            stats.gated_out += 1;
            push_failed("cb_fail", stats)?;
            return Ok(());
        }

        // Low-complexity cDNA (see crate::complexity): counted per cell and
        // tagged YL:i:1, or failed with --drop-low-complexity.
        let dominant = sh.complexity.check(&o.seq[r2s..r2e]);
        let low_complexity = dominant.is_some();
        if let Some(base) = dominant {
            stats.low_complexity += 1;
            stats.low_complexity_by_base.0[base] += 1;
            bump(sh.barcodes, &barcode_row, sample, 9);
            if sh.complexity.drop_low_complexity {
                push_failed("low_complexity", stats)?;
                return Ok(());
            }
        }

        let (r1_seq, r1_qual): (&[u8], &[u8]) = match &corrected_r1 {
            Some((s, q)) => (s, q),
            None => (&o.seq[r1s..r1e], &o.qual[r1s..r1e]),
        };
        if r1_seq.len() == 28 {
            stats.r1_eq_28 += 1;
        } else {
            stats.r1_not_28 += 1;
        }

        // Route to a target index: sample, Undetermined (= n_samples), or combined.
        // Reuses `early_outcome` (computed right after `extract()`, above) --
        // same inputs (`r1s`/`e.adapter_pos`/`o.seq` are unchanged since
        // then), so this is the identical result, not a fresh decision.
        let mut index_route = 5; // XI=5: demultiplexing was not requested.
        let mut sample_assigned = false;
        let target_idx = if sh.demuxer.is_some() {
            let outcome =
                early_outcome.expect("demuxer set and R1 located => early_outcome computed");
            index_route = match outcome.via {
                DemuxVia::Dual => 0,
                DemuxVia::I5Only => 1,
                DemuxVia::I7Only => 2,
                DemuxVia::Conflict => 3,
                DemuxVia::None => 4,
            };
            match outcome.via {
                DemuxVia::Dual => stats.via_dual += 1,
                DemuxVia::I5Only => stats.via_i5 += 1,
                DemuxVia::I7Only => stats.via_i7 += 1,
                DemuxVia::Conflict => stats.via_conflict += 1,
                DemuxVia::None => {}
            }
            match outcome.sample {
                Some(i) => {
                    stats.per_sample[i] += 1;
                    sample_assigned = true;
                    i
                }
                None => {
                    stats.undetermined += 1;
                    sh.n_samples
                }
            }
        } else {
            0
        };

        // R1 (CB+UMI) and R2 (cDNA) sit on the SAME strand of the canonical read,
        // but Illumina paired-end mates come from OPPOSITE strands. The first-strand
        // cDNA following polyT is the reverse complement of the mRNA sense, so the
        // canonical R2 slice maps antisense. Emit revcomp(R2) with reversed quality
        // to match the Illumina convention (R1 as-is, R2 flipped) — otherwise Cell
        // Ranger reports ~50% "Reads Mapped Antisense to Gene". The `--max-r2-len`
        // cap is applied here (at write, on the reverse-complemented read): the R2
        // primer reads from the adapter-proximal end, so keeping the first N bases
        // mirrors a sequencer stopping after N cycles. `push_fastq_r2_rc` writes the
        // revcomp straight into the buffer, touching only the emitted (≤N) bases.
        let buf = &mut bufs[target_idx];
        if !sh.minimap2 {
            push_fastq(&mut buf.0, id, r1_seq, r1_qual);
        }
        let raw_cb = &o.seq[r1s..r1s + sh.cb_len];
        let raw_umi = &o.seq[r1s + sh.cb_len..r1e];
        let raw_cb_qual = &o.qual[r1s..r1s + sh.cb_len];
        let mut header;
        let processing = ProcessingTags {
            orientation: match o.orientation {
                Orientation::P5Led => b'5',
                Orientation::P7Led => b'7',
                Orientation::Unknown => b'?',
            },
            barcode_route,
            index_route,
            // XA bits: 1=R1-primer tail, 2=poly-T anchor, 4=cDNA adapter
            // trimmed, 8=file-level expected/detected orientation applied.
            anchor_bits: (o.r1_tail_confirmed as u8)
                | (((e.polyt_run_len >= MIN_POLYT_FLANK) as u8) << 1)
                | ((e.adapter_trimmed as u8) << 2)
                | ((orient_cfg.expected.is_some() as u8) << 3),
            extracted_cdna_len: r2e - r2s,
        };
        let r2_id = if sh.minimap2 {
            let mut corrected_cb = r1_seq[..sh.cb_len].to_vec();
            corrected_cb.extend_from_slice(b"-1");
            header = minimap2_header(
                id,
                raw_cb,
                &corrected_cb,
                raw_umi,
                &r1_seq[sh.cb_len..],
                raw_cb_qual,
                &r1_qual[sh.cb_len..],
            );
            if sh.processing_tags {
                append_processing_tags(&mut header, processing);
            }
            if let Some(wl) = sh.probe_wl {
                if let Some(hit) = e.polyt_end.and_then(|pts| wl.locate(&o.seq, pts)) {
                    stats.probe_located += 1;
                    bump(sh.barcodes, &barcode_row, sample, 7);
                    if hit.id.is_some() {
                        stats.probe_hit += 1;
                        bump(sh.barcodes, &barcode_row, sample, 8);
                        if sample_assigned {
                            stats.probe_hit_sample_assigned += 1;
                        }
                    }
                    append_probe_tags(&mut header, &hit.raw, hit.id.as_deref());
                }
            }
            if low_complexity {
                header.extend_from_slice(crate::complexity::TAG);
            }
            &header
        } else if sh.processing_tags {
            header = id.to_vec();
            append_processing_tags(&mut header, processing);
            if low_complexity {
                header.extend_from_slice(crate::complexity::TAG);
            }
            &header
        } else {
            id
        };
        let emitted_r2 = push_fastq_r2_rc(
            &mut buf.1,
            r2_id,
            &o.seq,
            &o.qual,
            r2s,
            r2e,
            sh.extract_cfg.max_r2_len,
        );
        buf.2 += 1;
        if buf.0.len().max(buf.1.len()) >= FLUSH_BYTES {
            sh.targets[target_idx * sh.split + chunk].append(&buf.0, &buf.1, buf.2)?;
            buf.0.clear();
            buf.1.clear();
            buf.2 = 0;
        }

        stats.pairs_emitted += 1;
        bump(sh.barcodes, &barcode_row, sample, 3);
        if e.adapter_trimmed {
            stats.adapter_trimmed += 1;
        }
        let len = emitted_r2 as u64; // emitted (post-cap) R2 length
        stats.r2_len_sum += len;
        stats.r2_len_max = stats.r2_len_max.max(len);
    }
    Ok(())
}

pub(crate) fn classify_read(
    orientation: Orientation,
    has_r1: bool,
    has_r2: bool,
    cb_ok: bool,
) -> ReadClass {
    if orientation == Orientation::Unknown || !has_r1 {
        ReadClass::Unparseable
    } else if has_r2 && cb_ok {
        ReadClass::Full
    } else {
        ReadClass::Partial
    }
}
