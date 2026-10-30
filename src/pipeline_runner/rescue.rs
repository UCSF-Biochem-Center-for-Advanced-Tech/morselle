//! Cell-barcode rescue chain for a whitelist miss (moved out of
//! `process.rs` unchanged): Hamming-1 correction, then indel realign
//! (`--realign`), CB reframe (`--cb-reframe`), R1-tail re-lock
//! (`--r1-tail-relock`), CB insertion (`--cb-insertion`), lenient Hamming-2
//! (`--lenient`). Each rescue builds
//! a corrected 28 bp R1; the first that succeeds wins. Every corrected R1
//! takes its UMI right after the corrected CB as it sits in the read
//! ([`umi_start`]).

use crate::extract::Extracted;
use crate::orient::Oriented;
use crate::pipeline::CbStatus;

use super::relock::relock_r1;
use super::types::Stats;
use super::umi_anchor::umi_start;
use super::worker::Shared;

/// What the rescue chain decided for one read.
pub(crate) struct Rescue {
    /// Corrected R1 (CB + UMI) and qualities, when a rescue (or Hamming-1
    /// correction) succeeded.
    pub r1: Option<(Vec<u8>, Vec<u8>)>,
    /// The `XB` route.
    pub route: u8,
    /// Re-located CB start, when the R1-tail re-lock rescued the read.
    pub relocked_at: Option<usize>,
}

/// Corrected CB (+ its qualities) and the UMI after the CB's end in the read
/// (see [`umi_start`]); `None` if the read ends before the UMI does.
fn r1_with_umi(
    o: &Oriented,
    cb: &[u8],
    cb_q: &[u8],
    cb_start: usize,
    default_end: usize,
    umi_len: usize,
) -> Option<(Vec<u8>, Vec<u8>)> {
    let u = umi_start(&o.seq, cb, cb_start, default_end, umi_len);
    let end = u + umi_len;
    if end > o.seq.len() || end > o.qual.len() {
        return None;
    }
    let mut s = cb.to_vec();
    s.extend_from_slice(&o.seq[u..end]);
    let mut q = cb_q.to_vec();
    q.extend_from_slice(&o.qual[u..end]);
    Some((s, q))
}

/// Run the chain. `located_confident`: both CB flanks found (see
/// `process_one`); reframe and lenient require it.
pub(crate) fn rescue_barcode(
    sh: &Shared,
    o: &Oriented,
    e: &Extracted,
    cb_status: Option<CbStatus>,
    located_confident: bool,
    stats: &mut Stats,
) -> Rescue {
    let mut relocked_at = None;
    // CB rescue (synthesizes a corrected R1). M5 indel realign (--realign);
    // then CB reframe (--cb-reframe); then lenient Hamming-2 (--lenient).
    // Reframe and lenient are gated on confident location.
    // A unique Hamming-1 match is a correction, not merely a classification:
    // emit the canonical whitelist sequence in R1 so downstream CB tags and
    // fixed barcode templates refer to the same cell identifier.
    let mut corrected_r1: Option<(Vec<u8>, Vec<u8>)> = match (sh.whitelist, e.r1, cb_status) {
        (Some(w), Some((r1s, r1e)), Some(CbStatus::Within1)) => {
            w.correct_within1(&o.seq[r1s..r1s + sh.cb_len]).map(|cb| {
                let cb_q = &o.qual[r1s..r1s + sh.cb_len];
                r1_with_umi(o, &cb, cb_q, r1s, r1s + sh.cb_len, sh.umi_len).unwrap_or_else(|| {
                    let mut s = cb.to_vec();
                    s.extend_from_slice(&o.seq[r1s + sh.cb_len..r1e]);
                    (s, o.qual[r1s..r1e].to_vec())
                })
            })
        }
        _ => None,
    };

    // XB: how the emitted barcode passed or was rescued. Values 0/1 are
    // whitelist exact/within-1; 2/3/4 are indel/reframe/Hamming-2 rescue;
    // 5 means no whitelist and 6 means an unrescued miss emitted with gating off;
    // 7 is the R1-tail re-lock rescue, 8 the CB insertion rescue.
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
                            let mut q = o.qual[r1s..(r1s + sh.cb_len).min(o.qual.len())].to_vec();
                            q.resize(sh.cb_len, b'I');
                            let end = (r1s + sh.cb_len) as i64 + delta;
                            corrected_r1 = r1_with_umi(o, &cb, &q, r1s, end as usize, sh.umi_len);
                            if corrected_r1.is_some() {
                                barcode_route = 2;
                                stats.realigned += 1;
                            }
                        }
                    }
                }
            }
            // CB reframe rescue: a homopolymer indel shifted the 16 bp CB grab;
            // try nearby frames (+1 first) for a unique within-1 whitelist hit.
            // Gated on confident location. Off by default: against Illumina
            // molecules (CAT02 M2_GEX) only ~51-59% of these barcodes are right.
            if corrected_r1.is_none() && sh.cb_reframe && located_confident {
                if let Some(pts) = o.polyt_start {
                    if pts >= sh.umi_len {
                        if let Some((cb, off)) = w.reframe_cb(&o.seq, r1s, sh.reframe_window) {
                            let cb_s = (r1s as i64 + off).max(0) as usize;
                            let mut q = o.qual[cb_s..(cb_s + sh.cb_len).min(o.qual.len())].to_vec();
                            q.resize(sh.cb_len, b'I');
                            let end = cb_s + sh.cb_len;
                            corrected_r1 = r1_with_umi(o, &cb, &q, cb_s, end, sh.umi_len);
                            if corrected_r1.is_some() {
                                barcode_route = 3;
                                stats.reframed += 1;
                            }
                        }
                    }
                }
            }
            // R1-tail re-lock: an indel in the R1-primer tail can leave the
            // mismatch-only lock a base or two off; re-find the tail with
            // indels allowed and read the CB + UMI after it.
            if corrected_r1.is_none() && sh.r1_tail_relock {
                if let Some((s, q, cb_s)) =
                    relock_r1(w, &o.seq, &o.qual, r1s, sh.cb_len, sh.umi_len)
                {
                    let (cb, cb_q) = (&s[..sh.cb_len], &q[..sh.cb_len]);
                    let end = cb_s + sh.cb_len;
                    corrected_r1 = r1_with_umi(o, cb, cb_q, cb_s, end, sh.umi_len).or(Some((s, q)));
                    barcode_route = 7;
                    stats.relocked += 1;
                    relocked_at = Some(cb_s);
                }
            }
            // CB insertion rescue: one extra base inside the CB, at the
            // original lock; the UMI follows the 17 bp CB span (see umi_start).
            // Gated on confident location like reframe: on CAT02 the
            // unconfident ones landed on cell-like barcodes at 38.4% / 15.0%
            // (gated: 48.6% / 19.0%; reframe 37.8% / 19.5%).
            if corrected_r1.is_none() && sh.cb_insertion && located_confident {
                let umi_end = r1s + sh.cb_len + 1 + sh.umi_len;
                if umi_end <= o.seq.len() && umi_end <= o.qual.len() {
                    if let Some((cb, skip)) = w.insertion_cb(&o.seq[r1s..]) {
                        let mut q = o.qual[r1s..r1s + skip].to_vec();
                        q.extend_from_slice(&o.qual[r1s + skip + 1..r1s + sh.cb_len + 1]);
                        let end = r1s + sh.cb_len + 1;
                        corrected_r1 = r1_with_umi(o, &cb, &q, r1s, end, sh.umi_len);
                        if corrected_r1.is_some() {
                            barcode_route = 8;
                            stats.cb_insertion += 1;
                        }
                    }
                }
            }
            // Lenient within-2 substitution rescue: only when the CB is
            // confidently located (so the 2 mismatches are real base errors);
            // the UMI follows the CB as for every route (see umi_start).
            if corrected_r1.is_none() && sh.lenient && located_confident {
                if let Some(cb2) = w.correct_within2(&o.seq[r1s..r1s + sh.cb_len]) {
                    let cb_q = &o.qual[r1s..r1s + sh.cb_len];
                    let end = r1s + sh.cb_len;
                    corrected_r1 = r1_with_umi(o, &cb2, cb_q, r1s, end, sh.umi_len).or_else(|| {
                        let mut s = cb2.to_vec();
                        s.extend_from_slice(&o.seq[r1s + sh.cb_len..r1e]);
                        Some((s, o.qual[r1s..r1e].to_vec()))
                    });
                    barcode_route = 4;
                    stats.corrected2 += 1;
                }
            }
        }
    }
    Rescue {
        r1: corrected_r1,
        route: barcode_route,
        relocked_at,
    }
}
