use std::io::BufRead;
use std::path::PathBuf;

use anyhow::Result;
use rustc_hash::FxHashSet;

use super::types::{ClipClasses, HOMO, K};
pub use crate::anchors::revcomp;

/// All k-mers of `m` and its reverse complement — indel-robust motif detection.
pub fn kmers(m: &[u8]) -> Vec<[u8; K]> {
    let mut v = Vec::new();
    for s in [m.to_vec(), revcomp(m)] {
        for w in s.windows(K) {
            let mut a = [0u8; K];
            a.copy_from_slice(w);
            v.push(a);
        }
    }
    v.sort_unstable();
    v.dedup();
    v
}

pub fn has_kmer(clip: &[u8], set: &[[u8; K]]) -> bool {
    if clip.len() < K {
        return false;
    }
    clip.windows(K).any(|w| {
        set.binary_search(&{
            let mut a = [0u8; K];
            a.copy_from_slice(w);
            a
        })
        .is_ok()
    })
}

/// SplitMix64 — stable across machines, so the molecule subsample is identical
/// for every BAM (that is what makes the dumps joinable).
pub fn splitmix64(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E3779B97F4A7C15);
    let mut x = z;
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D049BB133111EB);
    x ^ (x >> 31)
}

/// Stable FNV-1a hash for variable-length molecule components. This is used only
/// for bounded-memory cardinality sampling, not identity reconstruction.
pub fn stable_hash(parts: &[&[u8]]) -> u64 {
    let mut value = 0xcbf29ce484222325u64;
    for part in parts {
        for &byte in *part {
            value ^= byte as u64;
            value = value.wrapping_mul(0x100000001b3);
        }
        value ^= 0xff;
        value = value.wrapping_mul(0x100000001b3);
    }
    value
}

/// 2-bit pack a barcode/UMI; None on any non-ACGT base.
pub fn pack(s: &[u8]) -> Option<u64> {
    let mut c = 0u64;
    for &b in s {
        c = (c << 2)
            | match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => return None,
            };
    }
    Some(c)
}

/// Decode a 2-bit packed code of `len` bases.
pub fn unpack(mut code: u64, len: usize) -> String {
    let mut s = vec![b'A'; len];
    for i in (0..len).rev() {
        s[i] = b"ACGT"[(code & 3) as usize];
        code >>= 2;
    }
    String::from_utf8(s).unwrap()
}

pub use crate::bam_helpers::{aux_char, aux_int, aux_str};

/// Classify a soft-clipped sequence: what did the aligner refuse to align?
pub fn classify(clip: &[u8], adapter_k: &[[u8; K]], tso_k: &[[u8; K]]) -> &'static str {
    let n = clip.len() as f64;
    let a = clip.iter().filter(|&&b| b == b'A').count() as f64 / n;
    let t = clip.iter().filter(|&&b| b == b'T').count() as f64 / n;
    if has_kmer(clip, adapter_k) {
        "adapter"
    } else if has_kmer(clip, tso_k) {
        "tso"
    } else if a >= HOMO {
        "polyA"
    } else if t >= HOMO {
        "polyT"
    } else {
        "unexplained"
    }
}

pub fn bump(c: &mut ClipClasses, kind: &str) {
    c.n_clips += 1;
    match kind {
        "adapter" => c.adapter += 1,
        "tso" => c.tso += 1,
        "polyA" => c.poly_a += 1,
        "polyT" => c.poly_t += 1,
        _ => c.unexplained += 1,
    }
}

pub fn load_called_cells(path: &PathBuf) -> Result<FxHashSet<String>> {
    let reader = crate::io::open_text_reader(path)?;
    let mut cells = FxHashSet::default();
    for line in reader.lines() {
        let line = line?;
        let barcode = line.split('\t').next().unwrap_or("").trim();
        if !barcode.is_empty() {
            cells.insert(barcode.to_string());
        }
    }
    Ok(cells)
}

pub fn nh_bucket(nh: i64) -> usize {
    match nh {
        ..=1 => 0,
        2 => 1,
        3..=5 => 2,
        6..=10 => 3,
        _ => 4,
    }
}
