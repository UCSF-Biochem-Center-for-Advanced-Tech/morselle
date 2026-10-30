//! `IndexLut` against the reference matcher (`HammingIndex`).
use super::*;
use crate::hamming_index::HammingIndex;

fn keys(n: usize, len: usize, seed: u32) -> Vec<Vec<u8>> {
    let mut s = seed;
    (0..n)
        .map(|_| {
            (0..len)
                .map(|_| {
                    s = s.wrapping_mul(1_103_515_245).wrapping_add(12345);
                    b"ACGT"[((s >> 16) & 3) as usize]
                })
                .collect()
        })
        .collect()
}

/// Same answers as `HammingIndex::find_windowed` for every frame of a
/// mixed haystack, dense and sparse, including ambiguous neighbours --
/// except that a shifted-frame mismatch hit abstains (`float_tie`)
/// exactly when a brute-force scan finds a competing sample: another
/// sample within `mm` mismatches in some frame, or one indel from the
/// `len ± 1` bp frame there.
#[test]
fn matches_hamming_index_exactly() {
    let mut ties = 0;
    for (n, len, mm) in [(14, 10, 1), (384, 10, 1), (96, 8, 2), (40, 12, 1)] {
        let ks = keys(n, len, 7 + n as u32);
        let mut old = HammingIndex::new();
        for (i, k) in ks.iter().enumerate() {
            old.insert(k, i, mm);
        }
        let lut = IndexLut::new(&ks, mm);
        let mut hay = Vec::new();
        for k in &ks {
            hay.extend_from_slice(b"GATTNACA");
            let mut v = k.clone();
            v[len / 2] = b'A';
            hay.extend_from_slice(&v);
        }
        for start in 0..hay.len() {
            for window in [0isize, 3, 12] {
                let got = lut.find(&hay, Some(start), window, true);
                let want = old.find_windowed(&hay, Some(start), len, window);
                let ctx = format!("n={n} len={len} start={start} window={window}");
                let Some(x) = want else {
                    assert_eq!(got.sample, None, "{ctx}");
                    continue;
                };
                let k = offset_order(window)
                    .find(|&k| old.find_windowed(&hay, frame(start, k), len, 0).is_some())
                    .unwrap();
                let at = frame(start, k).unwrap();
                let weak = k != 0 && hay[at..at + len] != ks[x][..];
                let tie = weak && competitor(&hay, start as isize, window, &ks, x, mm);
                assert_eq!(got.float_tie, tie, "{ctx}");
                assert_eq!(got.sample, (!tie).then_some(x), "{ctx}");
                ties += usize::from(tie);
            }
        }
    }
    assert!(ties > 0, "the haystack exercises no tie");
}

fn frame(start: usize, k: isize) -> Option<usize> {
    usize::try_from(start as isize + k).ok()
}

/// Brute force: does any frame (`len - 1`, `len` or `len + 1` bp) match
/// candidates that exclude `x` -- samples within `mm` mismatches at
/// `len`, or one indel away at `len ± 1`? A frame `x` shares with
/// others still counts for `x`.
fn competitor(hay: &[u8], start: isize, w: isize, ks: &[Vec<u8>], x: usize, mm: u32) -> bool {
    let len = ks[x].len();
    offset_order(w).any(|k| {
        [len - 1, len, len + 1].into_iter().any(|l| {
            let Some(f) = usize::try_from(start + k)
                .ok()
                .and_then(|p| hay.get(p..p + l))
                .filter(|f| f.iter().all(|b| b"ACGT".contains(b)))
            else {
                return false;
            };
            let hits: Vec<usize> = (0..ks.len())
                .filter(|&y| {
                    if l == len {
                        f.iter().zip(&ks[y]).filter(|(a, b)| a != b).count() as u32 <= mm
                    } else {
                        crate::demux::edit_distance(f, &ks[y]) == 1
                    }
                })
                .collect();
            !hits.is_empty() && !hits.contains(&x)
        })
    })
}

#[test]
fn non_acgt_frames_never_match() {
    let ks = vec![b"ACGTACGTAC".to_vec()];
    let lut = IndexLut::new(&ks, 1);
    assert_eq!(lut.get(b"ACGTACGTAC").sample, Some(0));
    assert_eq!(lut.get(b"ACGTNCGTAC").sample, None);
    assert_eq!(lut.get(b"acgtacgtac").sample, None);
}
