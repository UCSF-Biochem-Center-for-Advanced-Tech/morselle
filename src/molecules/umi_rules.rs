//! Opt-in UMI rules for SBX feature-barcoding libraries (CRISPR guides).
//!
//! Both rules are off by default; `scripts/run_flex_multi.sh` turns them on
//! for the CRISPR stream only. Why they exist (CAT06 L1 CRISPR, CP1 SBX vs
//! full-depth Illumina/Cell Ranger, 3,000 shared RZ211 cells, depth-matched
//! at ~10.1 M guide reads; `CAT06/read_depth_comparison_out/
//! diff_investigation/FINDINGS.md`):
//!
//! * SBX errors are mostly 1-base indels. Inside the 12-bp UMI an indel
//!   shifts every later base, so the error copy is many mismatches away from
//!   its parent and 1-mismatch correction never merges it. In a saturated
//!   CRISPR library (~450 reads per molecule) that left 6.1x Cell Ranger's
//!   molecule count; also merging one-indel neighbours cut it to 1.9x.
//! * The remaining excess was 1-read UMIs >= 3 edits from every parent:
//!   low-quality reads (17x more often >= 2 sub-Q20 UMI bases than reads of
//!   real molecules), not a misplaced UMI window. About 0.2% of reads, but
//!   at ~450 reads per molecule that is about one spurious molecule per real
//!   one. Folding 1-read UMIs into the family's top molecule when that top
//!   has >= 10 reads matched Cell Ranger: RZ211 CD45 call rate 59.4% vs
//!   60.2%, Neg 0.80% vs 0.93%, molecules 0.99x (any K in 5..50 within ~1
//!   point).
//!
//! Not for GEX by default: most GEX molecules have few reads, so a 1-read
//! UMI there is usually a real molecule, and folding is only safe where the
//! family is deep. The `K` floor keeps shallow families untouched, but GEX
//! impact has not been measured.
use std::collections::HashMap as StdHashMap;
use std::hash::BuildHasher;

use rustc_hash::FxHashMap as HashMap;

use super::types::{Dna, Evidence, RawKey};

/// Which opt-in UMI rules `correct_group` applies.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct UmiRules {
    /// Also treat same-length one-indel neighbours as correction candidates
    /// (same 2n-1 read-support rule as 1-mismatch neighbours).
    pub indels: bool,
    /// When > 0: after correction, fold every 1-read molecule into the
    /// (cell, sample, gene) family's top molecule if that top has at least
    /// this many reads.
    pub fold_singletons_min_top: u32,
}

impl UmiRules {
    pub fn is_default(self) -> bool {
        self == Self::default()
    }
}

fn seq_mask(len: usize) -> u64 {
    if len == 0 {
        0
    } else {
        (1_u64 << (2 * len)) - 1
    }
}

/// Same-length sequences one indel away, as read from a fixed-length UMI
/// window: a base deleted at `i` (the window then reads one base further,
/// modelled as any base appended) or a base inserted at `i` (the window's
/// last base falls off). Tag bits above the sequence (length, position tag)
/// are preserved, so neighbours stay in the same position class. May repeat
/// a sequence or yield a 1-mismatch neighbour; never yields `umi` itself.
pub fn one_indel_neighbors(umi: Dna) -> impl Iterator<Item = Dna> {
    let len = umi.len();
    let mask = seq_mask(len);
    let seq = umi.0 & mask;
    let tags = umi.0 & !mask;
    (0..len)
        .flat_map(move |i| {
            let tail_len = len - i; // bases from i to the end
            let prefix = seq >> (2 * tail_len);
            (0..4_u64).flat_map(move |b| {
                // delete base i, append b
                let after = seq & seq_mask(tail_len - 1);
                let deleted = (prefix << (2 * tail_len)) | (after << 2) | b;
                // insert b at i, drop the last base
                let shifted = (seq >> 2) & seq_mask(tail_len - 1);
                let inserted = (prefix << (2 * tail_len)) | (b << (2 * (tail_len - 1))) | shifted;
                [deleted, inserted]
            })
        })
        .filter(move |&s| s != seq)
        .map(move |s| Dna(s | tags))
}

/// Fold 1-read molecules into their family's top molecule when the top has
/// at least `min_top` reads. `keys` is one (cell, sample, gene) group and
/// `roots` already holds its corrected roots. Returns how many were folded.
pub fn fold_deep_singletons<S: BuildHasher>(
    keys: &[RawKey],
    evidence: &StdHashMap<RawKey, Evidence, S>,
    roots: &mut HashMap<RawKey, RawKey>,
    min_top: u32,
) -> u64 {
    if min_top == 0 || keys.len() < 2 {
        return 0;
    }
    let mut support: HashMap<RawKey, u64> = HashMap::default();
    for key in keys {
        *support.entry(roots[key]).or_default() += u64::from(evidence[key].reads);
    }
    // Explicit tie-break (support, then lowest UMI) keeps this deterministic.
    let Some((&top, &top_reads)) = support
        .iter()
        .max_by_key(|(root, reads)| (**reads, std::cmp::Reverse(root.umi)))
    else {
        return 0;
    };
    if top_reads < u64::from(min_top) {
        return 0;
    }
    let mut folded = 0;
    for key in keys {
        let root = roots[key];
        if root != top && support[&root] == 1 {
            roots.insert(*key, top);
            folded += 1;
        }
    }
    folded
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dna(s: &str) -> Dna {
        Dna::pack(s).unwrap()
    }

    /// Brute force: every same-length string from one deletion + append or
    /// one insertion + truncation.
    fn brute(s: &str) -> std::collections::BTreeSet<String> {
        let b = s.as_bytes();
        let mut out = std::collections::BTreeSet::new();
        for i in 0..b.len() {
            for c in b"ACGT" {
                let mut del: Vec<u8> = b[..i].iter().chain(&b[i + 1..]).copied().collect();
                del.push(*c);
                let mut ins: Vec<u8> = b[..i].to_vec();
                ins.push(*c);
                ins.extend_from_slice(&b[i..b.len() - 1]);
                out.insert(String::from_utf8(del).unwrap());
                out.insert(String::from_utf8(ins).unwrap());
            }
        }
        out.remove(s);
        out
    }

    #[test]
    fn indel_neighbors_match_brute_force() {
        for s in [
            "ACGTACGTACGT",
            "AAAAAAAAAAAA",
            "TCCTCAGATGTA",
            "GATTACA",
            "A",
        ] {
            let got: std::collections::BTreeSet<String> =
                one_indel_neighbors(dna(s)).map(|d| d.unpack()).collect();
            assert_eq!(got, brute(s), "{s}");
            assert!(one_indel_neighbors(dna(s)).all(|d| d.len() == s.len()));
        }
    }

    #[test]
    fn indel_neighbors_keep_position_tag() {
        let tagged = dna("ACGTACGTACGT").with_position_tag(77).unwrap();
        for n in one_indel_neighbors(tagged) {
            assert_eq!(n.0 & !seq_mask(12), tagged.0 & !seq_mask(12));
        }
    }

    #[test]
    fn shifted_umi_is_an_indel_neighbor() {
        // Real CAT06 shape: one inserted base early in the UMI shifts the rest.
        let parent = dna("TCCTCAGATGTA");
        let child = dna("TCGCTCAGATGT"); // G inserted at 2, last base lost
        assert!(parent.hamming(child).unwrap() > 1);
        assert!(one_indel_neighbors(parent).any(|n| n == child));
    }

    fn key(umi: &str) -> RawKey {
        RawKey {
            cell: Dna(1),
            sample: 0,
            gene: 0,
            umi: dna(umi),
        }
    }

    fn ev(reads: u32) -> Evidence {
        Evidence {
            reads,
            best_ordinal: reads as u64,
            best_mapq: 60,
            best_aligned: 90,
            min_read_hash: 0,
        }
    }

    #[test]
    fn folds_singletons_only_in_deep_families() {
        let keys = [
            key("AAAAAAAAAAAA"),
            key("CCCCCCCCCCCC"),
            key("GGGGGGGGGGGG"),
            key("TTTTTTTTTTTT"),
        ];
        let reads = [40, 2, 1, 1];
        let evidence: StdHashMap<RawKey, Evidence> =
            keys.iter().zip(reads).map(|(k, r)| (*k, ev(r))).collect();
        let fresh = || keys.iter().map(|k| (*k, *k)).collect::<HashMap<_, _>>();

        let mut roots = fresh();
        assert_eq!(fold_deep_singletons(&keys, &evidence, &mut roots, 10), 2);
        assert_eq!(roots[&keys[2]], keys[0]);
        assert_eq!(roots[&keys[3]], keys[0]);
        assert_eq!(roots[&keys[1]], keys[1], "2-read molecules are kept");

        let mut roots = fresh();
        assert_eq!(
            fold_deep_singletons(&keys, &evidence, &mut roots, 41),
            0,
            "top below K"
        );
        let mut roots = fresh();
        assert_eq!(
            fold_deep_singletons(&keys, &evidence, &mut roots, 0),
            0,
            "off"
        );
    }
}
