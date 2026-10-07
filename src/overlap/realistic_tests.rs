use super::alignment::overlap;
use super::types::{InsertClass, Merged, OverlapCfg};
use crate::anchors::{revcomp, R1_PRIMER};
use crate::extract::dt_run_end;
use crate::yc::decode;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

pub const CB: &[u8] = b"AAACCCAAGAAACACT";
pub const UMI: &[u8] = b"ACGTACGTACGT";
/// The dT block: RT-primer dT annealed to the polyA tail.
pub const DT: usize = 30;
/// CB(16) + UMI(12): where the dT run starts in R1 and in the merged read.
pub const BARCODE_LEN: usize = 28;
pub const ADAPTER: &[u8] = b"AGATCGGAAGAGC";

/// Deterministic complex sequence — no run longer than 3.
pub fn complex(len: usize, seed: u64) -> Vec<u8> {
    let bases = b"ACGT";
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let mut out: Vec<u8> = Vec::with_capacity(len);
    while out.len() < len {
        let base = bases[rng.gen_range(0..4)];
        let tail = out.len().saturating_sub(3);
        if out.len() >= 3 && out[tail..].iter().all(|&b| b == base) {
            continue;
        }
        out.push(base);
    }
    out
}

/// A cDNA core carrying four homopolymer runs of different bases. Ends on `C`
/// so `revcomp` does not begin with `T` and blur the dT boundary.
pub fn core_with_homopolymers() -> Vec<u8> {
    let mut core = Vec::new();
    core.extend_from_slice(&complex(24, 1));
    core.extend_from_slice(b"AAAAAAAAAA");
    core.extend_from_slice(&complex(22, 2));
    core.extend_from_slice(b"TTTTTTTTT");
    core.extend_from_slice(&complex(20, 3));
    core.extend_from_slice(b"GGGGGGGG");
    core.extend_from_slice(&complex(18, 4));
    core.extend_from_slice(b"CCCCCCC");
    core.extend_from_slice(&complex(18, 5));
    *core.last_mut().unwrap() = b'C';
    core
}

pub struct Pair {
    pub r1: Vec<u8>,
    pub r1q: Vec<u8>,
    pub r2: Vec<u8>,
    pub r2q: Vec<u8>,
    /// `CB | UMI | T(DT) | rc(core)` — what a complete merge yields.
    pub molecule: Vec<u8>,
    /// True when the mates reach each other at these read lengths.
    pub overlaps: bool,
}

/// Build a mate pair for a fragment with cDNA `core`, sequenced `r1_len` from
/// the barcode end and `r2_len` from the cDNA end.
pub fn build(core: &[u8], r1_len: usize, r2_len: usize) -> Pair {
    let mut molecule = Vec::new();
    molecule.extend_from_slice(CB);
    molecule.extend_from_slice(UMI);
    molecule.extend(std::iter::repeat_n(b'T', DT));
    molecule.extend_from_slice(&revcomp(core));

    let mut r1 = molecule.clone();
    r1.extend_from_slice(ADAPTER);  // read-through past a short fragment
    r1.truncate(r1_len);

    let mut r2 = revcomp(&molecule);  // core | A(DT) | rc(UMI) | rc(CB)
    r2.extend_from_slice(&revcomp(R1_PRIMER));  // read-through into the primer
    r2.truncate(r2_len);

    // In molecule coordinates R1 covers [0, a) and rc(R2) covers [M - b, M).
    let m = molecule.len();
    let a = r1_len.min(m);
    let b = r2_len.min(m);
    let overlaps = m.saturating_sub(b) < a;

    Pair {
        r1q: vec![b'I'; r1.len()],
        r2q: vec![b'I'; r2.len()],
        r1,
        r2,
        molecule,
        overlaps,
    }
}

fn merge(pair: &Pair, cfg: &OverlapCfg) -> Option<Merged> {
    // The raw mates, untrimmed. This is the whole point.
    overlap(&pair.r1, &pair.r1q, &pair.r2, &pair.r2q, cfg)
}

#[test]
fn merges_a_short_fragment_with_read_through_on_both_mates() {
    // 60 bp cDNA, 150 bp reads: R1 runs into adapter, R2 back into the primer.
    let core = complex(60, 11);
    let pair = build(&core, 150, 150);
    let merged = merge(&pair, &OverlapCfg::default()).expect("short fragment merges");

    assert_eq!(merged.seq, pair.molecule);
    assert_eq!(merged.class, InsertClass::Short);
    assert_eq!(merged.discordant, 0);
    // Both read-throughs are dropped purely by the overlap geometry.
    assert_eq!(merged.left_overhang, 0);
    assert_eq!(merged.right_overhang, 0);
}

#[test]
fn merges_a_medium_fragment_into_the_whole_molecule() {
    // 130 bp cDNA: the mates meet across the tail of the dT block.
    let core = complex(130, 12);
    let pair = build(&core, 150, 150);
    let merged = merge(&pair, &OverlapCfg::default()).expect("medium fragment merges");

    assert_eq!(merged.seq, pair.molecule);
    assert_eq!(merged.class, InsertClass::Medium);
    // Merging the raw mates gives an overhang on each side, unlike the
    // pre-trimmed geometry where R1 could never start first.
    assert!(merged.left_overhang > 0);
    assert!(merged.right_overhang > 0);
}

#[test]
fn a_long_fragment_does_not_merge() {
    let core = complex(400, 13);
    let pair = build(&core, 150, 150);
    assert!(!pair.overlaps, "fixture must not overlap");
    assert!(merge(&pair, &OverlapCfg::default()).is_none());
}

#[test]
fn merges_when_r1_is_sequenced_shorter_than_r2() {
    // R2 spans the whole molecule while R1 stops inside it, so R1 is
    // strictly interior to rc(R2) — overhangs on both sides. Probing
    // rc(R2)'s ends would find nothing; we index it and probe R1 instead.
    let core = complex(60, 23);
    let pair = build(&core, 100, 150);
    assert!(pair.r1.len() < pair.r2.len());

    let merged = merge(&pair, &OverlapCfg::default()).expect("an interior R1 must still merge");
    assert_eq!(merged.seq, pair.molecule);
    assert_eq!(merged.left_overhang, 0);
    assert!(merged.right_overhang > 0);
}

/// Does `haystack` contain `needle` as a contiguous run?
fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    needle.len() <= haystack.len() && haystack.windows(needle.len()).any(|w| w == needle)
}

#[test]
fn read_through_is_removed_by_geometry_not_by_matching() {
    // The R1 primer precedes the CB, so rc(R2)'s 5' read-through is primer
    // and R1's 3' read-through is adapter. Both sit outside the overlap.
    let core = complex(60, 26);
    let pair = build(&core, 150, 150);

    // The raw mates really do carry read-through, or this proves nothing.
    assert!(contains(&pair.r1, ADAPTER), "R1 must read into adapter");
    let rc2 = crate::anchors::revcomp(&pair.r2);
    assert!(
        contains(&rc2, &crate::anchors::R1_PRIMER[8..]),
        "rc(R2) must read back into the R1 primer"
    );

    let merged = merge(&pair, &OverlapCfg::default()).unwrap();
    assert_eq!(merged.seq, pair.molecule);
    assert!(!contains(&merged.seq, ADAPTER), "adapter is gone");
    assert!(
        !contains(&merged.seq, &crate::anchors::R1_PRIMER[8..]),
        "primer read-through is gone"
    );
    // Nothing was matched to achieve that — both overhangs were simply
    // outside the span the two mates agree on.
    assert_eq!(merged.left_overhang, 0);
    assert_eq!(merged.right_overhang, 0);
}

#[test]
fn read_through_is_removed_even_when_the_adapter_is_misread() {
    // Corrupt the adapter badly enough to defeat a Hamming-seeded trimmer.
    // The geometry never looks at it, so the merge is unaffected.
    let core = complex(60, 27);
    let mut pair = build(&core, 150, 150);
    let end = pair.molecule.len();
    for offset in [0usize, 2, 4, 6, 8] {
        if end + offset < pair.r1.len() {
            pair.r1[end + offset] = b'G';
        }
    }
    assert!(
        !contains(&pair.r1, ADAPTER),
        "adapter is no longer matchable"
    );

    let merged = merge(&pair, &OverlapCfg::default()).unwrap();
    assert_eq!(
        merged.seq, pair.molecule,
        "an unmatchable adapter is still trimmed, because geometry decided it"
    );
}

#[test]
fn short_r1_with_r2_reading_past_the_barcode() {
    // R1 100 bp, R2 150 bp reading 10 bp past the CB into the R1 primer.
    // Solving R2 = core | A(30) | rc(UMI)(12) | rc(CB)(16) | rc(primer)(10)
    // gives an 82 bp cDNA and a 140 bp molecule. R1 covers molecule [0,100)
    // and rc(R2) covers [-10,140): R1 is strictly interior, overhang on both
    // sides, and R1 never reaches the adapter.
    let core = complex(82, 28);
    let pair = build(&core, 100, 150);
    assert_eq!(pair.molecule.len(), 140);

    let merged = merge(&pair, &OverlapCfg::default()).expect("interior R1 merges");

    assert_eq!(merged.seq, pair.molecule);
    // Every base of R1 is inside the overlap: CB, UMI, dT and 42 bp of cDNA
    // all come out duplex-confirmed.
    assert_eq!(merged.overlap_len, 100);
    // The 10 bp of primer is dropped; R2 alone extends the molecule 40 bp.
    assert_eq!(merged.left_overhang, 0);
    assert_eq!(merged.right_overhang, 40);
    // R1 never saw adapter, but R2 read through — so this is a short insert.
    // A rule that only checks R1 would wrongly call it Medium.
    assert_eq!(merged.class, InsertClass::Short);
    assert_eq!(merged.yc.as_ref().unwrap().to_string(), "+100-40");

    // What Phase 3 will emit: strip CB+UMI+dT, keep the 82 bp cDNA, and
    // carry the tag across the trim and the flip into R2 orientation.
    let cdna_start = dt_run_end(&merged.seq, BARCODE_LEN);
    assert_eq!(cdna_start, 58);
    let mut tag = merged.yc.clone().unwrap();
    tag.trim(cdna_start as u64, 0);
    assert_eq!(tag.to_string(), "+42-40");
    tag.reverse_complement();
    // In sense orientation the 5'-most 40 bp are R2-only, the rest duplex.
    assert_eq!(tag.to_string(), "40-42+");
    assert_eq!(tag.sequence_length(), 82);

    let emitted = crate::anchors::revcomp(&merged.seq[cdna_start..]);
    assert_eq!(emitted, core, "the emitted read is the cDNA, sense strand");
}

#[test]
fn merges_across_four_homopolymer_runs() {
    let core = core_with_homopolymers();
    let pair = build(&core, 150, 150);
    let merged = merge(&pair, &OverlapCfg::default()).expect("homopolymers must not block");

    assert_eq!(merged.seq, pair.molecule);
    assert_eq!(merged.discordant, 0);

    // The tag is lossless with respect to the merged construct: it
    // reconstructs each mate's contribution, not the read-through the merge
    // discarded. R1 contributes a prefix, rc(R2) a suffix.
    let tag = merged.yc.as_ref().unwrap();
    let decoded = decode(tag, &merged.seq, &merged.qual).unwrap();
    assert!(
        pair.r1.starts_with(&decoded.r1_seq),
        "R1 contributes a prefix"
    );
    assert!(
        crate::anchors::revcomp(&pair.r2).ends_with(&decoded.r2_seq),
        "rc(R2) contributes a suffix"
    );
}

#[test]
fn the_dt_block_does_not_starve_the_seed_probe() {
    // The mates meet across the dT run, so the first usable seed on rc(R2)
    // is ~20 bp in. A fixed positional probe window would see only T.
    let core = complex(130, 24);
    let pair = build(&core, 150, 150);
    let rc2 = crate::anchors::revcomp(&pair.r2);
    assert!(
        rc2[..20].iter().all(|&b| b == b'T'),
        "fixture must start rc(R2) inside the dT block"
    );
    assert!(merge(&pair, &OverlapCfg::default()).is_some());
}

#[test]
fn the_merged_read_still_yields_cb_umi_and_cdna() {
    // What Phase 3 does with the merge: barcode at a fixed offset, cDNA
    // after the dT run — now backed by both strands.
    let core = complex(60, 25);
    let mut pair = build(&core, 150, 150);
    // Q20 inputs, so the summed duplex quality is visible below the Q40 cap.
    pair.r1q = vec![b'5'; pair.r1.len()];
    pair.r2q = vec![b'5'; pair.r2.len()];
    let merged = merge(&pair, &OverlapCfg::default()).unwrap();

    assert_eq!(&merged.seq[..16], CB);
    assert_eq!(&merged.seq[16..28], UMI);
    let cdna_start = dt_run_end(&merged.seq, BARCODE_LEN);
    assert_eq!(cdna_start, BARCODE_LEN + DT);
    assert_eq!(
        &merged.seq[cdna_start..],
        &crate::anchors::revcomp(&core)[..]
    );

    // The barcode sat inside the overlap, so both strands confirm it: Q20 +
    // Q20 saturates at the Q40 cap. On a short fragment, merging before
    // trimming buys duplex confirmation of the CB and UMI themselves.
    assert!(merged.qual[..28].iter().all(|&q| q == b'I'));
}

#[test]
fn tolerates_scattered_sequencing_mismatches() {
    let core = complex(130, 14);
    for every in [40usize, 20, 12] {
        let mut pair = build(&core, 150, 150);
        let mut corrupted = 0;
        for i in (30..pair.r2.len()).step_by(every) {
            pair.r2[i] = if pair.r2[i] == b'A' { b'C' } else { b'A' };
            corrupted += 1;
        }
        let merged = merge(&pair, &OverlapCfg::default())
            .unwrap_or_else(|| panic!("every {every} bp ({corrupted} errors) should merge"));
        assert_eq!(merged.seq.len(), pair.molecule.len());
        // Inside the overlap R1 is the higher-quality call and wins every
        // disagreement, so the molecule is recovered exactly. Beyond it R2 is
        // the only evidence there is and its errors are copied through —
        // which is precisely why merging is worth doing.
        let duplex_end = merged.left_overhang + merged.overlap_len;
        assert_eq!(
            merged.seq[..duplex_end],
            pair.molecule[..duplex_end],
            "every {every} bp: the overlap must be error-corrected"
        );
        assert!(merged.discordant > 0);
    }
}

#[test]
fn polya_tail_does_not_seed_a_false_merge() {
    // Far too long to overlap, but both mates carry long A/T runs.
    // Homopolymer seed masking is what keeps this from merging.
    let mut core = complex(200, 15);
    core.extend_from_slice(b"AAAAAAAAAAAAAAAAAAAA");
    core.extend_from_slice(&complex(200, 16));
    *core.last_mut().unwrap() = b'C';
    let pair = build(&core, 150, 150);

    assert!(!pair.overlaps);
    assert!(
        merge(&pair, &OverlapCfg::default()).is_none(),
        "poly-A/T must not be enough to anchor a merge"
    );
}

#[test]
fn an_overlap_that_is_only_homopolymer_does_not_merge() {
    // The mates meet, but only across a run of A. No complex seed exists in
    // the shared region, so we decline rather than merge on evidence that
    // means nothing. Built directly: the fixed 58 bp CB+UMI+dT offset cannot
    // place an overlap wholly inside an internal homopolymer.
    let mut r1 = complex(40, 21);
    r1.extend(std::iter::repeat_n(b'A', 40));
    let mut rc2 = vec![b'A'; 40];
    rc2.extend_from_slice(&complex(40, 22));
    let r2 = crate::anchors::revcomp(&rc2);

    assert!(overlap(
        &r1,
        &vec![b'I'; r1.len()],
        &r2,
        &vec![b'I'; r2.len()],
        &OverlapCfg::default()
    )
    .is_none());
}

#[test]
fn indel_between_mates_is_characterised() {
    // The merge is gapless by construction. An indel inside the overlap
    // frame-shifts everything after it, so the mismatch budget rejects the
    // alignment: we decline rather than emit a shifted chimera. The dT block
    // makes this worth watching — homopolymers are where Illumina indels live.
    let core = complex(130, 19);
    let mut pair = build(&core, 150, 150);
    pair.r2.insert(60, b'A');
    pair.r2.truncate(150);
    pair.r2q = vec![b'I'; pair.r2.len()];

    match merge(&pair, &OverlapCfg::default()) {
        None => {}  // preferred outcome
        Some(m) => assert!(
            m.seq.starts_with(&pair.molecule[..40]),
            "a partial merge must not corrupt the R1 side"
        ),
    }
}

#[test]
fn indel_outside_the_overlap_is_harmless() {
    let core = complex(130, 20);
    let mut pair = build(&core, 150, 150);
    pair.r2.insert(148, b'A');
    pair.r2.truncate(150);
    pair.r2q = vec![b'I'; pair.r2.len()];

    let merged = merge(&pair, &OverlapCfg::default()).expect("overlap is untouched");
    assert_eq!(merged.seq, pair.molecule);
}
