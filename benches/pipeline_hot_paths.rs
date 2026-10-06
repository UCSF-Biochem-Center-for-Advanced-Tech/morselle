//! Per-stage hot-path microbenchmarks (see `PROFILING.md`).
//!
//! Portable (no external data files — synthetic-but-realistic reads, same
//! shape as `orient.rs`'s own unit-test fixtures), so `cargo bench` runs
//! identically on any machine, independent of whether OS-level sampling
//! profiling (perf/samply) is available. Add one function here per stage's
//! hot path as new phases land; each `c.bench_function` is independently
//! runnable/filterable (`cargo bench -- extract_gex_read`).

use std::io::Write as _;

use criterion::{black_box, criterion_group, criterion_main, Criterion};

use sbx10x::anchors::{revcomp, P5, P7, R1_PRIMER};
use sbx10x::barcode::Whitelist;
use sbx10x::extract::{extract, ExtractCfg};
use sbx10x::hamming_index::HammingIndex;
use sbx10x::orient::{orient, FallbackMode, OrientCfg};
use sbx10x::probe_barcode::ProbeWhitelist;
use sbx10x::probe_set::{Probe, ProbeMatcher};

fn orient_cfg() -> OrientCfg {
    OrientCfg {
        cb_len: 16,
        umi_len: 12,
        min_homopolymer: 12,
        r1_tail_max_mm: 2,
        anchor_max_edits: 6,
        fallback: FallbackMode::HammingTail,
        expected: None,
    }
}

fn extract_cfg() -> ExtractCfg {
    ExtractCfg {
        min_cdna_len: 20,
        adapter_max_mm: 2,
        trim_tso: false,
        max_r2_len: 0,
        flex_r2_boundary: false,
        flex_min_cdna_len_no_anchor: 65,
    }
}

/// A canonical 3' GEX-shaped read: real dT-primer homopolymer, ~200bp cDNA.
fn gex_read() -> (Vec<u8>, Vec<u8>) {
    let mut seq = Vec::new();
    seq.extend_from_slice(P5);
    seq.extend_from_slice(b"ACGTACGTAC"); // i5
    seq.extend_from_slice(R1_PRIMER);
    seq.extend_from_slice(b"AAACCCAAGAAACCCA"); // CB (16)
    seq.extend_from_slice(b"ACGTACGTACGT"); // UMI (12)
    seq.extend(std::iter::repeat(b'T').take(30)); // real dT-primer run
    for _ in 0..10 {
        seq.extend_from_slice(b"GACTGACTGACTGACTGACTGACT"); // ~200bp cDNA
    }
    seq.extend_from_slice(b"AGATCGGAAGAGC"); // 3' adapter
    seq.extend_from_slice(b"TGCATGCATG"); // i7
    seq.extend_from_slice(P7);
    let qual = vec![b'I'; seq.len()];
    (seq, qual)
}

/// A Flex/CRISPR-shaped read, matching what FLEX_SUPPORT_PLAN.md's Phase 0
/// found empirically: no real homopolymer (orient()'s TailFallback path
/// fires, ~88% of real Flex reads per the orient_method_profile diagnostic),
/// then a probe barcode 10bp past the CB+UMI bracket.
fn flex_read(probe_barcode_fwd: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let mut seq = Vec::new();
    seq.extend_from_slice(P5);
    seq.extend_from_slice(b"ACGTACGTAC");
    seq.extend_from_slice(R1_PRIMER);
    seq.extend_from_slice(b"AAACCCAAGAAACCCA"); // CB
    seq.extend_from_slice(b"ACGTACGTACGT"); // UMI
    seq.extend_from_slice(b"GGGTT"); // short, non-qualifying spacer (no real run)
    seq.extend_from_slice(&revcomp(probe_barcode_fwd)); // probe barcode, revcomp'd
    seq.extend_from_slice(b"GACGACGACGACGACGACGACGACGACGAC"); // probe payload
    seq.extend_from_slice(b"AGATCGGAAGAGC");
    seq.extend_from_slice(b"TGCATGCATG");
    seq.extend_from_slice(P7);
    let qual = vec![b'I'; seq.len()];
    (seq, qual)
}

fn bench_orient_extract(c: &mut Criterion) {
    let (gex_seq, gex_qual) = gex_read();
    let (flex_seq, flex_qual) = flex_read(b"ACGTACGGGA");
    let ocfg = orient_cfg();
    let ecfg = extract_cfg();

    c.bench_function("orient_extract/gex_homopolymer_path", |b| {
        b.iter(|| {
            let o = orient(black_box(&gex_seq), black_box(&gex_qual), &ocfg);
            black_box(extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg, false))
        })
    });
    c.bench_function("orient_extract/flex_tail_fallback_path", |b| {
        b.iter(|| {
            let o = orient(black_box(&flex_seq), black_box(&flex_qual), &ocfg);
            black_box(extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg, false))
        })
    });
}

fn bench_hamming_index(c: &mut Criterion) {
    // 384 entries, same scale as the real Flex probe-barcode whitelist.
    let mut index: HammingIndex<u16> = HammingIndex::new();
    let mut rng_state = 0x1234_5678_u32;
    let mut next_seq = || -> [u8; 10] {
        let mut s = [0u8; 10];
        for b in &mut s {
            rng_state = rng_state.wrapping_mul(1_103_515_245).wrapping_add(12345);
            *b = b"ACGT"[((rng_state >> 16) & 3) as usize];
        }
        s
    };
    let mut probes = Vec::new();
    for i in 0..384u16 {
        let s = next_seq();
        index.insert(&s, i, 1);
        probes.push(s);
    }
    let query = probes[200];
    let mut haystack = vec![b'N'; 60];
    haystack[20..30].copy_from_slice(&query);

    c.bench_function("hamming_index/exact_get_384_entries_mm1", |b| {
        b.iter(|| black_box(index.get(black_box(&query))))
    });
    c.bench_function(
        "hamming_index/find_windowed_384_entries_mm1_window10",
        |b| b.iter(|| black_box(index.find_windowed(black_box(&haystack), Some(17), 10, 10))),
    );
}

fn bench_probe_barcode_locate(c: &mut Criterion) {
    let mut f = tempfile::NamedTempFile::new().unwrap();
    for i in 0..384u32 {
        let seq = std::str::from_utf8(&encode_dna(i)).unwrap().to_string();
        writeln!(f, "{seq}\t{seq}\tD-{i:03}").unwrap();
    }
    let wl = ProbeWhitelist::load(f.path(), 1, 10).unwrap();
    let (seq, _) = flex_read(&decode_dna(200));

    c.bench_function("probe_barcode/locate_windowed_384_entries", |b| {
        // polyt_end in this synthetic read = CB+UMI end (5bp short spacer, run_len 0).
        let polyt_end = P5.len() + 10 + R1_PRIMER.len() + 16 + 12;
        b.iter(|| black_box(wl.locate(black_box(&seq), polyt_end)))
    });
}

fn encode_dna(mut n: u32) -> [u8; 10] {
    let mut out = [b'A'; 10];
    for slot in out.iter_mut() {
        *slot = b"ACGT"[(n % 4) as usize];
        n /= 4;
    }
    out
}

fn decode_dna(n: u32) -> Vec<u8> {
    encode_dna(n).to_vec()
}

fn bench_cb_whitelist_classify(c: &mut Criterion) {
    // 6.8M-scale real whitelists are too slow to build per-bench-run; 200k is
    // enough to exercise the same 2-bit-packed FxHashSet code path realistically.
    let mut f = tempfile::NamedTempFile::new().unwrap();
    let mut rng_state = 0x9e37_79b9_u32;
    let mut cb = Vec::new();
    for _ in 0..200_000 {
        let mut s = [0u8; 16];
        for b in &mut s {
            rng_state = rng_state.wrapping_mul(1_103_515_245).wrapping_add(12345);
            *b = b"ACGT"[((rng_state >> 16) & 3) as usize];
        }
        writeln!(f, "{}", std::str::from_utf8(&s).unwrap()).unwrap();
        cb = s.to_vec();
    }
    let wl = Whitelist::load(f.path()).unwrap();
    let mut cb_1mm = cb.clone();
    cb_1mm[3] = if cb_1mm[3] == b'A' { b'C' } else { b'A' };

    c.bench_function("cb_whitelist/classify_exact_hit_200k", |b| {
        b.iter(|| black_box(wl.classify(black_box(&cb))))
    });
    c.bench_function("cb_whitelist/classify_within1_200k", |b| {
        b.iter(|| black_box(wl.classify(black_box(&cb_1mm))))
    });
}

fn bench_probeset_align_matcher(c: &mut Criterion) {
    let probes = vec![
        Probe {
            gene_id: "ENSG00000000001".into(),
            gene_name: "GENE1".into(),
            probe_id: "p1".into(),
            region: "spliced".into(),
            seq: b"AAAACCCCGGGGTTTTACGTACGTAAAACCCCGGGGTTTTACGTACGT01".to_vec(),
        },
        Probe {
            gene_id: "ENSG00000000002".into(),
            gene_name: "GENE2".into(),
            probe_id: "p2".into(),
            region: "spliced".into(),
            seq: b"TTTTGGGGCCCCAAAATGCATGCATTTTTGGGGCCCCAAAATGCATGC02".to_vec(),
        },
    ];
    let matcher = ProbeMatcher::build(&probes, 1, 15);
    let mut read = vec![b'N'; 15];
    read.extend_from_slice(b"AAAACCCCGGGGTTTTACGTACGTAAAACCCCGGGGTTTTACGTACGT01");
    read.extend_from_slice(b"NNNNNNNNNN");

    c.bench_function("probeset_align/probe_matcher_locate_exact", |b| {
        b.iter(|| black_box(matcher.locate(black_box(&read))))
    });
}

criterion_group!(
    benches,
    bench_orient_extract,
    bench_hamming_index,
    bench_probe_barcode_locate,
    bench_cb_whitelist_classify,
    bench_probeset_align_matcher,
);
criterion_main!(benches);
