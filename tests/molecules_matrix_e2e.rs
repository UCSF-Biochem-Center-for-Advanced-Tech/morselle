//! End-to-end regression test: annotated BAM -> `sbx10x-molecules` ->
//! `cr-bam-to-mtx`, run through the real binaries, checking the final
//! matrix counts under every collision policy and execution path.
//!
//! The fixture deliberately carries the tags real data has -- `UB` and `UR`
//! from extract, `ZG` from annotate on a multi-gene `GX` -- because the
//! 0.5.7 refactor passed every unit test while crashing on `UB` ("tag is
//! already present") and writing no `xf:25` representatives, so the matrix
//! came out empty. Nothing short of running both binaries catches that.
mod common;

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::Command;

use rust_htslib::bam::{
    self, header::HeaderRecord, record::Aux, record::Cigar, record::CigarString, Read as _,
};

const C1: &str = "AACCGGTTAACCGGTT-1";
const C2: &str = "TTGGCCAATTGGCCAA-1";
const U1: &str = "AAAAAAAAAAAA";
const U1_ERR: &str = "AAAAAAAAAAAC"; // 1 mismatch from U1
const U2: &str = "CCCCCCCCCCCC";
const U3: &str = "GGGGGGGGGGGG";
const U4: &str = "TTTTTTTTTTTT";

/// (name, pos, CB, GX, GN, UMI, ZG)
const READS: &[(&str, i64, &str, &str, &str, &str, Option<&str>)] = &[
    ("r1", 100, C1, "G1", "g1", U1, None),
    ("r2", 110, C1, "G1", "g1", U1, None),
    ("r3", 120, C1, "G1", "g1", U1, None),
    ("r4", 130, C1, "G1", "g1", U1_ERR, None), // corrects into U1 (single-cell)
    ("r5", 200, C1, "G2", "g2", U2, None),
    ("r6", 210, C1, "G2", "g2", U2, None),
    ("r7", 300, C1, "G1;G3", "g1;g3", U3, Some("G1")), // ZG resolves to G1
    ("r8", 400, C1, "G3", "g3", U1, None),             // collides with G1/U1
    ("r9", 500, C2, "G2", "g2", U4, None),
    ("r10", 510, C2, "G2", "g2", U4, None),
];

fn write_fixture(path: &Path) {
    let mut header = bam::Header::new();
    header.push_record(
        HeaderRecord::new(b"HD")
            .push_tag(b"VN", "1.6")
            .push_tag(b"SO", "coordinate"),
    );
    header.push_record(
        HeaderRecord::new(b"SQ")
            .push_tag(b"SN", "chr1")
            .push_tag(b"LN", 100_000),
    );
    let mut writer = bam::Writer::from_path(path, &header, bam::Format::Bam).unwrap();
    for &(name, pos, cb, gx, gn, umi, zg) in READS {
        let mut record = bam::Record::new();
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let cigar = CigarString(vec![Cigar::Match(seq.len() as u32)]);
        record.set(name.as_bytes(), Some(&cigar), seq, &[30; 40]);
        record.set_tid(0);
        record.set_pos(pos);
        record.set_mapq(60);
        record.unset_unmapped();
        record.push_aux(b"CB", Aux::String(cb)).unwrap();
        record.push_aux(b"UR", Aux::String(umi)).unwrap();
        record.push_aux(b"UB", Aux::String(umi)).unwrap();
        record.push_aux(b"GX", Aux::String(gx)).unwrap();
        record.push_aux(b"GN", Aux::String(gn)).unwrap();
        if let Some(zg) = zg {
            record.push_aux(b"ZG", Aux::String(zg)).unwrap();
        }
        writer.write(&record).unwrap();
    }
    drop(writer);
    bam::index::build(path, None, bam::index::Type::Bai, 1).unwrap();
}

fn run(bin: &str, args: &[&str]) {
    let output = Command::new(bin).args(args).output().unwrap();
    assert!(
        output.status.success(),
        "{bin} {args:?} failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

struct Run {
    /// (barcode, gene) -> count
    counts: BTreeMap<(String, String), f64>,
    is_real: bool,
    records: Vec<bam::Record>,
}

fn molecules_then_matrix(dir: &Path, name: &str, extra: &[&str]) -> Run {
    let input = dir.join("annotated.bam");
    if !input.exists() {
        write_fixture(&input);
    }
    let molecules = dir.join(format!("{name}.molecules.bam"));
    let outs = dir.join(format!("{name}.outs"));
    let mut args = vec![input.to_str().unwrap(), "-o", molecules.to_str().unwrap()];
    args.extend_from_slice(extra);
    run(
        &common::bin("sbx10x-molecules", env!("CARGO_BIN_EXE_sbx10x-molecules")).to_string_lossy(),
        &args,
    );
    run(
        &common::bin("cr-bam-to-mtx", env!("CARGO_BIN_EXE_cr-bam-to-mtx")).to_string_lossy(),
        &[
            molecules.to_str().unwrap(),
            "-o",
            outs.to_str().unwrap(),
            "--no-gzip",
        ],
    );
    assert!(
        !PathBuf::from(format!("{}.shards", molecules.display())).exists(),
        "shard directory must be cleaned up"
    );
    let mex = outs.join("raw_feature_bc_matrix");
    let lines = |f: &str| -> Vec<String> {
        std::fs::read_to_string(mex.join(f))
            .unwrap()
            .lines()
            .map(str::to_owned)
            .collect()
    };
    let features: Vec<String> = lines("features.tsv")
        .iter()
        .map(|l| l.split('\t').next().unwrap().to_owned())
        .collect();
    let barcodes = lines("barcodes.tsv");
    let mtx = lines("matrix.mtx");
    let is_real = mtx[0].contains(" real ");
    let mut counts = BTreeMap::new();
    for line in mtx.iter().filter(|l| !l.starts_with('%')).skip(1) {
        let f: Vec<&str> = line.split_whitespace().collect();
        let gene = features[f[0].parse::<usize>().unwrap() - 1].clone();
        let barcode = barcodes[f[1].parse::<usize>().unwrap() - 1].clone();
        counts.insert((barcode, gene), f[2].parse::<f64>().unwrap());
    }
    let records = bam::Reader::from_path(&molecules)
        .unwrap()
        .records()
        .map(Result::unwrap)
        .collect();
    Run {
        counts,
        is_real,
        records,
    }
}

fn expect(run: &Run, expected: &[(&str, &str, f64)]) {
    let want: BTreeMap<(String, String), f64> = expected
        .iter()
        .map(|&(bc, gene, n)| ((bc.to_owned(), gene.to_owned()), n))
        .collect();
    assert_eq!(
        run.counts.len(),
        want.len(),
        "got {:?}, want {want:?}",
        run.counts
    );
    for (k, &v) in &want {
        let got = run.counts.get(k).copied().unwrap_or(f64::NAN);
        assert!(
            (got - v).abs() < 1e-3,
            "{k:?}: got {got}, want {v} (all: {:?})",
            run.counts
        );
    }
}

fn tag_str(record: &bam::Record, tag: &[u8]) -> Option<String> {
    match record.aux(tag) {
        Ok(Aux::String(s)) => Some(s.to_owned()),
        _ => None,
    }
}

fn find<'a>(run: &'a Run, name: &str) -> &'a bam::Record {
    run.records
        .iter()
        .find(|r| r.qname() == name.as_bytes())
        .unwrap()
}

/// Invariants every mode must satisfy on the output BAM.
fn check_output_invariants(run: &Run) {
    assert_eq!(
        run.records.len(),
        READS.len(),
        "every input record is written"
    );
    for record in &run.records {
        assert!(matches!(record.aux(b"xf"), Ok(Aux::I32(17 | 25))));
        assert!(tag_str(record, b"UB").is_some());
        assert!(matches!(record.aux(b"XW"), Ok(Aux::Char(_))));
    }
    // annotate's ZG is preserved, and the multi-gene GX is materialized to
    // the ZG-selected gene so cr-bam-to-mtx accepts the representative.
    let r7 = find(run, "r7");
    assert_eq!(tag_str(r7, b"ZG").as_deref(), Some("G1"));
    assert_eq!(tag_str(r7, b"GX").as_deref(), Some("G1"));
    assert_eq!(tag_str(r7, b"GN").as_deref(), Some("g1"));
}

#[test]
fn single_cell_default_is_winner_take_all_with_umi_correction() {
    let dir = tempfile::tempdir().unwrap();
    for (name, extra) in [
        ("serial", vec![]),
        ("indexed", vec!["--evidence-workers", "2"]),
        ("sharded", vec!["--shards", "3", "--shard-concurrency", "2"]),
    ] {
        let run = molecules_then_matrix(dir.path(), name, &extra);
        check_output_invariants(&run);
        // G1: U1 (+U1_ERR corrected in, support 4) beats G3's 1 read; plus U3 via ZG.
        expect(&run, &[(C1, "G1", 2.0), (C1, "G2", 1.0), (C2, "G2", 1.0)]);
        assert!(
            !run.is_real,
            "{name}: winner-take-all must stay an integer matrix"
        );
        assert_eq!(
            tag_str(find(&run, "r4"), b"UB").as_deref(),
            Some(U1),
            "{name}"
        );
        assert!(
            matches!(find(&run, "r8").aux(b"XW"), Ok(Aux::Char(b'L'))),
            "{name}"
        );
        let representatives = run
            .records
            .iter()
            .filter(|r| matches!(r.aux(b"xf"), Ok(Aux::I32(25))))
            .count();
        assert_eq!(
            representatives, 4,
            "{name}: one representative per counted molecule"
        );
    }
}

#[test]
fn tagseq_fractional_opt_in_splits_without_umi_correction() {
    let dir = tempfile::tempdir().unwrap();
    for (name, extra) in [
        ("serial", vec!["--tagseq", "--fractional-collisions"]),
        (
            "sharded",
            vec!["--tagseq", "--fractional-collisions", "--shards", "2"],
        ),
    ] {
        let run = molecules_then_matrix(dir.path(), name, &extra);
        check_output_invariants(&run);
        // U1 collision: G1 has 3 reads, G3 has 1 -> 0.75 / 0.25.
        // U1_ERR is not corrected under --tagseq, so it is its own G1 molecule.
        expect(
            &run,
            &[
                (C1, "G1", 0.75 + 1.0 + 1.0),
                (C1, "G2", 1.0),
                (C1, "G3", 0.25),
                (C2, "G2", 1.0),
            ],
        );
        assert!(
            run.is_real,
            "{name}: fractional counts need a real-valued matrix"
        );
        let r8 = find(&run, "r8");
        assert!(matches!(r8.aux(b"XW"), Ok(Aux::Char(b'F'))));
        assert!(matches!(r8.aux(b"ZW"), Ok(Aux::Float(w)) if (w - 0.25).abs() < 1e-6));
        assert_eq!(tag_str(find(&run, "r4"), b"UB").as_deref(), Some(U1_ERR));
    }
}

#[test]
fn tagseq_without_fractional_counts_every_claimant_in_full() {
    let dir = tempfile::tempdir().unwrap();
    // Shared is the --tagseq default; the old opt-out flag is still accepted.
    for (name, args) in [
        ("default", vec!["--tagseq"]),
        (
            "legacy_flag",
            vec!["--tagseq", "--no-fractional-collisions"],
        ),
    ] {
        let run = molecules_then_matrix(dir.path(), name, &args);
        check_output_invariants(&run);
        expect(
            &run,
            &[
                (C1, "G1", 3.0),
                (C1, "G2", 1.0),
                (C1, "G3", 1.0),
                (C2, "G2", 1.0),
            ],
        );
        assert!(!run.is_real, "{name}");
        assert!(matches!(find(&run, "r8").aux(b"XW"), Ok(Aux::Char(b'S'))));
        assert!(
            find(&run, "r8").aux(b"ZW").is_err(),
            "ZW is only written for fractional molecules"
        );
    }
}

#[test]
fn fractional_single_cell_corrects_umis_then_splits() {
    let dir = tempfile::tempdir().unwrap();
    let run = molecules_then_matrix(dir.path(), "frac", &["--fractional-collisions"]);
    check_output_invariants(&run);
    // U1 group after correction: G1 support 4, G3 support 1 -> 0.8 / 0.2.
    expect(
        &run,
        &[
            (C1, "G1", 0.8 + 1.0),
            (C1, "G2", 1.0),
            (C1, "G3", 0.2),
            (C2, "G2", 1.0),
        ],
    );
    assert!(run.is_real);
}
