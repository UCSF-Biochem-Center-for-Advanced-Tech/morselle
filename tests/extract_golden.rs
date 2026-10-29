//! Golden-output regression test for `extract` (SBX arm).
//!
//! `tests/golden_outputs.rs` pins `molecules` + `matrix`; nothing pinned what
//! `extract` itself writes. This runs it on the real 100-read smoke fixture
//! with the real CE15204 sample sheet (no whitelist: it is not redistributed,
//! so CI does not have it) and compares a digest of every output FASTQ (records
//! sorted, since the pipelined path does not keep input order) plus a fixed set
//! of behaviour-defining counters against `tests/golden/extract.json`.
//!
//! The counter list is fixed on purpose: adding a new metric must not change
//! this golden, but changing how an existing read is classified, demultiplexed
//! or written must. New default-off options must leave it untouched. If a
//! change is *meant* to alter output, regenerate and explain why in the commit:
//!
//!     MORSELLE_UPDATE_GOLDEN=1 cargo test --release --test extract_golden
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use clap::Parser;
use morselle::config::Cli;
use serde_json::{json, Value};

mod common;

struct Case {
    name: &'static str,
    args: &'static [&'static str],
}

const CASES: &[Case] = &[
    Case {
        name: "smoke_sheet",
        args: &[],
    },
    Case {
        name: "smoke_sheet_tags_rigid_index",
        args: &["--processing-tags", "--index-window", "0"],
    },
    Case {
        name: "smoke_sheet_both_indexes",
        args: &["--processing-tags", "--demux-index", "both"],
    },
    // One-indel index matching (added with --index-indels; the cases above
    // were generated before it existed and must not change).
    Case {
        name: "smoke_sheet_index_indels",
        args: &["--processing-tags", "--index-indels"],
    },
];

/// Counters pinned by the golden (top level, then `demux.*`).
const TOP_COUNTERS: &[&str] = &[
    "total_reads",
    "oriented",
    "class_full",
    "class_partial",
    "class_unparseable",
    "pairs_emitted",
    "r1_eq_28",
    "r1_not_28",
    "adapter_trimmed",
    "r2_len_mean",
    "r2_len_max",
    "reads_failed",
];
const DEMUX_COUNTERS: &[&str] = &[
    "assigned",
    "undetermined",
    "via_dual",
    "via_i5_only",
    "via_i7_only",
    "via_conflict",
    "assigned_early",
    "undetermined_early",
    "index_unreachable",
    "per_sample",
    "per_sample_early",
];

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

fn golden_path() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/extract.json")
}

/// Digest of a FASTQ with its records sorted (order-independent).
fn fastq_digest(path: &Path) -> String {
    let text = common::read_text(path);
    let lines: Vec<&str> = text.lines().collect();
    let mut records: Vec<String> = lines.chunks(4).map(|r| r.join("\n")).collect();
    records.sort();
    common::digest(&records.join("\n"))
}

fn run_case(dir: &Path, case: &Case) -> Value {
    let out = dir.join(case.name);
    let failed = dir.join(format!("{}.failed.fastq.gz", case.name));
    let mut argv: Vec<String> = [
        "morselle",
        "--fastqs",
        fixture("smoke.fastq.gz").to_str().unwrap(),
        "--out-dir",
        out.to_str().unwrap(),
        "--bclconvert-samplesheet",
        fixture("CE15204_bclconvert_input.csv").to_str().unwrap(),
        "--threads",
        "2",
        "--failed-output",
        failed.to_str().unwrap(),
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    argv.extend(case.args.iter().map(|s| s.to_string()));
    morselle::run(&Cli::parse_from(&argv)).unwrap();

    let mut files = BTreeMap::new();
    let mut names: Vec<_> = std::fs::read_dir(&out)
        .unwrap()
        .map(|e| e.unwrap().file_name().into_string().unwrap())
        .filter(|n| n.ends_with(".fastq.gz"))
        .collect();
    names.sort();
    for name in names {
        let plain = out.join(name.trim_end_matches(".gz"));
        files.insert(name, json!(fastq_digest(&plain)));
    }
    files.insert(
        "failed.fastq.gz".to_owned(),
        json!(fastq_digest(
            &dir.join(format!("{}.failed.fastq", case.name))
        )),
    );

    let metrics: Value =
        serde_json::from_str(&common::read_text(&out.join("metrics_summary.json"))).unwrap();
    let mut counters = BTreeMap::new();
    for key in TOP_COUNTERS {
        counters.insert(key.to_string(), metrics[key].clone());
    }
    for key in DEMUX_COUNTERS {
        counters.insert(format!("demux.{key}"), metrics["demux"][key].clone());
    }
    json!({ "files": files, "counters": counters })
}

#[test]
fn extract_outputs_match_golden() {
    let update = std::env::var_os("MORSELLE_UPDATE_GOLDEN").is_some();
    let dir = tempfile::tempdir().unwrap();
    let results: BTreeMap<String, Value> = CASES
        .iter()
        .map(|c| (c.name.to_owned(), run_case(dir.path(), c)))
        .collect();

    if update {
        std::fs::write(
            golden_path(),
            serde_json::to_string_pretty(&results).unwrap() + "\n",
        )
        .unwrap();
        eprintln!("updated {}", golden_path().display());
        return;
    }
    let golden: BTreeMap<String, Value> = serde_json::from_str(
        &std::fs::read_to_string(golden_path()).expect("tests/golden/extract.json"),
    )
    .unwrap();
    let mut failures = Vec::new();
    for (name, got) in &results {
        match golden.get(name) {
            None => failures.push(format!("{name}: no golden entry")),
            Some(want) if want != got => failures.push(format!(
                "{name}:\n    golden: {}\n    got:    {}",
                serde_json::to_string(want).unwrap(),
                serde_json::to_string(got).unwrap()
            )),
            _ => {}
        }
    }
    assert!(
        failures.is_empty(),
        "extract output changed versus tests/golden/extract.json.\n\
         If this is unintended, it is a regression. If intended, regenerate with\n\
         MORSELLE_UPDATE_GOLDEN=1 cargo test --release --test extract_golden\n\
         and say why in the commit.\n\n{}",
        failures.join("\n")
    );
}

/// The early index-evidence counters partition the R1-located reads three
/// ways (full/partial, seen, identified) and agree with the via breakdown.
#[test]
fn index_evidence_counters_are_consistent() {
    let dir = tempfile::tempdir().unwrap();
    for case in CASES {
        let out = dir.path().join(case.name);
        run_case(dir.path(), case);
        let m: Value =
            serde_json::from_str(&common::read_text(&out.join("metrics_summary.json"))).unwrap();
        let d = &m["demux"];
        let ie = &d["index_evidence"];
        let n = |k: &str| ie[k].as_u64().unwrap();
        let located = m["total_reads"].as_u64().unwrap() - d["index_unreachable"].as_u64().unwrap();
        assert_eq!(n("reads"), located, "{}", case.name);
        assert_eq!(n("full") + n("partial"), n("reads"));
        assert_eq!(
            n("seen_both") + n("seen_i5_only") + n("seen_i7_only") + n("seen_none"),
            n("reads")
        );
        assert_eq!(
            n("identified_both")
                + n("identified_i5_only")
                + n("identified_i7_only")
                + n("identified_none"),
            n("reads")
        );
        // An identified index is always counted as seen; full needs both.
        assert!(n("full") <= n("identified_both"));
        assert!(n("identified_both") <= n("seen_both"));
        let indels = case.args.contains(&"--index-indels");
        assert_eq!(d["index_indels_i5"].as_bool().unwrap(), indels);
        assert_eq!(d["index_indels_i7"].as_bool().unwrap(), indels);
        if !indels {
            assert_eq!(
                n("identified_i5_via_indel") + n("identified_i7_via_indel"),
                0
            );
        }
    }
}
