//! Golden-output regression test for `molecules` + `matrix`.
//!
//! Runs both binaries on deterministic fixtures in every execution path
//! and collision policy, and compares a digest of every output record (all
//! tags, in order), the behaviour-defining metrics counters, and a digest
//! of the raw MEX matrix against `tests/golden/outputs.json`.
//!
//! A refactor or optimisation must leave this test passing untouched. If a
//! change is *meant* to alter output, regenerate and explain why in the
//! commit message:
//!
//!     SBXPRESS_UPDATE_GOLDEN=1 cargo test --release --test golden_outputs
//!
//! The non-fractional cases were seeded from, and are still reproduced by,
//! the 0.5.6 release binaries. To re-verify against any other build:
//!
//!     SBXPRESS_TEST_BIN_DIR=/path/to/old/target/release SBXPRESS_GOLDEN_LEGACY=1 \
//!         cargo test --release --test golden_outputs
//!
//! (`SBXPRESS_GOLDEN_LEGACY` maps flags to their pre-0.5.9 spelling and skips
//! the fractional cases, which 0.5.6 did not have.)
mod common;

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use common::{bam_digest, bin, mex_digest, molecule_counters, run, write_fixture, FixtureSpec};
use serde_json::{json, Value};

struct Case {
    name: &'static str,
    flex_fixture: bool,
    args: &'static [&'static str],
    /// Equivalent 0.5.6 flags, or `None` if 0.5.6 had no such mode.
    legacy_args: Option<&'static [&'static str]>,
}

const CASES: &[Case] = &[
    Case {
        name: "gex_serial",
        flex_fixture: false,
        args: &[],
        legacy_args: Some(&[]),
    },
    Case {
        name: "gex_indexed",
        flex_fixture: false,
        args: &["--evidence-workers", "2"],
        legacy_args: Some(&["--evidence-workers", "2"]),
    },
    Case {
        name: "gex_sharded",
        flex_fixture: false,
        args: &["--shards", "3", "--shard-concurrency", "2"],
        legacy_args: Some(&["--shards", "3", "--shard-concurrency", "2"]),
    },
    Case {
        name: "gex_saturation",
        flex_fixture: false,
        args: &[
            "--saturation",
            "{dir}/sat.tsv.gz",
            "--saturation-interval",
            "500",
        ],
        legacy_args: Some(&[
            "--saturation",
            "{dir}/sat.tsv.gz",
            "--saturation-interval",
            "500",
        ]),
    },
    Case {
        name: "tagseq_shared",
        flex_fixture: false,
        args: &["--tagseq", "--no-fractional-collisions"],
        legacy_args: Some(&["--tagseq"]),
    },
    Case {
        name: "tagseq_fractional",
        flex_fixture: false,
        args: &["--tagseq"],
        legacy_args: None,
    },
    Case {
        name: "fractional_single_cell",
        flex_fixture: false,
        args: &["--fractional-collisions"],
        legacy_args: None,
    },
    Case {
        name: "flex",
        flex_fixture: true,
        args: &["--flex"],
        legacy_args: Some(&["--flex"]),
    },
    Case {
        name: "flex_sharded",
        flex_fixture: true,
        args: &["--flex", "--shards", "2"],
        legacy_args: Some(&["--flex", "--shards", "2"]),
    },
];

fn golden_path() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/outputs.json")
}

fn run_case(dir: &Path, fixture: &Path, case: &Case, args: &[&str]) -> Value {
    let out = dir.join(format!("{}.bam", case.name));
    let args: Vec<String> = args
        .iter()
        .map(|a| a.replace("{dir}", dir.to_str().unwrap()))
        .collect();
    let mut all = vec![
        fixture.to_str().unwrap().to_owned(),
        "-o".into(),
        out.to_str().unwrap().to_owned(),
    ];
    all.extend(args);
    let all: Vec<&str> = all.iter().map(String::as_str).collect();
    run(
        &bin("sbx10x-molecules", env!("CARGO_BIN_EXE_sbx10x-molecules")),
        &all,
    );
    let outs = dir.join(format!("{}.outs", case.name));
    run(
        &bin("cr-bam-to-mtx", env!("CARGO_BIN_EXE_cr-bam-to-mtx")),
        &[
            out.to_str().unwrap(),
            "-o",
            outs.to_str().unwrap(),
            "--no-gzip",
        ],
    );
    let counters = molecule_counters(&PathBuf::from(format!("{}.metrics.json", out.display())));
    let mut result = json!({
        "bam": bam_digest(&out),
        "matrix": mex_digest(&outs.join("raw_feature_bc_matrix")),
        "counters": counters,
    });
    let saturation = dir.join("sat.tsv.gz");
    if saturation.exists() {
        result["saturation"] = json!(common::digest(&common::read_text(&dir.join("sat.tsv"))));
        std::fs::remove_file(saturation).unwrap();
    }
    result
}

#[test]
fn molecules_and_matrix_outputs_match_golden() {
    let update = std::env::var_os("SBXPRESS_UPDATE_GOLDEN").is_some();
    let legacy = std::env::var_os("SBXPRESS_GOLDEN_LEGACY").is_some();
    let dir = tempfile::tempdir().unwrap();
    let gex_dir = dir.path().join("gex");
    let flex_dir = dir.path().join("flex");
    std::fs::create_dir_all(&gex_dir).unwrap();
    std::fs::create_dir_all(&flex_dir).unwrap();
    // Small UMI pool: plenty of collisions for every policy to act on.
    let gex = write_fixture(
        &gex_dir,
        FixtureSpec {
            umi_pool: 300,
            ..FixtureSpec::default()
        },
    );
    let flex = write_fixture(
        &flex_dir,
        FixtureSpec {
            umi_pool: 300,
            probe_barcodes: Some(3),
            ..FixtureSpec::default()
        },
    );

    let mut results = BTreeMap::new();
    for case in CASES {
        let args = if legacy {
            match case.legacy_args {
                Some(args) => args,
                None => continue,
            }
        } else {
            case.args
        };
        let (dir, fixture) = if case.flex_fixture {
            (&flex_dir, &flex.bam)
        } else {
            (&gex_dir, &gex.bam)
        };
        results.insert(case.name.to_owned(), run_case(dir, fixture, case, args));
    }

    if update {
        assert!(
            !legacy,
            "update golden files from the current build, not a legacy one"
        );
        std::fs::create_dir_all(golden_path().parent().unwrap()).unwrap();
        std::fs::write(
            golden_path(),
            serde_json::to_string_pretty(&results).unwrap() + "\n",
        )
        .unwrap();
        eprintln!("updated {}", golden_path().display());
        return;
    }
    let golden: BTreeMap<String, Value> = serde_json::from_str(
        &std::fs::read_to_string(golden_path()).expect("tests/golden/outputs.json"),
    )
    .unwrap();
    let mut failures = Vec::new();
    for (name, got) in &results {
        let Some(want) = golden.get(name) else {
            failures.push(format!("{name}: no golden entry"));
            continue;
        };
        if got != want {
            failures.push(format!(
                "{name}:\n    golden: {}\n    got:    {}",
                serde_json::to_string(want).unwrap(),
                serde_json::to_string(got).unwrap()
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "molecules/matrix output changed versus tests/golden/outputs.json.\n\
         If this is unintended, it is a regression. If intended, regenerate with\n\
         SBXPRESS_UPDATE_GOLDEN=1 cargo test --release --test golden_outputs\n\
         and say why in the commit.\n\n{}",
        failures.join("\n")
    );
}
