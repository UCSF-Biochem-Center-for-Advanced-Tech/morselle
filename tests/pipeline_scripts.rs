//! Guards for the shell pipeline drivers in `scripts/`.
//!
//! `scripts/sbxpress-pipeline.sh` shipped in 0.5.4 with a bash syntax error
//! and with flags no `sbxpress` subcommand accepts; nothing ran it. These
//! tests parse every script and drive both pipeline scripts end to end
//! against stub `sbxpress`/`samtools` commands that validate each flag
//! against the real binary's `--help` (see tests/scripts/).
use std::path::Path;
use std::process::Command;

fn repo() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
}

#[test]
fn every_script_parses() {
    let mut failures = Vec::new();
    for entry in std::fs::read_dir(repo().join("scripts")).unwrap() {
        let path = entry.unwrap().path();
        if path.extension().is_some_and(|e| e == "sh") {
            let out = Command::new("bash").arg("-n").arg(&path).output().unwrap();
            if !out.status.success() {
                failures.push(format!("{}: {}", path.display(), String::from_utf8_lossy(&out.stderr)));
            }
        }
    }
    assert!(failures.is_empty(), "bash -n failed:\n{}", failures.join("\n"));
}

#[test]
fn pipeline_scripts_pass_shellcheck_when_available() {
    if Command::new("shellcheck").arg("--version").output().is_err() {
        eprintln!("shellcheck not installed; skipping");
        return;
    }
    let out = Command::new("shellcheck")
        .args(["-S", "warning"])
        .arg(repo().join("scripts/sbxpress-pipeline.sh"))
        .arg(repo().join("scripts/sbxpress-gex-pipeline.sh"))
        .arg(repo().join("scripts/run_flex_multi.sh"))
        .output()
        .unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stdout));
}

#[test]
fn pipeline_scripts_run_end_to_end_against_stubs() {
    let out = Command::new("bash")
        .arg(repo().join("tests/scripts/pipeline_smoke.sh"))
        .arg(env!("CARGO_BIN_EXE_sbxpress"))
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}\n{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
}
