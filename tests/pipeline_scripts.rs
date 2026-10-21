//! Guards for the shell pipeline drivers in `scripts/`.
//!
//! `scripts/morselle-pipeline.sh` shipped in 0.5.4 with a bash syntax error
//! and with flags no `morselle` subcommand accepts; nothing ran it. These
//! tests parse every script and drive both pipeline scripts end to end
//! against stub `morselle`/`samtools` commands that validate each flag
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
                failures.push(format!(
                    "{}: {}",
                    path.display(),
                    String::from_utf8_lossy(&out.stderr)
                ));
            }
        }
    }
    assert!(
        failures.is_empty(),
        "bash -n failed:\n{}",
        failures.join("\n")
    );
}

/// install.sh in 0.5.2–0.6.0 installed only 7 of the 11 private engines, so
/// an installed `morselle` could not run probeset-align, probe-demux, guides
/// or protospacer-call (found by the 0.6.0 server acceptance run on CAT06).
/// Every engine the dispatcher can name must be in install.sh's HELPERS.
#[test]
fn install_script_installs_every_dispatched_helper() {
    let dispatcher = std::fs::read_to_string(repo().join("src/bin/morselle.rs")).unwrap();
    let body = &dispatcher[dispatcher
        .find("fn engine(")
        .expect("engine() in morselle.rs")..];
    let body = &body[..body.find("\n}\n").unwrap()];
    let mut dispatched: Vec<String> = body
        .match_indices("=> Some(\"")
        .map(|(i, m)| {
            let rest = &body[i + m.len()..];
            rest[..rest.find('"').unwrap()].to_string()
        })
        .collect();
    let out = Command::new("bash")
        .arg("-c")
        .arg(r#"eval "$(sed -n '/^HELPERS=(/,/)/p' "$1")"; printf '%s\n' "${HELPERS[@]}""#)
        .arg("install-helpers")
        .arg(repo().join("install.sh"))
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let mut installed: Vec<String> = String::from_utf8_lossy(&out.stdout)
        .lines()
        .map(str::to_string)
        .collect();
    dispatched.sort();
    installed.sort();
    assert!(
        dispatched.len() >= 11,
        "parsed only {dispatched:?} from engine()"
    );
    assert_eq!(
        installed, dispatched,
        "install.sh HELPERS must match the engines morselle dispatches to"
    );
}

#[test]
fn pipeline_scripts_pass_shellcheck_when_available() {
    if Command::new("shellcheck")
        .arg("--version")
        .output()
        .is_err()
    {
        eprintln!("shellcheck not installed; skipping");
        return;
    }
    let out = Command::new("shellcheck")
        .args(["-S", "warning"])
        .arg(repo().join("scripts/morselle-pipeline.sh"))
        .arg(repo().join("scripts/morselle-gex-pipeline.sh"))
        .arg(repo().join("scripts/run_flex_multi.sh"))
        .arg(repo().join("scripts/compare_tagseq_modes.sh"))
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stdout)
    );
}

#[test]
fn pipeline_scripts_run_end_to_end_against_stubs() {
    let out = Command::new("bash")
        .arg(repo().join("tests/scripts/pipeline_smoke.sh"))
        .arg(env!("CARGO_BIN_EXE_morselle"))
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}\n{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
}

#[test]
fn pipeline_restart_and_reporting_python_regressions() {
    let out = Command::new("python3")
        .arg(repo().join("tests/scripts/python_regressions.py"))
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}\n{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
}
