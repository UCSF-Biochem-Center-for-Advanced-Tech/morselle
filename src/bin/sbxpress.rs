//! Unified public command for the SBX gene-expression workflow.
//!
//! The mature engines remain separate executables internally for now. Keeping
//! dispatch at the process boundary preserves their streaming and memory
//! behavior while presenting one stable user-facing CLI.

use std::env;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode};

const HELP: &str = "sbxpress — SBX single-cell gene-expression processing\n\n\
Usage: sbxpress <COMMAND> [ARGUMENTS]\n\n\
Commands:\n\
  extract    Extract and demultiplex tagged cDNA FASTQs from raw SBX reads\n\
  illumina   Convert paired Illumina 10x FASTQs to tagged R2-only FASTQ\n\
  minimap    Run minimap2 with transparent argument and stream forwarding\n\
  annotate   Add gene, transcript, region, and collision annotations to BAM\n\
  molecules  Correct UMIs, resolve gene collisions, and mark representatives\n\
  audit      Audit mapping, annotation, trimming, and molecule metrics\n\
  matrix     Create 10x-compatible raw feature-barcode matrices\n\
  cells      Call cells and create filtered feature-barcode matrices\n\
  probe-demux  Assign cell barcodes to samples by probe barcode (10x Flex)\n\
  guides     Assign Feature Barcoding reads (CRISPR guides, antibody/CMO tags) to features\n\
  protospacer-call  Call CRISPR guide identity per cell (adaptive UMI threshold)\n\
  probeset-align  Assign 10x Flex GEX reads to genes via the real probe-set reference\n\n\
Run 'sbxpress <COMMAND> --help' for command-specific options.\n";

fn helper_dir(executable: &Path) -> Option<PathBuf> {
    if let Some(value) = env::var_os("SBXPRESS_LIBEXEC") {
        return Some(PathBuf::from(value));
    }
    let bin = executable.parent()?;
    // Installed layout: <prefix>/bin/sbxpress and <prefix>/libexec/sbxpress/*.
    let installed = bin.parent()?.join("libexec/sbxpress");
    if installed.is_dir() {
        Some(installed)
    } else {
        // Cargo layout: target/{debug,release}/ contains every binary.
        Some(bin.to_path_buf())
    }
}

fn engine(command: &str) -> Option<&'static str> {
    match command {
        "extract" => Some("sbx10x"),
        "illumina" => Some("sbxpress-illumina"),
        "annotate" => Some("sbx10x-annotate"),
        "molecules" => Some("sbx10x-molecules"),
        "audit" => Some("bam-audit"),
        "matrix" => Some("cr-bam-to-mtx"),
        "cells" => Some("sbxpress-cells"),
        "probe-demux" => Some("sbxpress-probe-demux"),
        "guides" => Some("sbxpress-guides"),
        "protospacer-call" => Some("sbxpress-protospacer-call"),
        "probeset-align" => Some("sbxpress-probeset-align"),
        _ => None,
    }
}

fn run(mut command: Command) -> ExitCode {
    match command.status() {
        Ok(status) => ExitCode::from(status.code().unwrap_or(1) as u8),
        Err(error) => {
            eprintln!("sbxpress: {error}");
            ExitCode::FAILURE
        }
    }
}

fn main() -> ExitCode {
    let mut args = env::args_os();
    let argv0 = PathBuf::from(args.next().unwrap_or_else(|| "sbxpress".into()));
    let executable = env::current_exe().unwrap_or(argv0);
    let Some(subcommand) = args.next() else {
        print!("{HELP}");
        return ExitCode::SUCCESS;
    };
    let subcommand_text = subcommand.to_string_lossy();
    if matches!(subcommand_text.as_ref(), "-h" | "--help" | "help") {
        print!("{HELP}");
        return ExitCode::SUCCESS;
    }
    if matches!(subcommand_text.as_ref(), "-V" | "--version" | "version") {
        println!("sbxpress {}", env!("CARGO_PKG_VERSION"));
        return ExitCode::SUCCESS;
    }
    if subcommand_text == "minimap" {
        let mut child = Command::new("minimap2");
        child.args(args);
        return run(child);
    }
    let Some(helper) = engine(&subcommand_text) else {
        eprintln!("sbxpress: unknown command '{subcommand_text}'\n\n{HELP}");
        return ExitCode::from(2);
    };
    let Some(directory) = helper_dir(&executable) else {
        eprintln!("sbxpress: cannot determine helper executable directory");
        return ExitCode::FAILURE;
    };
    let helper = directory.join(helper);
    if !helper.is_file() {
        eprintln!(
            "sbxpress: helper not found: {} (set SBXPRESS_LIBEXEC to override)",
            helper.display()
        );
        return ExitCode::FAILURE;
    }
    let mut child = Command::new(helper);
    child.env("SBXPRESS_DISPATCH", "1");
    child.args(args);
    run(child)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_public_engine_has_a_private_helper() {
        assert_eq!(engine("extract"), Some("sbx10x"));
        assert_eq!(engine("illumina"), Some("sbxpress-illumina"));
        assert_eq!(engine("annotate"), Some("sbx10x-annotate"));
        assert_eq!(engine("molecules"), Some("sbx10x-molecules"));
        assert_eq!(engine("audit"), Some("bam-audit"));
        assert_eq!(engine("matrix"), Some("cr-bam-to-mtx"));
        assert_eq!(engine("cells"), Some("sbxpress-cells"));
        assert_eq!(engine("probe-demux"), Some("sbxpress-probe-demux"));
        assert_eq!(engine("guides"), Some("sbxpress-guides"));
        assert_eq!(
            engine("protospacer-call"),
            Some("sbxpress-protospacer-call")
        );
        assert_eq!(engine("probeset-align"), Some("sbxpress-probeset-align"));
        assert_eq!(engine("minimap"), None);
    }
}
