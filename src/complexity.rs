//! Low-complexity cDNA check at conversion (`illumina`, `extract`).
//!
//! DZEC05 (2026-10) had 28-56% of read pairs with a dead read 2: a two-colour
//! sequencer calls "no signal" as G, at fairly high quality (poly-G R2 reads
//! average Q25), so base quality does not catch them. Most did not map, but
//! 6.7M aligned to G runs inside LINC00486 and gave it 18-21% of all
//! molecules. The conversion stages already hold every cDNA read in memory,
//! so one pass over its bases is the cheapest place to see this. minimap2
//! masks low-complexity seeds internally but reports nothing per read.
//!
//! By default a flagged read is still emitted, with `YL:i:1` in its FASTQ
//! comment (minimap2 `-y` copies it into the BAM), and counted, so a normal
//! library's counts do not change and a bad one is visible in the metrics,
//! the per-barcode table and `molecules`' per-gene artifact report.
//! `--drop-low-complexity` sends flagged reads to the failed output instead.

/// FASTQ comment / BAM tag marking a low-complexity cDNA read.
pub const TAG: &[u8] = b"\tYL:i:1";

/// Thresholds for [`ComplexityCfg::is_low`] (shared by `illumina` and `extract`).
#[derive(clap::Args, Debug, Clone, Copy, PartialEq)]
pub struct ComplexityCfg {
    /// Flag a cDNA read as low-complexity when one base (A, C, G, T or N)
    /// makes up at least this fraction of it. 0 disables this test.
    #[arg(long, default_value_t = 0.8)]
    pub low_complexity_base_fraction: f64,

    /// Flag a cDNA read as low-complexity when its longest single-base run
    /// covers at least this fraction of it. 0 disables this test.
    #[arg(long, default_value_t = 0.5)]
    pub low_complexity_run_fraction: f64,

    /// Drop low-complexity reads (failure reason `low_complexity`) instead
    /// of emitting them tagged `YL:i:1`.
    #[arg(long, default_value_t = false)]
    pub drop_low_complexity: bool,
}

impl Default for ComplexityCfg {
    fn default() -> Self {
        Self {
            low_complexity_base_fraction: 0.8,
            low_complexity_run_fraction: 0.5,
            drop_low_complexity: false,
        }
    }
}

impl ComplexityCfg {
    /// Whether `seq` fails either test. Empty reads never do. One pass, no
    /// allocation.
    pub fn is_low(&self, seq: &[u8]) -> bool {
        if seq.is_empty() {
            return false;
        }
        let mut counts = [0usize; 5];
        let (mut run, mut longest) = (0usize, 0usize);
        let mut previous = 0u8;
        for &base in seq {
            let base = base.to_ascii_uppercase();
            counts[match base {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 4,
            }] += 1;
            run = if base == previous { run + 1 } else { 1 };
            longest = longest.max(run);
            previous = base;
        }
        let len = seq.len() as f64;
        let top = *counts.iter().max().unwrap_or(&0) as f64;
        let over = |count: f64, fraction: f64| fraction > 0.0 && count >= fraction * len;
        over(top, self.low_complexity_base_fraction)
            || over(longest as f64, self.low_complexity_run_fraction)
    }
}

/// Share of barcode-passing reads above which conversion warns. Tentative:
/// set from the 0.6.2 validation (normal libraries vs DZEC05).
pub const WARN_FRACTION: f64 = 0.05;

/// Print a warning when `flagged` is more than [`WARN_FRACTION`] of `total`.
pub fn warn_if_high(stage: &str, flagged: u64, total: u64) {
    if total > 0 && flagged as f64 > WARN_FRACTION * total as f64 {
        eprintln!(
            "[{stage}] WARNING: {flagged} of {total} reads ({:.1}%) have low-complexity cDNA \
             (e.g. poly-G from a dead read 2); check the library or sequencing run. \
             These reads are tagged YL:i:1 (--drop-low-complexity removes them).",
            100.0 * flagged as f64 / total as f64
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn poly_g_and_long_runs_are_low_complexity() {
        let cfg = ComplexityCfg::default();
        assert!(cfg.is_low(&[b'G'; 100]));
        // 80 G among 100 bases, scattered: base-fraction test.
        let scattered: Vec<u8> = (0..100)
            .map(|i| if i % 5 == 0 { b'A' } else { b'G' })
            .collect();
        assert!(cfg.is_low(&scattered));
        // A 50-base run of T in otherwise mixed sequence: run test.
        let mut run = b"ACGTACGTAC".repeat(5);
        run.extend_from_slice(&[b'T'; 50]);
        assert!(cfg.is_low(&run));
        // All N (no-call read).
        assert!(cfg.is_low(&[b'N'; 60]));
        // Lower-case input is treated like upper case.
        assert!(cfg.is_low(&[b'g'; 60]));
    }

    #[test]
    fn ordinary_cdna_is_not_low_complexity() {
        let cfg = ComplexityCfg::default();
        assert!(!cfg.is_low(b"GATTACACCATGGCTTGAACGTTAGCCATGCAAGTCCGATACGGATTTACAGGCATCAG"));
        // A real 3' read with a 20-base poly-A tail in 90 bases.
        let mut tail =
            b"GATTACACCATGGCTTGAACGTTAGCCATGCAAGTCCGATACGGATTTACAGGCATCAGTTCCGAA".to_vec();
        tail.extend_from_slice(&[b'A'; 20]);
        assert!(!cfg.is_low(&tail));
        assert!(!cfg.is_low(b""));
    }

    #[test]
    fn thresholds_of_zero_disable_each_test() {
        let off = ComplexityCfg {
            low_complexity_base_fraction: 0.0,
            low_complexity_run_fraction: 0.0,
            drop_low_complexity: false,
        };
        assert!(!off.is_low(&[b'G'; 100]));
        let run_only = ComplexityCfg {
            low_complexity_base_fraction: 0.0,
            ..ComplexityCfg::default()
        };
        // 80% G but no run reaches half the read: only the base test fires.
        let scattered: Vec<u8> = (0..100)
            .map(|i| if i % 5 == 0 { b'A' } else { b'G' })
            .collect();
        assert!(!run_only.is_low(&scattered));
    }
}
