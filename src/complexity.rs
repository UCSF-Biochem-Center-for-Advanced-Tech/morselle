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
//! By default a flagged read is dropped (failure reason `low_complexity`)
//! and counted: a dead read 2 is a library fault, not biology, and should
//! not be rescued by trimming. `--drop-low-complexity false` emits it
//! instead, with `YL:i:1` in its FASTQ comment (minimap2 `-y` copies it into
//! the BAM), so `molecules`' per-gene artifact report can name the genes it
//! lands on. Either way the metrics and per-barcode table count it.

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

    /// Drop low-complexity reads (failure reason `low_complexity`). On by
    /// default; `--drop-low-complexity false` emits them tagged `YL:i:1`
    /// instead.
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set,
          num_args = 0..=1, default_missing_value = "true")]
    pub drop_low_complexity: bool,
}

impl Default for ComplexityCfg {
    fn default() -> Self {
        Self {
            low_complexity_base_fraction: 0.8,
            low_complexity_run_fraction: 0.5,
            drop_low_complexity: true,
        }
    }
}

impl ComplexityCfg {
    /// Whether `seq` fails either test. Empty reads never do.
    pub fn is_low(&self, seq: &[u8]) -> bool {
        self.check(seq).is_some()
    }

    /// The read's most frequent base (index into [`BASES`]) when `seq` fails
    /// either test, else `None`. One pass, no allocation.
    pub fn check(&self, seq: &[u8]) -> Option<usize> {
        if seq.is_empty() {
            return None;
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
        let (top_index, &top) = counts
            .iter()
            .enumerate()
            .max_by_key(|&(index, &count)| (count, std::cmp::Reverse(index)))?;
        let over = |count: f64, fraction: f64| fraction > 0.0 && count >= fraction * len;
        (over(top as f64, self.low_complexity_base_fraction)
            || over(longest as f64, self.low_complexity_run_fraction))
        .then_some(top_index)
    }
}

/// Base names for [`ComplexityCfg::check`] indices and [`ByBase`] slots.
pub const BASES: [&str; 5] = ["A", "C", "G", "T", "N"];

/// Low-complexity reads by most frequent base (`low_complexity_by_base` in
/// the conversion metrics JSON).
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct ByBase(pub [u64; 5]);

impl ByBase {
    pub fn add(&mut self, other: &ByBase) {
        for (a, b) in self.0.iter_mut().zip(other.0) {
            *a += b;
        }
    }

    /// G (two-colour "no signal") plus N (no call): the dead-read signal.
    pub fn g_and_n(&self) -> u64 {
        self.0[2] + self.0[4]
    }
}

impl serde::Serialize for ByBase {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeMap;
        let mut map = serializer.serialize_map(Some(5))?;
        for (name, count) in BASES.iter().zip(self.0) {
            map.serialize_entry(name, &count)?;
        }
        map.end()
    }
}

/// Share of barcode-passing reads above which conversion warns. DZEC05 (dead
/// read 2) flags ~51%; normal libraries flag <= 0.5% once poly-A read-through
/// is trimmed (DZEC04 with --trim-r2-polya: 0.48%, 85% of it G).
pub const WARN_FRACTION: f64 = 0.05;

/// Print a warning when the flagged reads are more than [`WARN_FRACTION`]
/// of `total`. With `g_and_n_only`, only G- and N-dominated reads count:
/// poly-A/T is expected wherever R2 read-through is not trimmed
/// (`illumina` without `--trim-r2-polya`: DZEC04 flagged 4.6%, 90% A, with
/// no fault), so it is warned on only then.
/// `dropped` says which of the two the run did with them.
pub fn warn_if_high(stage: &str, by_base: &ByBase, total: u64, g_and_n_only: bool, dropped: bool) {
    if let Some(text) = warning(stage, by_base, total, g_and_n_only, dropped) {
        eprintln!("{text}");
    }
}

/// The text [`warn_if_high`] prints, if any.
pub fn warning(
    stage: &str,
    by_base: &ByBase,
    total: u64,
    g_and_n_only: bool,
    dropped: bool,
) -> Option<String> {
    let (flagged, what) = if g_and_n_only {
        (by_base.g_and_n(), "poly-G/N")
    } else {
        (by_base.0.iter().sum(), "low-complexity")
    };
    (total > 0 && flagged as f64 > WARN_FRACTION * total as f64).then(|| {
        format!(
            "[{stage}] WARNING: {flagged} of {total} reads ({:.1}%) have {what} cDNA \
             (A {}, C {}, G {}, T {}, N {}; poly-G/N suggests a dead read 2, poly-A \
             untrimmed read-through); check the library or sequencing run. {}",
            100.0 * flagged as f64 / total as f64,
            by_base.0[0],
            by_base.0[1],
            by_base.0[2],
            by_base.0[3],
            by_base.0[4],
            if dropped {
                "These reads were dropped (failure reason low_complexity)."
            } else {
                "These reads are kept, tagged YL:i:1 (--drop-low-complexity false)."
            }
        )
    })
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
    fn check_reports_the_most_frequent_base() {
        let cfg = ComplexityCfg::default();
        assert_eq!(cfg.check(&[b'G'; 100]), Some(2));
        assert_eq!(cfg.check(&[b'N'; 60]), Some(4));
        let mut run = b"ACGTACGTAC".repeat(5);
        run.extend_from_slice(&[b'T'; 50]);
        assert_eq!(cfg.check(&run), Some(3));
        assert_eq!(cfg.check(b"GATTACACCATGGCTTGAACGTTAGCCATGCAAGTCC"), None);
        let by_base = ByBase([5, 0, 3, 1, 2]);
        assert_eq!(by_base.g_and_n(), 5);
        assert_eq!(
            serde_json::to_string(&by_base).unwrap(),
            r#"{"A":5,"C":0,"G":3,"T":1,"N":2}"#
        );
    }

    #[test]
    fn warning_counts_poly_a_only_when_read_through_is_untrimmed() {
        // DZEC04 without --trim-r2-polya: ~4.6% flagged, mostly poly-A.
        let poly_a = ByBase([90, 0, 7, 3, 0]);
        assert!(warning("t", &poly_a, 1000, false, true).is_some());
        assert!(warning("t", &poly_a, 1000, true, true).is_none());
        // DZEC05: dead read 2, poly-G either way.
        let poly_g = ByBase([5, 0, 480, 5, 20]);
        let text = warning("t", &poly_g, 1000, true, true).unwrap();
        assert!(text.contains("500 of 1000 reads (50.0%) have poly-G/N cDNA"));
        assert!(text.contains("were dropped"));
        let kept = warning("t", &poly_g, 1000, false, false).unwrap();
        assert!(kept.contains("kept, tagged YL:i:1"));
        assert!(warning("t", &ByBase::default(), 0, false, true).is_none());
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
