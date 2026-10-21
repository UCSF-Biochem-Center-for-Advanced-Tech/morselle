//! Probe-barcode sample sheet parsing (the Flex design notes Phase 2).
//!
//! Two front-end formats feed one internal model, mirroring `sample_sheet.rs`'s
//! bclconvert/flat split:
//! - **Primary:** a real cellranger `multi` `config.csv` — only its `[samples]`
//!   section is read (`sample_id,probe_barcode_ids,...`; other columns and
//!   every other section are ignored), so a user can point `morselle`
//!   straight at an existing config file. `probe_barcode_ids` is `|`-delimited
//!   when a sample is pooled across more than one probe barcode (real
//!   example: `Pos,D-A01|D-A02`).
//! - **Fallback:** flat CSV, header `sample_id,probe_barcode_ids`, same
//!   `|`-delimited convention in that column.

use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::Path;

use anyhow::{bail, Context, Result};

/// One sample's declared probe-barcode identity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProbeSample {
    pub id: String,
    /// Declared probe-barcode ids (e.g. `["D-A01", "D-A02"]`) — a sample can
    /// be pooled across more than one probe barcode.
    pub probe_barcode_ids: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct ProbeSampleSheet {
    pub samples: Vec<ProbeSample>,
}

impl ProbeSampleSheet {
    /// Load from a path, auto-detecting format by the presence of a
    /// `[samples]` section header.
    pub fn load(path: &Path) -> Result<Self> {
        let text = fs::read_to_string(path)
            .with_context(|| format!("cannot read probe sample sheet {}", path.display()))?;
        let sheet = if text.lines().any(|l| l.trim() == "[samples]") {
            parse_config_csv(&text)
        } else {
            parse_flat(&text)
        }
        .with_context(|| format!("parsing probe sample sheet {}", path.display()))?;
        sheet.validate()?;
        Ok(sheet)
    }

    fn validate(&self) -> Result<()> {
        if self.samples.is_empty() {
            bail!("probe sample sheet has no samples");
        }
        let mut seen_ids = HashSet::new();
        let mut seen_probes: HashMap<String, String> = HashMap::new();
        for s in &self.samples {
            if !seen_ids.insert(s.id.clone()) {
                bail!("duplicate sample id '{}'", s.id);
            }
            if s.probe_barcode_ids.is_empty() {
                bail!("sample '{}' declares no probe_barcode_ids", s.id);
            }
            for pb in &s.probe_barcode_ids {
                if let Some(other) = seen_probes.insert(pb.clone(), s.id.clone()) {
                    bail!(
                        "probe barcode '{pb}' is declared by both '{other}' and '{}'",
                        s.id
                    );
                }
            }
        }
        Ok(())
    }
}

/// Parse a real cellranger `multi` `config.csv`, reading only `[samples]`.
fn parse_config_csv(text: &str) -> Result<ProbeSampleSheet> {
    let mut section = String::new();
    let mut header: Vec<String> = Vec::new();
    let mut samples = Vec::new();

    for raw in text.lines() {
        let line = raw.trim_end_matches('\r').trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('[') && line.ends_with(']') {
            section = line[1..line.len() - 1].to_ascii_lowercase();
            header.clear();
            continue;
        }
        if section != "samples" {
            continue;
        }
        let fields: Vec<&str> = line.split(',').map(str::trim).collect();
        if header.is_empty() {
            header = fields.iter().map(|s| s.to_ascii_lowercase()).collect();
            continue;
        }
        let row: HashMap<&str, &str> = header
            .iter()
            .map(String::as_str)
            .zip(fields.iter().copied())
            .collect();
        let id = row
            .get("sample_id")
            .copied()
            .context("[samples] row missing sample_id")?
            .to_string();
        let probes = row
            .get("probe_barcode_ids")
            .copied()
            .context("[samples] row missing probe_barcode_ids")?;
        samples.push(ProbeSample {
            id,
            probe_barcode_ids: split_probe_ids(probes),
        });
    }
    Ok(ProbeSampleSheet { samples })
}

/// Parse the flat fallback CSV: header `sample_id,probe_barcode_ids`.
fn parse_flat(text: &str) -> Result<ProbeSampleSheet> {
    let mut lines = text
        .lines()
        .map(|l| l.trim_end_matches('\r').trim())
        .filter(|l| !l.is_empty());
    let header = lines.next().context("empty probe sample sheet")?;
    let cols: Vec<String> = header
        .split(',')
        .map(|s| s.trim().to_ascii_lowercase())
        .collect();
    let col = |name: &str| cols.iter().position(|c| c == name);
    let (si, pi) = (
        col("sample_id").context("flat sheet missing 'sample_id' column")?,
        col("probe_barcode_ids").context("flat sheet missing 'probe_barcode_ids' column")?,
    );

    let mut samples = Vec::new();
    for line in lines {
        let f: Vec<&str> = line.split(',').map(str::trim).collect();
        let max = si.max(pi);
        if f.len() <= max {
            bail!("flat sheet row has too few columns: '{line}'");
        }
        samples.push(ProbeSample {
            id: f[si].to_string(),
            probe_barcode_ids: split_probe_ids(f[pi]),
        });
    }
    Ok(ProbeSampleSheet { samples })
}

fn split_probe_ids(field: &str) -> Vec<String> {
    field
        .split('|')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(String::from)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_temp(text: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(text.as_bytes()).unwrap();
        f
    }

    /// Trimmed shape of the real CAT06 config.csv (same [samples] layout,
    /// including a multi-probe-barcode pool and a single-probe-barcode
    /// sample), plus preceding sections that must be ignored.
    const CONFIG_CSV: &str = "\
[gene-expression]
reference,/some/ref
chemistry,auto

[libraries]
fastq_id,fastqs,lanes,feature_types
L1_GEX,/some/fastqs,any,Gene Expression

[samples]
sample_id,probe_barcode_ids,description,expect_cells,force_cells
Pos,D-A01|D-A02,,,
RZ213NR,D-A07,,,
Neg,D-B06|D-B07,,,
";

    #[test]
    fn parses_real_config_csv_shape() {
        let f = write_temp(CONFIG_CSV);
        let sheet = ProbeSampleSheet::load(f.path()).unwrap();
        assert_eq!(sheet.samples.len(), 3);
        assert_eq!(sheet.samples[0].id, "Pos");
        assert_eq!(sheet.samples[0].probe_barcode_ids, vec!["D-A01", "D-A02"]);
        assert_eq!(sheet.samples[1].id, "RZ213NR");
        assert_eq!(sheet.samples[1].probe_barcode_ids, vec!["D-A07"]);
    }

    #[test]
    fn parses_flat_fallback() {
        let f = write_temp("sample_id,probe_barcode_ids\nPos,D-A01|D-A02\nNeg,D-B06|D-B07\n");
        let sheet = ProbeSampleSheet::load(f.path()).unwrap();
        assert_eq!(sheet.samples.len(), 2);
        assert_eq!(sheet.samples[0].probe_barcode_ids, vec!["D-A01", "D-A02"]);
    }

    #[test]
    fn rejects_duplicate_sample_id() {
        let f = write_temp("sample_id,probe_barcode_ids\nA,D-A01\nA,D-A02\n");
        let err = ProbeSampleSheet::load(f.path()).unwrap_err();
        assert!(err.to_string().contains("duplicate sample id"));
    }

    #[test]
    fn rejects_probe_barcode_shared_across_samples() {
        let f = write_temp("sample_id,probe_barcode_ids\nA,D-A01\nB,D-A01\n");
        let err = ProbeSampleSheet::load(f.path()).unwrap_err();
        assert!(err.to_string().contains("declared by both"));
    }

    #[test]
    fn rejects_empty_probe_barcode_ids() {
        let f = write_temp("sample_id,probe_barcode_ids\nA,\n");
        let err = ProbeSampleSheet::load(f.path()).unwrap_err();
        assert!(err.to_string().contains("declares no probe_barcode_ids"));
    }
}
