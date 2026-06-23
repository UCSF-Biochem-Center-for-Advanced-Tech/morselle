//! Sample-sheet parsing into a normalized [`SampleSheet`] (plan §5).
//!
//! Two front-end formats feed one internal model:
//! - **Primary:** BCL Convert v2 sheet (the artifact used to demux the matching
//!   NovaSeq X run). We read `[Reads]` (cycle counts → CB/UMI/index lengths),
//!   `[BCLConvert_Settings]` (`Index2ColumnReverseComplement`, barcode
//!   mismatches), `[Header]` (`IndexOrientation`), and `[BCLConvert_Data]`
//!   (`Sample_ID,index(=i7),index2(=i5)`, collapsing duplicate per-lane rows).
//! - **Fallback:** flat CSV with header `sample,i7,i5`.
//!
//! i5 note: BCL Convert writes `index2` in the orientation implied by
//! `IndexOrientation`/`Index2ColumnReverseComplement`. We store it verbatim plus
//! the rev-comp flag and expose both orientations; M3 confirms empirically which
//! one the SBX-extracted i5 matches.

use std::collections::HashMap;
use std::fmt::Write as _;
use std::fs;
use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::anchors::revcomp;

/// Read-cycle layout from the `[Reads]` section.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Chemistry {
    /// Read 1 length (CB + UMI), e.g. 28.
    pub read1_cycles: usize,
    /// i7 length, e.g. 10.
    pub index1_cycles: usize,
    /// i5 length, e.g. 10.
    pub index2_cycles: usize,
    /// Read 2 (cDNA) length, e.g. 90.
    pub read2_cycles: usize,
}

impl Chemistry {
    /// UMI length given a cell-barcode length: `read1 - cb`.
    pub fn umi_len(&self, cb_len: usize) -> Result<usize> {
        self.read1_cycles
            .checked_sub(cb_len)
            .filter(|&u| u > 0)
            .with_context(|| {
                format!(
                    "Read1Cycles ({}) <= cb_len ({cb_len}); cannot derive UMI length",
                    self.read1_cycles
                )
            })
    }
}

/// One demultiplexing target.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SampleEntry {
    pub id: String,
    /// i7 (BCL Convert `index`), as written in the sheet.
    pub i7: Vec<u8>,
    /// i5 (BCL Convert `index2`), as written in the sheet.
    pub i5: Vec<u8>,
}

impl SampleEntry {
    /// i5 reverse-complemented (the other candidate orientation for matching).
    pub fn i5_revcomp(&self) -> Vec<u8> {
        revcomp(&self.i5)
    }
}

/// Which front-end format a sheet was parsed from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SourceFormat {
    BclConvert,
    Flat,
}

/// Normalized sample sheet shared by both parsers.
#[derive(Debug, Clone)]
pub struct SampleSheet {
    pub samples: Vec<SampleEntry>,
    /// The format this sheet was parsed from.
    pub source: SourceFormat,
    /// Present for BCL Convert sheets; `None` for the flat format.
    pub chemistry: Option<Chemistry>,
    /// `Index2ColumnReverseComplement` (true for NovaSeq X forward workflow).
    pub i5_reverse_complemented: bool,
    /// `IndexOrientation,Forward` → true.
    pub index_orientation_forward: bool,
    /// `BarcodeMismatchesIndex1`.
    pub i7_max_mm: Option<u32>,
    /// `BarcodeMismatchesIndex2`.
    pub i5_max_mm: Option<u32>,
}

impl SampleSheet {
    /// Load from a path, auto-detecting format by first non-empty line.
    pub fn load(path: &Path) -> Result<Self> {
        let text = fs::read_to_string(path)
            .with_context(|| format!("cannot read sample sheet {}", path.display()))?;
        let first = text
            .lines()
            .map(str::trim)
            .find(|l| !l.is_empty())
            .unwrap_or("");
        let sheet = if first.starts_with('[') {
            parse_bclconvert(&text)
        } else {
            parse_flat(&text)
        }
        .with_context(|| format!("parsing sample sheet {}", path.display()))?;
        sheet.validate()?;
        Ok(sheet)
    }

    /// Render the normalized sheet as flat `sample,i7,i5` CSV text.
    ///
    /// i7/i5 are emitted verbatim (as in the source sheet), so the output
    /// round-trips back through [`parse_flat`].
    pub fn to_flat_csv(&self) -> String {
        let mut out = String::from("sample,i7,i5\n");
        for s in &self.samples {
            let _ = writeln!(
                out,
                "{},{},{}",
                s.id,
                String::from_utf8_lossy(&s.i7),
                String::from_utf8_lossy(&s.i5),
            );
        }
        out
    }

    /// Write the flat `sample,i7,i5` CSV to `path` (creating parent dirs).
    pub fn write_flat_csv(&self, path: &Path) -> Result<()> {
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent)
                .with_context(|| format!("cannot create dir {}", parent.display()))?;
        }
        fs::write(path, self.to_flat_csv())
            .with_context(|| format!("cannot write flat sample sheet {}", path.display()))
    }

    /// Structural checks shared by both formats.
    fn validate(&self) -> Result<()> {
        if self.samples.is_empty() {
            bail!("sample sheet has no samples");
        }
        let mut seen_ids = HashMap::new();
        let mut seen_idx: HashMap<(Vec<u8>, Vec<u8>), String> = HashMap::new();
        for s in &self.samples {
            if let Some(prev) = seen_ids.insert(s.id.clone(), ()) {
                let _ = prev;
                bail!("duplicate sample id '{}'", s.id);
            }
            for (label, bc) in [("i7", &s.i7), ("i5", &s.i5)] {
                if bc.is_empty() || !bc.iter().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
                    bail!("sample '{}' has invalid {label} index '{}'", s.id, String::from_utf8_lossy(bc));
                }
            }
            if let Some(other) =
                seen_idx.insert((s.i7.clone(), s.i5.clone()), s.id.clone())
            {
                bail!(
                    "samples '{}' and '{}' share the same i7/i5 pair",
                    other, s.id
                );
            }
        }
        // If chemistry is known, index lengths must match the declared cycles.
        if let Some(chem) = self.chemistry {
            for s in &self.samples {
                if s.i7.len() != chem.index1_cycles {
                    bail!(
                        "sample '{}' i7 len {} != Index1Cycles {}",
                        s.id, s.i7.len(), chem.index1_cycles
                    );
                }
                if s.i5.len() != chem.index2_cycles {
                    bail!(
                        "sample '{}' i5 len {} != Index2Cycles {}",
                        s.id, s.i5.len(), chem.index2_cycles
                    );
                }
            }
        }
        Ok(())
    }
}

/// Parse the BCL Convert v2 multi-section format.
fn parse_bclconvert(text: &str) -> Result<SampleSheet> {
    let mut section = String::new();
    let mut reads: HashMap<String, usize> = HashMap::new();
    let mut settings: HashMap<String, String> = HashMap::new();
    let mut header: HashMap<String, String> = HashMap::new();
    let mut data_header: Vec<String> = Vec::new();
    // Preserve first-seen order while collapsing duplicate per-lane rows.
    let mut order: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    let mut by_key: HashMap<(Vec<u8>, Vec<u8>), SampleEntry> = HashMap::new();

    for raw in text.lines() {
        let line = raw.trim_end_matches('\r').trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('[') && line.ends_with(']') {
            section = line[1..line.len() - 1].to_ascii_lowercase();
            continue;
        }
        let fields: Vec<&str> = line.split(',').map(str::trim).collect();
        match section.as_str() {
            "header" => {
                if fields.len() >= 2 {
                    header.insert(fields[0].to_ascii_lowercase(), fields[1].to_string());
                }
            }
            "reads" => {
                if fields.len() >= 2 {
                    if let Ok(n) = fields[1].parse::<usize>() {
                        reads.insert(fields[0].to_ascii_lowercase(), n);
                    }
                }
            }
            "bclconvert_settings" => {
                if fields.len() >= 2 {
                    settings.insert(fields[0].to_ascii_lowercase(), fields[1].to_string());
                }
            }
            "bclconvert_data" => {
                if data_header.is_empty() {
                    data_header = fields.iter().map(|s| s.to_ascii_lowercase()).collect();
                    continue;
                }
                let row: HashMap<&str, &str> = data_header
                    .iter()
                    .map(String::as_str)
                    .zip(fields.iter().copied())
                    .collect();
                let id = row
                    .get("sample_id")
                    .copied()
                    .context("[BCLConvert_Data] row missing Sample_ID")?
                    .to_string();
                let i7 = row
                    .get("index")
                    .copied()
                    .context("[BCLConvert_Data] row missing index")?
                    .to_ascii_uppercase()
                    .into_bytes();
                let i5 = row
                    .get("index2")
                    .copied()
                    .context("[BCLConvert_Data] row missing index2")?
                    .to_ascii_uppercase()
                    .into_bytes();
                let key = (i7.clone(), i5.clone());
                match by_key.get(&key) {
                    Some(existing) if existing.id != id => bail!(
                        "i7/i5 pair shared by samples '{}' and '{}'",
                        existing.id, id
                    ),
                    Some(_) => {} // exact duplicate lane row — collapse
                    None => {
                        order.push(key.clone());
                        by_key.insert(key, SampleEntry { id, i7, i5 });
                    }
                }
            }
            _ => {} // ignore unknown sections (e.g. settings we don't use)
        }
    }

    let chemistry = match (
        reads.get("read1cycles"),
        reads.get("index1cycles"),
        reads.get("index2cycles"),
        reads.get("read2cycles"),
    ) {
        (Some(&r1), Some(&i1), Some(&i2), Some(&r2)) => Some(Chemistry {
            read1_cycles: r1,
            index1_cycles: i1,
            index2_cycles: i2,
            read2_cycles: r2,
        }),
        _ => None,
    };

    let i5_rc = settings
        .get("index2columnreversecomplement")
        .map(|v| is_truthy(v))
        .unwrap_or(false);
    let index_fwd = header
        .get("indexorientation")
        .map(|v| v.eq_ignore_ascii_case("forward"))
        .unwrap_or(true);

    let samples = order
        .into_iter()
        .map(|k| by_key.remove(&k).expect("key present"))
        .collect();

    Ok(SampleSheet {
        samples,
        source: SourceFormat::BclConvert,
        chemistry,
        i5_reverse_complemented: i5_rc,
        index_orientation_forward: index_fwd,
        i7_max_mm: settings
            .get("barcodemismatchesindex1")
            .and_then(|v| v.parse().ok()),
        i5_max_mm: settings
            .get("barcodemismatchesindex2")
            .and_then(|v| v.parse().ok()),
    })
}

/// Parse the flat fallback CSV: header `sample,i7,i5`.
fn parse_flat(text: &str) -> Result<SampleSheet> {
    let mut lines = text
        .lines()
        .map(|l| l.trim_end_matches('\r').trim())
        .filter(|l| !l.is_empty());
    let header = lines.next().context("empty flat sample sheet")?;
    let cols: Vec<String> = header.split(',').map(|s| s.trim().to_ascii_lowercase()).collect();
    let col = |name: &str| cols.iter().position(|c| c == name);
    let (si, i7i, i5i) = (
        col("sample").context("flat sheet missing 'sample' column")?,
        col("i7").context("flat sheet missing 'i7' column")?,
        col("i5").context("flat sheet missing 'i5' column")?,
    );

    let mut samples = Vec::new();
    for line in lines {
        let f: Vec<&str> = line.split(',').map(str::trim).collect();
        let max = si.max(i7i).max(i5i);
        if f.len() <= max {
            bail!("flat sheet row has too few columns: '{line}'");
        }
        samples.push(SampleEntry {
            id: f[si].to_string(),
            i7: f[i7i].to_ascii_uppercase().into_bytes(),
            i5: f[i5i].to_ascii_uppercase().into_bytes(),
        });
    }

    Ok(SampleSheet {
        samples,
        source: SourceFormat::Flat,
        chemistry: None,
        i5_reverse_complemented: false,
        index_orientation_forward: true,
        i7_max_mm: None,
        i5_max_mm: None,
    })
}

fn is_truthy(v: &str) -> bool {
    matches!(
        v.trim().to_ascii_lowercase().as_str(),
        "y" | "yes" | "true" | "1"
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn fixture(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures")
            .join(name)
    }

    #[test]
    fn parses_bclconvert_example() {
        let sheet = SampleSheet::load(&fixture("CE15204_bclconvert_input.csv")).unwrap();

        // 14 samples, lane duplicates collapsed.
        assert_eq!(sheet.samples.len(), 14);

        // Chemistry from [Reads].
        let chem = sheet.chemistry.unwrap();
        assert_eq!(chem.read1_cycles, 28);
        assert_eq!(chem.index1_cycles, 10);
        assert_eq!(chem.index2_cycles, 10);
        assert_eq!(chem.read2_cycles, 90);
        assert_eq!(chem.umi_len(16).unwrap(), 12);

        // Settings.
        assert!(sheet.i5_reverse_complemented);
        assert!(sheet.index_orientation_forward);
        assert_eq!(sheet.i7_max_mm, Some(1));
        assert_eq!(sheet.i5_max_mm, Some(1));

        // First and a known sample preserve sheet order + exact indices.
        assert_eq!(sheet.samples[0].id, "V2_GEX");
        assert_eq!(sheet.samples[0].i7, b"TTATTCGAGG");
        assert_eq!(sheet.samples[0].i5, b"AGCAGGACAG");

        let m2 = sheet.samples.iter().find(|s| s.id == "M2_GEX").unwrap();
        assert_eq!(m2.i7, b"TGTAGTCATT");
        assert_eq!(m2.i5, b"TACGATCAAG");
        // i5 revcomp helper (other candidate matching orientation).
        assert_eq!(m2.i5_revcomp(), b"CTTGATCGTA");

        assert_eq!(sheet.samples.last().unwrap().id, "M7_GEX");
    }

    #[test]
    fn bclconvert_to_flat_csv_roundtrips() {
        let sheet = SampleSheet::load(&fixture("CE15204_bclconvert_input.csv")).unwrap();
        assert_eq!(sheet.source, SourceFormat::BclConvert);

        let csv = sheet.to_flat_csv();
        assert!(csv.starts_with("sample,i7,i5\n"));
        assert!(csv.contains("M2_GEX,TGTAGTCATT,TACGATCAAG\n"));

        // Re-parsing the emitted flat CSV yields the same samples.
        let reparsed = parse_flat(&csv).unwrap();
        reparsed.validate().unwrap();
        assert_eq!(reparsed.source, SourceFormat::Flat);
        assert_eq!(reparsed.samples, sheet.samples);
    }

    #[test]
    fn parses_flat_fallback() {
        let sheet = SampleSheet::load(&fixture("flat_samplesheet.csv")).unwrap();
        assert_eq!(sheet.samples.len(), 3);
        assert!(sheet.chemistry.is_none());
        assert!(!sheet.i5_reverse_complemented);
        assert_eq!(sheet.samples[0].id, "V2_GEX");
        assert_eq!(sheet.samples[1].i7, b"TGTAGTCATT");
    }

    #[test]
    fn rejects_duplicate_index_pair() {
        let text = "sample,i7,i5\nA,ACGTACGTAC,TTTTTTTTTT\nB,ACGTACGTAC,TTTTTTTTTT\n";
        let err = parse_flat(text).unwrap().validate().unwrap_err();
        assert!(err.to_string().contains("share the same i7/i5"));
    }

    #[test]
    fn rejects_non_acgt_index() {
        let text = "sample,i7,i5\nA,ACGTACGTAN,TTTTTTTTTT\n";
        let err = parse_flat(text).unwrap().validate().unwrap_err();
        assert!(err.to_string().contains("invalid i7"));
    }
}
