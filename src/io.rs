//! FASTQ input discovery/reading and gzip-compressed output writers.

use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;

/// Expand the `--fastqs` arguments into a concrete, sorted list of FASTQ files.
///
/// Each argument may be a single file or a directory; directories are scanned
/// (non-recursively) for `*.fastq.gz` / `*.fq.gz` / `*.fastq` / `*.fq`.
pub fn discover_fastqs(inputs: &[PathBuf]) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    for input in inputs {
        let meta = fs::metadata(input)
            .with_context(|| format!("cannot stat input {}", input.display()))?;
        if meta.is_dir() {
            for entry in fs::read_dir(input)
                .with_context(|| format!("cannot read dir {}", input.display()))?
            {
                let path = entry?.path();
                if is_fastq(&path) {
                    files.push(path);
                }
            }
        } else if is_fastq(input) {
            files.push(input.clone());
        } else {
            anyhow::bail!("{} is not a recognized FASTQ file", input.display());
        }
    }
    files.sort();
    files.dedup();
    Ok(files)
}

fn is_fastq(path: &Path) -> bool {
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
    name.ends_with(".fastq.gz")
        || name.ends_with(".fq.gz")
        || name.ends_with(".fastq")
        || name.ends_with(".fq")
}

/// A buffered, gzip-compressed FASTQ writer.
pub struct FastqGzWriter {
    inner: GzEncoder<BufWriter<File>>,
    pub records: u64,
}

impl FastqGzWriter {
    pub fn create(path: &Path) -> Result<Self> {
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent)
                .with_context(|| format!("cannot create dir {}", parent.display()))?;
        }
        let file = File::create(path)
            .with_context(|| format!("cannot create output {}", path.display()))?;
        let inner = GzEncoder::new(BufWriter::new(file), Compression::default());
        Ok(Self { inner, records: 0 })
    }

    /// Write one FASTQ record. `id` excludes the leading '@'.
    pub fn write_record(&mut self, id: &[u8], seq: &[u8], qual: &[u8]) -> Result<()> {
        self.inner.write_all(b"@")?;
        self.inner.write_all(id)?;
        self.inner.write_all(b"\n")?;
        self.inner.write_all(seq)?;
        self.inner.write_all(b"\n+\n")?;
        self.inner.write_all(qual)?;
        self.inner.write_all(b"\n")?;
        self.records += 1;
        Ok(())
    }

    pub fn finish(self) -> Result<()> {
        self.inner.finish()?.flush()?;
        Ok(())
    }
}

/// A paired R1/R2 gzip writer for one sample (Cell Ranger naming:
/// `<prefix>_S<n>_L001_R{1,2}_001.fastq.gz`).
pub struct PairWriter {
    pub r1: FastqGzWriter,
    pub r2: FastqGzWriter,
}

impl PairWriter {
    pub fn create(dir: &Path, prefix: &str, s_num: usize) -> Result<Self> {
        let name = |read: u8| dir.join(format!("{prefix}_S{s_num}_L001_R{read}_001.fastq.gz"));
        Ok(PairWriter {
            r1: FastqGzWriter::create(&name(1))?,
            r2: FastqGzWriter::create(&name(2))?,
        })
    }

    pub fn write(
        &mut self,
        id: &[u8],
        r1_seq: &[u8],
        r1_qual: &[u8],
        r2_seq: &[u8],
        r2_qual: &[u8],
    ) -> Result<()> {
        self.r1.write_record(id, r1_seq, r1_qual)?;
        self.r2.write_record(id, r2_seq, r2_qual)?;
        Ok(())
    }

    pub fn finish(self) -> Result<()> {
        self.r1.finish()?;
        self.r2.finish()?;
        Ok(())
    }
}
