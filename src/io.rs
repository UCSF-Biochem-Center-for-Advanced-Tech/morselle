//! FASTQ input discovery/reading and gzip-compressed output writers.

use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read as IoRead, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

/// A shard output stream, gzip-compressed or plain. `None` for `level`
/// (rather than a separate bool) makes "uncompressed" a real absence of a
/// compression level, not a footgun default-level-0 gzip stream.
///
/// Exists because gzip encoding is real, measured CPU cost (~54% of
/// `extract`'s time in a profiled run — see docs/PROFILING.md) that's pure waste
/// when the output is about to be read once by a co-located process (a
/// mapper) rather than stored: `morselle extract --no-gzip` into a named
/// pipe lets `minimap2` read raw FASTQ directly with no compress/decompress
/// round trip on either end.
enum ShardSink {
    Gz(GzEncoder<BufWriter<File>>),
    Plain(BufWriter<File>),
}

impl ShardSink {
    fn create(path: &Path, level: Option<u32>) -> Result<Self> {
        let file = create_output(path)?;
        Ok(match level {
            Some(level) => ShardSink::Gz(GzEncoder::new(
                BufWriter::new(file),
                Compression::new(level),
            )),
            None => ShardSink::Plain(BufWriter::new(file)),
        })
    }

    fn write_all(&mut self, buf: &[u8]) -> io::Result<()> {
        match self {
            ShardSink::Gz(w) => w.write_all(buf),
            ShardSink::Plain(w) => w.write_all(buf),
        }
    }

    /// Finish the stream (gzip trailer, if any) and flush.
    fn finish(&mut self) -> Result<()> {
        match self {
            ShardSink::Gz(w) => {
                w.try_finish().context("gzip finish")?;
                w.get_mut().flush().context("flush")?;
            }
            ShardSink::Plain(w) => w.flush().context("flush")?,
        }
        Ok(())
    }
}

fn create_output(path: &Path) -> Result<File> {
    File::create(path).with_context(|| format!("cannot create output {}", path.display()))
}

/// FASTQ filename extension for a given compression level (`.fastq.gz` vs
/// `.fastq`) — content must match extension, or a downstream reader that
/// infers compression from the name (most of them) breaks.
fn fastq_ext(level: Option<u32>) -> &'static str {
    if level.is_some() {
        "fastq.gz"
    } else {
        "fastq"
    }
}

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

/// A paired R1/R2 gzip writer for one sample (10x reference naming:
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

/// A thread-safe append-only gzip FASTQ writer for a single stream.
///
/// Used for the `--failed-output` side channel, where every worker appends into
/// one file. Each worker compresses its own buffered batch into a complete gzip
/// member before taking the lock; only the write of the compressed bytes is
/// serialised. The file is a multi-member gzip, which gzip/zcat, Python's
/// `gzip` and `MultiGzDecoder` all read as one stream.
///
/// Compressing under the lock instead (one member, as `PairShardWriter` does)
/// made the whole extract single-threaded once failures were common: real
/// CAT02 3p full run, 53% of reads failed at ~630 bytes each, 4 of 64 threads
/// busy, ~65k reads/s against 1.3M reads/s without `--failed-output`.
pub struct ShardWriter {
    level: Option<u32>,
    inner: Mutex<FailedSink>,
}

struct FailedSink {
    out: BufWriter<File>,
    wrote: bool,
}

impl ShardWriter {
    pub fn create(path: &Path, level: Option<u32>) -> Result<Self> {
        if let Some(dir) = path.parent() {
            if !dir.as_os_str().is_empty() {
                fs::create_dir_all(dir)
                    .with_context(|| format!("cannot create dir {}", dir.display()))?;
            }
        }
        let file = create_output(path)?;
        Ok(ShardWriter {
            level,
            inner: Mutex::new(FailedSink {
                out: BufWriter::new(file),
                wrote: false,
            }),
        })
    }

    /// Append a pre-formatted batch of FASTQ bytes.
    pub fn append(&self, raw: &[u8]) -> Result<()> {
        if raw.is_empty() {
            return Ok(());
        }
        let member;
        let bytes = match self.level {
            Some(level) => {
                member = gzip_member(raw, level)?;
                &member[..]
            }
            None => raw,
        };
        let mut w = self.inner.lock().expect("shard writer lock poisoned");
        w.out.write_all(bytes)?;
        w.wrote = true;
        Ok(())
    }

    /// Flush. A gzip stream that got no records still gets one empty member,
    /// so the file is valid gzip rather than zero bytes.
    pub fn finish(&self) -> Result<()> {
        let mut w = self.inner.lock().expect("shard writer lock poisoned");
        if let (Some(level), false) = (self.level, w.wrote) {
            let empty = gzip_member(&[], level)?;
            w.out.write_all(&empty)?;
            w.wrote = true;
        }
        w.out.flush().context("flush")
    }
}

/// One complete gzip member holding `raw`.
fn gzip_member(raw: &[u8], level: u32) -> Result<Vec<u8>> {
    let mut enc = GzEncoder::new(Vec::with_capacity(raw.len() / 3), Compression::new(level));
    enc.write_all(raw)?;
    enc.finish().context("gzip finish")
}

/// A thread-safe, append-only paired writer used by the parallel pipeline.
///
/// Each output file is a SINGLE gzip member: callers pre-format a batch of FASTQ
/// records into raw bytes and `append` writes them straight into a per-file
/// `GzEncoder` under one lock per target. Both R1 and R2 are written under the
/// same lock so the Nth record in R1 stays the mate of the Nth in R2 (10x reference
/// rejects desynced pairs), and one member per file avoids the multi-member
/// boundary mis-framing some readers (10x reference) exhibit. Compression happens
/// under the lock; with a fast zlib backend it keeps up with the per-target data
/// rate, so the lock is rarely contended.
pub struct PairShardWriter {
    inner: Mutex<Pair>,
}

struct Pair {
    r1: Option<ShardSink>,
    r2: ShardSink,
}

impl PairShardWriter {
    /// `chunk` is the 1-based file-set number in the 10x reference name (`_001`…).
    /// `level: None` writes plain, uncompressed `.fastq` (see [`ShardSink`]).
    pub fn create(
        dir: &Path,
        prefix: &str,
        s_num: usize,
        chunk: usize,
        level: Option<u32>,
        emit_r1: bool,
    ) -> Result<Self> {
        fs::create_dir_all(dir).with_context(|| format!("cannot create dir {}", dir.display()))?;
        let ext = fastq_ext(level);
        let enc = |read: u8| -> Result<ShardSink> {
            let path = dir.join(format!("{prefix}_S{s_num}_L001_R{read}_{chunk:03}.{ext}"));
            ShardSink::create(&path, level)
        };
        Ok(PairShardWriter {
            inner: Mutex::new(Pair {
                r1: emit_r1.then(|| enc(1)).transpose()?,
                r2: enc(2)?,
            }),
        })
    }

    /// Append a pre-formatted batch of R1 and R2 FASTQ bytes. Both go in under one
    /// lock so the files stay mate-ordered. `_n` (record count) is unused.
    pub fn append(&self, r1_raw: &[u8], r2_raw: &[u8], _n: u64) -> Result<()> {
        if r1_raw.is_empty() && r2_raw.is_empty() {
            return Ok(());
        }
        let mut p = self.inner.lock().expect("pair writer lock poisoned");
        if let Some(r1) = &mut p.r1 {
            r1.write_all(r1_raw)?;
        }
        p.r2.write_all(r2_raw)?;
        Ok(())
    }

    /// Finish both streams (gzip trailer, if any) and flush.
    pub fn finish(&self) -> Result<()> {
        let mut p = self.inner.lock().expect("pair writer lock poisoned");
        if let Some(r1) = &mut p.r1 {
            r1.finish()?;
        }
        p.r2.finish()
    }
}

/// Open a text file for line-by-line reading, transparently decoding gzip when the extension is `.gz`.
pub fn open_text_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let input: Box<dyn IoRead> = if path.extension().is_some_and(|ext| ext == "gz") {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(Box::new(BufReader::new(input)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::read::MultiGzDecoder;

    #[test]
    fn no_gzip_writes_plain_readable_fastq_with_plain_extension() {
        let dir = tempfile::tempdir().unwrap();
        let w = PairShardWriter::create(dir.path(), "S", 1, 1, None, true).unwrap();
        w.append(b"@r1\nACGT\n+\nIIII\n", b"@r2\nTTTT\n+\nIIII\n", 1)
            .unwrap();
        w.finish().unwrap();

        let r1_path = dir.path().join("S_S1_L001_R1_001.fastq");
        let r2_path = dir.path().join("S_S1_L001_R2_001.fastq");
        assert!(
            r1_path.exists(),
            "expected a plain .fastq file, not .fastq.gz"
        );
        assert_eq!(fs::read(&r1_path).unwrap(), b"@r1\nACGT\n+\nIIII\n");
        assert_eq!(fs::read(&r2_path).unwrap(), b"@r2\nTTTT\n+\nIIII\n");
    }

    #[test]
    fn gzip_level_writes_valid_gzip_with_gz_extension() {
        let dir = tempfile::tempdir().unwrap();
        let w = PairShardWriter::create(dir.path(), "S", 1, 1, Some(4), true).unwrap();
        w.append(b"@r1\nACGT\n+\nIIII\n", b"@r2\nTTTT\n+\nIIII\n", 1)
            .unwrap();
        w.finish().unwrap();

        let r2_path = dir.path().join("S_S1_L001_R2_001.fastq.gz");
        assert!(r2_path.exists(), "expected a .fastq.gz file");
        let mut decoded = String::new();
        MultiGzDecoder::new(File::open(&r2_path).unwrap())
            .read_to_string(&mut decoded)
            .unwrap();
        assert_eq!(decoded, "@r2\nTTTT\n+\nIIII\n");
    }

    #[test]
    fn shard_writer_no_gzip_round_trips() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("failed.fastq");
        let w = ShardWriter::create(&path, None).unwrap();
        w.append(b"@r\nACGT\n+\nIIII\n").unwrap();
        w.finish().unwrap();
        assert_eq!(fs::read(&path).unwrap(), b"@r\nACGT\n+\nIIII\n");
    }

    #[test]
    fn shard_writer_gzip_keeps_every_batch_from_many_threads() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("failed.fastq.gz");
        let w = ShardWriter::create(&path, Some(6)).unwrap();
        std::thread::scope(|s| {
            for t in 0..8 {
                let w = &w;
                s.spawn(move || {
                    for i in 0..50 {
                        w.append(format!("@t{t}_{i}\nACGT\n+\nIIII\n").as_bytes())
                            .unwrap();
                    }
                });
            }
        });
        w.finish().unwrap();
        let mut decoded = String::new();
        MultiGzDecoder::new(File::open(&path).unwrap())
            .read_to_string(&mut decoded)
            .unwrap();
        let mut ids: Vec<&str> = decoded.lines().step_by(4).collect();
        ids.sort_unstable();
        let mut want: Vec<String> = (0..8)
            .flat_map(|t| (0..50).map(move |i| format!("@t{t}_{i}")))
            .collect();
        want.sort_unstable();
        assert_eq!(ids, want);
        assert_eq!(decoded.lines().count(), 8 * 50 * 4);
    }

    #[test]
    fn shard_writer_gzip_with_no_records_is_valid_empty_gzip() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("failed.fastq.gz");
        let w = ShardWriter::create(&path, Some(6)).unwrap();
        w.finish().unwrap();
        let mut decoded = String::new();
        MultiGzDecoder::new(File::open(&path).unwrap())
            .read_to_string(&mut decoded)
            .unwrap();
        assert!(decoded.is_empty());
        assert!(fs::metadata(&path).unwrap().len() > 0);
    }
}
