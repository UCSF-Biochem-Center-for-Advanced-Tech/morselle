//! FASTQ input discovery/reading and gzip-compressed output writers.

use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
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

/// A thread-safe, append-only gzip writer used by the parallel pipeline (M6).
///
/// Callers pre-format a batch of FASTQ records into a raw byte buffer; this type
/// gzip-compresses that batch into one independent member (off-lock) and appends
/// the compressed bytes under a brief I/O lock. Concatenated gzip members form a
/// valid gzip stream, so many threads can write one sample's file without
/// serializing on compression.
pub struct ShardedGzWriter {
    inner: Mutex<BufWriter<File>>,
    path: PathBuf,
    records: AtomicU64,
}

impl ShardedGzWriter {
    pub fn create(path: &Path) -> Result<Self> {
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent)
                .with_context(|| format!("cannot create dir {}", parent.display()))?;
        }
        let file = File::create(path)
            .with_context(|| format!("cannot create output {}", path.display()))?;
        Ok(ShardedGzWriter {
            inner: Mutex::new(BufWriter::new(file)),
            path: path.to_path_buf(),
            records: AtomicU64::new(0),
        })
    }

    /// Compress `raw` as one gzip member and append it; `n` = records within.
    pub fn append_block(&self, raw: &[u8], n: u64) -> Result<()> {
        if raw.is_empty() {
            return Ok(());
        }
        let bytes = gzip_member(raw)?;
        self.write_compressed(&bytes, n)
    }

    /// Write an already-compressed gzip member under this file's lock.
    fn write_compressed(&self, bytes: &[u8], n: u64) -> Result<()> {
        if bytes.is_empty() {
            return Ok(());
        }
        {
            let mut g = self.inner.lock().expect("writer lock poisoned");
            g.write_all(bytes)?;
        }
        self.records.fetch_add(n, Ordering::Relaxed);
        Ok(())
    }

    pub fn records(&self) -> u64 {
        self.records.load(Ordering::Relaxed)
    }

    pub fn finish(&self) -> Result<()> {
        self.inner.lock().expect("writer lock poisoned").flush()?;
        Ok(())
    }

    /// Rewrite this file as a SINGLE gzip member (the multi-member output is valid
    /// gzip, but some readers — e.g. Cell Ranger — mis-frame records at member
    /// boundaries). Decompress all members, recompress as one member into a temp
    /// file, fsync, then atomically rename over the original. Errors propagate, so
    /// a partial/failed recompress never replaces the good file. Call after
    /// `finish()`.
    pub fn finalize_single_member(&self, level: u32) -> Result<()> {
        let path = &self.path;
        let tmp = path.with_extension("gz.tmp");
        {
            let src = File::open(path)
                .with_context(|| format!("finalize: cannot open {}", path.display()))?;
            let mut dec = MultiGzDecoder::new(BufReader::new(src));
            let out = File::create(&tmp)
                .with_context(|| format!("finalize: cannot create {}", tmp.display()))?;
            let mut enc = GzEncoder::new(BufWriter::new(out), Compression::new(level));
            std::io::copy(&mut dec, &mut enc)
                .with_context(|| format!("finalize: recompress failed for {}", path.display()))?;
            let mut bw = enc.finish().context("finalize: gzip finish failed")?;
            bw.flush().context("finalize: flush failed")?;
            bw.get_ref().sync_all().context("finalize: fsync failed")?;
        }
        fs::rename(&tmp, path)
            .with_context(|| format!("finalize: rename {} -> {}", tmp.display(), path.display()))?;
        Ok(())
    }
}

/// A paired R1/R2 sharded writer for one output target (sample/Undetermined).
///
/// R1 and R2 must stay mate-ordered (the Nth record in R1 is the mate of the Nth
/// in R2) or Cell Ranger rejects the input. `gate` serializes each paired append
/// so the i-th block lands at the same position in both files — without it, two
/// threads can interleave R1 and R2 writes independently and desync the files.
pub struct PairShardWriter {
    pub r1: ShardedGzWriter,
    pub r2: ShardedGzWriter,
    gate: Mutex<()>,
}

impl PairShardWriter {
    pub fn create(dir: &Path, prefix: &str, s_num: usize) -> Result<Self> {
        let name = |read: u8| dir.join(format!("{prefix}_S{s_num}_L001_R{read}_001.fastq.gz"));
        Ok(PairShardWriter {
            r1: ShardedGzWriter::create(&name(1))?,
            r2: ShardedGzWriter::create(&name(2))?,
            gate: Mutex::new(()),
        })
    }

    /// Append a pre-formatted batch (`n` records) of R1 and R2 FASTQ bytes.
    /// Compression happens off-lock; the two writes happen under one lock so R1
    /// and R2 blocks are appended atomically and in the same order.
    pub fn append(&self, r1_raw: &[u8], r2_raw: &[u8], n: u64) -> Result<()> {
        if r1_raw.is_empty() && r2_raw.is_empty() {
            return Ok(());
        }
        let b1 = gzip_member(r1_raw)?;
        let b2 = gzip_member(r2_raw)?;
        let _g = self.gate.lock().expect("pair writer lock poisoned");
        self.r1.write_compressed(&b1, n)?;
        self.r2.write_compressed(&b2, n)?;
        Ok(())
    }

    pub fn finish(&self) -> Result<()> {
        self.r1.finish()?;
        self.r2.finish()
    }

    /// Rewrite both files as single gzip members (Cell Ranger compatibility).
    pub fn finalize_single_member(&self, level: u32) -> Result<()> {
        self.r1.finalize_single_member(level)?;
        self.r2.finalize_single_member(level)
    }
}

/// Compress `raw` into a single standalone gzip member (empty in → empty out).
fn gzip_member(raw: &[u8]) -> Result<Vec<u8>> {
    if raw.is_empty() {
        return Ok(Vec::new());
    }
    let mut enc = GzEncoder::new(Vec::with_capacity(raw.len() / 3 + 64), Compression::default());
    enc.write_all(raw)?;
    Ok(enc.finish()?)
}
