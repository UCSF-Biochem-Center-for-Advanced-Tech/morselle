//! Additive per-barcode counters shared by pipeline stages.
//!
//! The tables are opt-in (`--barcode-metrics`): a disabled table hands out a
//! scratch row and writes nothing, so a default run pays one branch per update.
//! Rows live in one flat map keyed `barcode\tsample`; a nested map per barcode
//! cost ~4x the memory on Flex data (1.79M barcode/sample rows for the CAT06
//! 5-chunk subset).
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use flate2::{write::GzEncoder, Compression};
use rustc_hash::FxHashMap;

pub const NO_BARCODE: &str = "__no_barcode__";
pub const UNRESOLVED_BARCODE: &str = "__unresolved_barcode__";

/// Keys up to this length are built on the stack for lookups.
const STACK_KEY: usize = 128;

#[derive(Debug, Clone)]
pub struct BarcodeTable<const N: usize> {
    enabled: bool,
    rows: FxHashMap<Box<str>, [u64; N]>,
    scratch: [u64; N],
}

impl<const N: usize> Default for BarcodeTable<N> {
    fn default() -> Self {
        Self::new(false)
    }
}

impl<const N: usize> BarcodeTable<N> {
    pub fn new(enabled: bool) -> Self {
        Self {
            enabled,
            rows: FxHashMap::default(),
            scratch: [0; N],
        }
    }

    pub fn is_enabled(&self) -> bool {
        self.enabled
    }

    pub fn len(&self) -> usize {
        self.rows.len()
    }

    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }

    pub fn row(&mut self, barcode: &str, sample: &str) -> &mut [u64; N] {
        if !self.enabled {
            return &mut self.scratch;
        }
        let len = barcode.len() + 1 + sample.len();
        let mut stack = [0u8; STACK_KEY];
        let owned;
        let key = if len <= STACK_KEY {
            stack[..barcode.len()].copy_from_slice(barcode.as_bytes());
            stack[barcode.len()] = b'\t';
            stack[barcode.len() + 1..len].copy_from_slice(sample.as_bytes());
            // Both halves are &str and the separator is ASCII.
            std::str::from_utf8(&stack[..len]).unwrap()
        } else {
            owned = format!("{barcode}\t{sample}");
            owned.as_str()
        };
        if !self.rows.contains_key(key) {
            self.rows.insert(key.into(), [0; N]);
        }
        self.rows.get_mut(key).unwrap()
    }

    pub fn merge(&mut self, other: Self) {
        if !other.enabled {
            return;
        }
        self.enabled = true;
        for (key, counts) in other.rows {
            let target = self.rows.entry(key).or_insert([0; N]);
            for (dst, src) in target.iter_mut().zip(counts) {
                *dst += src;
            }
        }
    }

    /// Moves the rows out, leaving an empty table with the same setting.
    pub fn take(&mut self) -> Self {
        std::mem::replace(self, Self::new(self.enabled))
    }

    /// Writes the table next to `metrics_path`; a disabled table writes nothing.
    pub fn write(&self, metrics_path: &Path, columns: &[&str; N]) -> Result<Option<PathBuf>> {
        if !self.enabled {
            return Ok(None);
        }
        let path = barcode_path(metrics_path);
        let file = File::create(&path).with_context(|| format!("creating {}", path.display()))?;
        let mut writer = GzEncoder::new(BufWriter::new(file), Compression::default());
        write!(writer, "barcode\tsample")?;
        for column in columns {
            write!(writer, "\t{column}")?;
        }
        writeln!(writer)?;
        // The tab sorts below every barcode character, so ordering whole keys
        // orders by barcode, then sample.
        let mut ordered: Vec<_> = self.rows.iter().collect();
        ordered.sort_unstable_by(|a, b| a.0.cmp(b.0));
        for (key, counts) in ordered {
            write!(writer, "{key}")?;
            for count in counts {
                write!(writer, "\t{count}")?;
            }
            writeln!(writer)?;
        }
        writer.finish()?.flush()?;
        Ok(Some(path))
    }
}

/// A table that many worker threads update in place. SBX `extract` sees
/// ~1.8M barcode/sample rows on Flex data; giving each of 64 workers its own
/// table cost 11.8 GB peak RSS (0.83 GB without), so rows live once, split
/// across locked shards to keep contention low.
pub struct SharedBarcodeTable<const N: usize> {
    shards: Vec<std::sync::Mutex<BarcodeTable<N>>>,
}

impl<const N: usize> SharedBarcodeTable<N> {
    pub fn new(shards: usize) -> Self {
        Self {
            shards: (0..shards.max(1))
                .map(|_| std::sync::Mutex::new(BarcodeTable::new(true)))
                .collect(),
        }
    }

    pub fn update(&self, barcode: &str, sample: &str, apply: impl FnOnce(&mut [u64; N])) {
        use std::hash::Hasher;
        let mut hasher = rustc_hash::FxHasher::default();
        hasher.write(barcode.as_bytes());
        hasher.write(sample.as_bytes());
        let shard = &self.shards[hasher.finish() as usize % self.shards.len()];
        apply(shard.lock().unwrap().row(barcode, sample));
    }

    pub fn into_table(self) -> BarcodeTable<N> {
        let mut table = BarcodeTable::new(true);
        for shard in self.shards {
            table.merge(shard.into_inner().unwrap());
        }
        table
    }
}

pub fn barcode_path(metrics_path: &Path) -> PathBuf {
    let path = metrics_path.to_string_lossy();
    let stem = path.strip_suffix(".json").unwrap_or(&path);
    PathBuf::from(format!("{stem}.barcodes.tsv.gz"))
}
