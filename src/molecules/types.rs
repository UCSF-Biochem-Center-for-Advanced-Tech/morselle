use serde::Serialize;

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq, Ord, PartialOrd)]
pub struct Dna(pub u64);

impl Dna {
    pub fn pack(value: &str) -> Option<Self> {
        let value = value.split('-').next().unwrap_or(value);
        if value.len() > 28 || value.is_empty() {
            return None;
        }
        let mut bits = 0_u64;
        for base in value.bytes() {
            bits = (bits << 2)
                | match base {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    _ => return None,
                };
        }
        Some(Self(bits | ((value.len() as u64) << 56)))
    }

    pub fn len(self) -> usize {
        (self.0 >> 56) as usize
    }

    pub fn is_empty(self) -> bool {
        self.len() == 0
    }

    pub fn hamming(self, other: Self) -> Option<u32> {
        if self.len() != other.len() {
            return None;
        }
        let mut xor = (self.0 ^ other.0) & ((1_u64 << (2 * self.len())) - 1);
        let mut distance = 0;
        while xor != 0 {
            distance += ((xor & 3) != 0) as u32;
            xor >>= 2;
        }
        Some(distance)
    }

    /// All sequences at Hamming distance exactly 1 (single-base substitution).
    pub fn one_mismatch_neighbors(self) -> impl Iterator<Item = Dna> {
        let len = self.len();
        let base = self.0;
        let len_tag = (len as u64) << 56;
        (0..len).flat_map(move |pos| {
            let shift = 2 * (len - 1 - pos);
            let orig = (base >> shift) & 3;
            (0..4u64).filter(move |&b| b != orig).map(move |b| {
                let cleared = base & !(3u64 << shift);
                Dna((cleared | (b << shift)) | len_tag)
            })
        })
    }

    /// This UMI tagged with a read-position class (see
    /// `helpers::position_tag`), for `--position-aware-umis`. The tag lives
    /// in bits 32..56, unused by sequences of <= 16 bases, so it takes part
    /// in equality/hashing (keys at different positions never meet) while
    /// `unpack`, `hamming` and `one_mismatch_neighbors` see only the bases
    /// -- 1-mismatch neighbours therefore stay at the same position.
    /// `None` for sequences longer than 16 bases.
    pub fn with_position_tag(self, tag: u32) -> Option<Self> {
        (self.len() <= 16).then(|| Self(self.0 | ((tag as u64 & POSITION_TAG_MASK) << 32)))
    }

    pub fn unpack(self) -> String {
        let len = self.len();
        let mut out = vec![b'A'; len];
        let mut bits = self.0;
        for base in out.iter_mut().rev() {
            *base = b"ACGT"[(bits & 3) as usize];
            bits >>= 2;
        }
        String::from_utf8(out).unwrap()
    }
}

/// Width of the position tag stored by [`Dna::with_position_tag`] (24 bits).
pub const POSITION_TAG_MASK: u64 = 0x00FF_FFFF;

/// No `PB` tag on this record (non-Flex data, or a probe barcode that
/// wasn't confidently located) -- `sample` is genuinely absent, not merely
/// unknown-but-present, so this must never collide with a real interned id.
pub const NO_SAMPLE: u32 = u32::MAX;

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct RawKey {
    pub cell: Dna,
    /// Interned `PB` (probe/sample barcode) tag, or `NO_SAMPLE` if absent.
    /// Required for Flex data: the same 16bp `CB` whitelist is shared
    /// across every sample in a pooled Flex run (confirmed on real CAT06
    /// data -- 99.8% of observed cell barcodes are shared by more than one
    /// sample), so `cell` alone is not a unique cell identity there. Every
    /// place that groups or looks up by `(cell, ...)` for correction or
    /// collision-resolution purposes must include `sample` too, or reads
    /// from different samples sharing a `CB` get silently merged.
    pub sample: u32,
    pub gene: u32,
    pub umi: Dna,
}

#[derive(Clone, Copy, Debug)]
pub struct Evidence {
    pub reads: u32,
    pub best_ordinal: u64,
    pub best_mapq: u8,
    pub best_aligned: u32,
    // A 32-bit rank has sub-read resolution even for this 1.6B-read run, and
    // keeps Evidence at its pre-saturation 24-byte layout.
    pub min_read_hash: u32,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Outcome {
    Unique,
    Won,
    Lost,
    Tie,
    /// --tagseq with --no-fractional-collisions: this (cell, UMI) was
    /// independently claimed by more than one gene and every claimant was
    /// counted in full, rather than picking a single winner.
    Shared,
    /// --fractional-collisions (default under --tagseq): every claimant is
    /// counted with weight `support / total support`; written as `ZW:f`.
    Fractional,
}

impl Outcome {
    /// Whether this molecule's representative read is counted (`xf:25`).
    pub fn is_counted(self) -> bool {
        matches!(
            self,
            Outcome::Unique | Outcome::Won | Outcome::Shared | Outcome::Fractional
        )
    }

    /// Single-character `XW:A` tag value.
    pub fn xw_code(self) -> u8 {
        match self {
            Outcome::Unique => b'U',
            Outcome::Won => b'W',
            Outcome::Lost => b'L',
            Outcome::Tie => b'T',
            Outcome::Shared => b'S',
            Outcome::Fractional => b'F',
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Decision {
    pub corrected: Dna,
    pub representative: u64,
    pub collision_degree: u32,
    pub outcome: Outcome,
    /// Matrix weight of this molecule: 1.0 for counted non-fractional
    /// outcomes, 0.0 for Lost/Tie, `support / total` for Fractional.
    pub weight: f32,
}

#[derive(Default, Serialize)]
pub struct Metrics {
    #[serde(skip)]
    pub barcode_rows: super::barcode_metrics::BarcodeRows,
    pub input_bam: String,
    pub output_bam: String,
    pub tagseq: bool,
    /// Policy flags this run used (see `sbx10x-molecules --help`).
    pub fractional_collisions: bool,
    pub position_aware_umis: bool,
    /// --min-mapq this run used (0 = no mapping filter).
    pub min_mapq: u8,
    pub records: u64,
    pub eligible_records: u64,
    /// Eligible records below --min-mapq counted because `annotate` tagged
    /// them `mm:i:1` (Cell Ranger's exonic rescue).
    pub rescued_records: u64,
    /// Otherwise-eligible records dropped below --min-mapq: no rescue, and
    /// `NH` absent or 1 (`low_mapq_records`) or `NH` > 1
    /// (`multimapper_records`). Not part of `eligible_records`.
    pub low_mapq_records: u64,
    pub multimapper_records: u64,
    /// Otherwise-eligible (mapped, valid CB+UMI+gene) records dropped for
    /// missing/empty `PB` under --flex. Never counted as a molecule:
    /// `probe-demux` requires the same non-empty `PB`, and letting them share
    /// one `NO_SAMPLE` bucket risks exactly the cross-sample UMI collision
    /// `RawKey::sample` exists to prevent.
    pub no_sample_records: u64,
    pub raw_cell_gene_umis: u64,
    pub corrected_umi_nodes: u64,
    pub molecules_before_gene_resolution: u64,
    pub molecules_counted: u64,
    pub umi_nodes_collapsed: u64,
    pub collision_groups: u64,
    pub collision_molecules_won: u64,
    pub collision_molecules_lost: u64,
    pub collision_ties: u64,
    /// --tagseq --no-fractional-collisions only: collision-group members
    /// counted independently instead of adjudicated.
    pub collision_molecules_shared: u64,
    /// --fractional-collisions only: collision-group members counted with a
    /// fractional weight. Each group adds exactly 1 to `molecules_counted`.
    pub collision_molecules_fractional: u64,
    pub representatives_tagged: u64,
    pub saturation: Option<SaturationMetrics>,
    pub elapsed_secs: f64,
}

/// Adds one shard's local counters into the run-wide accumulator. Used by
/// the concurrent shard path, where each shard builds its own `Metrics`
/// (starting from `Metrics::default()`) rather than sharing one `&mut
/// Metrics` across threads. `input_bam`/`output_bam`/`tagseq`/`saturation`/
/// `elapsed_secs` are set once at the top level, not per-shard, so they're
/// deliberately left alone here.
pub fn merge_metrics(target: &mut Metrics, source: &Metrics) {
    super::barcode_metrics::merge_rows(&mut target.barcode_rows, &source.barcode_rows, |s| s);
    target.records += source.records;
    target.eligible_records += source.eligible_records;
    target.rescued_records += source.rescued_records;
    target.low_mapq_records += source.low_mapq_records;
    target.multimapper_records += source.multimapper_records;
    target.no_sample_records += source.no_sample_records;
    target.raw_cell_gene_umis += source.raw_cell_gene_umis;
    target.corrected_umi_nodes += source.corrected_umi_nodes;
    target.molecules_before_gene_resolution += source.molecules_before_gene_resolution;
    target.molecules_counted += source.molecules_counted;
    target.umi_nodes_collapsed += source.umi_nodes_collapsed;
    target.collision_groups += source.collision_groups;
    target.collision_molecules_won += source.collision_molecules_won;
    target.collision_molecules_lost += source.collision_molecules_lost;
    target.collision_ties += source.collision_ties;
    target.collision_molecules_shared += source.collision_molecules_shared;
    target.collision_molecules_fractional += source.collision_molecules_fractional;
    target.representatives_tagged += source.representatives_tagged;
}

#[derive(Serialize)]
pub struct SaturationMetrics {
    pub output: String,
    pub interval_eligible_reads: u64,
    pub points: usize,
    pub observed_cells: usize,
}
