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

/// No `PB` tag on this record (non-Flex data, or a probe barcode that
/// wasn't confidently located) -- `sample` is genuinely absent, not merely
/// unknown-but-present, so this must never collide with a real interned id.
pub const NO_SAMPLE: u32 = u32::MAX;

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct RawKey {
    pub cell: Dna,
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
    pub min_read_hash: u32,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Outcome {
    Unique,
    Won,
    Lost,
    Tie,
    Shared,
    Fractional,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Decision {
    pub corrected: Dna,
    pub representative: u64,
    pub collision_degree: u32,
    pub outcome: Outcome,
    pub weight: f32,
}

#[derive(Default, Serialize)]
pub struct Metrics {
    pub input_bam: String,
    pub output_bam: String,
    pub tagseq: bool,
    pub records: u64,
    pub eligible_records: u64,
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
    pub collision_molecules_shared: u64,
    pub collision_molecules_fractional: u64,
    pub representatives_tagged: u64,
    pub saturation: Option<SaturationMetrics>,
    pub elapsed_secs: f64,
}

pub fn merge_metrics(target: &mut Metrics, source: &Metrics) {
    target.records += source.records;
    target.eligible_records += source.eligible_records;
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
