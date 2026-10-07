use std::path::PathBuf;
use clap::Parser;
use serde::Serialize;

pub const BIN: usize = 25;
pub const N_BINS: usize = 29;

pub const MIN_CLIP: usize = 5;
pub const K: usize = 8;
pub const HOMO: f64 = 0.8;

pub const ADAPTER: &[u8] = b"AGATCGGAAGAGC";
pub const TSO: &[u8] = b"AAGCAGTGGTATCAACGCAGAGTACATGGG";

#[derive(Parser, Debug)]
#[command(
    name = "bam-audit",
    version,
    about = "Single-pass audit of a 10x reference BAM"
)]
pub struct Cli {
    /// Input 10x reference BAM (possorted_genome_bam.bam).
    pub bam: PathBuf,

    /// Label for this BAM in the JSON (e.g. "ILLUMINA", "SBX p5").
    #[arg(long, default_value = "bam")]
    pub label: String,

    /// Output JSON summary.
    #[arg(short, long, default_value = "bam_audit.json")]
    pub output: PathBuf,

    /// BAM decompression threads.
    #[arg(long, default_value_t = 4)]
    pub threads: usize,

    /// Optional gzip TSV of hash-subsampled molecule representatives.
    #[arg(long)]
    pub molecules: Option<PathBuf>,

    /// Emit 1 in N molecules to `--molecules` (same N across BAMs!). 1 = all.
    #[arg(long, default_value_t = 64)]
    pub molecule_sample: u64,

    /// Optional called-cell barcode list.
    #[arg(long)]
    pub cells: Option<PathBuf>,
}

#[derive(Default)]
pub struct MolAccum {
    pub n_reads: u32,
    pub min_nh: u32,
    pub max_mapq: u8,
    pub region: u8,
    pub max_readlen: u32,
    pub counted: bool,
    pub eligible: bool,
    pub rep_gene: String,
}

#[derive(Default, Serialize)]
pub struct MappingDiagnostics {
    pub nh_histogram: [u64; 5],
    pub mapq_histogram: [u64; 4],
    pub nh_by_gene_degree: [[u64; 3]; 5],
    pub clip_fraction_histogram: [u64; 5],
    pub alignment_score_gap_histogram: [u64; 4],
    pub top_ambiguous_genes: Vec<NamedCount>,
}

#[derive(Serialize)]
pub struct NamedCount {
    pub name: String,
    pub reads: u64,
}

#[derive(Default, Serialize)]
pub struct MappingLenBin {
    pub len_lo: usize,
    pub len_hi: usize,
    pub reads: u64,
    pub unique_nh1: u64,
    pub gene_called: u64,
}

#[derive(Default, Serialize)]
pub struct MoleculeFunnel {
    pub reads_with_cb_umi: u64,
    pub raw_umi_nodes_est: u64,
    pub corrected_cb_umi_nodes_est: u64,
    pub gene_assigned_candidate_nodes_est: u64,
    pub eligible_nodes_est: u64,
    pub representatives: u64,
    pub representatives_in_called_cells: Option<u64>,
}

#[derive(Default, Serialize)]
pub struct SaturationPoint {
    pub read_fraction: f64,
    pub estimated_cb_umi_nodes: u64,
}

#[derive(Default, Serialize)]
pub struct ClipClasses {
    pub n_clips: u64,
    pub adapter: u64,
    pub tso: u64,
    pub poly_a: u64,
    pub poly_t: u64,
    pub unexplained: u64,
}

#[derive(Default, Serialize, Clone)]
pub struct LenBin {
    pub len_lo: usize,
    pub len_hi: usize,
    pub reads: u64,
    pub unique_nh1: u64,
    pub mapq255: u64,
    pub exonic: u64,
    pub intronic: u64,
    pub intergenic: u64,
    pub gene_called: u64,
    pub multi_gene: u64,
    pub mismatches: u64,
    pub aligned_bases: u64,
    pub clip5: u64,
    pub clip3: u64,
}

#[derive(Default, Serialize)]
pub struct MoleculesPerCell {
    pub mean: f64,
    pub median: f64,
    pub p10: f64,
    pub p90: f64,
    pub max: u64,
}

#[derive(Default, Serialize)]
pub struct CollisionOutcomeRecords {
    pub unique: u64,
    pub won: u64,
    pub lost: u64,
    pub tie: u64,
}

#[derive(Default, Serialize)]
pub struct Audit {
    pub label: String,
    pub bam: String,
    pub records: u64,
    pub primary: u64,
    pub mapped: u64,
    pub unmapped: u64,
    pub unique_nh1: u64,
    pub mapq255: u64,
    pub exonic: u64,
    pub intronic: u64,
    pub intergenic: u64,
    pub gene_called: u64,
    pub unique_gene_called: u64,
    pub multi_gene: u64,
    pub gene_degree: [u64; 5],
    pub transcript_degree: [u64; 5],
    pub terminal_collision_evaluated: u64,
    pub terminal_collision_resolved: u64,
    pub terminal_collision_unresolved: u64,
    pub terminal_collision_unsupported: u64,
    pub mismatches: u64,
    pub aligned_bases: u64,
    pub read_bases: u64,
    pub read_len_max: usize,
    pub softclip_bases: u64,
    pub heavily_clipped: u64,
    pub clip5: ClipClasses,
    pub clip3: ClipClasses,
    pub reads_with_cb_umi: u64,
    pub molecules: u64,
    pub xf_17_eligible: u64,
    pub xf_25_representative: u64,
    pub xf_confident_transcriptome: u64,
    pub xf_discarded_feature: u64,
    pub xf_discordant_genes: u64,
    pub xf_representative: u64,
    pub xf_single_feature: u64,
    pub xf_targeted_filtered: u64,
    pub umi_corrected_records: u64,
    pub raw_umi_nodes_sampled: u64,
    pub raw_umi_nodes_est: u64,
    pub collision_outcome_records: CollisionOutcomeRecords,
    pub mapping_diagnostics: MappingDiagnostics,
    pub aligned_length_bins: Vec<MappingLenBin>,
    pub molecule_funnel: MoleculeFunnel,
    pub read_support_histogram: [u64; 7],
    pub saturation_curve: Vec<SaturationPoint>,
    pub genes_detected: u64,
    pub cells_detected: u64,
    pub molecules_per_cell: MoleculesPerCell,
    pub molecules_emitted: u64,
    pub molecule_sample: u64,
    pub sampled_molecules: u64,
    pub sampled_reads_with_cb_umi: u64,
    pub distinct_molecules_est: u64,
    pub reads_per_molecule: f64,
    pub bins: Vec<LenBin>,
    pub elapsed_secs: f64,
    pub records_per_sec: u64,
}
