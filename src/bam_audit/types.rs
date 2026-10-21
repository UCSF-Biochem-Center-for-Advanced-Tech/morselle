use clap::Parser;
use serde::Serialize;
use std::path::PathBuf;

/// Longest read-length bin edge; everything above lands in the final bin.
pub const BIN: usize = 25;
pub const N_BINS: usize = 29; // 0..700+ in 25bp steps

pub const MIN_CLIP: usize = 5; // a clip must be this long to be classified
pub const K: usize = 8; // k-mer seed for indel-robust motif matching
pub const HOMO: f64 = 0.8; // fraction of one base to call a homopolymer clip

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

/// Per-molecule alignment evidence, accumulated across ALL of a molecule's reads
/// (not just the counted representative). This is what lets us adjudicate, for the
/// SAME molecule seen on two platforms, which one mapped it better: uniquely vs
/// multi-mapped (`min_nh`), into a gene body vs intergenic (`region`), and whether
/// it was confidently counted with a gene (`counted` / `rep_gene`).
#[derive(Default)]
pub struct MolAccum {
    pub n_reads: u32,
    pub min_nh: u32, // 1 = uniquely placed somewhere
    pub max_mapq: u8,
    pub region: u8, // bit2=exonic, bit1=intronic, bit0=intergenic (best kept)
    pub max_readlen: u32,
    pub counted: bool,    // any read is a counted representative (xf bit 8)
    pub eligible: bool,   // any read has 10x reference/morselle xf bit 16
    pub rep_gene: String, // gene of the counted representative
}

#[derive(Default, Serialize)]
pub struct MappingDiagnostics {
    /// NH buckets: 1, 2, 3-5, 6-10, 11+.
    pub nh_histogram: [u64; 5],
    /// MAPQ buckets: 0, 1-2, 3-254, 255.
    pub mapq_histogram: [u64; 4],
    /// Rows are NH buckets; columns are GX degree 0, 1, 2+.
    pub nh_by_gene_degree: [[u64; 3]; 5],
    /// Reads by soft-clip fraction: 0, (0,10], (10,25], (25,50], >50%.
    pub clip_fraction_histogram: [u64; 5],
    /// AS-XS score-gap buckets: missing either score, <=5, 6-20, >20.
    pub alignment_score_gap_histogram: [u64; 4],
    /// Gene IDs most often called on NH>1 reads (top 25, after the scan).
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
    // totals
    pub records: u64,
    pub primary: u64,
    pub mapped: u64,
    pub unmapped: u64,
    // mapping / annotation (primary mapped only)
    pub unique_nh1: u64,
    pub mapq255: u64,
    pub exonic: u64,
    pub intronic: u64,
    pub intergenic: u64,
    pub gene_called: u64,
    pub unique_gene_called: u64,
    pub multi_gene: u64,
    /// Compatibility degree histograms: buckets 0, 1, 2, 3, and 4+.
    pub gene_degree: [u64; 5],
    pub transcript_degree: [u64; 5],
    pub terminal_collision_evaluated: u64,
    pub terminal_collision_resolved: u64,
    pub terminal_collision_unresolved: u64,
    pub terminal_collision_unsupported: u64,
    // fidelity
    pub mismatches: u64,
    pub aligned_bases: u64,
    // lengths
    pub read_bases: u64,
    pub read_len_max: usize,
    // trimming
    pub softclip_bases: u64,
    pub heavily_clipped: u64,
    pub clip5: ClipClasses,
    pub clip3: ClipClasses,
    // molecules
    pub reads_with_cb_umi: u64,
    pub molecules: u64, // xf bit-8 representatives (10x reference's counted molecules)
    /// Exact 10x reference/morselle xf values most relevant to molecule counting.
    pub xf_17_eligible: u64,
    pub xf_25_representative: u64,
    /// Individual xf bit counts (a record may contribute to several fields).
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
    /// Sampled CB+UB read-support histogram: 1, 2, 3, 4, 5-9, 10-19, 20+ reads.
    pub read_support_histogram: [u64; 7],
    /// Expected observed CB+UB nodes after random read downsampling, estimated from
    /// the consistent molecule hash sample.
    pub saturation_curve: Vec<SaturationPoint>,
    pub genes_detected: u64,
    pub cells_detected: u64,
    pub molecules_per_cell: MoleculesPerCell,
    pub molecules_emitted: u64,
    /// Subsample rate used for both the dump and the distinct-molecule estimate.
    pub molecule_sample: u64,
    /// Distinct (CB,UMI) seen in the subsample, and the reads landing on them.
    /// The gate depends only on (CB,UMI), so every read of a sampled molecule is
    /// kept — making `reads_per_molecule` an unbiased estimate in bounded memory.
    pub sampled_molecules: u64,
    pub sampled_reads_with_cb_umi: u64,
    pub distinct_molecules_est: u64,
    pub reads_per_molecule: f64,
    // per length bin
    pub bins: Vec<LenBin>,
    // run info
    pub elapsed_secs: f64,
    pub records_per_sec: u64,
}
