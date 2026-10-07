use std::collections::HashMap;
use std::path::PathBuf;
use clap::Parser;
use serde::Serialize;

pub const BIN_SIZE: i64 = 16_384;

#[derive(Debug, Parser)]
#[command(
    name = "sbx10x-annotate",
    version,
    about = "Add 10x-compatible GX/GN/TX/AN/RE tags to a genomic BAM"
)]
pub struct Cli {
    /// Genomic BAM/SAM produced by minimap2, or `-` for standard input.
    #[arg(value_name = "BAM")]
    pub bam: PathBuf,

    /// GENCODE-compatible GTF, optionally gzip-compressed.
    #[arg(short, long, value_name = "GTF")]
    pub gtf: PathBuf,

    /// Tagged BAM output, or `-` for BAM on standard output.
    #[arg(short, long, value_name = "BAM")]
    pub output: PathBuf,

    /// Count same-strand gene-body overlaps as intronic assignments.
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    pub include_introns: bool,

    /// BAM decompression/compression threads.
    #[arg(short, long, default_value_t = 4)]
    pub threads: usize,

    /// Annotate records in parallel and permit output reordering within batches.
    /// Intended for pipelines that immediately coordinate-sort the output.
    #[arg(long)]
    pub unordered: bool,

    /// Records per parallel annotation batch when --unordered is enabled.
    #[arg(long, default_value_t = 8192)]
    pub batch_size: usize,

    /// Input is query-grouped primary/secondary minimap2 output. Add `NH` (reported
    /// hit count), `HI`, and an `mm`-like single-gene multimap rescue tag per group.
    /// Incompatible with --unordered because all alignments for a query stay together.
    #[arg(long)]
    pub multimap_tags: bool,

    /// Optional JSON metrics path (default: OUTPUT.metrics.json).
    #[arg(long)]
    pub metrics: Option<PathBuf>,

    /// Canonical gene feature manifest. Creates it from the filtered GTF when
    /// absent; when present, verifies byte-for-byte identity. Pass the same file
    /// to `sbxpress matrix --features` for stable rows across experiments.
    #[arg(long, value_name = "GENES.MTX")]
    pub genes_mtx: Option<PathBuf>,

    /// Apply the 10x GENCODE gene/transcript biotype, readthrough, and PAR rules.
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    pub reference_filter: bool,

    /// Remove numeric version suffixes from Ensembl gene/transcript IDs.
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    pub strip_id_versions: bool,

    /// Use terminal-exon support as a non-destructive preference signal for a 3' assay.
    #[arg(long)]
    pub three_prime: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Strand {
    Plus,
    Minus,
}

impl Strand {
    pub fn parse(value: &str) -> Option<Self> {
        match value {
            "+" => Some(Self::Plus),
            "-" => Some(Self::Minus),
            _ => None,
        }
    }

    pub fn matches_record(self, reverse: bool) -> bool {
        matches!((self, reverse), (Self::Plus, false) | (Self::Minus, true))
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Interval {
    pub start: i64,
    pub end: i64,
}

impl Interval {
    pub fn overlaps(self, other: Self) -> bool {
        self.start < other.end && other.start < self.end
    }

    pub fn contains(self, other: Self) -> bool {
        self.start <= other.start && self.end >= other.end
    }
}

#[derive(Debug)]
pub struct Gene {
    pub id: String,
    pub name: String,
    pub contig: usize,
    pub strand: Strand,
    pub span: Interval,
}

#[derive(Debug)]
pub struct Transcript {
    pub id: String,
    pub gene: usize,
    pub contig: usize,
    pub strand: Strand,
    pub exons: Vec<Interval>, // genomic ascending order
}

#[derive(Debug, Default)]
pub struct ContigIndex {
    pub exon_bins: Vec<Vec<usize>>, // transcript indexes
    pub gene_bins: Vec<Vec<usize>>,
}

#[derive(Debug, Default)]
pub struct Annotation {
    pub contig_names: Vec<String>,
    pub contig_ids: HashMap<String, usize>,
    pub genes: Vec<Gene>,
    pub transcripts: Vec<Transcript>,
    pub indexes: Vec<ContigIndex>,
    pub exons_seen: u64,
    pub exons_filtered_biotype: u64,
    pub exons_filtered_tag: u64,
}

#[derive(Debug)]
pub struct Call {
    pub region: u8,
    pub genes: Vec<usize>,
    pub tx: Vec<String>,
    pub antisense: Vec<String>,
    pub terminal_genes: Vec<usize>,
    pub terminal_tx: Vec<String>,
}

#[derive(Debug, Default, Serialize)]
pub struct Metrics {
    pub input_bam: String,
    pub output_bam: String,
    pub gtf: String,
    pub genes_mtx: Option<String>,
    pub include_introns: bool,
    pub reference_filter: bool,
    pub strip_id_versions: bool,
    pub three_prime: bool,
    pub unordered: bool,
    pub multimap_tags: bool,
    pub elapsed_secs: f64,
    pub records: u64,
    pub primary: u64,
    pub supplementary: u64,
    pub unmapped: u64,
    pub exonic: u64,
    pub intronic: u64,
    pub intergenic: u64,
    pub gene_assigned: u64,
    pub unique_gene: u64,
    pub multi_gene: u64,
    pub transcript_assigned: u64,
    pub antisense: u64,
    pub collision_records: u64,
    pub gene_degree_0: u64,
    pub gene_degree_1: u64,
    pub gene_degree_2: u64,
    pub gene_degree_3plus: u64,
    pub transcript_degree_0: u64,
    pub transcript_degree_1: u64,
    pub transcript_degree_2: u64,
    pub transcript_degree_3plus: u64,
    pub terminal_resolved_gene_collision: u64,
    pub terminal_unresolved_gene_collision: u64,
    pub terminal_unsupported_gene_collision: u64,
    pub genes_loaded: usize,
    pub transcripts_loaded: usize,
    pub gtf_exons_seen: u64,
    pub gtf_exons_filtered_biotype: u64,
    pub gtf_exons_filtered_tag: u64,
}

impl Metrics {
    pub fn add_record_metrics(&mut self, value: RecordMetrics) {
        self.records += value.records;
        self.primary += value.primary;
        self.supplementary += value.supplementary;
        self.unmapped += value.unmapped;
        self.exonic += value.exonic;
        self.intronic += value.intronic;
        self.intergenic += value.intergenic;
        self.gene_assigned += value.gene_assigned;
        self.unique_gene += value.unique_gene;
        self.multi_gene += value.multi_gene;
        self.transcript_assigned += value.transcript_assigned;
        self.antisense += value.antisense;
        self.collision_records += value.collision_records;
        self.gene_degree_0 += value.gene_degree_0;
        self.gene_degree_1 += value.gene_degree_1;
        self.gene_degree_2 += value.gene_degree_2;
        self.gene_degree_3plus += value.gene_degree_3plus;
        self.transcript_degree_0 += value.transcript_degree_0;
        self.transcript_degree_1 += value.transcript_degree_1;
        self.transcript_degree_2 += value.transcript_degree_2;
        self.transcript_degree_3plus += value.transcript_degree_3plus;
        self.terminal_resolved_gene_collision += value.terminal_resolved_gene_collision;
        self.terminal_unresolved_gene_collision += value.terminal_unresolved_gene_collision;
        self.terminal_unsupported_gene_collision += value.terminal_unsupported_gene_collision;
    }
}

#[derive(Debug, Default)]
pub struct RecordMetrics {
    pub records: u64,
    pub primary: u64,
    pub supplementary: u64,
    pub unmapped: u64,
    pub exonic: u64,
    pub intronic: u64,
    pub intergenic: u64,
    pub gene_assigned: u64,
    pub unique_gene: u64,
    pub multi_gene: u64,
    pub transcript_assigned: u64,
    pub antisense: u64,
    pub collision_records: u64,
    pub gene_degree_0: u64,
    pub gene_degree_1: u64,
    pub gene_degree_2: u64,
    pub gene_degree_3plus: u64,
    pub transcript_degree_0: u64,
    pub transcript_degree_1: u64,
    pub transcript_degree_2: u64,
    pub transcript_degree_3plus: u64,
    pub terminal_resolved_gene_collision: u64,
    pub terminal_unresolved_gene_collision: u64,
    pub terminal_unsupported_gene_collision: u64,
}
