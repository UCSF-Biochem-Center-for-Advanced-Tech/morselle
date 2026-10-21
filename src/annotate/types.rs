use clap::Parser;
use serde::Serialize;
use std::collections::HashMap;
use std::path::PathBuf;

pub const BIN_SIZE: i64 = 16_384;

#[derive(Debug, Parser)]
#[command(
    name = "morselle-annotate",
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

    /// Cell Ranger's confident-mapping MAPQ (60 for minimap2). A primary read
    /// below it whose genes come from compatible sense exonic transcripts of
    /// exactly one gene is tagged `mm:i:1` (rescued; MAPQ is kept), which
    /// `morselle molecules --min-mapq` counts. Intronic hits are never rescued.
    /// With --multimap-tags the rule spans the query group. 0 restores the
    /// pre-0.6.2 tags (no `mm`; the looser group rule with --multimap-tags).
    #[arg(long, default_value_t = 60)]
    pub confident_mapq: u8,

    /// Bases a read's first or last aligned block may extend past its exon and
    /// still be transcript-compatible (treated as soft-clipped). minimap2
    /// extends rather than splices a 1-3 base end anchor, which otherwise
    /// leaves the read exonic with no gene. Junctions must still match
    /// exactly. 0 requires exact exon ends (pre-0.6.2 behaviour).
    #[arg(long, default_value_t = 3)]
    pub max_end_overhang: u32,

    /// Optional JSON metrics path (default: OUTPUT.metrics.json).
    #[arg(long)]
    pub metrics: Option<PathBuf>,

    /// Also write per-barcode counters to `<metrics>.barcodes.tsv.gz` (MANUAL
    /// "Per-barcode read metrics"). Off by default; it adds run time and memory.
    #[arg(long)]
    pub barcode_metrics: bool,

    /// Canonical gene feature manifest. Creates it from the filtered GTF when
    /// absent; when present, verifies byte-for-byte identity. Pass the same file
    /// to `morselle matrix --features` for stable rows across experiments.
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

/// Per-record settings shared by every annotation path.
#[derive(Clone, Copy, Debug)]
pub struct RecordOptions {
    pub include_introns: bool,
    pub three_prime: bool,
    /// `Cli::confident_mapq`; 0 disables the per-record `mm` rescue.
    pub confident_mapq: u8,
    /// `Cli::max_end_overhang`; 0 requires exact exon ends.
    pub max_end_overhang: u32,
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
    /// `genes` come from splice-compatible sense exonic transcripts (not the
    /// intronic gene-body fallback).
    pub exonic_compatible: bool,
    /// Exonic-compatible only because of the end-overhang allowance.
    pub end_overhang: bool,
}

#[derive(Debug, Default, Serialize)]
pub struct Metrics {
    #[serde(skip)]
    pub per_barcode: crate::barcode_metrics::BarcodeTable<12>,
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
    pub confident_mapq: u8,
    pub max_end_overhang: u32,
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
    /// Primary reads below --confident-mapq tagged `mm:i:1`.
    pub rescued: u64,
    /// Records given genes only through --max-end-overhang.
    pub end_overhang: u64,
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
        self.rescued += value.rescued;
        self.end_overhang += value.end_overhang;
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
    pub rescued: u64,
    pub end_overhang: u64,
    /// The record's genes are exonic-compatible (`Call::exonic_compatible`);
    /// read by the query-group rescue, not a counter.
    pub exonic_compatible: bool,
}
