//! `bam-audit` — single-pass, O(1)-memory audit of a Cell Ranger BAM.
//!
//! Gathers everything needed for the R2 contamination audit and the long-read
//! mapping analysis, so a full (multi-hundred-GB) BAM can be summarized on the
//! server and only a small JSON moved back for analysis.
//!
//! Per BAM it reports, globally and per read-length bin:
//!   * mapping    — unique (NH==1), MAPQ 255, unmapped
//!   * annotation — RE tag region (Exonic / iNtronic / Intergenic), GX gene call,
//!                  multi-gene (ambiguous) GX
//!   * fidelity   — nM mismatches per 100 aligned bases
//!   * trimming   — 5'/3' soft clips, and WHAT they contain (adapter / TSO / polyA /
//!                  polyT / unexplained) using indel-robust 8-mer matching. Residual
//!                  adapter or polyT cannot align, so STAR soft-clips it; that is the
//!                  contamination signal.
//!   * molecules  — counted molecules (xf bit 8 representatives), reads/molecule,
//!                  and UMIs per read.
//!
//! Geometry note: sbx10x emits R2 reverse-complemented, so residual 3' read-through
//! adapter lands at the 5' end of R2 and residual polyT lands at the 3' end as polyA.
//! Both motifs are tested at BOTH ends so the result isn't biased by that expectation.
//!
//! Cross-BAM molecule concordance (does a long read call the same gene as a short
//! one?) would normally need a hash map over hundreds of millions of molecules.
//! Instead `--molecules` dumps a CONSISTENT hash-subsample: a molecule is emitted iff
//! `hash(CB,UMI) % --molecule-sample == 0`. The same molecules are therefore selected
//! in every BAM, so the dumps can simply be joined, at 1/N the size.
//!
//! Usage:
//!   bam-audit possorted_genome_bam.bam --label "SBX p5" -o p5_audit.json \
//!       --threads 8 --molecules p5_molecules.tsv.gz --molecule-sample 64

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::time::Instant;

use anyhow::{Context, Result};
use clap::{CommandFactory, FromArgMatches, Parser};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use rust_htslib::bam::{self, record::Cigar, Read};
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Serialize;

/// Longest read-length bin edge; everything above lands in the final bin.
const BIN: usize = 25;
const N_BINS: usize = 29; // 0..700+ in 25bp steps

const MIN_CLIP: usize = 5; // a clip must be this long to be classified
const K: usize = 8; // k-mer seed for indel-robust motif matching
const HOMO: f64 = 0.8; // fraction of one base to call a homopolymer clip

const ADAPTER: &[u8] = b"AGATCGGAAGAGC";
const TSO: &[u8] = b"AAGCAGTGGTATCAACGCAGAGTACATGGG";

#[derive(Parser, Debug)]
#[command(
    name = "bam-audit",
    version,
    about = "Single-pass audit of a Cell Ranger BAM"
)]
struct Cli {
    /// Input Cell Ranger BAM (possorted_genome_bam.bam).
    bam: PathBuf,

    /// Label for this BAM in the JSON (e.g. "ILLUMINA", "SBX p5").
    #[arg(long, default_value = "bam")]
    label: String,

    /// Output JSON summary.
    #[arg(short, long, default_value = "bam_audit.json")]
    output: PathBuf,

    /// BAM decompression threads.
    #[arg(long, default_value_t = 4)]
    threads: usize,

    /// Optional gzip TSV of hash-subsampled molecule representatives:
    /// `CB<TAB>UMI<TAB>GX<TAB>read_len`. Consistent across BAMs, so joinable.
    #[arg(long)]
    molecules: Option<PathBuf>,

    /// Emit 1 in N molecules to `--molecules` (same N across BAMs!). 1 = all.
    #[arg(long, default_value_t = 64)]
    molecule_sample: u64,

    /// Optional called-cell barcode list (Cell Ranger filtered barcodes.tsv.gz or
    /// sbxpress matrix barcodes.tsv.gz). Enables final-matrix molecule accounting.
    #[arg(long)]
    cells: Option<PathBuf>,
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

/// All k-mers of `m` and its reverse complement — indel-robust motif detection.
fn kmers(m: &[u8]) -> Vec<[u8; K]> {
    let mut v = Vec::new();
    for s in [m.to_vec(), revcomp(m)] {
        for w in s.windows(K) {
            let mut a = [0u8; K];
            a.copy_from_slice(w);
            v.push(a);
        }
    }
    v.sort_unstable();
    v.dedup();
    v
}

fn has_kmer(clip: &[u8], set: &[[u8; K]]) -> bool {
    if clip.len() < K {
        return false;
    }
    clip.windows(K).any(|w| {
        set.binary_search(&{
            let mut a = [0u8; K];
            a.copy_from_slice(w);
            a
        })
        .is_ok()
    })
}

/// SplitMix64 — stable across machines, so the molecule subsample is identical
/// for every BAM (that is what makes the dumps joinable).
fn splitmix64(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E3779B97F4A7C15);
    let mut x = z;
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D049BB133111EB);
    x ^ (x >> 31)
}

/// Stable FNV-1a hash for variable-length molecule components. This is used only
/// for bounded-memory cardinality sampling, not identity reconstruction.
fn stable_hash(parts: &[&[u8]]) -> u64 {
    let mut value = 0xcbf29ce484222325u64;
    for part in parts {
        for &byte in *part {
            value ^= byte as u64;
            value = value.wrapping_mul(0x100000001b3);
        }
        value ^= 0xff;
        value = value.wrapping_mul(0x100000001b3);
    }
    value
}

/// 2-bit pack a barcode/UMI; None on any non-ACGT base.
fn pack(s: &[u8]) -> Option<u64> {
    let mut c = 0u64;
    for &b in s {
        c = (c << 2)
            | match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => return None,
            };
    }
    Some(c)
}

/// Decode a 2-bit packed code of `len` bases.
fn unpack(mut code: u64, len: usize) -> String {
    let mut s = vec![b'A'; len];
    for i in (0..len).rev() {
        s[i] = b"ACGT"[(code & 3) as usize];
        code >>= 2;
    }
    String::from_utf8(s).unwrap()
}

/// Per-molecule alignment evidence, accumulated across ALL of a molecule's reads
/// (not just the counted representative). This is what lets us adjudicate, for the
/// SAME molecule seen on two platforms, which one mapped it better: uniquely vs
/// multi-mapped (`min_nh`), into a gene body vs intergenic (`region`), and whether
/// it was confidently counted with a gene (`counted` / `rep_gene`).
#[derive(Default)]
struct MolAccum {
    n_reads: u32,
    min_nh: u32, // 1 = uniquely placed somewhere
    max_mapq: u8,
    region: u8, // bit2=exonic, bit1=intronic, bit0=intergenic (best kept)
    max_readlen: u32,
    counted: bool,    // any read is a counted representative (xf bit 8)
    eligible: bool,   // any read has Cell Ranger/sbxpress xf bit 16
    rep_gene: String, // gene of the counted representative
}

#[derive(Default, Serialize)]
struct MappingDiagnostics {
    /// NH buckets: 1, 2, 3-5, 6-10, 11+.
    nh_histogram: [u64; 5],
    /// MAPQ buckets: 0, 1-2, 3-254, 255.
    mapq_histogram: [u64; 4],
    /// Rows are NH buckets; columns are GX degree 0, 1, 2+.
    nh_by_gene_degree: [[u64; 3]; 5],
    /// Reads by soft-clip fraction: 0, (0,10], (10,25], (25,50], >50%.
    clip_fraction_histogram: [u64; 5],
    /// AS-XS score-gap buckets: missing either score, <=5, 6-20, >20.
    alignment_score_gap_histogram: [u64; 4],
    /// Gene IDs most often called on NH>1 reads (top 25, after the scan).
    top_ambiguous_genes: Vec<NamedCount>,
}

#[derive(Serialize)]
struct NamedCount {
    name: String,
    reads: u64,
}

#[derive(Default, Serialize)]
struct MappingLenBin {
    len_lo: usize,
    len_hi: usize,
    reads: u64,
    unique_nh1: u64,
    gene_called: u64,
}

#[derive(Default, Serialize)]
struct MoleculeFunnel {
    reads_with_cb_umi: u64,
    raw_umi_nodes_est: u64,
    corrected_cb_umi_nodes_est: u64,
    gene_assigned_candidate_nodes_est: u64,
    eligible_nodes_est: u64,
    representatives: u64,
    representatives_in_called_cells: Option<u64>,
}

#[derive(Default, Serialize)]
struct SaturationPoint {
    read_fraction: f64,
    estimated_cb_umi_nodes: u64,
}

#[derive(Default, Serialize)]
struct ClipClasses {
    n_clips: u64,
    adapter: u64,
    tso: u64,
    poly_a: u64,
    poly_t: u64,
    unexplained: u64,
}

#[derive(Default, Serialize, Clone)]
struct LenBin {
    len_lo: usize,
    len_hi: usize,
    reads: u64,
    unique_nh1: u64,
    mapq255: u64,
    exonic: u64,
    intronic: u64,
    intergenic: u64,
    gene_called: u64,
    multi_gene: u64,
    mismatches: u64,
    aligned_bases: u64,
    clip5: u64,
    clip3: u64,
}

#[derive(Default, Serialize)]
struct MoleculesPerCell {
    mean: f64,
    median: f64,
    p10: f64,
    p90: f64,
    max: u64,
}

#[derive(Default, Serialize)]
struct CollisionOutcomeRecords {
    unique: u64,
    won: u64,
    lost: u64,
    tie: u64,
}

#[derive(Default, Serialize)]
struct Audit {
    label: String,
    bam: String,
    // totals
    records: u64,
    primary: u64,
    mapped: u64,
    unmapped: u64,
    // mapping / annotation (primary mapped only)
    unique_nh1: u64,
    mapq255: u64,
    exonic: u64,
    intronic: u64,
    intergenic: u64,
    gene_called: u64,
    unique_gene_called: u64,
    multi_gene: u64,
    /// Compatibility degree histograms: buckets 0, 1, 2, 3, and 4+.
    gene_degree: [u64; 5],
    transcript_degree: [u64; 5],
    terminal_collision_evaluated: u64,
    terminal_collision_resolved: u64,
    terminal_collision_unresolved: u64,
    terminal_collision_unsupported: u64,
    // fidelity
    mismatches: u64,
    aligned_bases: u64,
    // lengths
    read_bases: u64,
    read_len_max: usize,
    // trimming
    softclip_bases: u64,
    heavily_clipped: u64,
    clip5: ClipClasses,
    clip3: ClipClasses,
    // molecules
    reads_with_cb_umi: u64,
    molecules: u64, // xf bit-8 representatives (Cell Ranger's counted molecules)
    /// Exact Cell Ranger/sbxpress xf values most relevant to molecule counting.
    xf_17_eligible: u64,
    xf_25_representative: u64,
    /// Individual xf bit counts (a record may contribute to several fields).
    xf_confident_transcriptome: u64,
    xf_discarded_feature: u64,
    xf_discordant_genes: u64,
    xf_representative: u64,
    xf_single_feature: u64,
    xf_targeted_filtered: u64,
    umi_corrected_records: u64,
    raw_umi_nodes_sampled: u64,
    raw_umi_nodes_est: u64,
    collision_outcome_records: CollisionOutcomeRecords,
    mapping_diagnostics: MappingDiagnostics,
    aligned_length_bins: Vec<MappingLenBin>,
    molecule_funnel: MoleculeFunnel,
    /// Sampled CB+UB read-support histogram: 1, 2, 3, 4, 5-9, 10-19, 20+ reads.
    read_support_histogram: [u64; 7],
    /// Expected observed CB+UB nodes after random read downsampling, estimated from
    /// the consistent molecule hash sample.
    saturation_curve: Vec<SaturationPoint>,
    genes_detected: u64,
    cells_detected: u64,
    molecules_per_cell: MoleculesPerCell,
    molecules_emitted: u64,
    /// Subsample rate used for both the dump and the distinct-molecule estimate.
    molecule_sample: u64,
    /// Distinct (CB,UMI) seen in the subsample, and the reads landing on them.
    /// The gate depends only on (CB,UMI), so every read of a sampled molecule is
    /// kept — making `reads_per_molecule` an unbiased estimate in bounded memory.
    sampled_molecules: u64,
    sampled_reads_with_cb_umi: u64,
    distinct_molecules_est: u64,
    reads_per_molecule: f64,
    // per length bin
    bins: Vec<LenBin>,
    // run info
    elapsed_secs: f64,
    records_per_sec: u64,
}

fn aux_int(rec: &bam::Record, tag: &[u8]) -> Option<i64> {
    use rust_htslib::bam::record::Aux;
    match rec.aux(tag) {
        Ok(Aux::I8(v)) => Some(v as i64),
        Ok(Aux::U8(v)) => Some(v as i64),
        Ok(Aux::I16(v)) => Some(v as i64),
        Ok(Aux::U16(v)) => Some(v as i64),
        Ok(Aux::I32(v)) => Some(v as i64),
        Ok(Aux::U32(v)) => Some(v as i64),
        _ => None,
    }
}

fn aux_str<'a>(rec: &'a bam::Record, tag: &[u8]) -> Option<&'a str> {
    use rust_htslib::bam::record::Aux;
    match rec.aux(tag) {
        Ok(Aux::String(s)) => Some(s),
        _ => None,
    }
}

fn aux_char(rec: &bam::Record, tag: &[u8]) -> Option<u8> {
    use rust_htslib::bam::record::Aux;
    match rec.aux(tag) {
        Ok(Aux::Char(c)) => Some(c),
        _ => None,
    }
}

/// Classify a soft-clipped sequence: what did the aligner refuse to align?
fn classify(clip: &[u8], adapter_k: &[[u8; K]], tso_k: &[[u8; K]]) -> &'static str {
    let n = clip.len() as f64;
    let a = clip.iter().filter(|&&b| b == b'A').count() as f64 / n;
    let t = clip.iter().filter(|&&b| b == b'T').count() as f64 / n;
    if has_kmer(clip, adapter_k) {
        "adapter"
    } else if has_kmer(clip, tso_k) {
        "tso"
    } else if a >= HOMO {
        "polyA"
    } else if t >= HOMO {
        "polyT"
    } else {
        "unexplained"
    }
}

fn bump(c: &mut ClipClasses, kind: &str) {
    c.n_clips += 1;
    match kind {
        "adapter" => c.adapter += 1,
        "tso" => c.tso += 1,
        "polyA" => c.poly_a += 1,
        "polyT" => c.poly_t += 1,
        _ => c.unexplained += 1,
    }
}

fn load_called_cells(path: &PathBuf) -> Result<FxHashSet<String>> {
    let file =
        File::open(path).with_context(|| format!("cannot open cells file {}", path.display()))?;
    let reader: Box<dyn BufRead> = if path.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let mut cells = FxHashSet::default();
    for line in reader.lines() {
        let line = line?;
        let barcode = line.split('\t').next().unwrap_or("").trim();
        if !barcode.is_empty() {
            cells.insert(barcode.to_string());
        }
    }
    Ok(cells)
}

fn nh_bucket(nh: i64) -> usize {
    match nh {
        ..=1 => 0,
        2 => 1,
        3..=5 => 2,
        6..=10 => 3,
        _ => 4,
    }
}

fn main() -> Result<()> {
    let cli = if std::env::var_os("SBXPRESS_DISPATCH").is_some() {
        Cli::from_arg_matches(&Cli::command().bin_name("sbxpress audit").get_matches())?
    } else {
        Cli::parse()
    };
    let t0 = Instant::now();
    let adapter_k = kmers(ADAPTER);
    let tso_k = kmers(TSO);
    let called_cells = match &cli.cells {
        Some(path) => Some(load_called_cells(path)?),
        None => None,
    };

    let mut reader = bam::Reader::from_path(&cli.bam)
        .with_context(|| format!("cannot open BAM {}", cli.bam.display()))?;
    reader.set_threads(cli.threads.max(1)).ok();

    let mut mol_out = match &cli.molecules {
        Some(p) => Some(GzEncoder::new(
            BufWriter::new(
                File::create(p).with_context(|| format!("cannot create {}", p.display()))?,
            ),
            Compression::new(4),
        )),
        None => None,
    };

    let mut a = Audit {
        label: cli.label.clone(),
        bam: cli.bam.display().to_string(),
        molecule_sample: cli.molecule_sample.max(1),
        bins: (0..N_BINS)
            .map(|i| LenBin {
                len_lo: i * BIN,
                len_hi: i * BIN + BIN - 1,
                ..Default::default()
            })
            .collect(),
        aligned_length_bins: (0..N_BINS)
            .map(|i| MappingLenBin {
                len_lo: i * BIN,
                len_hi: i * BIN + BIN - 1,
                ..Default::default()
            })
            .collect(),
        ..Default::default()
    };
    // Per-molecule evidence for the subsample. Bounded: ~distinct_molecules / sample.
    let mut sampled: FxHashMap<u64, MolAccum> = FxHashMap::default();
    let mut counted_per_cell: FxHashMap<String, u64> = FxHashMap::default();
    let mut detected_genes: FxHashSet<String> = FxHashSet::default();
    let mut sampled_raw_umi_nodes: FxHashSet<u64> = FxHashSet::default();
    let mut sampled_gene_candidates: FxHashSet<u64> = FxHashSet::default();
    let mut representatives_in_called_cells = 0u64;
    let mut ambiguous_genes: FxHashMap<String, u64> = FxHashMap::default();

    let mut rec = bam::Record::new();
    while let Some(r) = reader.read(&mut rec) {
        r?;
        a.records += 1;
        if rec.is_secondary() || rec.is_supplementary() {
            continue;
        }
        a.primary += 1;
        if rec.is_unmapped() {
            a.unmapped += 1;
            continue;
        }
        a.mapped += 1;

        // --- read-orientation sequence + cigar (BAM stores both reference-oriented) ---
        let mut seq = rec.seq().as_bytes();
        let mut cig: Vec<Cigar> = rec.cigar().iter().copied().collect();
        if rec.is_reverse() {
            seq = revcomp(&seq);
            cig.reverse();
        }
        let rl = seq.len();
        if rl == 0 {
            continue;
        }
        let c5 = match cig.first() {
            Some(Cigar::SoftClip(n)) => *n as usize,
            _ => 0,
        };
        let c3 = match cig.last() {
            Some(Cigar::SoftClip(n)) => *n as usize,
            _ => 0,
        };
        let alen: usize = cig
            .iter()
            .map(|c| match c {
                Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => *n as usize,
                _ => 0,
            })
            .sum();

        // --- tags ---
        let nh = aux_int(&rec, b"NH").unwrap_or(1);
        let nm = aux_int(&rec, b"nM").unwrap_or(0).max(0) as u64;
        let re = aux_char(&rec, b"RE").unwrap_or(b'?');
        let gx = aux_str(&rec, b"GX").unwrap_or("");
        let xf = aux_int(&rec, b"xf").unwrap_or(0);

        let bi = (rl / BIN).min(N_BINS - 1);
        let b = &mut a.bins[bi];
        let abi = (alen / BIN).min(N_BINS - 1);
        let ab = &mut a.aligned_length_bins[abi];

        a.read_bases += rl as u64;
        a.read_len_max = a.read_len_max.max(rl);
        a.aligned_bases += alen as u64;
        a.mismatches += nm;
        a.softclip_bases += (c5 + c3) as u64;
        if (c5 + c3) as f64 > rl as f64 * 0.5 {
            a.heavily_clipped += 1;
        }
        b.reads += 1;
        b.aligned_bases += alen as u64;
        b.mismatches += nm;
        ab.reads += 1;

        if nh == 1 {
            a.unique_nh1 += 1;
            b.unique_nh1 += 1;
            ab.unique_nh1 += 1;
        }
        if rec.mapq() == 255 {
            a.mapq255 += 1;
            b.mapq255 += 1;
        }
        match re {
            b'E' => {
                a.exonic += 1;
                b.exonic += 1;
            }
            b'N' => {
                a.intronic += 1;
                b.intronic += 1;
            }
            b'I' => {
                a.intergenic += 1;
                b.intergenic += 1;
            }
            _ => {}
        }
        if !gx.is_empty() {
            a.gene_called += 1;
            b.gene_called += 1;
            ab.gene_called += 1;
            if gx.contains(';') {
                a.multi_gene += 1;
                b.multi_gene += 1;
            } else {
                a.unique_gene_called += 1;
            }
        }
        let gx_degree = if gx.is_empty() {
            0
        } else {
            gx.split(';').count()
        };
        let nh_bin = nh_bucket(nh);
        let gx_bin = gx_degree.min(2);
        a.mapping_diagnostics.nh_histogram[nh_bin] += 1;
        a.mapping_diagnostics.nh_by_gene_degree[nh_bin][gx_bin] += 1;
        let mapq_bin = match rec.mapq() {
            0 => 0,
            1..=2 => 1,
            3..=254 => 2,
            _ => 3,
        };
        a.mapping_diagnostics.mapq_histogram[mapq_bin] += 1;
        let clip_fraction_bin = match (c5 + c3) * 100 / rl {
            0 => 0,
            1..=10 => 1,
            11..=25 => 2,
            26..=50 => 3,
            _ => 4,
        };
        a.mapping_diagnostics.clip_fraction_histogram[clip_fraction_bin] += 1;
        let score_gap_bin = match (aux_int(&rec, b"AS"), aux_int(&rec, b"XS")) {
            (Some(ascore), Some(xscore)) if ascore - xscore <= 5 => 1,
            (Some(ascore), Some(xscore)) if ascore - xscore <= 20 => 2,
            (Some(_), Some(_)) => 3,
            _ => 0,
        };
        a.mapping_diagnostics.alignment_score_gap_histogram[score_gap_bin] += 1;
        if nh > 1 && !gx.is_empty() {
            for gene in gx.split(';').filter(|g| !g.is_empty()) {
                *ambiguous_genes.entry(gene.to_string()).or_default() += 1;
            }
        }
        let tx = aux_str(&rec, b"TX").unwrap_or("");
        let tx_degree = if tx.is_empty() {
            0
        } else {
            tx.split(';').count()
        };
        a.gene_degree[gx_degree.min(4)] += 1;
        a.transcript_degree[tx_degree.min(4)] += 1;
        if gx_degree > 1 {
            if let Some(preferred) = aux_int(&rec, b"XP") {
                a.terminal_collision_evaluated += 1;
                match preferred {
                    0 => a.terminal_collision_unsupported += 1,
                    1 => a.terminal_collision_resolved += 1,
                    _ => a.terminal_collision_unresolved += 1,
                }
            }
        }

        // --- soft-clip content (the contamination signal) ---
        if c5 >= MIN_CLIP {
            b.clip5 += 1;
            bump(&mut a.clip5, classify(&seq[..c5], &adapter_k, &tso_k));
        }
        if c3 >= MIN_CLIP {
            b.clip3 += 1;
            bump(&mut a.clip3, classify(&seq[rl - c3..], &adapter_k, &tso_k));
        }

        // --- molecules ---
        let cb_full = aux_str(&rec, b"CB");
        let cb = cb_full.map(|s| s.split('-').next().unwrap_or(s));
        let ub = aux_str(&rec, b"UB");
        let ur = aux_str(&rec, b"UR");
        if ur.zip(ub).is_some_and(|(raw, corrected)| raw != corrected) {
            a.umi_corrected_records += 1;
        }
        if xf == 17 {
            a.xf_17_eligible += 1;
        }
        if xf == 25 {
            a.xf_25_representative += 1;
        }
        a.xf_confident_transcriptome += ((xf & 1) != 0) as u64;
        a.xf_discarded_feature += ((xf & 2) != 0) as u64;
        a.xf_discordant_genes += ((xf & 4) != 0) as u64;
        a.xf_representative += ((xf & 8) != 0) as u64;
        a.xf_single_feature += ((xf & 16) != 0) as u64;
        a.xf_targeted_filtered += ((xf & 32) != 0) as u64;
        match aux_char(&rec, b"XW") {
            Some(b'U') => a.collision_outcome_records.unique += 1,
            Some(b'W') => a.collision_outcome_records.won += 1,
            Some(b'L') => a.collision_outcome_records.lost += 1,
            Some(b'T') => a.collision_outcome_records.tie += 1,
            _ => {}
        }
        if let (Some(cbs), Some(raw_umi)) = (cb, ur) {
            let key = stable_hash(&[cbs.as_bytes(), raw_umi.as_bytes()]);
            if cli.molecule_sample <= 1 || splitmix64(key) % cli.molecule_sample == 0 {
                sampled_raw_umi_nodes.insert(key);
            }
        }
        if let (Some(cbs), Some(ubs)) = (cb, ub) {
            a.reads_with_cb_umi += 1;
            if xf & 8 != 0 {
                a.molecules += 1;
                if cb_full.is_some_and(|full| {
                    called_cells
                        .as_ref()
                        .is_some_and(|cells| cells.contains(full))
                }) {
                    representatives_in_called_cells += 1;
                }
                *counted_per_cell.entry(cbs.to_string()).or_default() += 1;
                for gene in gx.split(';').filter(|gene| !gene.is_empty()) {
                    detected_genes.insert(gene.to_string());
                }
            }
            // Consistent hash-subsample keyed on (CB,UMI) only: the SAME molecules are
            // selected in every BAM (dumps join), and every read of a sampled molecule
            // is accumulated (so per-molecule evidence and reads/molecule are unbiased).
            if let Some((pcb, pub_)) = pack(cbs.as_bytes()).zip(pack(ubs.as_bytes())) {
                let key = (pcb << 24) | pub_; // lossless: CB in bits 24.., UMI in 0..23
                if cli.molecule_sample <= 1 || splitmix64(key) % cli.molecule_sample == 0 {
                    a.sampled_reads_with_cb_umi += 1;
                    let m = sampled.entry(key).or_default();
                    m.n_reads += 1;
                    m.min_nh = if m.n_reads == 1 {
                        nh as u32
                    } else {
                        m.min_nh.min(nh as u32)
                    };
                    m.max_mapq = m.max_mapq.max(rec.mapq());
                    m.max_readlen = m.max_readlen.max(rl as u32);
                    m.eligible |= xf & 16 != 0;
                    match re {
                        b'E' => m.region |= 4,
                        b'N' => m.region |= 2,
                        b'I' => m.region |= 1,
                        _ => {}
                    }
                    if xf & 8 != 0 {
                        m.counted = true;
                        if m.rep_gene.is_empty() && !gx.is_empty() {
                            m.rep_gene = gx.to_string();
                        }
                    }
                    if !gx.is_empty() {
                        let gene_key =
                            stable_hash(&[cbs.as_bytes(), ubs.as_bytes(), gx.as_bytes()]);
                        sampled_gene_candidates.insert(gene_key);
                    }
                }
            }
        }
    }

    // Dump one line per sampled molecule (all its reads aggregated). Columns:
    //   CB  UMI  counted  rep_gene  n_reads  min_nh  max_mapq  region  max_readlen
    // region: E>N>I>. (best hit). rep_gene is '-' when the molecule was not counted.
    if let Some(w) = mol_out.as_mut() {
        for (key, m) in &sampled {
            let cbs = unpack(key >> 24, 16);
            let ubs = unpack(key & 0xFF_FFFF, 12);
            let region = if m.region & 4 != 0 {
                'E'
            } else if m.region & 2 != 0 {
                'N'
            } else if m.region & 1 != 0 {
                'I'
            } else {
                '.'
            };
            let gene = if m.rep_gene.is_empty() {
                "-"
            } else {
                m.rep_gene.as_str()
            };
            writeln!(
                w,
                "{cbs}\t{ubs}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                m.counted as u8, gene, m.n_reads, m.min_nh, m.max_mapq, region, m.max_readlen
            )?;
        }
        a.molecules_emitted = sampled.len() as u64;
    }

    if let Some(w) = mol_out {
        w.finish()?.flush()?;
    }
    a.sampled_molecules = sampled.len() as u64;
    a.raw_umi_nodes_sampled = sampled_raw_umi_nodes.len() as u64;
    a.raw_umi_nodes_est = a.raw_umi_nodes_sampled * a.molecule_sample;
    a.molecule_funnel.reads_with_cb_umi = a.reads_with_cb_umi;
    a.molecule_funnel.raw_umi_nodes_est = a.raw_umi_nodes_est;
    a.molecule_funnel.corrected_cb_umi_nodes_est = a.sampled_molecules * a.molecule_sample;
    a.molecule_funnel.gene_assigned_candidate_nodes_est =
        sampled_gene_candidates.len() as u64 * a.molecule_sample;
    a.molecule_funnel.eligible_nodes_est =
        sampled.values().filter(|m| m.eligible).count() as u64 * a.molecule_sample;
    a.molecule_funnel.representatives = a.molecules;
    a.molecule_funnel.representatives_in_called_cells = called_cells
        .as_ref()
        .map(|_| representatives_in_called_cells);
    for m in sampled.values() {
        let bucket = match m.n_reads {
            1 => 0,
            2 => 1,
            3 => 2,
            4 => 3,
            5..=9 => 4,
            10..=19 => 5,
            _ => 6,
        };
        a.read_support_histogram[bucket] += 1;
    }
    a.saturation_curve = (1..=10)
        .map(|step| {
            let fraction = step as f64 / 10.0;
            let expected = sampled
                .values()
                .map(|m| 1.0 - (1.0 - fraction).powi(m.n_reads as i32))
                .sum::<f64>();
            SaturationPoint {
                read_fraction: fraction,
                estimated_cb_umi_nodes: (expected * a.molecule_sample as f64).round() as u64,
            }
        })
        .collect();
    let mut ambiguous_genes: Vec<_> = ambiguous_genes.into_iter().collect();
    ambiguous_genes.sort_unstable_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    a.mapping_diagnostics.top_ambiguous_genes = ambiguous_genes
        .into_iter()
        .take(25)
        .map(|(name, reads)| NamedCount { name, reads })
        .collect();
    a.genes_detected = detected_genes.len() as u64;
    a.cells_detected = counted_per_cell.len() as u64;
    let mut cell_counts: Vec<u64> = counted_per_cell.into_values().collect();
    cell_counts.sort_unstable();
    if !cell_counts.is_empty() {
        let quantile = |fraction: f64| -> f64 {
            let index = ((cell_counts.len() - 1) as f64 * fraction).round() as usize;
            cell_counts[index] as f64
        };
        a.molecules_per_cell = MoleculesPerCell {
            mean: a.molecules as f64 / cell_counts.len() as f64,
            median: quantile(0.5),
            p10: quantile(0.1),
            p90: quantile(0.9),
            max: *cell_counts.last().unwrap_or(&0),
        };
    }
    a.distinct_molecules_est = a.sampled_molecules * a.molecule_sample;
    a.reads_per_molecule = if a.sampled_molecules > 0 {
        a.sampled_reads_with_cb_umi as f64 / a.sampled_molecules as f64
    } else {
        0.0
    };
    a.elapsed_secs = t0.elapsed().as_secs_f64();
    a.records_per_sec = if a.elapsed_secs > 0.0 {
        (a.records as f64 / a.elapsed_secs) as u64
    } else {
        0
    };
    a.bins.retain(|b| b.reads > 0);
    a.aligned_length_bins.retain(|b| b.reads > 0);

    let json = serde_json::to_string_pretty(&a)?;
    std::fs::write(&cli.output, json)
        .with_context(|| format!("cannot write {}", cli.output.display()))?;

    let mapped = a.mapped;
    let pct = |x: u64, d: u64| {
        if d > 0 {
            100.0 * x as f64 / d as f64
        } else {
            0.0
        }
    };
    // Contamination is expressed as % of ALL mapped reads, not % of clips, so the
    // numbers are directly comparable across platforms.
    eprintln!(
        "RUN SUMMARY  {}\n  records {}  primary {}  unmapped {} ({:.2}%)\n  \
         unique(NH=1) {:.2}%  MAPQ255 {:.2}%  geneGX {:.2}%  multi-gene {:.3}%\n  \
         exonic {:.2}%  intronic {:.2}%  intergenic {:.2}%\n  \
         softclip {:.2}% of read bases   heavily-clipped {:.2}%   nM/100bp {:.3}\n  \
         contamination (% of mapped reads):  adapter {:.3}%  TSO {:.3}%  polyA(3') {:.3}%  polyT(5') {:.3}%\n  \
         counted molecules {}  genes {}  cells {}  median molecules/cell {:.1}\n  \
         reads/molecule {:.2}  UMIs/read {:.4}  UMI-corrected records {}\n  \
         subsample 1/{}  distinct molecules ~{}  dumped {}\n  \
         {:.1}s  {} rec/s  -> {}",
        a.label, a.records, a.primary, a.unmapped, pct(a.unmapped, a.primary),
        pct(a.unique_nh1, mapped), pct(a.mapq255, mapped), pct(a.gene_called, mapped), pct(a.multi_gene, mapped),
        pct(a.exonic, mapped), pct(a.intronic, mapped), pct(a.intergenic, mapped),
        pct(a.softclip_bases, a.read_bases), pct(a.heavily_clipped, mapped),
        if a.aligned_bases > 0 { 100.0 * a.mismatches as f64 / a.aligned_bases as f64 } else { 0.0 },
        pct(a.clip5.adapter + a.clip3.adapter, mapped),
        pct(a.clip5.tso + a.clip3.tso, mapped),
        pct(a.clip3.poly_a, mapped),
        pct(a.clip5.poly_t, mapped),
        a.molecules, a.genes_detected, a.cells_detected, a.molecules_per_cell.median,
        a.reads_per_molecule,
        if mapped > 0 { a.molecules as f64 / mapped as f64 } else { 0.0 },
        a.umi_corrected_records,
        a.molecule_sample, a.distinct_molecules_est, a.molecules_emitted,
        a.elapsed_secs, a.records_per_sec, cli.output.display()
    );
    Ok(())
}
