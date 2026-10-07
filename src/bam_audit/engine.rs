use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Instant;

use anyhow::{Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;
use rust_htslib::bam::{self, record::Cigar, Read};
use rustc_hash::{FxHashMap, FxHashSet};

use super::helpers::{
    aux_char, aux_int, aux_str, bump, classify, kmers, load_called_cells, nh_bucket, pack, revcomp,
    splitmix64, stable_hash, unpack,
};
use super::types::{
    Audit, Cli, LenBin, MappingLenBin, MolAccum, MoleculesPerCell, NamedCount, SaturationPoint,
    ADAPTER, BIN, MIN_CLIP, N_BINS, TSO,
};

pub fn run(cli: &Cli) -> Result<()> {
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
                let key = (pcb << 24) | pub_;  // lossless: CB in bits 24.., UMI in 0..23
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
