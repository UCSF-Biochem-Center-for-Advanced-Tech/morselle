//! Low-complexity cDNA (`YL:i:1`, crate::complexity) through its real
//! neighbours: conversion tags it, probe alignment and `annotate` carry it,
//! and `molecules` reports the gene in `artifacts` without dropping reads.
mod common;

use common::{barcode_table, bin, read_bam, run};
use rust_htslib::bam::record::Aux;
use std::path::Path;

fn binary(name: &str) -> std::path::PathBuf {
    match name {
        "morselle-illumina" => bin(name, env!("CARGO_BIN_EXE_morselle-illumina")),
        "morselle-probeset-align" => bin(name, env!("CARGO_BIN_EXE_morselle-probeset-align")),
        "morselle-annotate" => bin(name, env!("CARGO_BIN_EXE_morselle-annotate")),
        "morselle-molecules" => bin(name, env!("CARGO_BIN_EXE_morselle-molecules")),
        _ => unreachable!(),
    }
}

fn json(path: &Path) -> serde_json::Value {
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

/// 12-base UMI number `i`, distinct for every i < 4^12.
fn umi(mut i: usize) -> String {
    (0..12)
        .map(|_| {
            let base = b"ACGT"[i % 4] as char;
            i /= 4;
            base
        })
        .collect()
}

fn molecules(input: &Path, root: &Path) -> serde_json::Value {
    let output = root.join("molecules.bam");
    let metrics = root.join("molecules.json");
    run(
        &binary("morselle-molecules"),
        &[
            input.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
            "--metrics",
            metrics.to_str().unwrap(),
        ],
    );
    json(&metrics)
}

#[test]
fn illumina_flags_low_complexity_and_molecules_names_the_gene() {
    let dir = tempfile::tempdir().unwrap();
    let root = dir.path();
    let barcode = "AAACCCAAGAAACACT";
    // A G-rich probe (90% G, so its reads are low-complexity) and an ordinary one.
    let poly_g = format!("{}ACGTA", "G".repeat(45));
    let normal = "GATTACACCATGGCTTGAACGTTAGCCATGCAAGTCCGATACGGATTTAC";
    let whitelist = root.join("whitelist.txt");
    std::fs::write(&whitelist, format!("{barcode}\n")).unwrap();
    let (mut r1, mut r2) = (String::new(), String::new());
    for i in 0..130 {
        let seq = if i < 110 { poly_g.as_str() } else { normal };
        let r1_seq = format!("{barcode}{}", umi(i));
        r1.push_str(&format!("@r{i}/1\n{r1_seq}\n+\n{}\n", "I".repeat(28)));
        r2.push_str(&format!("@r{i}/2\n{seq}\n+\n{}\n", "I".repeat(seq.len())));
    }
    std::fs::write(root.join("R1.fastq"), r1).unwrap();
    std::fs::write(root.join("R2.fastq"), r2).unwrap();
    let tagged = root.join("tagged.fastq.gz");
    let conversion = root.join("conversion.json");
    run(
        &binary("morselle-illumina"),
        &[
            "--r1",
            root.join("R1.fastq").to_str().unwrap(),
            "--r2",
            root.join("R2.fastq").to_str().unwrap(),
            "--whitelist",
            whitelist.to_str().unwrap(),
            "--output",
            tagged.to_str().unwrap(),
            "--metrics",
            conversion.to_str().unwrap(),
            "--barcode-metrics",
            // Keep and tag (the default drops them) so molecules can name the gene.
            "--drop-low-complexity",
            "false",
        ],
    );
    let metrics = json(&conversion);
    assert_eq!(metrics["records_emitted"], 130);
    assert_eq!(metrics["low_complexity"], 110);
    let table = barcode_table(&root.join("conversion.barcodes.tsv.gz"));
    assert_eq!(
        table[&(format!("{barcode}-1"), ".".to_owned())]["low_complexity"],
        110
    );

    let probe_csv = root.join("probe_set.csv");
    std::fs::write(&probe_csv, format!("#probe_set_file_format=3.0\n#panel_name=test\ngene_id,probe_seq,probe_id,included,region,gene_name\nENSG00000000001,{poly_g},ENSG00000000001|GRICH|p1,TRUE,spliced,GRICH\nENSG00000000002,{normal},ENSG00000000002|NORMAL|p1,TRUE,spliced,NORMAL\n")).unwrap();
    let aligned = root.join("aligned.bam");
    run(
        &binary("morselle-probeset-align"),
        &[
            tagged.to_str().unwrap(),
            "--probe-set",
            probe_csv.to_str().unwrap(),
            "--output",
            aligned.to_str().unwrap(),
        ],
    );
    let records = read_bam(&aligned);
    let tagged_records = records
        .iter()
        .filter(|r| matches!(r.aux(b"YL"), Ok(Aux::I32(1))))
        .count();
    assert_eq!(tagged_records, 110);

    let metrics = molecules(&aligned, root);
    // Flagged, not dropped: every read is still eligible.
    assert_eq!(metrics["eligible_records"], 130);
    let artifacts = &metrics["artifacts"];
    assert_eq!(artifacts["gene_assigned_reads"], 130);
    assert_eq!(artifacts["low_complexity_reads"], 110);
    assert_eq!(artifacts["low_mapq_reads"], 0);
    let genes = artifacts["genes"].as_array().unwrap();
    assert_eq!(genes.len(), 1, "{genes:?}");
    assert_eq!(genes[0]["gene_name"], "GRICH");
    assert_eq!(genes[0]["low_complexity_reads"], 110);
    assert_eq!(genes[0]["flagged_fraction"], 1.0);
}

/// Gene A (+): exon 101-200, intron to 500.
const GTF: &str = "chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id \"GENEA\"; transcript_id \"TXA\"; gene_name \"A\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\";
chr1\ttest\texon\t501\t600\t.\t+\t.\tgene_id \"GENEA\"; transcript_id \"TXA\"; gene_name \"A\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\";
";

#[test]
fn annotate_keeps_yl_and_low_mapq_intronic_gene_is_listed() {
    let dir = tempfile::tempdir().unwrap();
    let root = dir.path();
    let gtf = root.join("genes.gtf");
    std::fs::write(&gtf, GTF).unwrap();
    // 120 low-MAPQ intronic poly-G reads in gene A (the DZEC05 LINC00486
    // pattern), and 10 confident exonic reads.
    let mut sam = String::from("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n");
    for i in 0..130 {
        let (pos, mapq, seq, yl) = if i < 120 {
            (301, 3, "G".repeat(50), "\tYL:i:1")
        } else {
            (111, 60, "ACGT".repeat(12) + "AC", "")
        };
        let u = umi(i);
        sam.push_str(&format!(
            "r{i}\t0\tchr1\t{pos}\t{mapq}\t50M\t*\t0\t0\t{seq}\t{}\tCB:Z:AAACCCAAGAAACACT-1\tUR:Z:{u}\tUB:Z:{u}{yl}\n",
            "I".repeat(50)
        ));
    }
    let mapped = root.join("mapped.sam");
    std::fs::write(&mapped, sam).unwrap();
    let annotated = root.join("annotated.bam");
    run(
        &binary("morselle-annotate"),
        &[
            mapped.to_str().unwrap(),
            "--gtf",
            gtf.to_str().unwrap(),
            "--output",
            annotated.to_str().unwrap(),
            "--include-introns",
            "true",
        ],
    );
    let records = read_bam(&annotated);
    let r0 = records.iter().find(|r| r.qname() == b"r0").unwrap();
    assert!(matches!(
        r0.aux(b"YL"),
        Ok(Aux::I32(1)) | Ok(Aux::U8(1)) | Ok(Aux::I8(1))
    ));

    let metrics = molecules(&annotated, root);
    // The default --min-mapq 60 drops the intronic reads; the report still
    // attributes them to the gene, whatever --min-mapq was.
    assert_eq!(metrics["eligible_records"], 10);
    assert_eq!(metrics["low_mapq_records"], 120);
    let artifacts = &metrics["artifacts"];
    assert_eq!(artifacts["gene_assigned_reads"], 130);
    assert_eq!(artifacts["low_complexity_reads"], 120);
    assert_eq!(artifacts["low_mapq_reads"], 120);
    let genes = artifacts["genes"].as_array().unwrap();
    assert_eq!(genes.len(), 1);
    assert_eq!(genes[0]["gene_id"], "GENEA");
    assert_eq!(genes[0]["reads"], 130);
}
