//! End-to-end integration tests for Illumina + Flex and Illumina + CAT Tagseq modes.

use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
use tempfile::tempdir;

use rust_htslib::bam::{self, Read as _};
use sbx10x::cell_calling;
use sbx10x::probe_set::{self, Cli as ProbeSetCli};

fn write_file(path: &PathBuf, content: &str) {
    let mut f = File::create(path).expect("failed to create temp file");
    f.write_all(content.as_bytes())
        .expect("failed to write temp file");
}

#[test]
fn test_illumina_flex_probeset_align_integration() {
    let dir = tempdir().expect("failed to create temp dir");
    let probe_set_path = dir.path().join("probe_set.csv");
    let r2_fastq_path = dir.path().join("sample_R2.fastq");
    let output_bam_path = dir.path().join("output.bam");
    let genes_out_path = dir.path().join("genes.mtx");
    let metrics_path = dir.path().join("metrics.json");

    // 1. Write a 10x Flex probe set CSV fixture (50bp probes)
    let probe_seq_1 = "AAAACCCCGGGGTTTTACGTACGTAAAACCCCGGGGTTTTACGTACGT01";
    let probe_seq_2 = "TTTTGGGGCCCCAAAATGCATGCATTTTTGGGGCCCCAAAATGCATGC02";
    let probe_csv_content = format!(
        "\
#probe_set_file_format=3.0
#panel_name=test_flex
gene_id,probe_seq,probe_id,included,region,gene_name
ENSG00000000001,{probe_seq_1},ENSG00000000001|GENE1|p1,TRUE,spliced,GENE1
ENSG00000000002,{probe_seq_2},ENSG00000000002|GENE2|p2,TRUE,spliced,GENE2
"
    );
    write_file(&probe_set_path, &probe_csv_content);

    // 2. Write an Illumina-style R2 FASTQ carrying 10x tags in comments
    let fastq_content = format!(
        "\
@read_flex_1\tCR:Z:AAACCCAAGAAACACT\tCB:Z:AAACCCAAGAAACACT-1\tUR:Z:ACGTACGTACGT\tUB:Z:ACGTACGTACGT
{probe_seq_1}
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read_flex_2\tCR:Z:TTTGGGCCCAAAACCC\tCB:Z:TTTGGGCCCAAAACCC-1\tUR:Z:TGCATGCATGCA\tUB:Z:TGCATGCATGCA
{probe_seq_2}
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"
    );
    write_file(&r2_fastq_path, &fastq_content);

    // 3. Run sbxpress probeset-align engine
    let cli = ProbeSetCli {
        r2_fastq: r2_fastq_path,
        probe_set: probe_set_path,
        output: output_bam_path.clone(),
        genes_out: Some(genes_out_path.clone()),
        metrics: Some(metrics_path.clone()),
        max_mm: 1,
        window: 15,
        threads: 2,
    };
    probe_set::run(&cli).expect("probeset-align execution failed");

    // 4. Verify BAM outputs
    assert!(output_bam_path.exists(), "output BAM must exist");
    assert!(genes_out_path.exists(), "gene manifest must exist");
    assert!(metrics_path.exists(), "metrics JSON must exist");

    let mut bam_reader =
        bam::Reader::from_path(&output_bam_path).expect("failed to open output BAM");
    let header = bam_reader.header();
    assert_eq!(header.target_count(), 2, "must contain 2 target genes");

    let mut record_count = 0;
    for result in bam_reader.records() {
        let record = result.expect("failed to read BAM record");
        record_count += 1;
        assert!(!record.is_unmapped(), "record must be mapped to probe");

        // Verify tags (GX, GN, CB, UB)
        let gx = record.aux(b"GX").expect("must carry GX gene ID tag");
        let cb = record.aux(b"CB").expect("must carry CB barcode tag");
        if let bam::record::Aux::String(gx_str) = gx {
            assert!(
                gx_str.starts_with("ENSG0000000000"),
                "GX tag must match gene ID"
            );
        } else {
            panic!("GX tag is not a string");
        }
        if let bam::record::Aux::String(cb_str) = cb {
            assert!(
                cb_str.contains("-1"),
                "CB tag must carry gem group suffix"
            );
        } else {
            panic!("CB tag is not a string");
        }
    }
    assert_eq!(record_count, 2, "both reads must be mapped");
}

#[test]
fn test_illumina_cat_tagseq_cell_calling_integration() {
    let dir = tempdir().expect("failed to create temp dir");
    let barcodes_path = dir.path().join("tagseq_barcodes.txt");

    // Write a whitelist file with 2 specific cell barcodes
    let barcodes_content = "AAACCCAAGAAACACT-1\nTTTGGGCCCAAAACCC-1\n";
    write_file(&barcodes_path, barcodes_content);

    let raw_barcodes = vec![
        "AAACCCAAGAAACACT-1".to_string(),
        "UNWANTED_BARCODE-1".to_string(),
        "TTTGGGCCCAAAACCC-1".to_string(),
    ];

    // Resolve tagseq barcodes against the raw matrix barcodes
    let resolved =
        cell_calling::resolve_tagseq_barcodes(&barcodes_path, &raw_barcodes).unwrap();

    assert_eq!(resolved.len(), 2, "must resolve exactly 2 barcodes");
    assert_eq!(resolved[0], 0, "first barcode maps to index 0");
    assert_eq!(resolved[1], 2, "second barcode maps to index 2");
}
