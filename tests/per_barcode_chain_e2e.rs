//! Real Illumina conversion -> Flex probe alignment -> molecules -> matrix.
mod common;

use common::{barcode_table, bin, read_mex, run};
use std::io::Read as _;

fn binary(name: &str) -> std::path::PathBuf {
    match name {
        "sbxpress-illumina" => bin(name, env!("CARGO_BIN_EXE_sbxpress-illumina")),
        "sbxpress" => bin(name, env!("CARGO_BIN_EXE_sbxpress")),
        "sbxpress-probeset-align" => bin(name, env!("CARGO_BIN_EXE_sbxpress-probeset-align")),
        "sbx10x-molecules" => bin(name, env!("CARGO_BIN_EXE_sbx10x-molecules")),
        "cr-bam-to-mtx" => bin(name, env!("CARGO_BIN_EXE_cr-bam-to-mtx")),
        _ => unreachable!(),
    }
}

#[test]
fn sbx_extract_probe_molecules_matrix_barcode_funnel() {
    let dir = tempfile::tempdir().unwrap();
    let root = dir.path();
    let source =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/smoke.fastq.gz");
    let extracted = root.join("extracted");
    run(
        &binary("sbxpress"),
        &[
            "extract",
            "--fastqs",
            source.to_str().unwrap(),
            "--out-dir",
            extracted.to_str().unwrap(),
            "--threads",
            "1",
            "--minimap2",
            "--barcode-metrics",
        ],
    );
    let r2 = extracted.join("sbxpress_S1_L001_R2_001.fastq.gz");
    let mut fastq = String::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(&r2).unwrap())
        .read_to_string(&mut fastq)
        .unwrap();
    let mut lines = fastq.lines();
    let header = lines.next().unwrap();
    let sequence = lines.next().unwrap();
    let cb = header
        .split('\t')
        .find_map(|field| field.strip_prefix("CB:Z:"))
        .unwrap();
    let probe = &sequence[..50];
    let probe_csv = root.join("probe_set.csv");
    std::fs::write(&probe_csv, format!("#probe_set_file_format=3.0\n#panel_name=test\ngene_id,probe_seq,probe_id,included,region,gene_name\nENSG00000000001,{probe},ENSG00000000001|GENE1|p1,TRUE,spliced,GENE1\n")).unwrap();
    let aligned = root.join("aligned.bam");
    run(
        &binary("sbxpress-probeset-align"),
        &[
            r2.to_str().unwrap(),
            "--probe-set",
            probe_csv.to_str().unwrap(),
            "--output",
            aligned.to_str().unwrap(),
        ],
    );
    let molecules = root.join("molecules.bam");
    run(
        &binary("sbx10x-molecules"),
        &[
            aligned.to_str().unwrap(),
            "--output",
            molecules.to_str().unwrap(),
            "--barcode-metrics",
        ],
    );
    let extract_rows = barcode_table(&extracted.join("metrics_summary.barcodes.tsv.gz"));
    let molecule_rows = barcode_table(&root.join("molecules.bam.metrics.barcodes.tsv.gz"));
    let key = (cb.to_owned(), ".".to_owned());
    let e = &extract_rows[&key];
    let m = &molecule_rows[&key];
    assert!(
        e["reads"] >= e["emitted"]
            && e["emitted"] >= m["eligible_reads"]
            && m["eligible_reads"] >= m["molecules_counted"]
    );
    assert!(m["molecules_counted"] > 0);
    let outs = root.join("outs");
    run(
        &binary("cr-bam-to-mtx"),
        &[
            molecules.to_str().unwrap(),
            "--output",
            outs.to_str().unwrap(),
            "--no-gzip",
        ],
    );
    let (counts, _) = read_mex(&outs.join("raw_feature_bc_matrix"));
    let matrix_cell: f64 = counts
        .iter()
        .filter(|((barcode, _), _)| barcode == cb)
        .map(|(_, count)| count)
        .sum();
    assert_eq!(matrix_cell, m["molecules_counted"] as f64);
}

#[test]
fn illumina_probe_molecules_matrix_barcode_funnel() {
    let dir = tempfile::tempdir().unwrap();
    let root = dir.path();
    let barcode = "AAACCCAAGAAACACT";
    let umi = "ACGTACGTACGT";
    let probe = "AAAACCCCGGGGTTTTACGTACGTAAAACCCCGGGGTTTTACGTACGT01";
    let whitelist = root.join("whitelist.txt");
    std::fs::write(&whitelist, format!("{barcode}\n")).unwrap();
    let r1 = root.join("R1.fastq");
    let r2 = root.join("R2.fastq");
    let mut r1_text = String::new();
    let mut r2_text = String::new();
    for (name, cb) in [
        ("first", barcode),
        ("duplicate", barcode),
        ("rejected", "GGGGGGGGGGGGGGGG"),
    ] {
        let r1_seq = format!("{cb}{umi}");
        r1_text.push_str(&format!(
            "@{name}/1\n{r1_seq}\n+\n{}\n",
            "I".repeat(r1_seq.len())
        ));
        r2_text.push_str(&format!(
            "@{name}/2\n{probe}\n+\n{}\n",
            "I".repeat(probe.len())
        ));
    }
    std::fs::write(&r1, r1_text).unwrap();
    std::fs::write(&r2, r2_text).unwrap();
    let tagged = root.join("tagged.fastq.gz");
    let conversion_metrics = root.join("conversion.json");
    run(
        &binary("sbxpress-illumina"),
        &[
            "--r1",
            r1.to_str().unwrap(),
            "--r2",
            r2.to_str().unwrap(),
            "--whitelist",
            whitelist.to_str().unwrap(),
            "--output",
            tagged.to_str().unwrap(),
            "--metrics",
            conversion_metrics.to_str().unwrap(),
            "--gem-group",
            "2",
            "--barcode-metrics",
        ],
    );
    let conversion = barcode_table(&root.join("conversion.barcodes.tsv.gz"));
    let cell = (format!("{barcode}-2"), ".".to_owned());
    assert_eq!(conversion[&cell]["reads"], 2);
    assert_eq!(conversion[&cell]["emitted"], 2);
    assert_eq!(
        conversion[&("__unresolved_barcode__".to_owned(), ".".to_owned())]["failed_cb_unmatched"],
        1
    );

    let probe_csv = root.join("probe_set.csv");
    std::fs::write(&probe_csv, format!("#probe_set_file_format=3.0\n#panel_name=test\ngene_id,probe_seq,probe_id,included,region,gene_name\nENSG00000000001,{probe},ENSG00000000001|GENE1|p1,TRUE,spliced,GENE1\n")).unwrap();
    let aligned = root.join("aligned.bam");
    run(
        &binary("sbxpress-probeset-align"),
        &[
            tagged.to_str().unwrap(),
            "--probe-set",
            probe_csv.to_str().unwrap(),
            "--output",
            aligned.to_str().unwrap(),
        ],
    );
    let molecules = root.join("molecules.bam");
    let molecule_metrics = root.join("molecules.json");
    run(
        &binary("sbx10x-molecules"),
        &[
            aligned.to_str().unwrap(),
            "--output",
            molecules.to_str().unwrap(),
            "--metrics",
            molecule_metrics.to_str().unwrap(),
            "--barcode-metrics",
        ],
    );
    let table = barcode_table(&root.join("molecules.barcodes.tsv.gz"));
    assert_eq!(table[&cell]["eligible_reads"], 2);
    assert_eq!(table[&cell]["raw_cell_gene_umis"], 1);
    assert_eq!(table[&cell]["molecules_counted"], 1);
    let outs = root.join("outs");
    run(
        &binary("cr-bam-to-mtx"),
        &[
            molecules.to_str().unwrap(),
            "--output",
            outs.to_str().unwrap(),
            "--no-gzip",
        ],
    );
    let (counts, _) = read_mex(&outs.join("raw_feature_bc_matrix"));
    let matrix_cell: f64 = counts
        .iter()
        .filter(|((cb, _), _)| cb == &cell.0)
        .map(|(_, n)| n)
        .sum();
    assert_eq!(matrix_cell, 1.0);
}
