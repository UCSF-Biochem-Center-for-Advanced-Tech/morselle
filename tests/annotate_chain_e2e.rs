//! Real Illumina conversion -> genomic SAM -> annotate -> molecules, checking
//! the per-barcode tables of every stage against each other and their pools.
//!
//! Tests cannot depend on an installed minimap2, so the SAM is written the way
//! `minimap2 -a -y` writes it: one primary record per read (`--secondary=no`),
//! FASTQ comment tags copied verbatim, unplaced reads as unmapped records.
mod common;

use common::{barcode_table, bin, run, Rng};
use std::collections::BTreeMap;
use std::io::Read as _;

type Table = BTreeMap<(String, String), BTreeMap<String, u64>>;

fn binary(name: &str) -> std::path::PathBuf {
    match name {
        "sbxpress-illumina" => bin(name, env!("CARGO_BIN_EXE_sbxpress-illumina")),
        "sbx10x-annotate" => bin(name, env!("CARGO_BIN_EXE_sbx10x-annotate")),
        "sbx10x-molecules" => bin(name, env!("CARGO_BIN_EXE_sbx10x-molecules")),
        _ => unreachable!(),
    }
}

fn json(path: &std::path::Path) -> serde_json::Value {
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

fn assert_sums_to_pool(table: &Table, pool: &serde_json::Value, columns: &[(&str, &str)]) {
    for (column, field) in columns {
        let sum: u64 = table.values().map(|row| row[*column]).sum();
        assert_eq!(sum, pool[*field].as_u64().unwrap(), "{column} vs {field}");
    }
}

/// Minimap2 `-a -y` output for the tagged reads: exact placements on `genome`.
fn genomic_sam(tagged: &std::path::Path, genome: &str) -> String {
    let mut fastq = String::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(tagged).unwrap())
        .read_to_string(&mut fastq)
        .unwrap();
    let lines: Vec<&str> = fastq.lines().collect();
    let mut placed = Vec::new();
    let mut unmapped = Vec::new();
    for record in lines.chunks(4) {
        let mut header = record[0].trim_start_matches('@').split('\t');
        let name = header.next().unwrap().split(' ').next().unwrap();
        let tags = header.collect::<Vec<_>>().join("\t");
        let (seq, qual) = (record[1], record[3]);
        match genome.find(seq) {
            Some(offset) => placed.push((
                offset,
                format!(
                    "{name}\t0\tchr1\t{}\t60\t{}M\t*\t0\t0\t{seq}\t{qual}\t{tags}",
                    offset + 1,
                    seq.len()
                ),
            )),
            None => unmapped.push(format!(
                "{name}\t4\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{qual}\t{tags}"
            )),
        }
    }
    placed.sort();
    let mut sam = format!(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:{}\n",
        genome.len()
    );
    for line in placed.into_iter().map(|(_, line)| line).chain(unmapped) {
        sam.push_str(&line);
        sam.push('\n');
    }
    sam
}

#[test]
fn illumina_annotate_molecules_barcode_funnel() {
    let dir = tempfile::tempdir().unwrap();
    let root = dir.path();
    let genome = Rng::new(61).dna(1200);
    let (cell_a, cell_b) = ("AAACCCAAGAAACACT", "TTTGTTGAGTCGTTAC");
    let whitelist = root.join("whitelist.txt");
    std::fs::write(&whitelist, format!("{cell_a}\n{cell_b}\n")).unwrap();
    // Exon 101-700 (+). Reads: a PCR duplicate pair and a second molecule in
    // cell A; an intergenic and an unplaceable read in cell B; one read whose
    // barcode is off the whitelist.
    let reads = [
        ("a1", cell_a, "ACGTACGTACGT", genome[150..210].to_owned()),
        ("a2", cell_a, "ACGTACGTACGT", genome[150..210].to_owned()),
        ("a3", cell_a, "TTTTGGGGCCCC", genome[300..360].to_owned()),
        ("b1", cell_b, "ACGTACGTACGT", genome[900..960].to_owned()),
        ("b2", cell_b, "CATGCATGCATG", Rng::new(5).dna(60)),
        (
            "rejected",
            "GGGGGGGGGGGGGGGG",
            "ACGTACGTACGT",
            genome[400..460].to_owned(),
        ),
    ];
    let (r1, r2) = (root.join("R1.fastq"), root.join("R2.fastq"));
    let (mut r1_text, mut r2_text) = (String::new(), String::new());
    for (name, cb, umi, cdna) in &reads {
        r1_text.push_str(&format!("@{name}/1\n{cb}{umi}\n+\n{}\n", "I".repeat(28)));
        r2_text.push_str(&format!(
            "@{name}/2\n{cdna}\n+\n{}\n",
            "I".repeat(cdna.len())
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
            "--barcode-metrics",
        ],
    );
    let mapped = root.join("mapped.sam");
    std::fs::write(&mapped, genomic_sam(&tagged, &genome)).unwrap();
    let gtf = root.join("genes.gtf");
    std::fs::write(
        &gtf,
        "chr1\ttest\texon\t101\t700\t.\t+\t.\tgene_id \"ENSG000001.1\"; transcript_id \"ENST000001.1\"; gene_name \"GENE1\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\";\n",
    )
    .unwrap();
    let (annotated, annotation_metrics) =
        (root.join("annotated.bam"), root.join("annotation.json"));
    run(
        &binary("sbx10x-annotate"),
        &[
            mapped.to_str().unwrap(),
            "--gtf",
            gtf.to_str().unwrap(),
            "--output",
            annotated.to_str().unwrap(),
            "--metrics",
            annotation_metrics.to_str().unwrap(),
            "--three-prime",
            "--barcode-metrics",
        ],
    );
    let (molecules, molecule_metrics) = (root.join("molecules.bam"), root.join("molecules.json"));
    run(
        &binary("sbx10x-molecules"),
        &[
            annotated.to_str().unwrap(),
            "--output",
            molecules.to_str().unwrap(),
            "--metrics",
            molecule_metrics.to_str().unwrap(),
            "--barcode-metrics",
        ],
    );

    let conversion = barcode_table(&root.join("conversion.barcodes.tsv.gz"));
    let annotation = barcode_table(&root.join("annotation.barcodes.tsv.gz"));
    let molecule = barcode_table(&root.join("molecules.barcodes.tsv.gz"));
    assert_sums_to_pool(
        &conversion,
        &json(&conversion_metrics),
        &[("reads", "records_read"), ("emitted", "records_emitted")],
    );
    let columns = [
        "records",
        "primary",
        "unmapped",
        "exonic",
        "intronic",
        "intergenic",
        "antisense",
        "gene_assigned",
        "unique_gene",
        "multi_gene",
    ];
    assert_sums_to_pool(
        &annotation,
        &json(&annotation_metrics),
        &columns.map(|c| (c, c)),
    );
    assert_sums_to_pool(
        &molecule,
        &json(&molecule_metrics),
        &[
            ("eligible_reads", "eligible_records"),
            ("molecules_counted", "molecules_counted"),
        ],
    );

    let a = (format!("{cell_a}-1"), ".".to_owned());
    let b = (format!("{cell_b}-1"), ".".to_owned());
    assert_eq!(conversion[&a]["emitted"], 3);
    assert_eq!(conversion[&b]["emitted"], 2);
    assert_eq!(
        conversion[&("__unresolved_barcode__".to_owned(), ".".to_owned())]["failed_cb_unmatched"],
        1
    );
    assert_eq!(annotation[&a]["gene_assigned"], 3);
    assert_eq!(annotation[&b]["unmapped"], 1);
    assert_eq!(annotation[&b]["gene_assigned"], 0);
    assert_eq!(molecule[&a]["eligible_reads"], 3);
    assert_eq!(molecule[&a]["molecules_counted"], 2);
    assert!(molecule
        .get(&b)
        .is_none_or(|row| row["eligible_reads"] == 0));

    // Default --secondary=no route: each stage can only lose a cell's reads.
    for (cell, ann) in &annotation {
        let emitted = conversion[cell]["emitted"];
        let (eligible, counted) = molecule.get(cell).map_or((0, 0), |row| {
            (row["eligible_reads"], row["molecules_counted"])
        });
        assert_eq!(ann["records"], emitted, "{cell:?}");
        assert!(
            ann["records"] >= ann["primary"] && ann["primary"] >= eligible && eligible >= counted,
            "{cell:?}"
        );
    }
}
