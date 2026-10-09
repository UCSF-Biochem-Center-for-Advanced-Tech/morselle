//! Layer B (Cell Ranger's counting rule) through its real neighbours:
//! minimap2-style SAM -> `annotate` (`mm:i:1` rescue, `NH` with
//! `--multimap-tags`) -> `molecules --min-mapq` -> `matrix`.
use super::*;

fn annotate() -> std::path::PathBuf {
    bin("sbx10x-annotate", env!("CARGO_BIN_EXE_sbx10x-annotate"))
}

const CELL: &str = "AAACCCAAGAAACACT-1";

/// Gene A (+): exons 101-200 and 301-400; gene B (+): exon 601-700.
const GTF: &str = "chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id \"GENEA\"; transcript_id \"TXA\"; gene_name \"A\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\";
chr1\ttest\texon\t301\t400\t.\t+\t.\tgene_id \"GENEA\"; transcript_id \"TXA\"; gene_name \"A\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\";
chr1\ttest\texon\t601\t700\t.\t+\t.\tgene_id \"GENEB\"; transcript_id \"TXB\"; gene_name \"B\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\";
";

/// (name, UMI, [(flag, 1-based pos, MAPQ)]): primary first, then `-N`
/// secondaries, as minimap2 writes a query group.
/// (flag, 1-based position, MAPQ) of one reported hit.
type Hit = (u16, u32, u8);

const READS: [(&str, &str, &[Hit]); 5] = [
    // Unique exonic hit: counted.
    ("confident", "AAAAAAAAAAAA", &[(0, 111, 60)]),
    // Low MAPQ; the only other hit is intergenic: rescued as gene A.
    ("rescued", "CCCCCCCCCCCC", &[(0, 121, 5), (256, 801, 0)]),
    // Low-MAPQ intronic hit in gene A: intronic hits are never rescued.
    ("intronic", "GGGGGGGGGGGG", &[(0, 221, 3)]),
    // Exonic in both A and B: two genes, not rescued.
    ("two_genes", "TTTTTTTTTTTT", &[(0, 131, 0), (256, 611, 0)]),
    // Confidently mapped despite a weaker secondary: counted.
    (
        "confident_nh2",
        "ACACACACACAC",
        &[(0, 141, 60), (256, 621, 0)],
    ),
];

fn sam(secondaries: bool) -> String {
    let mut text = String::from("@HD\tVN:1.6\tSO:unsorted\tGO:query\n@SQ\tSN:chr1\tLN:1000\n");
    for (name, umi, hits) in READS {
        for &(flag, pos, mapq) in hits.iter().filter(|h| secondaries || h.0 == 0) {
            let (seq, qual) = if flag & 256 == 0 {
                ("A".repeat(50), "I".repeat(50))
            } else {
                ("*".to_owned(), "*".to_owned())
            };
            text.push_str(&format!(
                "{name}\t{flag}\tchr1\t{pos}\t{mapq}\t50M\t*\t0\t0\t{seq}\t{qual}\tCB:Z:{CELL}\tUR:Z:{umi}\tUB:Z:{umi}\n"
            ));
        }
    }
    text
}

fn json(path: &Path) -> serde_json::Value {
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

/// Runs the chain; returns (annotated records, molecule metrics, matrix total).
fn chain(
    dir: &Path,
    name: &str,
    annotate_args: &[&str],
    molecule_args: &[&str],
) -> (Vec<bam::Record>, serde_json::Value, f64) {
    let gtf = dir.join("genes.gtf");
    std::fs::write(&gtf, GTF).unwrap();
    let secondaries = annotate_args.contains(&"--multimap-tags");
    let mapped = dir.join(format!("{name}.sam"));
    std::fs::write(&mapped, sam(secondaries)).unwrap();
    let annotated = dir.join(format!("{name}.annotated.bam"));
    let mut args = vec![
        mapped.to_str().unwrap(),
        "--gtf",
        gtf.to_str().unwrap(),
        "--output",
        annotated.to_str().unwrap(),
    ];
    args.extend(annotate_args);
    run(&annotate(), &args);
    let output = dir.join(format!("{name}.molecules.bam"));
    let mut flags = vec![
        annotated.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--barcode-metrics",
    ];
    flags.extend(molecule_args);
    run(&molecules(), &flags);
    let metrics_path = std::path::PathBuf::from(format!("{}.metrics.json", output.display()));
    let metrics = json(&metrics_path);
    // The per-barcode gate columns reconcile with the pool counters.
    let table = barcode_table(&sbx10x::barcode_metrics::barcode_path(&metrics_path));
    for (column, field) in [
        ("rescued_reads", "rescued_records"),
        ("low_mapq_reads", "low_mapq_records"),
        ("multimapper_reads", "multimapper_records"),
    ] {
        let sum: u64 = table.values().map(|row| row[column]).sum();
        assert_eq!(sum, metrics[field].as_u64().unwrap(), "{name}: {column}");
    }
    let outs = dir.join(format!("{name}.outs"));
    run(
        &matrix(),
        &[
            output.to_str().unwrap(),
            "-o",
            outs.to_str().unwrap(),
            "--no-gzip",
        ],
    );
    let (counts, _) = read_mex(&outs.join("raw_feature_bc_matrix"));
    let total = counts.values().sum();
    (read_bam(&annotated), metrics, total)
}

fn int_tag(records: &[bam::Record], name: &str, tag: &[u8]) -> Option<i64> {
    let record = records
        .iter()
        .find(|r| r.qname() == name.as_bytes() && !r.is_secondary())
        .unwrap();
    match record.aux(tag) {
        Ok(Aux::I32(v)) => Some(v.into()),
        Ok(Aux::U8(v)) => Some(v.into()),
        Ok(Aux::I8(v)) => Some(v.into()),
        _ => None,
    }
}

#[test]
fn grouped_route_counts_confident_and_rescued_reads_only() {
    let dir = tempfile::tempdir().unwrap();
    let (records, metrics, total) = chain(dir.path(), "grouped", &["--multimap-tags"], &[]);
    assert_eq!(int_tag(&records, "rescued", b"mm"), Some(1));
    assert_eq!(int_tag(&records, "rescued", b"NH"), Some(2));
    assert_eq!(int_tag(&records, "two_genes", b"mm"), None);
    assert_eq!(int_tag(&records, "two_genes", b"NH"), Some(2));
    assert_eq!(int_tag(&records, "intronic", b"mm"), None);
    assert_eq!(int_tag(&records, "confident_nh2", b"mm"), None);
    assert_eq!(metrics["min_mapq"], 60);
    assert_eq!(metrics["rescued_records"], 1);
    assert_eq!(metrics["low_mapq_records"], 1);
    assert_eq!(metrics["multimapper_records"], 1);
    assert_eq!(metrics["eligible_records"], 3);
    // confident, rescued, confident_nh2.
    assert_eq!(total, 3.0);
}

#[test]
fn default_route_rescues_single_gene_exonic_reads() {
    let dir = tempfile::tempdir().unwrap();
    let (records, metrics, total) = chain(dir.path(), "default", &[], &[]);
    // Without secondaries a low-MAPQ exonic read naming one gene is rescued,
    // as Cell Ranger does with one reported alignment; NH is never written.
    assert_eq!(int_tag(&records, "rescued", b"mm"), Some(1));
    assert_eq!(int_tag(&records, "two_genes", b"mm"), Some(1));
    assert_eq!(int_tag(&records, "intronic", b"mm"), None);
    assert_eq!(int_tag(&records, "rescued", b"NH"), None);
    assert_eq!(metrics["rescued_records"], 2);
    assert_eq!(metrics["low_mapq_records"], 1);
    assert_eq!(metrics["multimapper_records"], 0);
    assert_eq!(total, 4.0);
}

#[test]
fn min_mapq_zero_and_confident_mapq_zero_restore_unfiltered_counts() {
    let dir = tempfile::tempdir().unwrap();
    let (_, metrics, total) = chain(
        dir.path(),
        "unfiltered",
        &["--multimap-tags"],
        &["--min-mapq", "0"],
    );
    assert_eq!(metrics["min_mapq"], 0);
    assert_eq!(metrics["low_mapq_records"], 0);
    assert_eq!(metrics["multimapper_records"], 0);
    assert_eq!(total, 5.0);
    // annotate --confident-mapq 0: no rescue tag, so only MAPQ-60 reads count.
    let (records, metrics, total) = chain(dir.path(), "no_rescue", &["--confident-mapq", "0"], &[]);
    assert_eq!(int_tag(&records, "rescued", b"mm"), None);
    assert_eq!(metrics["rescued_records"], 0);
    assert_eq!(metrics["low_mapq_records"], 3);
    assert_eq!(total, 2.0);
}
