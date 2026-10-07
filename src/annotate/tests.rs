use std::fs;
use rust_htslib::bam::{self, header::HeaderRecord, record::{Aux, Cigar, CigarString}, Read};

use super::call::{compatible_position, transcript_cigar};
use super::pipeline::run;
use super::types::{Cli, Interval, Strand, Transcript};

#[test]
fn spliced_transcript_compatibility_and_position() {
    let tx = Transcript {
        id: "tx".into(),
        gene: 0,
        contig: 0,
        strand: Strand::Plus,
        exons: vec![
            Interval {
                start: 100,
                end: 150,
            },
            Interval {
                start: 200,
                end: 260,
            },
        ],
    };
    assert_eq!(
        compatible_position(
            &tx,
            &[
                Interval {
                    start: 120,
                    end: 150
                },
                Interval {
                    start: 200,
                    end: 230
                }
            ]
        ),
        Some(20)
    );
    assert_eq!(
        compatible_position(
            &tx,
            &[
                Interval {
                    start: 120,
                    end: 140
                },
                Interval {
                    start: 200,
                    end: 230
                }
            ]
        ),
        None
    );
}

#[test]
fn transcript_cigar_drops_introns_and_reverses_minus() {
    assert_eq!(
        transcript_cigar("5S20M100N30M2I4M", Strand::Plus),
        "5S50M2I4M"
    );
    assert_eq!(
        transcript_cigar("5S20M100N30M2I4M", Strand::Minus),
        "4M2I50M5S"
    );
}

#[test]
fn minus_position_uses_high_genomic_block() {
    let tx = Transcript {
        id: "tx".into(),
        gene: 0,
        contig: 0,
        strand: Strand::Minus,
        exons: vec![
            Interval {
                start: 100,
                end: 150,
            },
            Interval {
                start: 200,
                end: 260,
            },
        ],
    };
    assert_eq!(
        compatible_position(
            &tx,
            &[
                Interval {
                    start: 120,
                    end: 150
                },
                Interval {
                    start: 200,
                    end: 230
                },
            ],
        ),
        Some(30),
    );
}

#[test]
fn unordered_annotation_matches_ordered_tags_and_counts() {
    let temp = tempfile::tempdir().unwrap();
    let gtf = temp.path().join("genes.gtf");
    fs::write(
        &gtf,
        "chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id \"ENSG000001.1\"; transcript_id \"ENST000001.1\"; gene_name \"GENE1\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\";\n",
    )
    .unwrap();

    let input = temp.path().join("input.bam");
    let mut header = bam::Header::new();
    header.push_record(
        HeaderRecord::new(b"SQ")
            .push_tag(b"SN", "chr1")
            .push_tag(b"LN", 1000),
    );
    let mut writer = bam::Writer::from_path(&input, &header, bam::Format::Bam).unwrap();
    let cigar = CigarString(vec![Cigar::Match(20)]);
    for index in 0..128 {
        let mut record = bam::Record::new();
        let name = format!("read-{index:03}");
        record.set(
            name.as_bytes(),
            Some(&cigar),
            b"ACGTACGTACGTACGTACGT",
            &[30; 20],
        );
        record.set_tid(0);
        record.set_pos(if index % 2 == 0 { 110 } else { 300 });
        writer.write(&record).unwrap();
    }
    drop(writer);

    let run_mode = |unordered: bool, stem: &str| {
        let output = temp.path().join(format!("{stem}.bam"));
        let metrics = temp.path().join(format!("{stem}.json"));
        run(&Cli {
            bam: input.clone(),
            gtf: gtf.clone(),
            output: output.clone(),
            include_introns: true,
            threads: 4,
            unordered,
            batch_size: 17,
            multimap_tags: false,
            metrics: Some(metrics.clone()),
            genes_mtx: None,
            reference_filter: true,
            strip_id_versions: true,
            three_prime: true,
        })
        .unwrap();
        let mut reader = bam::Reader::from_path(output).unwrap();
        let mut tags = reader
            .records()
            .map(|record| {
                let record = record.unwrap();
                let gx = match record.aux(b"GX").ok() {
                    Some(Aux::String(value)) => value.to_owned(),
                    _ => String::new(),
                };
                let re = match record.aux(b"RE").ok() {
                    Some(Aux::Char(value)) => value,
                    _ => b'.',
                };
                (record.qname().to_vec(), gx, re)
            })
            .collect::<Vec<_>>();
        tags.sort_unstable();
        let metrics: serde_json::Value =
            serde_json::from_slice(&fs::read(metrics).unwrap()).unwrap();
        (tags, metrics)
    };

    let (ordered_tags, ordered_metrics) = run_mode(false, "ordered");
    let (unordered_tags, unordered_metrics) = run_mode(true, "unordered");
    assert_eq!(ordered_tags, unordered_tags);
    for field in [
        "records",
        "primary",
        "exonic",
        "intergenic",
        "gene_assigned",
        "unique_gene",
        "collision_records",
    ] {
        assert_eq!(ordered_metrics[field], unordered_metrics[field], "{field}");
    }
    assert_eq!(unordered_metrics["unordered"], true);
}
