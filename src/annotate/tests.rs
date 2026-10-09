use rust_htslib::bam::{
    self,
    header::HeaderRecord,
    record::{Aux, Cigar, CigarString},
    Read,
};
use std::fs;
use std::io::Read as IoRead;

use super::call::{compatible_position, compatible_position_with, transcript_cigar};
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

fn iv(start: i64, end: i64) -> Interval {
    Interval { start, end }
}

fn two_exon_tx(strand: Strand) -> Transcript {
    Transcript {
        id: "tx".into(),
        gene: 0,
        contig: 0,
        strand,
        exons: vec![iv(100, 150), iv(200, 260)],
    }
}

#[test]
fn read_end_overhang_is_soft_clipped_up_to_the_limit() {
    let plus = two_exon_tx(Strand::Plus);
    // Exact ends need no allowance.
    assert_eq!(
        compatible_position_with(&plus, &[iv(110, 140)], 3),
        Some((10, false))
    );
    // 2 bp before exon 1: compatible only with the allowance, position clamped.
    assert_eq!(compatible_position_with(&plus, &[iv(98, 140)], 0), None);
    assert_eq!(
        compatible_position_with(&plus, &[iv(98, 140)], 3),
        Some((0, true))
    );
    // 3 bp is the limit; 4 bp is not compatible.
    assert_eq!(
        compatible_position_with(&plus, &[iv(97, 140)], 3),
        Some((0, true))
    );
    assert_eq!(compatible_position_with(&plus, &[iv(96, 140)], 3), None);
    // Spliced read running 3 bp past the end of exon 2.
    assert_eq!(
        compatible_position_with(&plus, &[iv(120, 150), iv(200, 263)], 3),
        Some((20, true))
    );
    // Internal junctions stay exact: exon 1 end overshot by 1 bp.
    assert_eq!(
        compatible_position_with(&plus, &[iv(120, 151), iv(200, 230)], 3),
        None
    );
    // Minus strand: position counts from the clamped high end.
    let minus = two_exon_tx(Strand::Minus);
    assert_eq!(
        compatible_position_with(&minus, &[iv(120, 150), iv(200, 262)], 3),
        Some((0, true))
    );
    assert_eq!(
        compatible_position(&minus, &[iv(120, 150), iv(200, 262)]),
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
        record.unset_unmapped();
        if index < 127 {
            let cb = if index % 2 == 0 {
                "AAAAAAAAAAAAAAAA-1"
            } else {
                "CCCCCCCCCCCCCCCC-1"
            };
            record.push_aux(b"CB", Aux::String(cb)).unwrap();
        }
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
            confident_mapq: 60,
            max_end_overhang: 3,
            metrics: Some(metrics.clone()),
            barcode_metrics: true,
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
        let mut table = String::new();
        flate2::read::MultiGzDecoder::new(
            fs::File::open(crate::barcode_metrics::barcode_path(&metrics)).unwrap(),
        )
        .read_to_string(&mut table)
        .unwrap();
        let metrics: serde_json::Value =
            serde_json::from_slice(&fs::read(metrics).unwrap()).unwrap();
        (tags, metrics, table)
    };

    let (ordered_tags, ordered_metrics, ordered_table) = run_mode(false, "ordered");
    let (unordered_tags, unordered_metrics, unordered_table) = run_mode(true, "unordered");
    assert_eq!(ordered_tags, unordered_tags);
    assert_eq!(ordered_table, unordered_table);
    assert!(ordered_table.contains("__no_barcode__\t.\t1\t1\t0"));
    let mut lines = ordered_table.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let rows: Vec<Vec<&str>> = lines.map(|line| line.split('\t').collect()).collect();
    for field in [
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
    ] {
        let index = header.iter().position(|&column| column == field).unwrap();
        let sum: u64 = rows
            .iter()
            .map(|row| row[index].parse::<u64>().unwrap())
            .sum();
        assert_eq!(sum, ordered_metrics[field].as_u64().unwrap(), "{field}");
    }
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
