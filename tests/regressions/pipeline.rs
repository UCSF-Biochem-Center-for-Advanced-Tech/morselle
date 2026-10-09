//! Real fixture records carried through molecules, matrix, and cells.
use super::*;
use rust_htslib::bam::Read as _;

fn subset(dir: &Path, name: &str, pb: bool, duplicate: bool) -> std::path::PathBuf {
    let source = write_fixture(dir, FixtureSpec::default());
    let reader = bam::Reader::from_path(&source.bam).unwrap();
    let header = bam::Header::from_template(reader.header());
    let mut record = read_bam(&source.bam)
        .into_iter()
        .find(|r| {
            !r.is_unmapped()
                && !r.is_secondary()
                && !r.is_supplementary()
                && s(r, b"GX").is_some_and(|g| !g.contains(';'))
                && s(r, b"CB").is_some()
                && s(r, b"UR").is_some()
        })
        .unwrap();
    if pb {
        let _ = record.remove_aux(b"PB");
        record.push_aux(b"PB", Aux::String("BC001")).unwrap();
    }
    let path = dir.join(format!("{name}.bam"));
    let mut writer = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
    writer.write(&record).unwrap();
    if duplicate {
        writer.write(&record).unwrap();
    }
    drop(writer);
    bam::index::build(&path, None, bam::index::Type::Bai, 1).unwrap();
    path
}

fn count_chain(dir: &Path, input: &Path, name: &str, args: &[&str]) -> BTreeMap<String, u64> {
    let output = dir.join(format!("{name}.molecules.bam"));
    let mut flags = vec![
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--barcode-metrics",
    ];
    flags.extend(args);
    run(&molecules(), &flags);
    let counters = molecule_counters(&std::path::PathBuf::from(format!(
        "{}.metrics.json",
        output.display()
    )));
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
    check_matrix_total(&counts, &counters);
    assert_eq!(counters["representatives_tagged"], 1);
    assert_eq!(counts.values().sum::<f64>(), 1.0);
    run(
        &cells(),
        &[
            outs.join("raw_feature_bc_matrix.h5").to_str().unwrap(),
            "-o",
            outs.to_str().unwrap(),
            "--force-cells",
            "1",
        ],
    );
    counters
}

#[test]
fn probe_barcode_records_count_outside_flex_in_every_execution_path() {
    let dir = tempfile::tempdir().unwrap();
    let input = subset(dir.path(), "pb", true, false);
    let mut baseline = None;
    for (name, args) in [
        ("serial", vec![]),
        ("indexed", vec!["--evidence-workers", "2"]),
        ("indexed4", vec!["--evidence-workers", "4"]),
        ("sharded", vec!["--shards", "3"]),
        ("sharded4", vec!["--shards", "4"]),
        ("flex", vec!["--flex"]),
    ] {
        count_chain(dir.path(), &input, name, &args);
        let rows = barcode_table(
            &dir.path()
                .join(format!("{name}.molecules.bam.metrics.barcodes.tsv.gz")),
        );
        if let Some(expected) = &baseline {
            assert_eq!(&rows, expected, "{name} per-barcode counters");
        } else {
            baseline = Some(rows);
        }
    }
}

#[test]
fn repeated_record_identity_counts_once_in_serial_indexed_and_sharded_chains() {
    let dir = tempfile::tempdir().unwrap();
    let input = subset(dir.path(), "duplicates", false, true);
    for (name, args) in [
        ("serial", vec![]),
        ("indexed", vec!["--evidence-workers", "2"]),
        ("sharded", vec!["--shards", "3"]),
    ] {
        let counters = count_chain(dir.path(), &input, name, &args);
        assert_eq!(counters["eligible_records"], 2);
        let rows = barcode_table(
            &dir.path()
                .join(format!("{name}.molecules.bam.metrics.barcodes.tsv.gz")),
        );
        let only = rows.values().next().unwrap();
        assert_eq!(rows.len(), 1);
        assert_eq!(
            (
                only["eligible_reads"],
                only["raw_cell_gene_umis"],
                only["molecules_counted"],
                only["collision_groups"],
                only["representatives_tagged"]
            ),
            (2, 1, 1, 0, 1)
        );
    }
}

#[test]
fn mismatched_reference_dictionaries_are_rejected_before_output() {
    let dir = tempfile::tempdir().unwrap();
    let source = write_fixture(dir.path(), FixtureSpec::default());
    let reader = bam::Reader::from_path(&source.bam).unwrap();
    let bytes = String::from_utf8(reader.header().as_bytes().to_vec())
        .unwrap()
        .replace("chr", "other_chr");
    let header = bam::Header::from_template(&bam::HeaderView::from_bytes(bytes.as_bytes()));
    let other = dir.path().join("other.bam");
    let mut writer = bam::Writer::from_path(&other, &header, bam::Format::Bam).unwrap();
    for record in read_bam(&source.bam) {
        writer.write(&record).unwrap();
    }
    drop(writer);
    let output = dir.path().join("out.bam");
    let result = std::process::Command::new(molecules())
        .args([
            source.bam.to_str().unwrap(),
            other.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--shards",
            "2",
        ])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(
        String::from_utf8_lossy(&result.stderr).contains("incompatible BAM reference dictionary")
    );
    assert!(!output.exists());
}

#[test]
fn fractional_saturation_final_depth_matches_molecule_and_matrix_totals() {
    let dir = tempfile::tempdir().unwrap();
    let source = write_fixture(
        dir.path(),
        FixtureSpec {
            umi_pool: 20,
            ..FixtureSpec::default()
        },
    );
    let output = dir.path().join("fractional.bam");
    let saturation = dir.path().join("saturation.tsv.gz");
    run(
        &molecules(),
        &[
            source.bam.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--tagseq",
            "--fractional-collisions",
            "--saturation",
            saturation.to_str().unwrap(),
            "--saturation-interval",
            "1000",
        ],
    );
    let counters = molecule_counters(&std::path::PathBuf::from(format!(
        "{}.metrics.json",
        output.display()
    )));
    let text = common::read_text(&dir.path().join("saturation.tsv"));
    let final_run = text.lines().rfind(|l| l.starts_with("run\t")).unwrap();
    let mass: u64 = final_run.split('\t').nth(4).unwrap().parse().unwrap();
    assert_eq!(mass, counters["molecules_counted"]);
    let outs = dir.path().join("outs");
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
    check_matrix_total(&counts, &counters);
}

#[test]
fn legacy_annotation_module_retains_cli_and_run_entry_points() {
    let _: fn(&sbx10x::annotate_bam::Cli) -> anyhow::Result<()> = sbx10x::annotate_bam::run;
}

#[test]
fn highly_unequal_collisions_keep_tiny_weights_through_matrix_and_cells() {
    let dir = tempfile::tempdir().unwrap();
    let input = subset(dir.path(), "source", false, false);
    let reader = bam::Reader::from_path(&input).unwrap();
    let header = bam::Header::from_template(reader.header());
    let record = read_bam(&input).remove(0);
    let skewed = dir.path().join("skewed.bam");
    let mut writer = bam::Writer::from_path(&skewed, &header, bam::Format::Bam).unwrap();
    for _ in 0..100_000 {
        writer.write(&record).unwrap();
    }
    let mut minority = record.clone();
    minority.remove_aux(b"GX").unwrap();
    let _ = minority.remove_aux(b"ZG");
    minority.push_aux(b"GX", Aux::String("MINORITY")).unwrap();
    writer.write(&minority).unwrap();
    drop(writer);
    let output = dir.path().join("weighted.bam");
    run(
        &molecules(),
        &[
            skewed.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--tagseq",
            "--fractional-collisions",
        ],
    );
    let outs = dir.path().join("outs");
    run(
        &matrix(),
        &[
            output.to_str().unwrap(),
            "-o",
            outs.to_str().unwrap(),
            "--no-gzip",
        ],
    );
    let (mex, is_real) = read_mex(&outs.join("raw_feature_bc_matrix"));
    assert!(is_real);
    let minority = mex
        .iter()
        .find(|((_, gene), _)| gene == "MINORITY")
        .unwrap()
        .1;
    assert!(*minority > 0.0 && *minority < 0.00001);
    let raw = outs.join("raw_feature_bc_matrix.h5");
    let before = hdf5::File::open(&raw)
        .unwrap()
        .dataset("matrix/data")
        .unwrap()
        .read_raw::<f64>()
        .unwrap();
    assert!(before.contains(minority));
    run(
        &cells(),
        &[
            raw.to_str().unwrap(),
            "-o",
            outs.to_str().unwrap(),
            "--force-cells",
            "1",
        ],
    );
    let after = hdf5::File::open(outs.join("filtered_feature_bc_matrix.h5"))
        .unwrap()
        .dataset("matrix/data")
        .unwrap()
        .read_raw::<f64>()
        .unwrap();
    assert_eq!(before, after);
}

#[test]
fn indexed_representative_is_the_best_physical_occurrence() {
    let dir = tempfile::tempdir().unwrap();
    let source = subset(dir.path(), "source", false, false);
    let reader = bam::Reader::from_path(&source).unwrap();
    let header = bam::Header::from_template(reader.header());
    let mut record = read_bam(&source).remove(0);
    let input = dir.path().join("qualities.bam");
    let mut writer = bam::Writer::from_path(&input, &header, bam::Format::Bam).unwrap();
    record.set_mapq(20);
    writer.write(&record).unwrap();
    record.set_mapq(60);
    writer.write(&record).unwrap();
    // Identical best-quality duplicate must not produce a second xf:25.
    writer.write(&record).unwrap();
    drop(writer);
    bam::index::build(&input, None, bam::index::Type::Bai, 1).unwrap();
    count_chain(dir.path(), &input, "indexed", &["--evidence-workers", "2"]);
    let records = read_bam(&dir.path().join("indexed.molecules.bam"));
    assert_eq!(
        records.iter().map(xf).collect::<Vec<_>>(),
        vec![Some(17), Some(25), Some(17)]
    );
    assert_eq!(records[1].mapq(), 60);
}

#[test]
fn sharded_output_header_declares_unsorted_records() {
    let dir = tempfile::tempdir().unwrap();
    let input = subset(dir.path(), "source", false, false);
    let out = dir.path().join("sharded.bam");
    run(
        &molecules(),
        &[
            input.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--shards",
            "2",
        ],
    );
    let reader = bam::Reader::from_path(out).unwrap();
    assert!(String::from_utf8_lossy(reader.header().as_bytes()).contains("SO:unsorted"));
}

#[test]
fn fractional_counts_are_rejected_by_uncalibrated_automatic_cell_calling() {
    let dir = tempfile::tempdir().unwrap();
    let source = write_fixture(
        dir.path(),
        FixtureSpec {
            umi_pool: 20,
            ..FixtureSpec::default()
        },
    );
    let output = dir.path().join("molecules.bam");
    run(
        &molecules(),
        &[
            source.bam.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--tagseq",
            "--fractional-collisions",
        ],
    );
    let outs = dir.path().join("outs");
    run(
        &matrix(),
        &[output.to_str().unwrap(), "-o", outs.to_str().unwrap()],
    );
    let result = std::process::Command::new(cells())
        .args([
            outs.join("raw_feature_bc_matrix.h5").to_str().unwrap(),
            "-o",
            outs.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(
        String::from_utf8_lossy(&result.stderr).contains("not calibrated for fractional counts")
    );
    assert!(!outs.join("filtered_feature_bc_matrix.h5").exists());
}

#[test]
fn fractional_molecule_weights_survive_protospacer_calls_without_rounding() {
    let dir = tempfile::tempdir().unwrap();
    let input = subset(dir.path(), "source", false, false);
    let reader = bam::Reader::from_path(&input).unwrap();
    let header = bam::Header::from_template(reader.header());
    let record = read_bam(&input).remove(0);
    let weighted = dir.path().join("guides.bam");
    let mut writer = bam::Writer::from_path(&weighted, &header, bam::Format::Bam).unwrap();
    for umi in [
        "AAAAAAAAAAAA",
        "CCCCCCCCCCCC",
        "GGGGGGGGGGGG",
        "TTTTTTTTTTTT",
    ] {
        for (gene, support) in [("GUIDE", 9), ("OTHER", 1)] {
            let mut r = record.clone();
            for (tag, value) in [(b"UR", umi), (b"UB", umi), (b"GX", gene), (b"GN", gene)] {
                let _ = r.remove_aux(tag);
                r.push_aux(tag, Aux::String(value)).unwrap();
            }
            for _ in 0..support {
                writer.write(&r).unwrap();
            }
        }
    }
    drop(writer);
    let output = dir.path().join("weighted.bam");
    run(
        &molecules(),
        &[
            weighted.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--tagseq",
            "--fractional-collisions",
        ],
    );
    let outs = dir.path().join("outs");
    run(
        &matrix(),
        &[output.to_str().unwrap(), "-o", outs.to_str().unwrap()],
    );
    run(
        &cells(),
        &[
            outs.join("raw_feature_bc_matrix.h5").to_str().unwrap(),
            "-o",
            outs.to_str().unwrap(),
            "--force-cells",
            "1",
        ],
    );
    let protospacer = bin(
        "sbxpress-protospacer-call",
        env!("CARGO_BIN_EXE_sbxpress-protospacer-call"),
    );
    run(
        &protospacer,
        &[
            outs.join("filtered_feature_bc_matrix.h5").to_str().unwrap(),
            "-o",
            outs.to_str().unwrap(),
        ],
    );
    let text = std::fs::read_to_string(outs.join("protospacer_calls_per_cell.csv")).unwrap();
    let fields: Vec<_> = text.lines().nth(1).unwrap().split(',').collect();
    assert_eq!(fields[2], "GUIDE");
    let count: f64 = fields[3].parse().unwrap();
    assert!((count - 3.6).abs() < 0.00001);
    assert_ne!(count, count.round());
}

#[test]
fn matrix_rejects_invalid_weights_on_real_molecule_output() {
    let dir = tempfile::tempdir().unwrap();
    let input = subset(dir.path(), "source", false, false);
    let output = dir.path().join("molecule.bam");
    run(
        &molecules(),
        &[input.to_str().unwrap(), "-o", output.to_str().unwrap()],
    );
    let reader = bam::Reader::from_path(&output).unwrap();
    let header = bam::Header::from_template(reader.header());
    let record = read_bam(&output).remove(0);
    for (name, weight) in [
        ("negative", -0.25),
        ("nan", f32::NAN),
        ("inf", f32::INFINITY),
    ] {
        let path = dir.path().join(format!("{name}.bam"));
        let mut altered = record.clone();
        altered.push_aux(b"ZW", Aux::Float(weight)).unwrap();
        let mut writer = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
        writer.write(&altered).unwrap();
        drop(writer);
        let result = std::process::Command::new(matrix())
            .args([
                path.to_str().unwrap(),
                "-o",
                dir.path().join(name).to_str().unwrap(),
            ])
            .output()
            .unwrap();
        assert!(!result.status.success());
        assert!(
            String::from_utf8_lossy(&result.stderr).contains("ZW must be finite and non-negative")
        );
    }
}
