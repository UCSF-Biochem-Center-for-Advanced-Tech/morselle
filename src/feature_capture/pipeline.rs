use std::fs::{self, File};
use std::io::{BufWriter, Write as _};
use std::path::Path;
use std::time::Instant;

use anyhow::{ensure, Context, Result};
use needletail::parse_fastx_file;
use rust_htslib::bam::{
    self, header::HeaderRecord, record::Aux, record::Cigar, record::CigarString,
};
use tracing::info;

use super::matcher::FeatureMatcher;
use super::types::{Cli, Feature, FeatureCount, FeatureReference, Metrics, SYNTHETIC_CONTIG_LEN};

pub fn parse_tag<'a>(comment: &'a str, tag: &str) -> Option<&'a str> {
    let start = comment.find(tag)? + tag.len();
    let rest = &comment[start..];
    let end = rest.find('\t').unwrap_or(rest.len());
    Some(&rest[..end])
}

fn features_tsv_bytes(features: &[Feature]) -> Vec<u8> {
    let mut out = Vec::new();
    for f in features {
        out.extend_from_slice(f.id.as_bytes());
        out.push(b'\t');
        out.extend_from_slice(f.name.as_bytes());
        out.push(b'\t');
        out.extend_from_slice(f.feature_type.as_bytes());
        out.push(b'\n');
    }
    out
}

pub fn establish_features_tsv(path: &Path, features: &[Feature]) -> Result<()> {
    let expected = features_tsv_bytes(features);
    if path.exists() {
        let observed = fs::read(path).with_context(|| format!("reading {}", path.display()))?;
        ensure!(
            observed == expected,
            "existing feature manifest {} does not match this feature reference; \
             use a matching CSV or choose a new path",
            path.display()
        );
        return Ok(());
    }
    if let Some(parent) = path.parent().filter(|p| !p.as_os_str().is_empty()) {
        fs::create_dir_all(parent)?;
    }
    fs::write(path, &expected).with_context(|| format!("writing {}", path.display()))
}

pub fn run(cli: &Cli) -> Result<()> {
    let started = Instant::now();
    let reference = FeatureReference::load(&cli.feature_reference)?;
    info!(features = reference.features.len(), path = %cli.feature_reference.display(), "loaded feature reference");
    let matcher = FeatureMatcher::build(
        &reference.features,
        cli.anchor_max_mm,
        cli.capture_max_mm,
        cli.anchor_window,
    );

    if let Some(path) = &cli.features_out {
        establish_features_tsv(path, &reference.features)?;
    }

    let mut header = bam::Header::new();
    for f in &reference.features {
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", &f.id)
                .push_tag(b"LN", SYNTHETIC_CONTIG_LEN),
        );
    }
    if let Some(parent) = cli.output.parent().filter(|p| !p.as_os_str().is_empty()) {
        fs::create_dir_all(parent)?;
    }
    let mut writer = bam::Writer::from_path(&cli.output, &header, bam::Format::Bam)
        .with_context(|| format!("creating {}", cli.output.display()))?;

    let mut reads = 0u64;
    let mut matched = 0u64;
    let mut per_feature = vec![0u64; reference.features.len()];

    let mut fastq = parse_fastx_file(&cli.r2_fastq)
        .with_context(|| format!("opening {}", cli.r2_fastq.display()))?;
    while let Some(rec) = fastq.next() {
        let rec = rec.context("reading R2 FASTQ record")?;
        reads += 1;
        let header_line = String::from_utf8_lossy(rec.id());
        let (qname, comment) = header_line.split_once('\t').unwrap_or((&header_line, ""));
        let seq = rec.seq();
        let qual = rec.qual().unwrap_or(&[]);

        let mut record = bam::Record::new();
        let hit = matcher.locate(&seq);
        if let Some(feature_idx) = hit {
            matched += 1;
            per_feature[feature_idx] += 1;
            let cigar = CigarString(vec![Cigar::Match(seq.len() as u32)]);
            record.set(qname.as_bytes(), Some(&cigar), &seq, qual);
            record.set_tid(feature_idx as i32);
            record.set_pos(0);
            record.set_mapq(60);
            record.unset_unmapped();
            let f = &reference.features[feature_idx];
            record.push_aux(b"GX", Aux::String(&f.id))?;
            record.push_aux(b"GN", Aux::String(&f.name))?;
        } else {
            record.set(qname.as_bytes(), None, &seq, qual);
        }
        for tag in ["CR", "CB", "UR", "UB", "PR", "PB"] {
            if let Some(value) = parse_tag(comment, &format!("{tag}:Z:")) {
                if !value.is_empty() {
                    record.push_aux(tag.as_bytes(), Aux::String(value))?;
                }
            }
        }
        writer.write(&record)?;
    }
    drop(writer);

    let metrics = Metrics {
        reads,
        matched,
        hit_rate_pct: if reads > 0 {
            100.0 * matched as f64 / reads as f64
        } else {
            0.0
        },
        per_feature: reference
            .features
            .iter()
            .zip(per_feature)
            .map(|(f, calls)| FeatureCount {
                id: f.id.clone(),
                calls,
            })
            .collect(),
    };
    if let Some(path) = &cli.metrics {
        let mut w = BufWriter::new(File::create(path)?);
        w.write_all(serde_json::to_string_pretty(&metrics)?.as_bytes())?;
    }

    info!(
        reads = metrics.reads,
        matched = metrics.matched,
        hit_rate_pct = format!("{:.1}", metrics.hit_rate_pct),
        elapsed = format!("{:.1}s", started.elapsed().as_secs_f64()),
        out = %cli.output.display(),
        "RUN SUMMARY"
    );
    Ok(())
}
