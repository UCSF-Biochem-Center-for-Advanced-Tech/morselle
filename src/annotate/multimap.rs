//! `--multimap-tags`: annotate one minimap2 query group (primary plus `-N`
//! secondaries, kept together) and tag the group's mapping ambiguity.
use std::collections::HashSet;

use anyhow::Result;
use rust_htslib::bam::{self, record::Aux};

use super::barcode_metrics::observe;
use super::pipeline::{annotate_record_with, replace_i32_tag};
use super::types::{Annotation, Metrics, RecordMetrics, RecordOptions};

/// Annotate one consecutive minimap2 query group and add mapper/annotation
/// ambiguity tags without a name sort or an intermediate BAM. `NH` is the number
/// of *reported* primary+secondary hits, bounded by minimap2's `-N` setting.
pub(super) fn write_multimap_group(
    group: Vec<bam::Record>,
    annotation: &Annotation,
    contig_names: &[Vec<u8>],
    opts: &RecordOptions,
    writer: &mut bam::Writer,
    metrics: &mut Metrics,
) -> Result<()> {
    // The group, not each record, decides the rescue.
    let per_record = RecordOptions {
        confident_mapq: 0,
        ..*opts
    };
    let mut records = Vec::with_capacity(group.len());
    for record in group {
        records.push(annotate_record_with(
            record,
            annotation,
            contig_names,
            &per_record,
        )?);
    }
    let hits = records
        .iter()
        .enumerate()
        .filter_map(|(index, (record, _))| {
            (!record.is_unmapped() && !record.is_supplementary()).then_some(index)
        })
        .collect::<Vec<_>>();
    let nh = hits.len() as i32;
    let rescued = if opts.confident_mapq == 0 {
        loose_rescue(&records, &hits)
    } else {
        cell_ranger_rescue(&records, &hits, opts.confident_mapq)
    };
    for (hi, &index) in hits.iter().enumerate() {
        let (record, record_metrics) = &mut records[index];
        replace_i32_tag(record, b"NH", Some(nh))?;
        replace_i32_tag(record, b"HI", Some((hi + 1) as i32))?;
        let mm = rescued.contains(&index);
        replace_i32_tag(record, b"mm", mm.then_some(1))?;
        record_metrics.rescued = (mm && opts.confident_mapq > 0) as u64;
    }
    for (record, record_metrics) in records {
        observe(metrics, &record, &record_metrics);
        metrics.add_record_metrics(record_metrics);
        writer.write(&record)?;
    }
    Ok(())
}

fn genes(record: &bam::Record) -> impl Iterator<Item = String> + '_ {
    let gx = match record.aux(b"GX") {
        Ok(Aux::String(gx)) => gx,
        _ => "",
    };
    gx.split(';')
        .filter(|gene| !gene.is_empty())
        .map(str::to_owned)
}

/// Pre-0.6.2 rule, kept for `--confident-mapq 0`.
fn loose_rescue(records: &[(bam::Record, RecordMetrics)], hits: &[usize]) -> Vec<usize> {
    let mut union = HashSet::new();
    for &index in hits {
        union.extend(genes(&records[index].0));
    }
    // The reference-pipeline `mm` convention additionally records a MAPQ promotion. We deliberately
    // retain minimap2's native MAPQ; here `mm=1` means the compatible-gene rescue
    // condition alone: multiple reported hits and exactly one compatible gene.
    if hits.len() > 1 && union.len() == 1 {
        hits.to_vec()
    } else {
        Vec::new()
    }
}

/// Cell Ranger's `rescue_alignments_se`: nothing to rescue if any hit is
/// confidently mapped; otherwise, if the exonic-compatible hits name exactly
/// one gene, that read is rescued. Cell Ranger promotes the first such hit
/// to primary; we keep minimap2's flags and MAPQ, so only a primary hit can
/// carry `mm:1` here (a rescue that needs a secondary promoted is left
/// uncounted until the grouped route writes one record per read).
fn cell_ranger_rescue(
    records: &[(bam::Record, RecordMetrics)],
    hits: &[usize],
    confident_mapq: u8,
) -> Vec<usize> {
    if hits
        .iter()
        .any(|&index| records[index].0.mapq() >= confident_mapq)
    {
        return Vec::new();
    }
    let compatible = hits
        .iter()
        .copied()
        .filter(|&index| records[index].1.exonic_compatible)
        .collect::<Vec<_>>();
    let union = compatible
        .iter()
        .flat_map(|&index| genes(&records[index].0))
        .collect::<HashSet<_>>();
    if union.len() != 1 {
        return Vec::new();
    }
    compatible
        .into_iter()
        .filter(|&index| !records[index].0.is_secondary())
        .collect()
}
