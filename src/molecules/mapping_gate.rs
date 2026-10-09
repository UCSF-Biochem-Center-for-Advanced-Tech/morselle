//! `--min-mapq` (Part 2 Layer B): Cell Ranger's rule that only confidently
//! mapped reads, or reads `annotate` rescued with `mm:i:1`, are counted.
use rust_htslib::bam;

use super::helpers::MoleculeTags;
use super::types::Metrics;

/// How a record fares under `molecules --min-mapq` (Layer B, Cell Ranger's
/// counting rule): counted if MAPQ reaches the cutoff or `annotate` rescued
/// it (`mm:i:1`). `NH` only names the reason a low-MAPQ read is dropped:
/// Cell Ranger counts a confidently mapped read whatever its other hits.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MappingGate {
    Confident,
    Rescued,
    LowMapq,
    Multimapper,
}

pub fn mapping_gate(record: &bam::Record, tags: &MoleculeTags, min_mapq: u8) -> MappingGate {
    if min_mapq == 0 || record.mapq() >= min_mapq {
        MappingGate::Confident
    } else if tags.mm == Some(1) {
        MappingGate::Rescued
    } else if tags.nh.is_some_and(|nh| nh > 1) {
        MappingGate::Multimapper
    } else {
        MappingGate::LowMapq
    }
}

/// Gate outcomes of one evidence pass (serial run, shard or region).
#[derive(Clone, Copy, Debug, Default)]
pub struct GateCounts {
    pub rescued: u64,
    pub low_mapq: u64,
    pub multimapper: u64,
}

impl GateCounts {
    /// Counts `gate`; false when the record must be skipped.
    pub fn admit(&mut self, gate: MappingGate) -> bool {
        match gate {
            MappingGate::Confident => true,
            MappingGate::Rescued => {
                self.rescued += 1;
                true
            }
            MappingGate::LowMapq => {
                self.low_mapq += 1;
                false
            }
            MappingGate::Multimapper => {
                self.multimapper += 1;
                false
            }
        }
    }

    pub fn add_to(self, metrics: &mut Metrics) {
        metrics.rescued_records += self.rescued;
        metrics.low_mapq_records += self.low_mapq;
        metrics.multimapper_records += self.multimapper;
    }
}
