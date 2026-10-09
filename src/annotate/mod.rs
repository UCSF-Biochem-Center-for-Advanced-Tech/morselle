mod barcode_metrics;
pub mod call;
pub mod index;
mod multimap;
pub mod pipeline;
pub mod types;

#[cfg(test)]
mod tests;

pub use pipeline::{annotate_record, annotate_record_with, run};
pub use types::{
    Annotation, Call, Cli, ContigIndex, Gene, Interval, Metrics, RecordMetrics, RecordOptions,
    Strand, Transcript,
};
