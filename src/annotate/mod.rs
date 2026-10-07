pub mod call;
pub mod index;
pub mod pipeline;
pub mod types;

#[cfg(test)]
mod tests;

pub use pipeline::{annotate_record, run};
pub use types::{
    Annotation, Call, Cli, ContigIndex, Gene, Interval, Metrics, RecordMetrics, Strand, Transcript,
};
