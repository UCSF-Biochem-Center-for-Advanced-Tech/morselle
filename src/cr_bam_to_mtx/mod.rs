//! Reconstruct a 10x reference raw gene-expression matrix from its barcoded BAM.
mod io_helpers;
mod matrix_ops;
mod types;

#[cfg(test)]
mod tests;

pub use io_helpers::*;
pub use matrix_ops::*;
pub use types::*;

#[cfg(test)]
mod precision_tests;
