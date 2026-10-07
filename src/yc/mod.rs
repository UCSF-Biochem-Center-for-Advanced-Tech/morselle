pub mod codec;
pub mod tag;
pub mod types;

#[cfg(test)]
mod tests;

pub use codec::{decode, encode_overlap};
pub use types::{
    code_bases, complement_code, pair_code, parse_u64, BaseType, DecodedPair, DuplexOp, EncodeStats,
    YcError, YcTag, CONSENSUS_GAP_QUAL, GAP,
};
