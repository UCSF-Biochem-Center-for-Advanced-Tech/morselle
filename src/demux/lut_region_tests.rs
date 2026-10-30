//! `IndexLut`: packed-region path vs per-frame path.
use super::*;

/// The packed-region path and the per-frame path agree on every start,
/// including starts near both read ends, N bases, and the indel pass.
#[test]
fn region_path_matches_per_frame_path() {
    let ks: Vec<Vec<u8>> = vec![
        b"ACGTTGCATG".to_vec(),
        b"TTGACCAGGA".to_vec(),
        b"CAGTCAGTCC".to_vec(),
    ];
    let mut lut = IndexLut::new(&ks, 1);
    lut.add_indels();
    let mut hay = Vec::new();
    for k in &ks {
        hay.extend_from_slice(b"GANTACA");
        let mut d = k.clone();
        d.remove(3);
        hay.extend_from_slice(&d);
        hay.extend_from_slice(b"CC");
        hay.extend_from_slice(k);
    }
    for start in 0..hay.len() + 8 {
        for window in [0isize, 1, 3] {
            for use_indels in [false, true] {
                assert_eq!(
                    lut.find(&hay, Some(start), window, use_indels),
                    lut.find_slow(&hay, start as isize, window, use_indels),
                    "start={start} window={window} indels={use_indels}"
                );
            }
        }
    }
}
