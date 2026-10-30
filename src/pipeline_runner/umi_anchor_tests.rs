use super::*;

const CB: &[u8] = b"ACGTACGGTCAGCATC";
const UMI_A: &[u8] = b"GACCATGAGCGA"; // does not end in T
const UMI_T: &[u8] = b"GACCATGAGCGT"; // ends in T
const UMI_TGTA: &[u8] = b"GACCAGCTGTGA"; // T-rich end, absorbed by a tolerant run

/// `prefix` + barcode as read + UMI + poly-T + cDNA; returns (read, CB start).
fn read(cb_read: &[u8], umi: &[u8]) -> (Vec<u8>, usize) {
    let mut s = b"CTTCCGATCT".to_vec();
    let start = s.len();
    s.extend_from_slice(cb_read);
    s.extend_from_slice(umi);
    s.extend_from_slice(&[b'T'; 30]);
    s.extend_from_slice(b"GATTACAGATTACA");
    (s, start)
}

fn umi_of(cb_read: &[u8], umi: &[u8], default_shift: isize) -> Vec<u8> {
    let (s, st) = read(cb_read, umi);
    let def = (st as isize + 16 + default_shift) as usize;
    let u = umi_start(&s, CB, st, def, 12);
    s[u..u + 12].to_vec()
}

fn with_insert(at: usize, b: u8) -> Vec<u8> {
    let mut v = CB.to_vec();
    v.insert(at, b);
    v
}

fn with_delete(at: usize) -> Vec<u8> {
    let mut v = CB.to_vec();
    v.remove(at);
    v
}

#[test]
fn exact_barcode_keeps_the_positional_umi_even_when_it_ends_in_t() {
    for umi in [UMI_A, UMI_T, UMI_TGTA] {
        assert_eq!(umi_of(CB, umi, 0), umi);
    }
}

#[test]
fn inner_insertion_or_deletion_places_the_umi_after_the_real_barcode_end() {
    for at in 1..15 {
        for umi in [UMI_A, UMI_T, UMI_TGTA] {
            // Route default already accounts for the indel (+1 / -1).
            assert_eq!(umi_of(&with_insert(at, b'G'), umi, 1), umi, "ins {at}");
            assert_eq!(umi_of(&with_delete(at), umi, -1), umi, "del {at}");
            // Route default ignores it (fixed frame): alignment still finds it
            // when the indel is the only edit.
            assert_eq!(
                umi_of(&with_insert(at, b'G'), umi, 0),
                umi,
                "ins {at} fixed"
            );
        }
    }
}

#[test]
fn last_base_indel_that_looks_like_a_substitution_is_resolved_by_poly_t() {
    // Insertion before the last CB base: Hamming-1 at the fixed frame, and
    // the poly-T one base late says the barcode is 17 bp.
    for umi in [UMI_A, UMI_TGTA] {
        assert_eq!(umi_of(&with_insert(15, b'G'), umi, 0), umi);
    }
    // Deletion of the last CB base, UMI not ending in T: Hamming-1 at the
    // fixed frame (UMI's first base read as the CB's last), poly-T one base
    // early.
    assert_eq!(umi_of(&with_delete(15), UMI_A, 0), UMI_A);
}

#[test]
fn a_real_substitution_keeps_the_fixed_frame() {
    let mut cb = CB.to_vec();
    cb[7] = b'A';
    for umi in [UMI_A, UMI_T, UMI_TGTA] {
        assert_eq!(umi_of(&cb, umi, 0), umi);
    }
}

#[test]
fn no_poly_t_falls_back_to_the_route_default() {
    let mut s = b"CTTCCGATCT".to_vec();
    s.extend_from_slice(CB);
    s.extend_from_slice(UMI_A);
    s.extend_from_slice(b"GATTACAGATTACAGATTACA");
    assert_eq!(umi_start(&s, CB, 10, 26, 12), 26);
    // Read too short for any alternative end.
    assert_eq!(umi_start(&s[..12], CB, 10, 26, 12), 26);
}

/// The case the read cannot settle: an insertion before the last CB base
/// followed by a UMI ending in T reads exactly like a last-base substitution
/// followed by a UMI that does not. The exact poly-T fit picks the
/// substitution reading; either way the CB is right and the UMI is one of
/// the two readings.
#[test]
fn last_base_insertion_with_a_t_ending_umi_is_ambiguous() {
    let read_cb = with_insert(15, b'G');
    let got = umi_of(&read_cb, UMI_T, 0);
    let mut sub_reading = vec![CB[15]];
    sub_reading.extend_from_slice(&UMI_T[..11]);
    assert_eq!(got, sub_reading);
}
