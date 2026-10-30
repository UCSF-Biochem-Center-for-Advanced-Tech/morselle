//! `Whitelist` lookups and rescues.
use super::*;

#[test]
fn encode_roundtrip_and_reject() {
    assert_eq!(encode16(b"AAAAAAAAAAAAAAAA"), Some(0));
    assert_eq!(encode16(b"TTTTTTTTTTTTTTTT"), Some(u32::MAX));
    // distinct barcodes give distinct codes
    assert_ne!(encode16(b"ACGTACGTACGTACGT"), encode16(b"TGCATGCATGCATGCA"));
    assert_eq!(encode16(b"AAAAAAAAAAAAAAAN"), None); // N
    assert_eq!(encode16(b"AAAA"), None); // wrong length
}

fn wl(barcodes: &[&[u8]]) -> Whitelist {
    let codes: FxHashSet<u32> = barcodes.iter().map(|b| encode16(b).unwrap()).collect();
    Whitelist { codes }
}

#[test]
fn classify_exact_within1_fail() {
    let w = wl(&[b"AAACCCAAGAAACACT", b"TTTTTTTTTTTTTTTT"]);
    assert_eq!(w.classify(b"AAACCCAAGAAACACT"), CbStatus::OnWhitelist);
    // one substitution from the first member
    assert_eq!(w.classify(b"AAACCCAAGAAACACA"), CbStatus::Within1);
    // two substitutions -> fail
    assert_eq!(w.classify(b"AAACCCAAGAAACAGG"), CbStatus::Fail);
    // N -> fail
    assert_eq!(w.classify(b"NAACCCAAGAAACACT"), CbStatus::Fail);
}

#[test]
fn decode_roundtrip() {
    for bc in [
        b"AAACCCAAGAAACACT",
        b"ACGTACGTACGTACGT",
        b"TTTTTTTTTTTTTTTT",
    ] {
        assert_eq!(&decode16(encode16(bc).unwrap()), bc);
    }
}

#[test]
fn realign_insertion() {
    let truth = b"AAACCCAAGAAACACT";
    let w = wl(&[truth]);
    // region = CB with one inserted base ('G' after position 4) + UMI tail
    // truth: AAACC|CAAGAAACACT -> insert G -> AAACCGCAAGAAACACT (17) + UMI...
    let region = b"AAACCGCAAGAAACACTNNNNNNNNNNNN";
    assert_eq!(w.realign_cb(region, 1), Some(*truth));
}

#[test]
fn realign_deletion() {
    let truth = b"AAACCCAAGAAACACT";
    let w = wl(&[truth]);
    // region = CB missing its 6th base (the 2nd 'C' group) -> 15 bases + UMI
    // truth AAACC C AAGAAACACT -> delete one C -> AAACCAAGAAACACT (15)
    let region = b"AAACCAAGAAACACTNNNNNNNNNNNN";
    assert_eq!(w.realign_cb(region, -1), Some(*truth));
}

#[test]
fn within2_unique_and_ambiguous() {
    let truth = b"AAACCCAAGAAACACT";
    let w = wl(&[truth, b"TTTTTTTTTTTTTTTT"]);
    // exact
    assert_eq!(w.correct_within2(truth), Some(*truth));
    // two substitutions from truth -> unique
    assert_eq!(w.correct_within2(b"AAACCCAAGAAACAGG"), Some(*truth));
    // far from everything
    assert_eq!(w.correct_within2(b"GGGGGGGGGGGGGGGG"), None);
    // ambiguous: equidistant (1 mm) from two whitelist members
    let a = b"AAAAAAAAAAAAAAAA";
    let b = b"AAAAAAAAAAAAAAAC";
    let w2 = wl(&[a, b]);
    // "AAAAAAAAAAAAAAAG" is 1mm from both a and b -> ambiguous
    assert_eq!(w2.correct_within2(b"AAAAAAAAAAAAAAAG"), None);
}

#[test]
fn within1_unique_and_ambiguous() {
    let truth = b"AAACCCAAGAAACACT";
    let w = wl(&[truth, b"TTTTTTTTTTTTTTTT"]);
    assert_eq!(w.correct_within1(truth), Some(*truth)); // exact
    assert_eq!(w.correct_within1(b"AAACCCAAGAAACACA"), Some(*truth)); // 1 mm
    assert_eq!(w.correct_within1(b"AAACCCAAGAAACAGG"), None); // 2 mm -> none
                                                              // ambiguous: 1 mm from two whitelist members
    let w2 = wl(&[b"AAAAAAAAAAAAAAAA", b"AAAAAAAAAAAAAAAC"]);
    assert_eq!(w2.correct_within1(b"AAAAAAAAAAAAAAAG"), None);
    assert_eq!(w2.classify(b"AAAAAAAAAAAAAAAG"), CbStatus::Fail);
}

#[test]
fn reframe_recovers_shifted_cb() {
    let truth = b"AAACCCAAGAAACACT";
    let w = wl(&[truth, b"TTTTTTTTTTTTTTTT"]);
    // Rigid cb_start grabs one base too early (a deletion put the CB one base
    // 3'): prepend a junk base so truth sits at offset +1. Within-1 noise on
    // the CB ('T' at the last pos) must still resolve uniquely.
    let mut seq = vec![b'G']; // [0] junk, rigid frame would grab seq[0..16]
    seq.extend_from_slice(b"AAACCCAAGAAACACT"); // truth at [1..17]
    seq.extend_from_slice(b"NNNNNNNNNNNN");
    let (cb, off) = w.reframe_cb(&seq, 0, 2).expect("unique +1 hit");
    assert_eq!(&cb, truth);
    assert_eq!(off, 1);
    // no nearby unique hit -> None
    let junk: Vec<u8> = b"GGGGGGGGGGGGGGGGGGGG".to_vec();
    assert_eq!(w.reframe_cb(&junk, 1, 2), None);
}

#[test]
fn realign_ambiguous_or_none() {
    let w = wl(&[b"AAACCCAAGAAACACT", b"TTTTTTTTTTTTTTTT"]);
    // garbage region -> no whitelist hit
    assert_eq!(w.realign_cb(b"GGGGGGGGGGGGGGGGGGGGGGGGGGGG", 1), None);
    // delta 0 not handled
    assert_eq!(w.realign_cb(b"AAACCCAAGAAACACTNNNNNNNNNNNN", 0), None);
}

#[test]
fn insertion_cb_drops_one_inner_base_and_requires_uniqueness() {
    let cb = b"AAACCCAAGAAACACT";
    let w = wl(&[cb]);
    for at in 1..16 {
        let mut r = cb[..at].to_vec();
        r.push(b'G');
        r.extend_from_slice(&cb[at..]);
        r.extend_from_slice(b"TTTT");
        let (got, skip) = w.insertion_cb(&r).unwrap_or_else(|| panic!("at {at}"));
        assert_eq!(&got, cb);
        // A homopolymer insertion can be dropped at any of its bases.
        assert!(r[skip] == b'G' && skip >= 1 && skip <= 15, "at {at}");
    }
    // An exact CB (no insertion) is not an insertion hit. (The caller
    // only asks for whitelist misses; here the next base differs from the
    // CB's last, or dropping base 15 would re-create the exact CB.)
    let mut exact = cb.to_vec();
    exact.extend_from_slice(b"GTTT");
    assert_eq!(w.insertion_cb(&exact).map(|x| x.0), None);
    // A one-deletion reading naming another barcode makes it ambiguous.
    let mut r = cb[..8].to_vec();
    r.push(b'G');
    r.extend_from_slice(&cb[8..]);
    r.push(b'T');
    let mut del: Vec<u8> = r[..15].to_vec();
    del.insert(5, b'C');
    let w2 = wl(&[cb, &del[..16]]);
    assert!(w2.insertion_cb(&r).is_none());
    assert!(w.insertion_cb(&r[..16]).is_none(), "too short");
}
