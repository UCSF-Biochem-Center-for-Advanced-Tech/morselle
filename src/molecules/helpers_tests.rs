//! The single-pass / allocation-free tag helpers must match the original
//! per-tag implementations exactly.
use rust_htslib::bam::{self, record::Aux};

use super::helpers::{choose_gene_parts, MoleculeTags};

/// The pre-optimization `choose_gene`, kept verbatim as the reference.
fn reference_choose_gene(gx: &str, gn: Option<&str>, zg: Option<&str>) -> Option<(String, String)> {
    let gn = gn.map(str::to_owned).unwrap_or_else(|| gx.to_owned());
    let genes: Vec<&str> = gx.split(';').filter(|v| !v.is_empty()).collect();
    let names: Vec<&str> = gn.split(';').collect();
    if genes.len() == 1 {
        return Some((
            genes[0].to_owned(),
            names.first().copied().unwrap_or(genes[0]).to_owned(),
        ));
    }
    let zg = zg?.to_owned();
    if zg.contains(';') || zg.is_empty() {
        return None;
    }
    let index = genes.iter().position(|gene| *gene == zg)?;
    Some((
        zg,
        names.get(index).copied().unwrap_or(genes[index]).to_owned(),
    ))
}

#[test]
fn choose_gene_parts_matches_reference_on_random_inputs() {
    let atoms = ["", "G1", "G2", "G3", ";", "g1", "g2"];
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    for _ in 0..20_000 {
        let mut field = || -> String {
            let len = next() % 6;
            (0..len)
                .map(|_| atoms[(next() % atoms.len() as u64) as usize])
                .collect::<Vec<_>>()
                .join("")
        };
        let (a, b, c) = (field(), field(), field());
        let gn = (next() % 3 != 0).then_some(b.as_str());
        let zg = (next() % 3 != 0).then_some(c.as_str());
        let got = choose_gene_parts(&a, gn, zg).map(|(g, n)| (g.to_owned(), n.to_owned()));
        assert_eq!(
            got,
            reference_choose_gene(&a, gn, zg),
            "gx={a:?} gn={gn:?} zg={zg:?}"
        );
    }
}

fn record_with(tags: Vec<(&[u8], Aux)>) -> bam::Record {
    let mut record = bam::Record::new();
    record.set(b"q", None, b"ACGT", &[30; 4]);
    for (tag, value) in tags {
        record.push_aux(tag, value).unwrap();
    }
    record
}

#[test]
fn molecule_tags_match_per_tag_lookups() {
    let record = record_with(vec![
        (b"NM", Aux::I32(1)),
        (b"CB", Aux::String("AAAC-1")),
        (b"UR", Aux::String("ACGT")),
        (b"UB", Aux::String("ACGA")),
        (b"GX", Aux::String("G1;G2")),
        (b"GN", Aux::String("g1;g2")),
        (b"ZG", Aux::String("G2")),
        (b"PB", Aux::Char(b'X')), // wrong type: ignored, like aux_string
    ]);
    let tags = MoleculeTags::read(&record);
    assert_eq!(tags.cb, Some("AAAC-1"));
    assert_eq!(tags.umi(), Some("ACGT"));
    assert_eq!(tags.gene(), Some(("G2", "g2")));
    assert_eq!(tags.pb, None);
    assert!(!tags.has_output_tags);

    let rerun = record_with(vec![(b"UB", Aux::String("ACGA")), (b"xf", Aux::I32(25))]);
    let tags = MoleculeTags::read(&rerun);
    assert_eq!(tags.umi(), Some("ACGA"));
    assert!(tags.has_output_tags);
}

/// The raw-byte scanner must agree with `record.aux(tag)` (htslib lookup +
/// rust-htslib decode) for every wanted tag, on records mixing every aux
/// type, duplicate tags, wrong-typed tags, and invalid UTF-8.
#[test]
fn raw_scan_matches_per_tag_aux_lookups() {
    use rust_htslib::bam::record::AuxArray;
    let wanted: [&[u8; 2]; 7] = [b"CB", b"UR", b"UB", b"GX", b"GN", b"ZG", b"PB"];
    let pool: [&[u8; 2]; 11] = [
        b"CB", b"UR", b"UB", b"GX", b"GN", b"ZG", b"PB", b"NM", b"CY", b"xf", b"zz",
    ];
    let mut state: u64 = 0xD1B5_4A32_D192_ED03;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    for case in 0..3_000 {
        let mut record = bam::Record::new();
        let seq_len = (next() % 5) as usize;
        record.set(b"read", None, &b"ACGTA"[..seq_len], &vec![30; seq_len]);
        let n_tags = next() % 9;
        for _ in 0..n_tags {
            let tag = pool[(next() % pool.len() as u64) as usize];
            let ints = [1i16, -2, 3];
            // push_aux refuses duplicates, so append raw for repeats.
            let value = match next() % 8 {
                0 => Aux::String("AAAC-1"),
                1 => Aux::String("G1;G2"),
                2 => Aux::I32(7),
                3 => Aux::Char(b'X'),
                4 => Aux::Float(0.5),
                5 => Aux::ArrayI16(AuxArray::from(&ints[..])),
                6 => Aux::HexByteArray("1AE3"),
                _ => Aux::String(""),
            };
            if record.aux(tag).is_err() {
                record.push_aux(tag, value).unwrap();
            } else if let Aux::String(s) = value {
                // duplicate occurrence appended after the first
                let mut data = crate::bam_helpers::raw_data(&record).to_vec();
                data.extend_from_slice(tag);
                data.push(b'Z');
                data.extend_from_slice(s.as_bytes());
                data.push(0);
                record.set_data(&data);
            }
        }
        if next() % 5 == 0 {
            // a non-wanted Z tag holding invalid UTF-8
            let mut data = crate::bam_helpers::raw_data(&record).to_vec();
            data.extend_from_slice(b"QQZ\xff\xfe\0");
            record.set_data(&data);
        }
        let tags = MoleculeTags::read(&record);
        let got = [
            tags.cb, tags.ur, tags.ub, tags.gx, tags.gn, tags.zg, tags.pb,
        ];
        for (tag, value) in wanted.iter().zip(got) {
            let expected = match record.aux(*tag) {
                Ok(Aux::String(s)) => Some(s),
                _ => None,
            };
            assert_eq!(
                value,
                expected,
                "case {case}, tag {}",
                String::from_utf8_lossy(*tag)
            );
        }
        assert_eq!(
            tags.has_output_tags,
            record.aux(b"xf").is_ok(),
            "case {case}"
        );
    }
}
