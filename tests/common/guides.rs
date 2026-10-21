//! A Flex CRISPR guide-capture BAM shaped like real `morselle guides`
//! output (CAT06 `aligned/TSeq100_L1_CRI.bam`): one synthetic contig per
//! guide (`CD45`, `inactive`), records mapped at its start with MAPQ 60,
//! `GX`/`GN` = the guide, then `CR`/`CB`/`UR`/`UB`/`PR`/`PB` in that order.
//!
//! It carries what made real SBX CRISPR overcount (see
//! `src/molecules/umi_rules.rs`): saturated families (20-60 reads per
//! molecule) whose copies include 1-base UMI indels and occasional garbled
//! 1-read UMIs, next to shallow families of genuine 1-3-read molecules.
//! Every error is placed deterministically so each rule's effect is exact:
//!
//! * each deep molecule gets one 2-read indel variant and one 1-read indel
//!   variant of its UMI, plus a 1-mismatch copy (merged by default);
//! * each deep family gets `GARBLED_PER_DEEP_FAMILY` random 1-read UMIs;
//! * shallow families have no errors and a top below 10 reads.
use std::path::{Path, PathBuf};

use rust_htslib::bam::{
    self, header::HeaderRecord, record::Aux, record::Cigar, record::CigarString,
};

use super::Rng;

pub const GARBLED_PER_DEEP_FAMILY: usize = 2;

pub struct GuideFixture {
    pub bam: PathBuf,
    /// Distinct source molecules (the truth a correct count recovers).
    pub molecules: usize,
    /// Deep molecules (each carries one 2-read indel variant).
    pub deep_molecules: usize,
    pub deep_families: usize,
    /// Molecules sequenced by exactly one read (all in shallow families).
    pub true_singletons: usize,
}

struct GuideRead {
    tid: i32,
    name: String,
    cell: String,
    umi: String,
    pb: String,
}

/// One base inserted at `at` (the 12-bp window then loses its last base).
fn insertion(umi: &str, at: usize, base: u8) -> String {
    let mut b = umi.as_bytes().to_vec();
    b.insert(at, base);
    b.truncate(umi.len());
    String::from_utf8(b).unwrap()
}

/// Base at `at` deleted (the window reads one base further: `base`).
fn deletion(umi: &str, at: usize, base: u8) -> String {
    let mut b = umi.as_bytes().to_vec();
    b.remove(at);
    b.push(base);
    String::from_utf8(b).unwrap()
}

fn hamming(a: &str, b: &str) -> usize {
    a.bytes().zip(b.bytes()).filter(|(x, y)| x != y).count()
}

fn mismatch(umi: &str, at: usize) -> String {
    let mut b = umi.as_bytes().to_vec();
    b[at] = match b[at] {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    };
    String::from_utf8(b).unwrap()
}

pub fn write_guide_fixture(dir: &Path, seed: u64) -> GuideFixture {
    let mut rng = Rng::new(seed);
    let guides = ["CD45", "inactive"];
    let pbs = ["BC001", "BC002"];
    let cells: Vec<String> = (0..8).map(|_| rng.dna(16)).collect();
    let mut reads = Vec::new();
    let mut fx = GuideFixture {
        bam: dir.join("guides.bam"),
        molecules: 0,
        deep_molecules: 0,
        deep_families: 0,
        true_singletons: 0,
    };
    let mut m = 0_usize;
    let push = |reads: &mut Vec<GuideRead>,
                tid: i32,
                m: usize,
                copy: String,
                cell: &str,
                umi: String,
                pb: &str| {
        reads.push(GuideRead {
            tid,
            name: format!("SBX:{m}:{copy}"),
            cell: cell.to_owned(),
            umi,
            pb: pb.to_owned(),
        });
    };
    for (c, cell) in cells.iter().enumerate() {
        let pb = pbs[c % 2];
        // Cells 0-3: deep CD45 families (the guide they carry). Cells 4-7:
        // shallow background of the other guide.
        let (tid, deep) = if c < 4 { (0, true) } else { (1, false) };
        let n_mol = 3 + rng.below(4) as usize;
        for _ in 0..n_mol {
            let umi = rng.dna(12);
            if deep {
                let reads_n = 20 + rng.below(41) as usize;
                for k in 0..reads_n {
                    push(&mut reads, tid, m, k.to_string(), cell, umi.clone(), pb);
                }
                // An indel inside a homopolymer can land 1 mismatch from
                // the parent, which plain 1-mismatch correction already
                // merges. Keep only the shifting kind (>= 2 mismatches),
                // the kind the indel rule exists for.
                let (at, two, one) = loop {
                    let at = 1 + rng.below(9) as usize;
                    let two = insertion(&umi, at, b"ACGT"[rng.below(4) as usize]);
                    let one = deletion(&umi, (at + 3) % 12, b"ACGT"[rng.below(4) as usize]);
                    if hamming(&two, &umi) >= 2
                        && hamming(&one, &umi) >= 2
                        && hamming(&one, &two) >= 2
                    {
                        break (at, two, one);
                    }
                };
                push(&mut reads, tid, m, "i0".into(), cell, two.clone(), pb);
                push(&mut reads, tid, m, "i1".into(), cell, two, pb);
                push(&mut reads, tid, m, "d0".into(), cell, one, pb);
                push(
                    &mut reads,
                    tid,
                    m,
                    "s0".into(),
                    cell,
                    mismatch(&umi, (at + 6) % 12),
                    pb,
                );
                fx.deep_molecules += 1;
            } else {
                let reads_n = 1 + rng.below(3) as usize;
                for k in 0..reads_n {
                    push(&mut reads, tid, m, k.to_string(), cell, umi.clone(), pb);
                }
                fx.true_singletons += (reads_n == 1) as usize;
            }
            m += 1;
        }
        if deep {
            fx.deep_families += 1;
            for g in 0..GARBLED_PER_DEEP_FAMILY {
                // Low-quality read: several errors, nowhere near a parent.
                // Its own molecule id, so it never counts as a truth molecule.
                let id = 1_000_000 + c * 10 + g;
                push(&mut reads, tid, id, "0".into(), cell, rng.dna(12), pb);
            }
        }
    }
    fx.molecules = m;
    reads.sort_by(|a, b| (a.tid, &a.name).cmp(&(b.tid, &b.name)));

    let mut header = bam::Header::new();
    header.push_record(
        HeaderRecord::new(b"HD")
            .push_tag(b"VN", "1.6")
            .push_tag(b"SO", "coordinate"),
    );
    for name in guides {
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", name)
                .push_tag(b"LN", 1_000_000),
        );
    }
    let mut writer = bam::Writer::from_path(&fx.bam, &header, bam::Format::Bam).unwrap();
    let seq = b"TTAGCTCTTAAACCTGTTACCCTAAGAACAAACCCATATAAGAAAACCTGAATACGCGGTTGAATCCGCTCCGGTCCTAGCAAACCCATA";
    for r in &reads {
        let mut record = bam::Record::new();
        let cigar = CigarString(vec![Cigar::Match(seq.len() as u32)]);
        record.set(r.name.as_bytes(), Some(&cigar), seq, &vec![20; seq.len()]);
        record.set_tid(r.tid);
        record.set_pos(0);
        record.set_mapq(60);
        record.unset_unmapped();
        let guide = guides[r.tid as usize];
        record.push_aux(b"GX", Aux::String(guide)).unwrap();
        record.push_aux(b"GN", Aux::String(guide)).unwrap();
        record.push_aux(b"CR", Aux::String(&r.cell)).unwrap();
        record
            .push_aux(b"CB", Aux::String(&format!("{}-1", r.cell)))
            .unwrap();
        record.push_aux(b"UR", Aux::String(&r.umi)).unwrap();
        record.push_aux(b"UB", Aux::String(&r.umi)).unwrap();
        record.push_aux(b"PR", Aux::String(&r.pb)).unwrap();
        record.push_aux(b"PB", Aux::String(&r.pb)).unwrap();
        writer.write(&record).unwrap();
    }
    drop(writer);
    bam::index::build(&fx.bam, None, bam::index::Type::Bai, 1).unwrap();
    fx
}
