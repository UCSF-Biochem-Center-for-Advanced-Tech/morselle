//! `extract --whitelist` end to end, including the `--r1-tail-relock` rescue
//! (`XB:i:7`), on generated reads.
//!
//! The extract golden runs without a whitelist (10x whitelists are external
//! reference data and are not committed), so the CB gate and every whitelist
//! rescue were covered only by unit tests. Here a synthetic whitelist and
//! synthetic SBX reads shaped like the smoke fixture's (leader, P5, i5, R1
//! primer, CB, UMI, polyT, cDNA, rc(R2), i7, rc(P7); V2_GEX indexes from the
//! real CE15204 sheet) carry every one-base deletion and insertion in the R1
//! primer's last 15 bp, and (for `--cb-insertion`, `XB:i:8`) one extra base
//! at every inner CB position. Each read's name records its planted CB and
//! UMI, so a rescue can be checked against the truth rather than a digest.
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use clap::Parser;
use morselle::config::Cli;
use serde_json::Value;

const P5: &str = "AATGATACGGCGACCACCGAGAT";
const R1_PRIMER: &str = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT";
const R2_PRIMER: &str = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT";
const P7: &str = "ATCTCGTATGCCGTCTTCTGCTTG";
/// V2_GEX in CE15204 (sheet orientation).
const I7: &str = "TTATTCGAGG";
const I5: &str = "AGCAGGACAG";
/// Reads per planted primer variant.
const PER_VARIANT: usize = 6;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

fn rc(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            _ => 'A',
        })
        .collect()
}

struct Rng(u64);

impl Rng {
    fn seq(&mut self, n: usize) -> String {
        (0..n)
            .map(|_| {
                self.0 = self
                    .0
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                ['A', 'C', 'G', 'T'][(self.0 >> 33) as usize & 3]
            })
            .collect()
    }
}

/// The R1 primer intact, then every one-base deletion and insertion in its
/// last 15 bp (the tail the re-lock searches), named by kind and position.
fn primer_variants() -> Vec<(String, String)> {
    let p = R1_PRIMER;
    let mut out = vec![("exact".to_owned(), p.to_owned())];
    for i in p.len() - 15..p.len() {
        out.push((format!("del{i}"), format!("{}{}", &p[..i], &p[i + 1..])));
        for b in ["A", "C", "G", "T"] {
            out.push((format!("ins{i}{b}"), format!("{}{b}{}", &p[..i], &p[i..])));
        }
    }
    out
}

/// One planted read: name prefix, R1 primer as read, and how the CB and UMI
/// appear in the read (`None`: as planted).
struct Variant {
    name: String,
    primer: String,
    /// Insert a `G` at this CB position, and end the UMI in `T` (so polyT
    /// seems to start one base early and the indel realign, which needs the
    /// polyT one base late, leaves the read to `--cb-insertion`).
    cb_insert_at: Option<usize>,
    /// Other CB edits, with a random UMI (about a quarter end in `T`).
    cb_edit: CbEdit,
}

#[derive(Clone, Copy)]
enum CbEdit {
    None,
    /// Insert a `G` at this CB position.
    Insert(usize),
    /// Delete the base at this CB position.
    Delete(usize),
}

fn variants() -> Vec<Variant> {
    let mut out: Vec<Variant> = primer_variants()
        .into_iter()
        .map(|(name, primer)| Variant {
            name,
            primer,
            cb_insert_at: None,
            cb_edit: CbEdit::None,
        })
        .collect();
    for at in 1..16 {
        out.push(Variant {
            name: format!("cbins{at}"),
            primer: R1_PRIMER.to_owned(),
            cb_insert_at: Some(at),
            cb_edit: CbEdit::None,
        });
    }
    for at in 1..16 {
        for (name, edit) in [
            (format!("cbinsr{at}"), CbEdit::Insert(at)),
            (format!("cbdel{at}"), CbEdit::Delete(at)),
        ] {
            out.push(Variant {
                name,
                primer: R1_PRIMER.to_owned(),
                cb_insert_at: None,
                cb_edit: edit,
            });
        }
    }
    out
}

/// Writes `reads.fastq` and `whitelist.txt` into `dir`.
fn write_inputs(dir: &Path) -> (PathBuf, PathBuf) {
    let mut rng = Rng(0x5eed);
    let whitelist: Vec<String> = (0..400).map(|_| rng.seq(16)).collect();
    let mut fastq = String::new();
    let mut n = 0;
    for v in variants() {
        for k in 0..PER_VARIANT {
            let cb = &whitelist[n % whitelist.len()];
            n += 1;
            let (cb_read, umi) = match (v.cb_insert_at, v.cb_edit) {
                (Some(at), _) => (format!("{}G{}", &cb[..at], &cb[at..]), rng.seq(11) + "T"),
                (None, CbEdit::None) => (cb.clone(), rng.seq(12)),
                (None, CbEdit::Insert(at)) => (format!("{}G{}", &cb[..at], &cb[at..]), rng.seq(12)),
                (None, CbEdit::Delete(at)) => {
                    (format!("{}{}", &cb[..at], &cb[at + 1..]), rng.seq(12))
                }
            };
            let seq = format!(
                "AACAATA{P5}CTACAC{}{}{cb_read}{umi}{}{}{}{I7}{}A",
                rc(I5),
                v.primer,
                "T".repeat(30),
                rng.seq(250),
                rc(R2_PRIMER),
                rc(P7),
            );
            let qual = "I".repeat(seq.len());
            let name = &v.name;
            fastq.push_str(&format!("@{name}_{k}_{cb}{umi}\n{seq}\n+\n{qual}\n"));
        }
    }
    let reads = dir.join("reads.fastq");
    let wl = dir.join("whitelist.txt");
    std::fs::write(&reads, fastq).unwrap();
    std::fs::write(&wl, whitelist.join("\n") + "\n").unwrap();
    (reads, wl)
}

/// One emitted pair: R1 (CB + UMI), R2 header tags, R2 sequence, and the
/// sample file it was written to.
struct Pair {
    r1: String,
    xb: String,
    r2: String,
    sample: String,
}

/// Run extract; return emitted pairs by read name, and the metrics.
fn run(dir: &Path, name: &str, extra: &[&str]) -> (BTreeMap<String, Pair>, Value) {
    let (reads, wl) = write_inputs(dir);
    let out = dir.join(name);
    let mut argv: Vec<String> = [
        "morselle",
        "--fastqs",
        reads.to_str().unwrap(),
        "--out-dir",
        out.to_str().unwrap(),
        "--bclconvert-samplesheet",
        fixture("CE15204_bclconvert_input.csv").to_str().unwrap(),
        "--whitelist",
        wl.to_str().unwrap(),
        "--orientation",
        "p5-led",
        "--threads",
        "2",
        "--no-gzip",
        "--processing-tags",
        // Written when reframe was on by default (off since 2026-10-30);
        // these cases exercise every rescue route, reframe included.
        "--cb-reframe",
        "true",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    argv.extend(extra.iter().map(|s| s.to_string()));
    morselle::run(&Cli::parse_from(&argv)).unwrap();

    let mut pairs = BTreeMap::new();
    for entry in std::fs::read_dir(&out).unwrap() {
        let path = entry.unwrap().path();
        let fname = path.file_name().unwrap().to_str().unwrap().to_owned();
        if !fname.ends_with("_R1_001.fastq") {
            continue;
        }
        let r1 = std::fs::read_to_string(&path).unwrap();
        let r2 = std::fs::read_to_string(out.join(fname.replace("_R1_", "_R2_"))).unwrap();
        let (r1, r2): (Vec<&str>, Vec<&str>) = (r1.lines().collect(), r2.lines().collect());
        for (a, b) in r1.chunks(4).zip(r2.chunks(4)) {
            let id = a[0][1..].split_whitespace().next().unwrap().to_owned();
            let xb = b[0]
                .split('\t')
                .find_map(|t| t.strip_prefix("XB:i:"))
                .unwrap_or("?")
                .to_owned();
            pairs.insert(
                id,
                Pair {
                    r1: a[1].to_owned(),
                    xb,
                    r2: b[1].to_owned(),
                    sample: fname.split("_S").next().unwrap().to_owned(),
                },
            );
        }
    }
    let metrics: Value =
        serde_json::from_str(&std::fs::read_to_string(out.join("metrics_summary.json")).unwrap())
            .unwrap();
    (pairs, metrics)
}

/// Position of a planted CB insertion or deletion, from the read name.
fn cb_edit_at(id: &str) -> Option<usize> {
    let case = id.split('_').next()?;
    let digits = ["cbinsr", "cbins", "cbdel"]
        .iter()
        .find_map(|p| case.strip_prefix(p))?;
    digits.parse().ok()
}

/// The planted CB + UMI recorded in a read name (`case_k_<cb><umi>`).
fn planted(id: &str) -> &str {
    id.rsplit('_').next().unwrap()
}

#[test]
fn relock_rescues_tail_indels_with_the_planted_cb_and_umi() {
    let dir = tempfile::tempdir().unwrap();
    let (off, m_off) = run(dir.path(), "off", &[]);
    let (on, m_on) = run(dir.path(), "on", &["--r1-tail-relock"]);

    // Off: the option is truly off.
    assert_eq!(m_off["whitelist"]["relocked"], 0);
    assert!(off.values().all(|p| p.xb != "7"));
    // Intact primer: exact route, the planted CB and UMI, either way.
    for (id, p) in off.iter().filter(|(id, _)| id.starts_with("exact_")) {
        assert_eq!((p.xb.as_str(), p.r1.as_str()), ("0", planted(id)), "{id}");
    }
    assert_eq!(
        off.keys().filter(|id| id.starts_with("exact_")).count(),
        PER_VARIANT
    );

    // On: every read emitted without the rescue is emitted unchanged.
    for (id, p) in &off {
        let q = on
            .get(id)
            .unwrap_or_else(|| panic!("{id} lost with relock on"));
        assert_eq!(
            (&q.r1, &q.xb, &q.r2, &q.sample),
            (&p.r1, &p.xb, &p.r2, &p.sample),
            "{id}"
        );
    }
    // Every re-locked read carries exactly the planted CB and UMI, comes
    // only from an edited primer, and was not emitted without the rescue.
    let relocked: Vec<_> = on.iter().filter(|(_, p)| p.xb == "7").collect();
    for (id, p) in &relocked {
        assert_eq!(p.r1, planted(id), "{id}");
        assert!(!id.starts_with("exact_") && !off.contains_key(*id), "{id}");
    }
    assert_eq!(m_on["whitelist"]["relocked"], relocked.len());
    // Measured when this test was written: 45 of 450 edited reads, deletions
    // and insertions both (the earlier rescues take most of the others).
    assert!(relocked.len() >= 30, "only {} re-locked", relocked.len());
    assert!(relocked.iter().any(|(id, _)| id.starts_with("del")));
    assert!(relocked.iter().any(|(id, _)| id.starts_with("ins")));
    // Demux re-anchored on the re-located CB: every re-locked read reaches
    // its sample (before re-anchoring, one here went Undetermined).
    for (id, p) in &relocked {
        assert_eq!(p.sample, "V2_GEX", "{id}");
    }
}

/// Every rescue route that emits a CB on these reads emits the planted one
/// (the whitelist is random, so a wrong CB would have to be a coincidence).
#[test]
fn every_emitted_cb_and_umi_is_the_planted_one() {
    let dir = tempfile::tempdir().unwrap();
    let (on, _) = run(dir.path(), "on", &["--r1-tail-relock"]);
    let mut routes = BTreeMap::<String, usize>::new();
    for (id, p) in &on {
        assert_eq!(&p.r1[..16], &planted(id)[..16], "{id} XB {}", p.xb);
        *routes.entry(p.xb.clone()).or_default() += 1;
    }
    // Exact, realign, reframe and re-lock all occur here.
    for xb in ["0", "2", "3", "7"] {
        assert!(routes.contains_key(xb), "no XB {xb} read: {routes:?}");
    }
}

#[test]
fn cb_insertion_rescues_an_extra_base_in_the_barcode() {
    let dir = tempfile::tempdir().unwrap();
    let (off, m_off) = run(dir.path(), "off", &["--r1-tail-relock"]);
    let (on, m_on) = run(dir.path(), "on", &["--r1-tail-relock", "--cb-insertion"]);
    assert_eq!(m_off["whitelist"]["cb_insertion"], 0);
    // Nothing emitted without the rescue changes with it.
    for (id, p) in &off {
        let q = on
            .get(id)
            .unwrap_or_else(|| panic!("{id} lost with --cb-insertion"));
        assert_eq!((&q.r1, &q.xb, &q.sample), (&p.r1, &p.xb, &p.sample), "{id}");
    }
    let rescued: Vec<_> = on.iter().filter(|(_, p)| p.xb == "8").collect();
    assert_eq!(m_on["whitelist"]["cb_insertion"], rescued.len());
    for (id, p) in &rescued {
        // The planted CB and UMI, from a planted insertion only.
        assert_eq!(p.r1, planted(id), "{id}");
        assert!(id.starts_with("cbins") && !off.contains_key(*id), "{id}");
        assert_eq!(p.sample, "V2_GEX", "{id}");
    }
    // Measured when written: 72 of the 90 planted (the earlier rescues take
    // some of the rest).
    assert!(rescued.len() >= 60, "only {} rescued", rescued.len());
}

/// With every rescue on, each emitted R1 is the planted CB and UMI, whatever
/// the route: in particular a UMI that ends in `T` (the first base polyT
/// would claim) is kept, not shifted by one base.
#[test]
fn every_rescued_umi_is_the_planted_one() {
    let dir = tempfile::tempdir().unwrap();
    let (on, _) = run(dir.path(), "on", &["--r1-tail-relock", "--cb-insertion"]);
    let mut routes = BTreeMap::<String, usize>::new();
    let (mut wrong, mut ambiguous) = (Vec::new(), 0);
    for (id, p) in &on {
        *routes.entry(p.xb.clone()).or_default() += 1;
        assert_eq!(&p.r1[..16], &planted(id)[..16], "{id} XB {}", p.xb);
        if p.r1 == planted(id) {
            continue;
        }
        // Two cases the read cannot settle: an edit at the CB's last bases
        // that leaves it on the whitelist (XB 0, which keeps its fixed
        // frame), and one near the end followed by a UMI ending in T (reads
        // like a substitution plus a UMI that does not).
        let near_end = cb_edit_at(id).is_some_and(|at| at >= 13);
        if near_end && (p.xb == "0" || planted(id).ends_with('T')) {
            ambiguous += 1;
        } else {
            wrong.push(format!("{id} XB {} got {}", p.xb, p.r1));
        }
    }
    assert!(
        wrong.is_empty(),
        "{} wrong of {}: {wrong:#?}",
        wrong.len(),
        on.len()
    );
    // Measured when written: 13 of 689 (75 wrong before the UMI anchor).
    assert!(ambiguous <= 20, "{ambiguous} ambiguous");
    for xb in ["0", "2", "3", "7", "8"] {
        assert!(routes.contains_key(xb), "no XB {xb} read: {routes:?}");
    }
    // Planted UMIs ending in T reach the rescue routes, too.
    assert!(on
        .iter()
        .any(|(id, p)| ["2", "3"].contains(&p.xb.as_str()) && planted(id).ends_with('T')));
}
