//! Shared helpers for the end-to-end / golden tests: a deterministic
//! generator for annotated BAMs shaped like real `annotate` output, a
//! binary runner, and stable digests of BAM / MEX outputs.
//!
//! The fixture deliberately carries everything real data has and earlier
//! unit-test fixtures did not: minimap2 + extract + annotate tags in their
//! real order (so `UB` and `ZG` already exist), multi-gene `GX` resolved by
//! `ZG`, 1-mismatch UMI errors, cross-gene UMI collisions, PCR duplicates,
//! unmapped and untagged reads, and optional Flex `PB` probe barcodes.
#![allow(dead_code)]

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::io::Read as _;
use std::path::{Path, PathBuf};
use std::process::Command;

use rust_htslib::bam::{
    self, header::HeaderRecord, record::Aux, record::Cigar, record::CigarString, Read as _,
};

pub mod guides;

/// Read a stage's gzip TSV into named integer counters without hiding sample rows.
pub fn barcode_table(path: &Path) -> BTreeMap<(String, String), BTreeMap<String, u64>> {
    let mut text = String::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(path).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    assert_eq!(&header[..2], &["barcode", "sample"]);
    lines
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            assert_eq!(fields.len(), header.len());
            let counters = header[2..]
                .iter()
                .zip(&fields[2..])
                .map(|(&name, value)| (name.to_owned(), value.parse().unwrap()))
                .collect();
            ((fields[0].to_owned(), fields[1].to_owned()), counters)
        })
        .collect()
}

/// Small xorshift64 RNG: deterministic across platforms and Rust versions.
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Self(seed.max(1))
    }
    pub fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    pub fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    pub fn chance(&mut self, per_mille: u64) -> bool {
        self.below(1000) < per_mille
    }
    pub fn dna(&mut self, len: usize) -> String {
        (0..len)
            .map(|_| b"ACGT"[self.below(4) as usize] as char)
            .collect()
    }
}

#[derive(Clone, Copy)]
pub struct FixtureSpec {
    pub seed: u64,
    /// Molecules to draw (each becomes 1..=5 reads).
    pub molecules: usize,
    /// Cell barcodes; the first `real_cells` get ~20x more molecules.
    pub cells: usize,
    pub real_cells: usize,
    pub genes: usize,
    /// Size of the UMI pool. Small pools force cross-gene collisions.
    pub umi_pool: usize,
    /// Tag reads with Flex probe barcodes `BC001`.. (`None`: no `PB` tag).
    pub probe_barcodes: Option<usize>,
    /// PCR duplicates of a molecule share its start, soft clip and strand
    /// (as real duplicates do), and molecules use both strands. Off keeps
    /// the original per-read random positions the golden fixtures use.
    pub shared_start: bool,
}

impl Default for FixtureSpec {
    fn default() -> Self {
        Self {
            seed: 0x5EED_0001,
            molecules: 4_000,
            cells: 60,
            real_cells: 12,
            genes: 40,
            umi_pool: 600,
            probe_barcodes: None,
            shared_start: false,
        }
    }
}

pub struct Fixture {
    pub bam: PathBuf,
    /// Distinct cell barcodes, without the `-1` suffix (whitelist format).
    pub whitelist: PathBuf,
    pub cells: Vec<String>,
    pub probe_barcodes: Vec<String>,
}

struct Read {
    tid: i32,
    pos: i64,
    name: String,
    cb: Option<String>,
    umi: String,
    gx: String,
    gn: String,
    zg: Option<String>,
    pb: Option<String>,
    mapq: u8,
    unmapped: bool,
    clip: u32,
    reverse: bool,
}

/// Write a coordinate-sorted, indexed annotated BAM plus its whitelist.
pub fn write_fixture(dir: &Path, spec: FixtureSpec) -> Fixture {
    let mut rng = Rng::new(spec.seed);
    let cells: Vec<String> = (0..spec.cells).map(|_| rng.dna(16)).collect();
    let umis: Vec<String> = (0..spec.umi_pool).map(|_| rng.dna(12)).collect();
    let genes: Vec<(String, String, i32, i64)> = (0..spec.genes)
        .map(|g| {
            let (tid, pos) = (rng.below(2) as i32, 1_000 + rng.below(900_000) as i64);
            (format!("ENSG{g:011}"), format!("GENE{g}"), tid, pos)
        })
        .collect();
    let probe_barcodes: Vec<String> = (1..=spec.probe_barcodes.unwrap_or(0))
        .map(|i| format!("BC{i:03}"))
        .collect();
    // Each cell is assigned one probe barcode (two for every 5th cell, the
    // composite-GEM case probe-demux has to handle).
    let cell_pbs: Vec<Vec<String>> = (0..spec.cells)
        .map(|c| {
            if probe_barcodes.is_empty() {
                return vec![];
            }
            let first = probe_barcodes[c % probe_barcodes.len()].clone();
            if c % 5 == 0 && probe_barcodes.len() > 1 {
                vec![
                    first,
                    probe_barcodes[(c + 1) % probe_barcodes.len()].clone(),
                ]
            } else {
                vec![first]
            }
        })
        .collect();

    let mut reads = Vec::new();
    for m in 0..spec.molecules {
        let cell = if rng.chance(800) {
            rng.below(spec.real_cells as u64) as usize
        } else {
            rng.below(spec.cells as u64) as usize
        };
        let gene = rng.below(spec.genes as u64) as usize;
        let umi = umis[rng.below(spec.umi_pool as u64) as usize].clone();
        let pb = (!cell_pbs[cell].is_empty())
            .then(|| cell_pbs[cell][rng.below(cell_pbs[cell].len() as u64) as usize].clone());
        let multi_gene = rng.chance(50);
        let other = (gene + 1 + rng.below(spec.genes as u64 - 1) as usize) % spec.genes;
        let copies = 1 + rng.below(5);
        // Drawn only in shared_start mode, so default fixtures are unchanged.
        let start = spec.shared_start.then(|| {
            (
                genes[gene].3 + rng.below(400) as i64,
                rng.below(3) as u32 * 10,
                rng.chance(500),
            )
        });
        for copy in 0..copies {
            let mut read_umi = umi.clone();
            if copy > 0 && rng.chance(40) {
                // 1-mismatch sequencing error on a duplicate
                let p = rng.below(12) as usize;
                let mut b = read_umi.into_bytes();
                b[p] = b"ACGT"[((b"ACGT".iter().position(|&x| x == b[p]).unwrap() as u64
                    + 1
                    + rng.below(3))
                    % 4) as usize];
                read_umi = String::from_utf8(b).unwrap();
            }
            let (gx, gn, zg) = if multi_gene {
                (
                    format!("{};{}", genes[gene].0, genes[other].0),
                    format!("{};{}", genes[gene].1, genes[other].1),
                    Some(genes[gene].0.clone()),
                )
            } else {
                (genes[gene].0.clone(), genes[gene].1.clone(), None)
            };
            reads.push(Read {
                tid: genes[gene].2,
                pos: genes[gene].3 + rng.below(400) as i64,
                name: format!("SBX:{m}:{copy}"),
                cb: (!rng.chance(10)).then(|| format!("{}-1", cells[cell])),
                umi: read_umi,
                gx,
                gn,
                zg,
                // ~2% of Flex reads lack a confident probe barcode.
                pb: pb.clone().filter(|_| !rng.chance(20)),
                mapq: [60, 60, 60, 3, 255][rng.below(5) as usize],
                unmapped: rng.chance(15),
                clip: rng.below(3) as u32 * 10,
                reverse: false,
            });
            if let (Some((pos, clip, reverse)), Some(read)) = (start, reads.last_mut()) {
                (read.pos, read.clip, read.reverse) = (pos, clip, reverse);
            }
        }
    }
    reads.sort_by(|a, b| (a.tid, a.pos, &a.name).cmp(&(b.tid, b.pos, &b.name)));

    let path = dir.join("annotated.bam");
    let mut header = bam::Header::new();
    header.push_record(
        HeaderRecord::new(b"HD")
            .push_tag(b"VN", "1.6")
            .push_tag(b"SO", "coordinate"),
    );
    for name in ["chr1", "chr2"] {
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", name)
                .push_tag(b"LN", 1_000_000),
        );
    }
    let mut writer = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
    let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    for r in &reads {
        let mut record = bam::Record::new();
        let cigar = CigarString(if r.clip > 0 {
            vec![
                Cigar::Match(seq.len() as u32 - r.clip),
                Cigar::SoftClip(r.clip),
            ]
        } else {
            vec![Cigar::Match(seq.len() as u32)]
        });
        record.set(r.name.as_bytes(), Some(&cigar), seq, &[30; 60]);
        record.set_tid(r.tid);
        record.set_pos(r.pos);
        record.set_mapq(r.mapq);
        if r.reverse {
            record.set_flags(record.flags() | 0x10);
        }
        if r.unmapped {
            record.set_unmapped();
        } else {
            record.unset_unmapped();
        }
        // Tag order as minimap2 -y + annotate produce it.
        record.push_aux(b"NM", Aux::I32(1)).unwrap();
        record.push_aux(b"AS", Aux::I32(100)).unwrap();
        if let Some(cb) = &r.cb {
            record.push_aux(b"CR", Aux::String(&cb[..16])).unwrap();
            record.push_aux(b"CB", Aux::String(cb)).unwrap();
            record
                .push_aux(b"CY", Aux::String("IIIIIIIIIIIIIIII"))
                .unwrap();
        }
        record.push_aux(b"UR", Aux::String(&r.umi)).unwrap();
        record.push_aux(b"UB", Aux::String(&r.umi)).unwrap();
        record.push_aux(b"UY", Aux::String("IIIIIIIIIIII")).unwrap();
        if let Some(pb) = &r.pb {
            record.push_aux(b"PR", Aux::String(pb)).unwrap();
            record.push_aux(b"PB", Aux::String(pb)).unwrap();
        }
        record.push_aux(b"RE", Aux::Char(b'E')).unwrap();
        record.push_aux(b"GX", Aux::String(&r.gx)).unwrap();
        record.push_aux(b"GN", Aux::String(&r.gn)).unwrap();
        if let Some(zg) = &r.zg {
            record.push_aux(b"ZG", Aux::String(zg)).unwrap();
        }
        writer.write(&record).unwrap();
    }
    drop(writer);
    bam::index::build(&path, None, bam::index::Type::Bai, 1).unwrap();

    let whitelist = dir.join("whitelist.txt");
    std::fs::write(&whitelist, cells.join("\n") + "\n").unwrap();
    Fixture {
        bam: path,
        whitelist,
        cells,
        probe_barcodes,
    }
}

/// Path of a crate binary. `MORSELLE_TEST_BIN_DIR` points the tests at a
/// different build (e.g. an older release) to compare behaviour.
pub fn bin(name: &str, default: &str) -> PathBuf {
    match std::env::var_os("MORSELLE_TEST_BIN_DIR") {
        Some(dir) => Path::new(&dir).join(name),
        None => PathBuf::from(default),
    }
}

pub fn run(bin: &Path, args: &[&str]) {
    let out = Command::new(bin).args(args).output().unwrap();
    assert!(
        out.status.success(),
        "{} {args:?} failed:\n{}",
        bin.display(),
        String::from_utf8_lossy(&out.stderr)
    );
}

pub fn read_bam(path: &Path) -> Vec<bam::Record> {
    bam::Reader::from_path(path)
        .unwrap()
        .records()
        .map(Result::unwrap)
        .collect()
}

/// Stable textual form of one aux value (independent of rust-htslib's
/// `Debug` output, so library upgrades can't churn golden digests).
fn aux_text(value: &Aux) -> String {
    match value {
        Aux::Char(v) => format!("A:{}", *v as char),
        Aux::I8(v) => format!("i:{v}"),
        Aux::U8(v) => format!("i:{v}"),
        Aux::I16(v) => format!("i:{v}"),
        Aux::U16(v) => format!("i:{v}"),
        Aux::I32(v) => format!("i:{v}"),
        Aux::U32(v) => format!("i:{v}"),
        Aux::Float(v) => format!("f:{:08x}", v.to_bits()),
        Aux::Double(v) => format!("d:{:016x}", v.to_bits()),
        Aux::String(v) => format!("Z:{v}"),
        Aux::HexByteArray(v) => format!("H:{v}"),
        other => format!("B:{other:?}"),
    }
}

/// One line per record: name, flag, position, and every aux tag in order.
pub fn record_text(record: &bam::Record) -> String {
    let mut line = format!(
        "{}\t{}\t{}\t{}",
        String::from_utf8_lossy(record.qname()),
        record.flags(),
        record.tid(),
        record.pos()
    );
    for (tag, value) in record.aux_iter().map(Result::unwrap) {
        let _ = write!(
            line,
            "\t{}:{}",
            String::from_utf8_lossy(tag),
            aux_text(&value)
        );
    }
    line
}

/// FNV-1a 64-bit digest, hex.
pub fn digest(text: &str) -> String {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for &b in text.as_bytes() {
        h ^= b as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    format!("{h:016x}")
}

pub fn bam_digest(path: &Path) -> String {
    let text: String = read_bam(path)
        .iter()
        .map(|r| record_text(r) + "\n")
        .collect();
    digest(&text)
}

/// Read `path`, or `path.gz` if only the gzipped file exists.
pub fn read_text(path: &Path) -> String {
    use std::io::Read as _;
    if path.exists() {
        return std::fs::read_to_string(path).unwrap();
    }
    let gz = PathBuf::from(format!("{}.gz", path.display()));
    let mut text = String::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(&gz).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    text
}

/// Raw MEX matrix as (barcode, feature id) -> value, read from a
/// `raw_feature_bc_matrix` / `filtered_feature_bc_matrix` dir (plain or gzip).
pub fn read_mex(dir: &Path) -> (BTreeMap<(String, String), f64>, bool) {
    let lines =
        |f: &str| -> Vec<String> { read_text(&dir.join(f)).lines().map(str::to_owned).collect() };
    let features: Vec<String> = lines("features.tsv")
        .iter()
        .map(|l| l.split('\t').next().unwrap().to_owned())
        .collect();
    let barcodes = lines("barcodes.tsv");
    let mtx = lines("matrix.mtx");
    let is_real = mtx[0].contains(" real ");
    let mut counts = BTreeMap::new();
    for line in mtx.iter().filter(|l| !l.starts_with('%')).skip(1) {
        let f: Vec<&str> = line.split_whitespace().collect();
        let feature = features[f[0].parse::<usize>().unwrap() - 1].clone();
        let barcode = barcodes[f[1].parse::<usize>().unwrap() - 1].clone();
        counts.insert((barcode, feature), f[2].parse::<f64>().unwrap());
    }
    (counts, is_real)
}

/// Digest of a MEX directory's content (version metadata line excluded).
pub fn mex_digest(dir: &Path) -> String {
    let mut text = String::new();
    for f in ["matrix.mtx", "barcodes.tsv", "features.tsv"] {
        for line in read_text(&dir.join(f)).lines() {
            if !line.contains("software_version") {
                text.push_str(line);
                text.push('\n');
            }
        }
    }
    digest(&text)
}

/// The run-wide counters from a molecules metrics JSON that define its
/// behaviour (paths and timings excluded).
pub fn molecule_counters(metrics_path: &Path) -> BTreeMap<String, u64> {
    let json: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(metrics_path).unwrap()).unwrap();
    [
        "records",
        "eligible_records",
        "no_sample_records",
        "raw_cell_gene_umis",
        "corrected_umi_nodes",
        "molecules_before_gene_resolution",
        "molecules_counted",
        "umi_nodes_collapsed",
        "collision_groups",
        "collision_molecules_won",
        "collision_molecules_lost",
        "collision_ties",
        "collision_molecules_shared",
        "collision_molecules_fractional",
        "representatives_tagged",
    ]
    .iter()
    .map(|k| (k.to_string(), json[*k].as_u64().unwrap_or(0)))
    .collect()
}

/// Source molecule of a fixture read (its name is `SBX:<molecule>:<copy>`).
pub fn molecule_id(record: &bam::Record) -> u64 {
    let name = String::from_utf8_lossy(record.qname()).to_string();
    name.split(':').nth(1).unwrap().parse().unwrap()
}
