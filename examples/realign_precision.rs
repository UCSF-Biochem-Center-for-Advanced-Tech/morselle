//! Dev tool: precision + recall of the single-indel CB realignment (`--realign`).
//!
//! Builds synthetic CB regions from known whitelist barcodes carrying a single
//! insertion or deletion (+ optional extra substitution), followed by a random
//! UMI, and runs `Whitelist::realign_cb`. Classifies: correct / wrong / abstain.
//! A genuine single-indel barcode should be correct-or-abstain (never wrong); the
//! leak case (indel + substitution) tests false assignment when the read isn't a
//! clean single-indel.
//!
//! Usage: cargo run --release --example realign_precision -- <whitelist> [samples]

use morselle::barcode::{decode16, Whitelist};

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn upto(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}
const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

fn main() -> anyhow::Result<()> {
    let mut a = std::env::args().skip(1);
    let wlp = a
        .next()
        .expect("usage: realign_precision <whitelist> [samples]");
    let n: usize = a.next().map(|s| s.parse().unwrap()).unwrap_or(300_000);
    let wl = Whitelist::load(std::path::Path::new(&wlp))?;
    let codes = wl.codes_sorted();
    let mut rng = Rng(0x5EED_1234_ABCD);

    let umi = |rng: &mut Rng| -> Vec<u8> { (0..12).map(|_| BASES[rng.upto(4)]).collect() };

    println!(
        "whitelist {} barcodes; {} samples per scenario\n",
        codes.len(),
        n
    );
    println!(
        "{:<34} {:>9} {:>9} {:>9}",
        "scenario (delta passed)", "correct%", "wrong%", "abstain%"
    );
    println!("{}", "-".repeat(64));

    let mut run = |label: &str, delta: i32, mutate: &dyn Fn(&mut Rng, &[u8; 16]) -> Vec<u8>| {
        let (mut correct, mut wrong, mut abstain) = (0u64, 0u64, 0u64);
        for _ in 0..n {
            let truth = decode16(codes[rng.upto(codes.len())]);
            let mut region = mutate(&mut rng, &truth);
            region.extend(umi(&mut rng)); // random UMI after the CB region
            match wl.realign_cb(&region, delta) {
                Some(res) if res == truth => correct += 1,
                Some(_) => wrong += 1,
                None => abstain += 1,
            }
        }
        let t = n as f64;
        println!(
            "{:<34} {:>8.2}% {:>8.2}% {:>8.2}%",
            label,
            100.0 * correct as f64 / t,
            100.0 * wrong as f64 / t,
            100.0 * abstain as f64 / t
        );
    };

    // 1 insertion (delta +1): insert a random base at a random position → 17 bp CB.
    run("1 insertion            (Δ=+1)", 1, &|rng, truth| {
        let pos = rng.upto(17);
        let mut v = Vec::with_capacity(17);
        v.extend_from_slice(&truth[..pos.min(16)]);
        v.insert(pos.min(v.len()), BASES[rng.upto(4)]);
        // rebuild cleanly: insert into full 16
        let mut out = truth.to_vec();
        out.insert(rng.upto(17).min(16), BASES[rng.upto(4)]);
        let _ = v;
        out
    });
    // 1 deletion (delta -1): drop one base → 15 bp CB.
    run("1 deletion             (Δ=-1)", -1, &|rng, truth| {
        let mut out = truth.to_vec();
        out.remove(rng.upto(16));
        out
    });
    // LEAK: 1 insertion + 1 substitution (still span Δ=+1, but not a clean indel).
    run("1 insertion + 1 subst  (Δ=+1)", 1, &|rng, truth| {
        let mut out = truth.to_vec();
        let sp = rng.upto(16);
        let cur = out[sp];
        loop {
            let b = BASES[rng.upto(4)];
            if b != cur {
                out[sp] = b;
                break;
            }
        }
        out.insert(rng.upto(17).min(16), BASES[rng.upto(4)]);
        out
    });
    // LEAK: 1 deletion + 1 substitution (span Δ=-1).
    run("1 deletion + 1 subst   (Δ=-1)", -1, &|rng, truth| {
        let mut out = truth.to_vec();
        let sp = rng.upto(16);
        let cur = out[sp];
        loop {
            let b = BASES[rng.upto(4)];
            if b != cur {
                out[sp] = b;
                break;
            }
        }
        out.remove(rng.upto(15));
        out
    });
    Ok(())
}
