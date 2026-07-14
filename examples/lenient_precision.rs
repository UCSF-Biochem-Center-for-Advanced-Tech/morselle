//! Dev tool: precision + recall of the lenient within-2 CB rescue.
//!
//! Simulates known whitelist barcodes carrying *exactly* N substitution errors,
//! runs `Whitelist::correct_within2`, and classifies the outcome:
//!   correct  — returned the true source barcode (recall)
//!   wrong    — returned a DIFFERENT whitelist barcode (precision leak)
//!   abstain  — None (ambiguous / no unique within-2 hit)
//! A genuine ≤2-error barcode can only be correct-or-abstain (the true barcode is
//! always within 2), so `wrong` is ~0 for N≤2. `wrong` at N≥3 is the real leak: a
//! >2-error (or off-whitelist) CB masquerading as within-2 of the wrong entry.
//! Weight those against the confident distance distribution (cb_distance) to get
//! the practical false-assignment rate of enabling `--lenient`.
//!
//! Usage: cargo run --release --example lenient_precision -- <whitelist> [samples_per_level]

use sbx10x::barcode::{decode16, Whitelist};

/// SplitMix64 deterministic PRNG.
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

fn main() -> anyhow::Result<()> {
    let mut a = std::env::args().skip(1);
    let wlp = a
        .next()
        .expect("usage: lenient_precision <whitelist> [samples]");
    let n: usize = a.next().map(|s| s.parse().unwrap()).unwrap_or(300_000);
    let wl = Whitelist::load(std::path::Path::new(&wlp))?;
    let codes = wl.codes_sorted();
    let bases = b"ACGT";
    let mut rng = Rng(0xC0FFEE_1234_5678);

    println!(
        "whitelist {} barcodes; {} samples per error level\n",
        codes.len(),
        n
    );
    println!(
        "{:>7} {:>10} {:>10} {:>10}   interpretation",
        "N errs", "correct%", "wrong%", "abstain%"
    );
    println!("{}", "-".repeat(70));
    let mut wrong_by_n = [0f64; 6];
    let mut correct_by_n = [0f64; 6];
    for nerr in 1..=4usize {
        let (mut correct, mut wrong, mut abstain) = (0u64, 0u64, 0u64);
        for _ in 0..n {
            let truth = decode16(codes[rng.upto(codes.len())]);
            let mut cb = truth;
            let mut used = [false; 16];
            let mut applied = 0;
            while applied < nerr {
                let p = rng.upto(16);
                if used[p] {
                    continue;
                }
                used[p] = true;
                let cur = cb[p];
                loop {
                    let b = bases[rng.upto(4)];
                    if b != cur {
                        cb[p] = b;
                        break;
                    }
                }
                applied += 1;
            }
            match wl.correct_within2(&cb) {
                Some(res) if res == truth => correct += 1,
                Some(_) => wrong += 1,
                None => abstain += 1,
            }
        }
        let tot = n as f64;
        let note = match nerr {
            1 | 2 => "true CB within-2 → recover or abstain",
            _ => "true CB NOT within-2 → any hit is WRONG (leak)",
        };
        println!(
            "{:>7} {:>9.2}% {:>9.2}% {:>9.2}%   {}",
            nerr,
            100.0 * correct as f64 / tot,
            100.0 * wrong as f64 / tot,
            100.0 * abstain as f64 / tot,
            note
        );
        wrong_by_n[nerr] = wrong as f64 / tot;
        correct_by_n[nerr] = correct as f64 / tot;
    }

    // Practical precision estimate: weight the per-N correct/wrong RATES by the
    // observed confident-location distance distribution (cb_distance on run_5p).
    // Lenient only fires on d>=2 (d0/d1 already assigned). A read at distance d is a
    // real barcode with ~d errors, so it recovers/leaks at the N=d simulated rate
    // (d>=4 capped at N=4). This is the honest precision: it does NOT assume d2 is
    // fully recovered — on this dense whitelist within-2 mostly abstains.
    let conf: [(usize, f64); 4] = [(2, 0.041), (3, 0.039), (4, 0.025), (5, 0.002)];
    let eligible: f64 = conf.iter().map(|(_, f)| f).sum();
    let mut good = 0.0;
    let mut bad = 0.0;
    for (d, f) in conf {
        let n = d.min(4);
        good += f * correct_by_n[n];
        bad += f * wrong_by_n[n];
    }
    println!(
        "\nApplied to run_5p confident CB-fails (lenient fires on d>=2 = {:.1}% of confident CBs):",
        100.0 * eligible
    );
    println!(
        "  correct assignments (recall)   = {:.3}% of confident CBs",
        100.0 * good
    );
    println!(
        "  WRONG assignments  (leak)      = {:.3}% of confident CBs",
        100.0 * bad
    );
    println!(
        "  => rescue precision correct/(correct+wrong) ~ {:.1}%",
        100.0 * good / (good + bad).max(1e-9)
    );
    Ok(())
}
