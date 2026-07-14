//! Microbenchmark: whitelist membership structure for CB classification.
//!
//! Compares the two candidate backends for `Whitelist::contains_code` — a sorted
//! `Vec<u32>` + binary search (the old design) vs an `FxHashSet<u32>` (the current
//! design) — across whitelist sizes from 96 to ~2M, on two query streams:
//!   * HIT-heavy  — queries are on the whitelist (classify returns after 1 probe);
//!   * FAIL-heavy — queries are absent (classify runs the full 48 within-1 probes).
//! The classification logic (exact + 48 single-substitution neighbor probes) is
//! identical to `barcode::Whitelist::classify`, so this isolates the effect of the
//! membership structure alone.
//!
//! Run: `cargo run --release --example bench_whitelist [-- N1 N2 ...]`
//! (defaults to 96, 10000, 736320, 2000000).

use std::hint::black_box;
use std::time::Instant;

use rustc_hash::FxHashSet;

/// SplitMix64 — cheap, high-quality deterministic PRNG for reproducible codes.
struct SplitMix(u64);
impl SplitMix {
    fn next_u32(&mut self) -> u32 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        ((z ^ (z >> 31)) & 0xFFFF_FFFF) as u32
    }
}

/// Classify a code exactly as `Whitelist::classify` does: exact hit → 0, within
/// Hamming-1 → 1, else fail → 2. `contains` is the membership test under study.
#[inline]
fn classify<F: Fn(u32) -> bool>(code: u32, contains: &F) -> u8 {
    if contains(code) {
        return 0;
    }
    for pos in 0..16 {
        let shift = pos * 2;
        let cur = (code >> shift) & 0b11;
        for alt in 0u32..4 {
            if alt != cur {
                let neighbor = (code & !(0b11 << shift)) | (alt << shift);
                if contains(neighbor) {
                    return 1;
                }
            }
        }
    }
    2
}

/// Time `classify` over `queries` with membership `contains`; returns ns/query.
/// Reports the min over a few passes to suppress scheduler noise.
fn time_stream<F: Fn(u32) -> bool>(queries: &[u32], contains: &F) -> f64 {
    const PASSES: usize = 5;
    // Warm up caches/branch predictors.
    let mut warm = 0u64;
    for &q in queries {
        warm += classify(q, contains) as u64;
    }
    black_box(warm);

    let mut best = f64::INFINITY;
    for _ in 0..PASSES {
        let t = Instant::now();
        let mut acc = 0u64;
        for &q in queries {
            acc += classify(q, contains) as u64;
        }
        black_box(acc);
        let ns = t.elapsed().as_nanos() as f64 / queries.len() as f64;
        best = best.min(ns);
    }
    best
}

fn main() {
    let sizes: Vec<usize> = {
        let args: Vec<String> = std::env::args().skip(1).collect();
        if args.is_empty() {
            vec![96, 10_000, 736_320, 2_000_000]
        } else {
            args.iter()
                .map(|s| s.parse().expect("N must be an integer"))
                .collect()
        }
    };
    const QUERIES: usize = 300_000;

    println!(
        "{:>10}  {:>6}  {:<10}  {:>12}  {:>12}  {:>8}",
        "N", "stream", "backend", "ns/query", "Mquery/s", "speedup"
    );
    println!("{}", "-".repeat(74));

    for &n in &sizes {
        // Build the whitelist: n distinct random codes.
        let mut rng = SplitMix(0xABBA_1234_5678_9ABC ^ n as u64);
        let mut set: FxHashSet<u32> = FxHashSet::default();
        set.reserve(n);
        while set.len() < n {
            set.insert(rng.next_u32());
        }
        let mut vec: Vec<u32> = set.iter().copied().collect();
        vec.sort_unstable();

        // HIT-heavy queries: sample existing members.
        let hits: Vec<u32> = (0..QUERIES)
            .map(|i| vec[(i * 2654435761) % vec.len()])
            .collect();
        // FAIL-heavy queries: codes guaranteed absent from the whitelist.
        let mut fails: Vec<u32> = Vec::with_capacity(QUERIES);
        while fails.len() < QUERIES {
            let c = rng.next_u32();
            if !set.contains(&c) {
                fails.push(c);
            }
        }

        let vec_contains = |c: u32| vec.binary_search(&c).is_ok();
        let set_contains = |c: u32| set.contains(&c);

        let vec_bytes = vec.len() * std::mem::size_of::<u32>();
        // hashbrown: ~1 control byte + 4 key bytes per slot at ~7/8 load.
        let set_bytes = set.capacity() * (std::mem::size_of::<u32>() + 1);
        println!(
            "N={n}  (Vec {:.1} MB  |  FxHashSet ~{:.1} MB, cap {})",
            vec_bytes as f64 / 1e6,
            set_bytes as f64 / 1e6,
            set.capacity()
        );

        for (label, queries) in [("hit", &hits), ("fail", &fails)] {
            let v = time_stream(queries, &vec_contains);
            let s = time_stream(queries, &set_contains);
            let row = |backend: &str, ns: f64, base: f64| {
                println!(
                    "{:>10}  {:>6}  {:<10}  {:>12.2}  {:>12.2}  {:>7.2}x",
                    n,
                    label,
                    backend,
                    ns,
                    1000.0 / ns,
                    base / ns
                );
            };
            row("sorted-vec", v, v);
            row("fxhashset", s, v); // speedup vs the sorted-vec baseline
        }
        println!();
    }
}
