use statrs::distribution::{Continuous, Normal};
use super::types::{Component, BOUNDARY_SCAN_STEPS, INIT_PERCENTILES, MIN_STD};

pub fn fit_gmm(x: &[f64]) -> (Component, Component) {
    let (km_low, km_high) = kmeans1d_init(x);
    std::iter::once(em_once_from(x, km_low, km_high))
        .chain(
            INIT_PERCENTILES
                .iter()
                .map(|&(lo_pct, hi_pct)| em_once(x, lo_pct, hi_pct)),
        )
        .max_by(|a, b| a.2.partial_cmp(&b.2).unwrap())
        .map(|(low, high, _ll)| {
            if low.mean <= high.mean {
                (low, high)
            } else {
                (high, low)
            }
        })
        .unwrap()
}

pub fn kmeans1d_init(x: &[f64]) -> (f64, f64) {
    let (mut c0, mut c1) = (
        x.iter().cloned().fold(f64::INFINITY, f64::min),
        x.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
    );
    for _ in 0..20 {
        let (mut sum0, mut n0, mut sum1, mut n1) = (0.0, 0.0, 0.0, 0.0);
        for &xi in x {
            if (xi - c0).abs() <= (xi - c1).abs() {
                sum0 += xi;
                n0 += 1.0;
            } else {
                sum1 += xi;
                n1 += 1.0;
            }
        }
        if n0 > 0.0 {
            c0 = sum0 / n0;
        }
        if n1 > 0.0 {
            c1 = sum1 / n1;
        }
    }
    (c0.min(c1), c0.max(c1))
}

pub fn em_once(x: &[f64], low_pct: f64, high_pct: f64) -> (Component, Component, f64) {
    let mut sorted = x.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let pct = |p: f64| sorted[((sorted.len() - 1) as f64 * p).round() as usize];
    em_once_from(x, pct(low_pct), pct(high_pct))
}

pub fn em_once_from(x: &[f64], init_low_mean: f64, init_high_mean: f64) -> (Component, Component, f64) {
    let n = x.len() as f64;
    let mean_all = x.iter().sum::<f64>() / n;
    let std_all = (x.iter().map(|v| (v - mean_all).powi(2)).sum::<f64>() / n)
        .sqrt()
        .max(MIN_STD);

    let mut low = Component {
        weight: 0.5,
        mean: init_low_mean,
        std: std_all.max(MIN_STD),
    };
    let mut high = Component {
        weight: 0.5,
        mean: init_high_mean,
        std: std_all.max(MIN_STD),
    };
    if (high.mean - low.mean).abs() < MIN_STD {
        high.mean = low.mean + 1.0;
    }

    let mut prev_ll = f64::NEG_INFINITY;
    let mut ll = prev_ll;
    for _ in 0..200 {
        let dist_low = Normal::new(low.mean, low.std.max(MIN_STD)).unwrap();
        let dist_high = Normal::new(high.mean, high.std.max(MIN_STD)).unwrap();

        let mut r_high = Vec::with_capacity(x.len());
        ll = 0.0;
        for &xi in x {
            let pl = low.weight * dist_low.pdf(xi);
            let ph = high.weight * dist_high.pdf(xi);
            let total = (pl + ph).max(f64::MIN_POSITIVE);
            r_high.push(ph / total);
            ll += total.ln();
        }
        if (ll - prev_ll).abs() < 1e-8 {
            break;
        }
        prev_ll = ll;

        let sum_rh: f64 = r_high.iter().sum();
        let sum_rl: f64 = n - sum_rh;
        if sum_rh < 1e-6 || sum_rl < 1e-6 {
            break;
        }
        let mean_h = x.iter().zip(&r_high).map(|(xi, r)| xi * r).sum::<f64>() / sum_rh;
        let mean_l = x
            .iter()
            .zip(&r_high)
            .map(|(xi, r)| xi * (1.0 - r))
            .sum::<f64>()
            / sum_rl;
        let var_h = x
            .iter()
            .zip(&r_high)
            .map(|(xi, r)| r * (xi - mean_h).powi(2))
            .sum::<f64>()
            / sum_rh;
        let var_l = x
            .iter()
            .zip(&r_high)
            .map(|(xi, r)| (1.0 - r) * (xi - mean_l).powi(2))
            .sum::<f64>()
            / sum_rl;

        low = Component {
            weight: sum_rl / n,
            mean: mean_l,
            std: var_l.sqrt().max(MIN_STD),
        };
        high = Component {
            weight: sum_rh / n,
            mean: mean_h,
            std: var_h.sqrt().max(MIN_STD),
        };
    }

    (low, high, ll)
}

pub fn decision_boundary(low: &Component, high: &Component) -> f64 {
    let midpoint = (low.mean + high.mean) / 2.0;
    if high.mean <= low.mean {
        return midpoint;
    }
    let dist_low = Normal::new(low.mean, low.std.max(MIN_STD)).unwrap();
    let dist_high = Normal::new(high.mean, high.std.max(MIN_STD)).unwrap();
    let diff_at =
        |x: f64| high.weight.ln() + dist_high.ln_pdf(x) - low.weight.ln() - dist_low.ln_pdf(x);

    let step = (high.mean - low.mean) / BOUNDARY_SCAN_STEPS as f64;
    let mut prev_x = low.mean;
    let mut prev_diff = diff_at(prev_x);
    for i in 1..=BOUNDARY_SCAN_STEPS {
        let x = low.mean + step * i as f64;
        let diff = diff_at(x);
        if prev_diff.signum() != diff.signum() {
            return (prev_x + x) / 2.0;
        }
        prev_x = x;
        prev_diff = diff;
    }
    midpoint
}

pub fn fit_threshold(counts: &[u64], min_umi: u64) -> Option<u64> {
    let nonzero: Vec<f64> = counts
        .iter()
        .filter(|&&c| c > 0)
        .map(|&c| (c as f64).ln_1p())
        .collect();
    if nonzero.is_empty() {
        return None;
    }
    if nonzero.len() < 2 || nonzero.iter().all(|&v| (v - nonzero[0]).abs() < 1e-9) {
        return Some(min_umi);
    }
    let (low, high) = fit_gmm(&nonzero);
    let boundary = decision_boundary(&low, &high);
    let fitted = boundary.exp_m1().ceil().max(0.0) as u64;
    Some(fitted.max(min_umi))
}
