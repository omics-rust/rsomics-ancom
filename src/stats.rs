//! One-way ANOVA p-values (scipy `f_oneway`) and the row-wise multiple-testing
//! corrections (statsmodels `holm` / `fdr_bh`) that ANCOM runs over the pairwise
//! log-ratio p-value matrix.

/// scipy `f_oneway` p-value for a single test: F = MSb/MSw, then the
/// F-distribution survival function `I_{dfw/(dfw+dfb*F)}(dfw/2, dfb/2)`.
/// Zero within-group variance (constant log-ratios) gives `F = nan`, matching
/// scipy; the caller treats a nan p-value as a non-rejection.
pub fn f_oneway_pvalue(groups: &[&[f64]]) -> f64 {
    let k = groups.len();
    let mut n_total = 0usize;
    let mut grand = 0.0_f64;
    for g in groups {
        n_total += g.len();
        grand += g.iter().sum::<f64>();
    }
    grand /= n_total as f64;

    let mut ss_between = 0.0_f64;
    let mut ss_within = 0.0_f64;
    for g in groups {
        let m = g.iter().sum::<f64>() / g.len() as f64;
        ss_between += g.len() as f64 * (m - grand).powi(2);
        for &x in *g {
            ss_within += (x - m).powi(2);
        }
    }

    let df_between = (k - 1) as f64;
    let df_within = (n_total - k) as f64;
    let f = (ss_between / df_between) / (ss_within / df_within);
    if !f.is_finite() {
        return f64::NAN;
    }
    let x = df_within / (df_within + df_between * f);
    betainc(df_within / 2.0, df_between / 2.0, x)
}

/// statsmodels Holm-Bonferroni step-down over a row, returned in original order.
/// nan p-values are sorted last (numpy argsort) and stay nan.
pub fn holm(pvals: &[f64]) -> Vec<f64> {
    let n = pvals.len();
    let order = argsort(pvals);
    let mut acc = f64::NEG_INFINITY;
    let mut sorted_corr = vec![0.0_f64; n];
    for (rank, &idx) in order.iter().enumerate() {
        let raw = pvals[idx] * (n - rank) as f64;
        acc = nanmax(acc, raw);
        sorted_corr[rank] = clip1(acc);
    }
    unsort(&order, &sorted_corr)
}

/// statsmodels Benjamini-Hochberg over a row, returned in original order.
pub fn fdr_bh(pvals: &[f64]) -> Vec<f64> {
    let n = pvals.len();
    let order = argsort(pvals);
    let mut raw = vec![0.0_f64; n];
    for (rank, &idx) in order.iter().enumerate() {
        raw[rank] = pvals[idx] * n as f64 / (rank + 1) as f64;
    }
    let mut acc = f64::INFINITY;
    for rank in (0..n).rev() {
        acc = nanmin(acc, raw[rank]);
        raw[rank] = clip1(acc);
    }
    unsort(&order, &raw)
}

fn argsort(v: &[f64]) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..v.len()).collect();
    idx.sort_by(|&a, &b| match (v[a].is_nan(), v[b].is_nan()) {
        (true, true) => a.cmp(&b),
        (true, false) => std::cmp::Ordering::Greater,
        (false, true) => std::cmp::Ordering::Less,
        (false, false) => v[a].partial_cmp(&v[b]).unwrap().then(a.cmp(&b)),
    });
    idx
}

fn unsort(order: &[usize], sorted_vals: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0_f64; order.len()];
    for (rank, &idx) in order.iter().enumerate() {
        out[idx] = sorted_vals[rank];
    }
    out
}

/// `np.clip(x, 0, 1)` — NaN passes through unchanged.
fn clip1(x: f64) -> f64 {
    if x.is_nan() { x } else { x.clamp(0.0, 1.0) }
}

fn nanmax(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        a.max(b)
    }
}

fn nanmin(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        a.min(b)
    }
}

/// Lanczos approximation, g=7, n=9 — same coefficients SciPy's Cephes uses, so
/// the `betainc` it feeds matches to full f64 precision.
fn ln_gamma(x: f64) -> f64 {
    const G: f64 = 7.0;
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        std::f64::consts::PI.ln() - (std::f64::consts::PI * x).sin().ln() - ln_gamma(1.0 - x)
    } else {
        let x = x - 1.0;
        let mut a = C[0];
        let t = x + G + 0.5;
        for (i, &c) in C.iter().enumerate().skip(1) {
            a += c / (x + i as f64);
        }
        0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
    }
}

/// Regularized incomplete beta `I_x(a, b)` via the Lentz continued fraction
/// (Numerical Recipes `betai`), accurate to ~1e-15 — well past the 1e-6 compat bar.
fn betainc(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let bt = (ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (1.0 - x).ln()).exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        bt * betacf(a, b, x) / a
    } else {
        1.0 - bt * betacf(b, a, 1.0 - x) / b
    }
}

fn betacf(a: f64, b: f64, x: f64) -> f64 {
    const TINY: f64 = 1e-30;
    const EPS: f64 = 3e-16;
    let qab = a + b;
    let qap = a + 1.0;
    let qam = a - 1.0;
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < TINY {
        d = TINY;
    }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..=300 {
        let m = m as f64;
        let m2 = 2.0 * m;
        let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if d.abs() < TINY {
            d = TINY;
        }
        c = 1.0 + aa / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        h *= d * c;
        let aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < TINY {
            d = TINY;
        }
        c = 1.0 + aa / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < EPS {
            break;
        }
    }
    h
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn f_oneway_matches_scipy() {
        // f_oneway([1,2,3,2.5], [5,6,4.5], [2,2.1,1.9,2.2,2.0]) -> p = 0.0001297431658002538
        let a = [1.0, 2.0, 3.0, 2.5];
        let b = [5.0, 6.0, 4.5];
        let c = [2.0, 2.1, 1.9, 2.2, 2.0];
        let p = f_oneway_pvalue(&[&a, &b, &c]);
        assert!((p - 0.000_129_743_165_800_253_8).abs() < 1e-12, "p = {p}");
    }

    #[test]
    fn f_oneway_constant_is_nan() {
        let a = [1.0, 1.0];
        let b = [1.0, 1.0];
        assert!(f_oneway_pvalue(&[&a, &b]).is_nan());
    }

    #[test]
    fn holm_matches_statsmodels() {
        let p = holm(&[0.01, 0.04, 0.03, 0.005]);
        let want = [0.03, 0.06, 0.06, 0.02];
        for (g, w) in p.iter().zip(want) {
            assert!((g - w).abs() < 1e-12, "{g} vs {w}");
        }
    }

    #[test]
    fn holm_preserves_nan() {
        let p = holm(&[0.01, f64::NAN, 0.04, 0.2]);
        let want = [0.04, f64::NAN, 0.12, 0.4];
        for (g, w) in p.iter().zip(want) {
            if w.is_nan() {
                assert!(g.is_nan());
            } else {
                assert!((g - w).abs() < 1e-12, "{g} vs {w}");
            }
        }
    }

    #[test]
    fn holm_diagonal_zero_harmless() {
        let p = holm(&[0.0, 0.02, 0.5, 0.01]);
        let want = [0.0, 0.04, 0.5, 0.03];
        for (g, w) in p.iter().zip(want) {
            assert!((g - w).abs() < 1e-12, "{g} vs {w}");
        }
    }

    #[test]
    fn bh_matches_statsmodels() {
        let p = fdr_bh(&[0.0, 0.02, 0.5, 0.01]);
        let want = [0.0, 0.026_666_666_666_666_67, 0.5, 0.02];
        for (g, w) in p.iter().zip(want) {
            assert!((g - w).abs() < 1e-12, "{g} vs {w}");
        }
    }
}
