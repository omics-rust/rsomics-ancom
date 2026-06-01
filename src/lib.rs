use std::fmt::Write as _;
use std::io::{BufRead, Write};

use rayon::prelude::*;
use rsomics_common::{Result, RsomicsError};

mod stats;
mod table;

pub use table::{Table, read_grouping};

#[derive(Clone, Copy)]
pub enum Correction {
    Holm,
    BenjaminiHochberg,
    None,
}

impl Correction {
    pub fn parse(s: &str) -> Result<Self> {
        match s.to_lowercase().as_str() {
            "holm" | "holm-bonferroni" => Ok(Self::Holm),
            "bh" | "fdr_bh" | "benjamini-hochberg" => Ok(Self::BenjaminiHochberg),
            "none" => Ok(Self::None),
            other => Err(RsomicsError::InvalidInput(format!(
                "unknown correction '{other}' (holm | benjamini-hochberg | none)"
            ))),
        }
    }
}

pub struct AncomResult {
    pub features: Vec<String>,
    pub w: Vec<usize>,
    pub signif: Vec<bool>,
}

/// ANCOM (Mandal et al. 2015), matching `skbio.stats.composition.ancom`.
/// For each feature, W counts how many of the pairwise log-ratio tests against
/// the other features reject the null at `alpha` (after the row-wise multiple-
/// testing correction); features whose W clears the tau/theta detection cutoff
/// are flagged significant.
///
/// # Errors
/// Non-positive counts (log undefined), out-of-range alpha/tau/theta, or a
/// grouping that is all-distinct or single-group — all rejected as in skbio.
pub fn ancom(
    table: &Table,
    grouping: &[String],
    alpha: f64,
    tau: f64,
    theta: f64,
    correction: Correction,
) -> Result<AncomResult> {
    if grouping.len() != table.n_samples() {
        return Err(RsomicsError::InvalidInput(format!(
            "grouping has {} entries but the table has {} samples",
            grouping.len(),
            table.n_samples()
        )));
    }
    for &v in &table.data {
        if v <= 0.0 {
            return Err(RsomicsError::InvalidInput(
                "table contains zero or negative values — add a pseudocount first (log undefined)"
                    .into(),
            ));
        }
    }
    for (name, val) in [("alpha", alpha), ("tau", tau), ("theta", theta)] {
        if !(0.0 < val && val < 1.0) {
            return Err(RsomicsError::InvalidInput(format!(
                "`{name}`={val} is not within 0 and 1"
            )));
        }
    }

    let labels = factorize(grouping);
    let n_groups = labels.iter().copied().max().map_or(0, |m| m + 1);
    if n_groups == grouping.len() {
        return Err(RsomicsError::InvalidInput(
            "all grouping values are unique — no within-group variance".into(),
        ));
    }
    if n_groups == 1 {
        return Err(RsomicsError::InvalidInput(
            "all grouping values are the same — no between-group variance".into(),
        ));
    }

    let m = table.n_features();

    // log matrix grouped by label: log_by_group[g] is a flat (n_g × m) block,
    // so each pair's per-group log-ratio vectors are gathered with stride-m
    // reads instead of chasing the row-major table.
    let group_sizes: Vec<usize> = (0..n_groups)
        .map(|g| labels.iter().filter(|&&l| l == g).count())
        .collect();
    let mut log_by_group: Vec<Vec<f64>> = group_sizes
        .iter()
        .map(|&sz| Vec::with_capacity(sz * m))
        .collect();
    for (i, &g) in labels.iter().enumerate() {
        log_by_group[g].extend(table.row(i).iter().map(|&x| x.ln()));
    }

    // pairwise p-value matrix, row-major m×m, diagonal 0 (skbio convention)
    let mut pmat = vec![0.0_f64; m * m];
    let upper: Vec<(usize, usize)> = (0..m)
        .flat_map(|i| ((i + 1)..m).map(move |j| (i, j)))
        .collect();

    let pvals: Vec<f64> = upper
        .par_iter()
        .map(|&(i, j)| {
            let mut per_group: Vec<Vec<f64>> = Vec::with_capacity(n_groups);
            for g in 0..n_groups {
                let block = &log_by_group[g];
                let ng = group_sizes[g];
                let mut ratios = Vec::with_capacity(ng);
                for r in 0..ng {
                    ratios.push(block[r * m + i] - block[r * m + j]);
                }
                per_group.push(ratios);
            }
            let refs: Vec<&[f64]> = per_group.iter().map(|v| v.as_slice()).collect();
            stats::f_oneway_pvalue(&refs)
        })
        .collect();

    for (&(i, j), &p) in upper.iter().zip(&pvals) {
        pmat[i * m + j] = p;
        pmat[j * m + i] = p;
    }

    // correct each row (the 0 diagonal participates, as in skbio), then 1-diag
    match correction {
        Correction::Holm => apply_rows(&mut pmat, m, stats::holm),
        Correction::BenjaminiHochberg => apply_rows(&mut pmat, m, stats::fdr_bh),
        Correction::None => {}
    }
    for i in 0..m {
        pmat[i * m + i] = 1.0;
    }

    let w: Vec<usize> = (0..m)
        .map(|i| (0..m).filter(|&j| pmat[i * m + j] < alpha).count())
        .collect();

    let signif = detect(&w, m, tau, theta);

    Ok(AncomResult {
        features: table.features.clone(),
        w,
        signif,
    })
}

fn apply_rows(pmat: &mut [f64], m: usize, f: fn(&[f64]) -> Vec<f64>) {
    for i in 0..m {
        let corrected = f(&pmat[i * m..(i + 1) * m]);
        pmat[i * m..(i + 1) * m].copy_from_slice(&corrected);
    }
}

/// The ANCOM detection rule (skbio): pick a W cutoff from the staircase of
/// proportions exceeding `c_start - {0.05..0.25}`, governed by tau, then flag
/// features whose W clears it. Below theta nothing is detected.
fn detect(w: &[usize], n_feats: usize, tau: f64, theta: f64) -> Vec<bool> {
    let w_max = *w.iter().max().unwrap() as f64;
    let c_start = w_max / n_feats as f64;
    if c_start < theta {
        return vec![false; w.len()];
    }
    let cutoff: [f64; 5] = std::array::from_fn(|k| c_start - (0.05 + 0.05 * k as f64));

    // prop_cut[k] = fraction of features with W > n_feats * cutoff[k]
    let prop_cut: [f64; 5] = std::array::from_fn(|k| {
        let thr = n_feats as f64 * cutoff[k];
        w.iter().filter(|&&wi| wi as f64 > thr).count() as f64 / w.len() as f64
    });

    // dels[k] = |prop_cut[k] - prop_cut[k+1]| (np.roll(-1)); dels[-1] forced 0
    let mut dels = [0.0_f64; 5];
    for k in 0..4 {
        dels[k] = (prop_cut[k] - prop_cut[k + 1]).abs();
    }

    let nu = if dels[0] < tau && dels[1] < tau && dels[2] < tau {
        cutoff[1]
    } else if dels[0] >= tau && dels[1] < tau && dels[2] < tau {
        cutoff[2]
    } else if dels[1] >= tau && dels[2] < tau && dels[3] < tau {
        cutoff[3]
    } else {
        cutoff[4]
    };

    let thr = nu * n_feats as f64;
    w.iter().map(|&wi| wi as f64 >= thr).collect()
}

/// Dense integer codes in `np.unique` order (sorted unique labels). Group order
/// is irrelevant to the symmetric f_oneway p-values; it would only matter for
/// percentile columns, which this tool does not emit.
fn factorize(grouping: &[String]) -> Vec<usize> {
    let mut order: Vec<&String> = grouping.iter().collect();
    order.sort_unstable();
    order.dedup();
    grouping
        .iter()
        .map(|g| order.binary_search(&g).unwrap())
        .collect()
}

impl AncomResult {
    /// Write the skbio result table: a `feature\tW\tSignif` header then one row
    /// per feature. W is integer-exact; Signif is `True`/`False` (skbio column).
    ///
    /// # Errors
    /// Propagates write errors.
    pub fn write_tsv<W: Write>(&self, mut out: W) -> Result<()> {
        writeln!(out, "feature\tW\tSignif").map_err(RsomicsError::Io)?;
        let mut line = String::new();
        for ((feat, &w), &s) in self.features.iter().zip(&self.w).zip(&self.signif) {
            line.clear();
            line.push_str(feat);
            line.push('\t');
            let _ = write!(line, "{w}");
            line.push('\t');
            line.push_str(if s { "True" } else { "False" });
            writeln!(out, "{line}").map_err(RsomicsError::Io)?;
        }
        Ok(())
    }
}

/// # Errors
/// Propagates parse, grouping, compute, and write errors.
#[allow(clippy::too_many_arguments)]
pub fn run<W: Write>(
    table_reader: impl BufRead,
    grouping_reader: impl BufRead,
    out: W,
    delim: char,
    alpha: f64,
    tau: f64,
    theta: f64,
    correction: Correction,
) -> Result<()> {
    let table = Table::parse(table_reader, delim)?;
    let grouping = read_grouping(grouping_reader, &table.samples, delim)?;
    let res = ancom(&table, &grouping, alpha, tau, theta, correction)?;
    res.write_tsv(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn doc_table() -> Table {
        // skbio docstring example: 6 samples, 7 features
        let txt = "\tb1\tb2\tb3\tb4\tb5\tb6\tb7\n\
            s1\t12\t11\t10\t10\t10\t10\t10\n\
            s2\t9\t11\t12\t10\t10\t10\t10\n\
            s3\t1\t11\t10\t11\t10\t5\t9\n\
            s4\t22\t21\t9\t10\t10\t10\t10\n\
            s5\t20\t22\t10\t10\t13\t10\t10\n\
            s6\t23\t21\t14\t10\t10\t10\t10\n";
        Table::parse(txt.as_bytes(), '\t').unwrap()
    }

    fn doc_grouping() -> Vec<String> {
        [
            "treatment",
            "treatment",
            "treatment",
            "placebo",
            "placebo",
            "placebo",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect()
    }

    #[test]
    fn doc_example_w_and_signif() {
        let res = ancom(
            &doc_table(),
            &doc_grouping(),
            0.05,
            0.02,
            0.1,
            Correction::Holm,
        )
        .unwrap();
        assert_eq!(res.w, vec![0, 4, 0, 1, 1, 0, 1]);
        assert_eq!(
            res.signif,
            vec![false, true, false, false, false, false, false]
        );
    }

    #[test]
    fn zero_count_errors() {
        let txt = "\tf1\tf2\ns1\t0\t5\ns2\t3\t2\n";
        let t = Table::parse(txt.as_bytes(), '\t').unwrap();
        let g = vec!["a".to_string(), "b".to_string()];
        assert!(ancom(&t, &g, 0.05, 0.02, 0.1, Correction::Holm).is_err());
    }

    #[test]
    fn single_group_errors() {
        let g = vec!["x".to_string(); 6];
        assert!(ancom(&doc_table(), &g, 0.05, 0.02, 0.1, Correction::Holm).is_err());
    }
}
