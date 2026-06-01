# rsomics-ancom

ANCOM — Analysis of Composition of Microbiomes (Mandal et al. 2015) — as a single
fast CLI. Equivalent to `skbio.stats.composition.ancom`.

For each feature, ANCOM performs a one-way ANOVA on the log-ratio against every
other feature and counts how many of those tests reject the null after a row-wise
multiple-comparisons correction. That count is the **W** statistic. Features
whose W clears the tau/theta detection cutoff are flagged significant.

```
rsomics-ancom table.tsv --grouping groups.tsv [--alpha 0.05] \
    [--multiple-comparisons-correction holm-bonferroni] [-o result.tsv]
```

- `table.tsv` — feature table: header row of feature IDs (corner cell ignored),
  then one `sample_id  count...` line per sample. All counts must be positive
  (add a pseudocount first — the log of zero is undefined).
- `--grouping` — `sample_id<TAB>group` per sample.

Output is `feature<TAB>W<TAB>Signif`, one row per feature.

## Origin

This crate is an independent Rust reimplementation of `skbio.stats.composition.ancom`
based on:

- Mandal et al., "Analysis of composition of microbiomes: a novel method for
  studying microbial composition", *Microbial Ecology in Health & Disease* (2015),
  26. DOI: 10.3402/mehd.v26.27663
- The scikit-bio implementation (Modified BSD License), read and cited: the
  default significance test is SciPy's one-way ANOVA `f_oneway`, the default
  multiple-comparisons correction is statsmodels Holm-Bonferroni, and the W
  detection cutoff follows the tau/theta staircase rule in `ancom()`.

The pairwise one-way ANOVA p-value is `I_x(dfw/2, dfb/2)` with
`x = dfw / (dfw + dfb·F)` (scipy `fdtrc`); the Holm / Benjamini-Hochberg
corrections reproduce statsmodels `multipletests` step-down/step-up exactly,
including its NaN handling for constant-ratio (zero-variance) pairs. W is
integer-exact and the `Signif` decisions are identical to scikit-bio; corrected
p-values match to ~1e-9.

License: MIT OR Apache-2.0.
Upstream credit: scikit-bio https://scikit-bio.org/ (Modified BSD License).
