use std::collections::HashMap;
use std::path::PathBuf;
use std::process::Command;

fn ours_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-ancom"))
}

fn golden(name: &str) -> String {
    format!("{}/tests/golden/{}", env!("CARGO_MANIFEST_DIR"), name)
}

fn oracle_script() -> String {
    format!("{}/tests/oracle_skbio.py", env!("CARGO_MANIFEST_DIR"))
}

/// feature -> (W, Signif)
fn parse_result(text: &str) -> HashMap<String, (u64, bool)> {
    text.lines()
        .skip(1)
        .filter_map(|l| {
            let mut f = l.split('\t');
            let feat = f.next()?.to_string();
            let w: u64 = f.next()?.trim().parse().ok()?;
            let s = f.next()?.trim() == "True";
            Some((feat, (w, s)))
        })
        .collect()
}

fn ours(table: &str, groups: &str, correction: &str) -> HashMap<String, (u64, bool)> {
    let out = Command::new(ours_bin())
        .arg(golden(table))
        .args(["--grouping", &golden(groups)])
        .args(["--multiple-comparisons-correction", correction])
        .output()
        .expect("run rsomics-ancom");
    assert!(
        out.status.success(),
        "ours failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    parse_result(&String::from_utf8(out.stdout).unwrap())
}

/// Committed skbio-captured W + Signif. Runs with no scikit-bio present — the
/// always-on regression gate.
fn check_committed(table: &str, groups: &str, correction: &str, expected_file: &str) {
    let expected = parse_result(&std::fs::read_to_string(golden(expected_file)).unwrap());
    let got = ours(table, groups, correction);
    assert_eq!(got.len(), expected.len(), "{table} feature count");
    for (feat, &(w, s)) in &expected {
        let &(gw, gs) = got
            .get(feat)
            .unwrap_or_else(|| panic!("missing feature {feat}"));
        assert_eq!(gw, w, "{table} {feat} W");
        assert_eq!(gs, s, "{table} {feat} Signif");
    }
}

#[test]
fn committed_small_holm() {
    check_committed(
        "small_table.tsv",
        "small_groups.tsv",
        "holm-bonferroni",
        "small_expected_holm.tsv",
    );
}

#[test]
fn committed_small_bh() {
    check_committed(
        "small_table.tsv",
        "small_groups.tsv",
        "benjamini-hochberg",
        "small_expected_bh.tsv",
    );
}

#[test]
fn committed_med_holm() {
    check_committed(
        "med_table.tsv",
        "med_groups.tsv",
        "holm-bonferroni",
        "med_expected_holm.tsv",
    );
}

/// scikit-bio is the named oracle; loud-skip if it (or python) is unavailable.
/// `RSOMICS_SKBIO_PYTHON` overrides the interpreter (e.g. an isolated venv).
fn skbio_python() -> Option<String> {
    let mut candidates = Vec::new();
    if let Ok(p) = std::env::var("RSOMICS_SKBIO_PYTHON") {
        candidates.push(p);
    }
    candidates.push("python3".into());
    candidates.push("python".into());
    for py in candidates {
        let probe = Command::new(&py)
            .args(["-c", "import skbio.stats.composition"])
            .output();
        if let Ok(out) = probe
            && out.status.success()
        {
            return Some(py);
        }
    }
    eprintln!("SKIP: scikit-bio not importable — install `scikit-bio` to run the differential");
    None
}

fn oracle(py: &str, table: &str, groups: &str, correction: &str) -> HashMap<String, (u64, bool)> {
    let out = Command::new(py)
        .arg(oracle_script())
        .arg(golden(table))
        .arg(golden(groups))
        .arg(correction)
        .output()
        .expect("run scikit-bio oracle");
    assert!(
        out.status.success(),
        "oracle failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    parse_result(&String::from_utf8(out.stdout).unwrap())
}

/// Live differential: W integer-exact and Signif decision-exact vs skbio.
fn differential(table: &str, groups: &str, our_corr: &str, skbio_corr: &str) {
    let Some(py) = skbio_python() else { return };
    let o = ours(table, groups, our_corr);
    let t = oracle(&py, table, groups, skbio_corr);
    assert_eq!(o.len(), t.len(), "{table} feature count");
    for (feat, &(tw, ts)) in &t {
        let &(ow, os) = o
            .get(feat)
            .unwrap_or_else(|| panic!("missing feature {feat}"));
        assert_eq!(ow, tw, "{table} {feat} W: ours {ow} vs skbio {tw}");
        assert_eq!(os, ts, "{table} {feat} Signif: ours {os} vs skbio {ts}");
    }
}

#[test]
fn matches_skbio_small_holm() {
    differential(
        "small_table.tsv",
        "small_groups.tsv",
        "holm-bonferroni",
        "holm",
    );
}

#[test]
fn matches_skbio_small_bh() {
    differential(
        "small_table.tsv",
        "small_groups.tsv",
        "benjamini-hochberg",
        "fdr_bh",
    );
}

#[test]
fn matches_skbio_med_holm() {
    differential("med_table.tsv", "med_groups.tsv", "holm-bonferroni", "holm");
}
