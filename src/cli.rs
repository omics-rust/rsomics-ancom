use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_ancom::{Correction, run};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-ancom", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// Feature table TSV (samples × features counts, corner cell ignored); reads stdin when "-" or omitted.
    #[arg(default_value = "-")]
    input: PathBuf,

    /// Grouping TSV: one `id<TAB>group` line per sample.
    #[arg(long)]
    grouping: PathBuf,

    /// Significance level for each pairwise test.
    #[arg(long, default_value_t = 0.05)]
    alpha: f64,

    /// Cutoff-selection constant; smaller is more conservative.
    #[arg(long, default_value_t = 0.02)]
    tau: f64,

    /// Lower bound on the W proportion for any feature to be detected.
    #[arg(long, default_value_t = 0.1)]
    theta: f64,

    /// Multiple-comparisons correction: holm-bonferroni | benjamini-hochberg | none.
    #[arg(long, default_value = "holm-bonferroni")]
    multiple_comparisons_correction: String,

    /// Parse inputs as comma-separated instead of tab-separated.
    #[arg(long, default_value_t = false)]
    csv: bool,

    /// Output path; writes stdout when "-".
    #[arg(short = 'o', long, default_value = "-")]
    output: String,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        self.common.install_rayon_pool()?;
        let delim = if self.csv { ',' } else { '\t' };
        let correction = Correction::parse(&self.multiple_comparisons_correction)?;

        let table_reader: Box<dyn std::io::BufRead> = if self.input.as_os_str() == "-" {
            Box::new(BufReader::new(std::io::stdin().lock()))
        } else {
            Box::new(BufReader::new(File::open(&self.input).map_err(|e| {
                RsomicsError::InvalidInput(format!("{}: {e}", self.input.display()))
            })?))
        };
        let grouping_reader = BufReader::new(File::open(&self.grouping).map_err(|e| {
            RsomicsError::InvalidInput(format!("{}: {e}", self.grouping.display()))
        })?);
        let mut out: Box<dyn Write> = if self.output == "-" && self.common.json {
            Box::new(std::io::sink())
        } else if self.output == "-" {
            Box::new(BufWriter::new(std::io::stdout().lock()))
        } else {
            Box::new(BufWriter::new(
                File::create(&self.output).map_err(RsomicsError::Io)?,
            ))
        };

        run(
            table_reader,
            grouping_reader,
            &mut out,
            delim,
            self.alpha,
            self.tau,
            self.theta,
            correction,
        )?;
        out.flush().map_err(RsomicsError::Io)
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "ANCOM differential abundance test from a feature table and grouping.",
    origin: Some(Origin {
        upstream: "scikit-bio skbio.stats.composition.ancom",
        upstream_license: "BSD-3-Clause",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.3402/mehd.v26.27663"),
    }),
    usage_lines: &[
        "[table.tsv] --grouping groups.tsv [--alpha 0.05] [--multiple-comparisons-correction holm-bonferroni] [-o result.tsv]",
    ],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: None,
                long: "grouping",
                aliases: &[],
                value: Some("<path>"),
                type_hint: None,
                required: true,
                default: None,
                description: "Grouping TSV (id<TAB>group per sample).",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "alpha",
                aliases: &[],
                value: Some("<float>"),
                type_hint: Some("f64"),
                required: false,
                default: Some("0.05"),
                description: "Significance level for each pairwise test.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "tau",
                aliases: &[],
                value: Some("<float>"),
                type_hint: Some("f64"),
                required: false,
                default: Some("0.02"),
                description: "Cutoff-selection constant (smaller = conservative).",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "theta",
                aliases: &[],
                value: Some("<float>"),
                type_hint: Some("f64"),
                required: false,
                default: Some("0.1"),
                description: "Lower bound on the W proportion for detection.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "multiple-comparisons-correction",
                aliases: &[],
                value: Some("<method>"),
                type_hint: None,
                required: false,
                default: Some("holm-bonferroni"),
                description: "holm-bonferroni | benjamini-hochberg | none.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "csv",
                aliases: &[],
                value: None,
                type_hint: None,
                required: false,
                default: Some("false"),
                description: "Parse inputs as comma-separated.",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("String"),
                required: false,
                default: Some("-"),
                description: "Output path (- for stdout).",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "ANCOM with default Holm correction",
            command: "rsomics-ancom table.tsv --grouping groups.tsv",
        },
        Example {
            description: "Benjamini-Hochberg correction, write to a file",
            command: "rsomics-ancom table.tsv --grouping groups.tsv --multiple-comparisons-correction benjamini-hochberg -o ancom.tsv",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
