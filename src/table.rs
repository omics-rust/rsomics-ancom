use std::io::BufRead;

use rsomics_common::{Result, RsomicsError};

/// A samples × features count matrix with row (sample) and column (feature) IDs,
/// stored row-major. Mirrors what skbio ingests from a DataFrame.
pub struct Table {
    pub samples: Vec<String>,
    pub features: Vec<String>,
    pub data: Vec<f64>,
}

impl Table {
    pub fn n_samples(&self) -> usize {
        self.samples.len()
    }

    pub fn n_features(&self) -> usize {
        self.features.len()
    }

    pub fn row(&self, i: usize) -> &[f64] {
        let m = self.n_features();
        &self.data[i * m..(i + 1) * m]
    }

    /// Parse a TSV/CSV feature table: header line is feature IDs (the first
    /// cell, the corner, is ignored), each following line is a sample ID
    /// followed by its counts.
    pub fn parse(reader: impl BufRead, delim: char) -> Result<Self> {
        let mut lines = reader.lines();
        let header = lines
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput("empty table".into()))?
            .map_err(RsomicsError::Io)?;
        let features: Vec<String> = header
            .split(delim)
            .skip(1)
            .map(|s| s.trim().to_string())
            .collect();
        if features.is_empty() {
            return Err(RsomicsError::InvalidInput(
                "table header has no feature columns".into(),
            ));
        }
        let m = features.len();

        let mut samples = Vec::new();
        let mut data = Vec::new();
        for line in lines {
            let line = line.map_err(RsomicsError::Io)?;
            if line.trim().is_empty() {
                continue;
            }
            let mut cells = line.split(delim);
            let id = cells
                .next()
                .ok_or_else(|| RsomicsError::InvalidInput("table row has no cells".into()))?;
            samples.push(id.trim().to_string());
            let before = data.len();
            for cell in cells {
                let v: f64 = cell.trim().parse().map_err(|_| {
                    RsomicsError::InvalidInput(format!("non-numeric table value: '{cell}'"))
                })?;
                data.push(v);
            }
            if data.len() - before != m {
                return Err(RsomicsError::InvalidInput(format!(
                    "sample '{}' has {} values, expected {m}",
                    samples.last().unwrap(),
                    data.len() - before
                )));
            }
        }
        if samples.is_empty() {
            return Err(RsomicsError::InvalidInput("table has no samples".into()));
        }
        Ok(Self {
            samples,
            features,
            data,
        })
    }
}

/// Read `id<TAB>group` grouping lines and return one label per table sample, in
/// table order. Comment/blank lines are skipped; a sample missing from the file
/// is an error (skbio's superset rule).
pub fn read_grouping(reader: impl BufRead, samples: &[String], delim: char) -> Result<Vec<String>> {
    use std::collections::HashMap;
    let mut map: HashMap<String, String> = HashMap::new();
    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let mut it = line.splitn(2, delim);
        let id = it.next().unwrap().trim();
        let grp = it
            .next()
            .ok_or_else(|| {
                RsomicsError::InvalidInput(format!("grouping line lacks a group column: '{line}'"))
            })?
            .trim();
        map.insert(id.to_string(), grp.to_string());
    }
    samples
        .iter()
        .map(|id| {
            map.get(id).cloned().ok_or_else(|| {
                RsomicsError::InvalidInput(format!("sample '{id}' has no entry in the grouping"))
            })
        })
        .collect()
}
