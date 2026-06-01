use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_ancom::{Correction, Table, ancom};

fn synth(n_samples: usize, n_features: usize) -> (Table, Vec<String>) {
    let mut header = String::new();
    for j in 0..n_features {
        header.push('\t');
        header.push_str(&format!("f{j}"));
    }
    let mut txt = header;
    txt.push('\n');
    let mut state = 0x2545_F491_4F6C_DD1D_u64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state % 100 + 1) as f64
    };
    let mut grouping = Vec::with_capacity(n_samples);
    for s in 0..n_samples {
        txt.push_str(&format!("s{s}"));
        for _ in 0..n_features {
            txt.push('\t');
            txt.push_str(&format!("{}", next() as u64));
        }
        txt.push('\n');
        grouping.push(format!("g{}", s % 2));
    }
    (Table::parse(txt.as_bytes(), '\t').unwrap(), grouping)
}

fn bench(c: &mut Criterion) {
    let (table, grouping) = synth(40, 300);
    c.bench_function("ancom_40x300", |b| {
        b.iter(|| {
            ancom(
                black_box(&table),
                black_box(&grouping),
                0.05,
                0.02,
                0.1,
                Correction::Holm,
            )
            .unwrap()
        })
    });
}

criterion_group!(benches, bench);
criterion_main!(benches);
