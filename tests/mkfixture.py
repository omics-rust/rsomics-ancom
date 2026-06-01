#!/usr/bin/env python3
"""Generate deterministic ANCOM fixtures (feature table + grouping).

Usage: mkfixture.py n_samples n_features seed out_table.tsv out_groups.tsv [n_groups]

Counts are positive integers (ANCOM requires positive values). A handful of
features are spiked to differ across groups so the W statistics are non-trivial.
"""
import sys

import numpy as np

n_samples = int(sys.argv[1])
n_features = int(sys.argv[2])
seed = int(sys.argv[3])
out_table = sys.argv[4]
out_groups = sys.argv[5]
n_groups = int(sys.argv[6]) if len(sys.argv) > 6 else 2

rng = np.random.default_rng(seed)

labels = np.array([i % n_groups for i in range(n_samples)])

base = rng.integers(10, 100, size=(n_samples, n_features)).astype(float)
# spike ~5% of features to be group-dependent
n_spike = max(1, n_features // 20)
spike_idx = rng.choice(n_features, size=n_spike, replace=False)
for f in spike_idx:
    shift = rng.uniform(2.0, 6.0, size=n_groups)
    for s in range(n_samples):
        base[s, f] *= shift[labels[s]]
counts = np.rint(base).astype(int)
counts = np.maximum(counts, 1)

samples = [f"s{i}" for i in range(n_samples)]
features = [f"f{j}" for j in range(n_features)]

with open(out_table, "w") as fh:
    fh.write("\t" + "\t".join(features) + "\n")
    for i, s in enumerate(samples):
        fh.write(s + "\t" + "\t".join(str(int(c)) for c in counts[i]) + "\n")

with open(out_groups, "w") as fh:
    for i, s in enumerate(samples):
        fh.write(f"{s}\tg{labels[i]}\n")
