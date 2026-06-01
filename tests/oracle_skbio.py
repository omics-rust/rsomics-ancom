#!/usr/bin/env python3
"""skbio.stats.composition.ancom oracle.

Usage: oracle_skbio.py table.tsv grouping.tsv [correction] > result.tsv
Emits `feature<TAB>W<TAB>Signif` matching the rsomics-ancom output.
"""
import sys

import pandas as pd
from skbio.stats.composition import ancom

table_path, grouping_path = sys.argv[1], sys.argv[2]
correction = sys.argv[3] if len(sys.argv) > 3 else "holm"
if correction == "none":
    correction = None

table = pd.read_csv(table_path, sep="\t", index_col=0)
grouping = pd.read_csv(grouping_path, sep="\t", header=None, index_col=0).iloc[:, 0]
grouping = grouping.loc[table.index]

res, _ = ancom(table, grouping, p_adjust=correction)

print("feature\tW\tSignif")
for feat in table.columns:
    w = int(res.loc[feat, "W"])
    signif = bool(res.loc[feat, "Signif"])
    print(f"{feat}\t{w}\t{signif}")
