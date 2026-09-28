#!/usr/bin/env python
"""Create full and alternative-event AltAnalyze count matrices."""

from __future__ import print_function

import os
import sys

import numpy as np
import pandas as pd


if len(sys.argv) != 2:
    raise SystemExit("Usage: prune-altanalyze.py ALTANALYZE_OUTPUT")

output = os.path.abspath(sys.argv[1])
expression = os.path.join(output, "ExpressionInput")
counts_path = os.path.join(expression, "counts.original.txt")
annotation_path = os.path.join(
    output,
    "AltResults",
    "AlternativeOutput",
    "Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt",
)

counts = pd.read_csv(counts_path, sep="\t", index_col=0)
counts.index = [item.split("=")[0] for item in counts.index]
counts = counts.loc[np.logical_not(counts.index.duplicated()).tolist(), :]
counts.to_csv(os.path.join(expression, "counts.original.full.txt"), sep="\t")

annotations = pd.read_csv(annotation_path, sep="\t")
events = set(":".join(item.split("|")[0].split(":")[1:]) for item in annotations["UID"])
counts.loc[counts.index.isin(events), :].to_csv(
    os.path.join(expression, "counts.original.pruned.txt"), sep="\t"
)
