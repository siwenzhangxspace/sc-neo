#!/usr/bin/env python3
"""Write a validated manifest subset for one cell type."""
from __future__ import annotations
import argparse
from scneo_utils import read_manifest, write_tsv

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--manifest", required=True)
parser.add_argument("--cell-type", required=True)
parser.add_argument("--out", required=True)
args = parser.parse_args()
rows = [row for row in read_manifest(args.manifest) if row["cell_type"] == args.cell_type]
if not rows:
    raise ValueError(f"manifest has no {args.cell_type} cells")
write_tsv(args.out, rows, ["patient_id", "cell_id", "cell_type", "fastq_1", "fastq_2", "read_layout"])
