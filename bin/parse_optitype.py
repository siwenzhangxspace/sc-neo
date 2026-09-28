#!/usr/bin/env python3
"""Parse one per-patient OptiType result into a normalized HLA row.

A patient's tumor cells are merged into a single-end FASTQ pool before running
OptiType once per patient, so this parser emits one patient-level row including
the number of contributing cells and the reads/objective diagnostics from the
OptiType result TSV (v1.5.0 column names `nof_reads`/`obj`, with legacy
`Reads`/`Objective` accepted as fallback).
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

from scneo_utils import read_tsv, stringify, uniq, write_tsv


def normalize_allele(value: str) -> str:
    value = (value or "").strip()
    if not value or value in {"-", "NA", "None"}:
        return ""
    value = value.replace("HLA-", "")
    if not value.startswith(("A*", "B*", "C*")):
        if value.startswith(("A", "B", "C")) and "*" not in value:
            value = value[0] + "*" + value[1:]
    return "HLA-" + value


def parse_optitype_tsv(path: Path) -> tuple[List[str], Dict[str, str]]:
    rows = read_tsv(path)
    if not rows:
        return [], {}
    row = rows[0]
    alleles = []
    for col in ["A1", "A2", "B1", "B2", "C1", "C2", "HLA-A1", "HLA-A2", "HLA-B1", "HLA-B2", "HLA-C1", "HLA-C2"]:
        if col in row:
            value = normalize_allele(row.get(col, ""))
            if value:
                alleles.append(value)
    diagnostics = {}
    for new_col, legacy_col in [("nof_reads", "Reads"), ("obj", "Objective")]:
        value = row.get(new_col, row.get(legacy_col, ""))
        if value not in ("", None):
            diagnostics[new_col] = stringify(value)
    return uniq(alleles), diagnostics


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patient-id", required=True)
    parser.add_argument("--n-cells", type=int, default=1, help="number of tumor cells merged into this patient pool")
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    candidates = sorted(Path(args.input_dir).rglob("*.tsv"))
    alleles: List[str] = []
    diagnostics: Dict[str, str] = {}
    for candidate in candidates:
        if candidate.name == Path(args.out).name:
            continue
        parsed, diag = parse_optitype_tsv(candidate)
        if parsed:
            alleles = parsed
            diagnostics = diag
            break
    row = {
        "patient_id": args.patient_id,
        "n_cells": str(args.n_cells),
        "hla_alleles": ",".join(alleles),
        "hla_a": ",".join([a for a in alleles if a.startswith("HLA-A")]),
        "hla_b": ",".join([a for a in alleles if a.startswith("HLA-B")]),
        "hla_c": ",".join([a for a in alleles if a.startswith("HLA-C")]),
        "nof_reads": diagnostics.get("nof_reads", ""),
        "obj": diagnostics.get("obj", ""),
        "optitype_status": "ok" if alleles else "no_hla_result",
    }
    write_tsv(
        args.out,
        [row],
        ["patient_id", "n_cells", "hla_alleles", "hla_a", "hla_b", "hla_c", "nof_reads", "obj", "optitype_status"],
    )


if __name__ == "__main__":
    main()
