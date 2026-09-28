#!/usr/bin/env python3
"""Summarize per-patient OptiType HLA calls into the patient-level allele table.

Each patient is typed once from a merged tumor-cell FASTQ pool, so the per-patient
files already carry the final allele set. This script merges them into one table,
keeps the coverage diagnostics (nof_reads/obj), and derives a completeness status:
`ok` (6 alleles), `partial` (1-5), `no_hla` (no alleles), or `missing` (no file).
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

from scneo_utils import read_manifest, read_tsv, write_tsv

FIELDS = ["patient_id", "n_cells", "hla_alleles", "hla_a", "hla_b", "hla_c", "nof_reads", "obj", "hla_status"]


def status_for(alleles: List[str]) -> str:
    n = len(alleles)
    if n >= 6:
        return "ok"
    if n >= 1:
        return "partial"
    return "no_hla"


def summarize(manifest_path: str, per_patient: List[str]) -> List[Dict[str, str]]:
    manifest = read_manifest(manifest_path)
    patients = sorted({row["patient_id"] for row in manifest})

    rows_by_patient: Dict[str, Dict[str, str]] = {}
    for path in per_patient:
        rows = read_tsv(path)
        if not rows:
            continue
        rows_by_patient[rows[0].get("patient_id", "")] = rows[0]

    patient_rows = []
    for patient_id in patients:
        row = rows_by_patient.get(patient_id)
        if row is None:
            patient_rows.append({key: "" for key in FIELDS} | {"patient_id": patient_id, "hla_status": "missing"})
            continue
        alleles = [a for a in (row.get("hla_alleles", "") or "").split(",") if a]
        patient_rows.append(
            {
                "patient_id": patient_id,
                "n_cells": row.get("n_cells", ""),
                "hla_alleles": ",".join(sorted(alleles)),
                "hla_a": ",".join(sorted(a for a in alleles if a.startswith("HLA-A"))),
                "hla_b": ",".join(sorted(a for a in alleles if a.startswith("HLA-B"))),
                "hla_c": ",".join(sorted(a for a in alleles if a.startswith("HLA-C"))),
                "nof_reads": row.get("nof_reads", ""),
                "obj": row.get("obj", ""),
                "hla_status": status_for(alleles),
            }
        )
    return patient_rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--per-patient", nargs="*", default=[])
    parser.add_argument("--out-patient", required=True)
    args = parser.parse_args()

    write_tsv(args.out_patient, summarize(args.manifest, args.per_patient), FIELDS)


if __name__ == "__main__":
    main()
