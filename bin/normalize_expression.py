#!/usr/bin/env python3
"""Normalize per-cell gene counts as log(1 + raw/total*10000) and write matrices."""

from __future__ import annotations

import argparse
import math
from collections import defaultdict

from scneo_utils import read_tsv, safe_float, write_tsv


LONG_FIELDS = ["patient_id", "cell_id", "gene_id", "gene_name", "raw_count", "normalized_expression"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", nargs="+", required=True)
    parser.add_argument("--out-long", required=True)
    parser.add_argument("--out-raw-matrix", required=True)
    parser.add_argument("--out-normalized-matrix", required=True)
    args = parser.parse_args()

    rows = []
    for path in args.counts:
        rows.extend(read_tsv(path))
    normalized = normalize_rows(rows)
    write_tsv(args.out_long, normalized, LONG_FIELDS)
    write_matrix(args.out_raw_matrix, normalized, "raw_count")
    write_matrix(args.out_normalized_matrix, normalized, "normalized_expression")


def normalize_rows(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    totals: dict[tuple[str, str], float] = defaultdict(float)
    for row in rows:
        totals[(row.get("patient_id", ""), row.get("cell_id", ""))] += (
            safe_float(row.get("raw_count"), 0.0) or 0.0
        )

    out = []
    for row in rows:
        raw = safe_float(row.get("raw_count"), 0.0) or 0.0
        total = totals[(row.get("patient_id", ""), row.get("cell_id", ""))]
        normalized = math.log1p(raw / total * 10000.0) if total > 0 else 0.0
        out.append(
            {
                "patient_id": row.get("patient_id", ""),
                "cell_id": row.get("cell_id", ""),
                "gene_id": row.get("gene_id", ""),
                "gene_name": row.get("gene_name", ""),
                "raw_count": int(raw) if raw.is_integer() else raw,
                "normalized_expression": f"{normalized:.10g}",
            }
        )
    return out


def write_matrix(path: str, rows: list[dict[str, object]], value_field: str) -> None:
    cells = sorted({(str(row["patient_id"]), str(row["cell_id"])) for row in rows})
    cell_columns = [f"{patient}|{cell}" for patient, cell in cells]
    values: dict[tuple[str, str], dict[str, object]] = defaultdict(dict)
    names: dict[tuple[str, str], str] = {}
    for row in rows:
        key = (str(row["gene_id"]), str(row["gene_name"]))
        names[key] = str(row["gene_name"])
        values[key][f"{row['patient_id']}|{row['cell_id']}"] = row.get(value_field, "")

    matrix_rows = []
    for gene_id, gene_name in sorted(values):
        matrix = {"gene_id": gene_id, "gene_name": gene_name}
        for cell in cell_columns:
            matrix[cell] = values[(gene_id, gene_name)].get(cell, 0)
        matrix_rows.append(matrix)
    write_tsv(path, matrix_rows, ["gene_id", "gene_name", *cell_columns])


if __name__ == "__main__":
    main()
