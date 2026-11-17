#!/usr/bin/env python3
"""Produce gene TPM from featureCounts lengths or import measured per-cell TPM."""

from __future__ import annotations

import argparse
import math
from collections import defaultdict
from scneo_utils import read_tsv, write_tsv

LONG_FIELDS = ["patient_id", "cell_id", "gene_id", "gene_name", "raw_count",
               "gene_length", "normalized_expression", "tpm", "expression_unit", "expression_source"]


def normalize_rows(rows):
    totals = defaultdict(float)
    rate_totals = defaultdict(float)
    parsed = []
    seen = set()
    for row in rows:
        patient, cell = row["patient_id"], row["cell_id"]
        key = (patient, cell, row["gene_id"])
        if key in seen:
            raise ValueError(f"Duplicate count measurement {key}")
        seen.add(key)
        raw = float(row["raw_count"])
        if not row.get("gene_length"):
            raise ValueError("Counts lack gene_length. Regenerate with parse_featurecounts.py.")
        length = float(row["gene_length"])
        if not math.isfinite(raw) or raw < 0 or not math.isfinite(length) or length <= 0:
            raise ValueError(f"Invalid count or gene length for {key}")
        totals[(patient, cell)] += raw
        rate_totals[(patient, cell)] += raw / length
        parsed.append((row, raw, length))
    result = []
    for row, raw, length in parsed:
        key = (row["patient_id"], row["cell_id"])
        total, rate_total = totals[key], rate_totals[key]
        result.append({
            **row, "raw_count": raw, "gene_length": length,
            "normalized_expression": format(math.log1p(raw / total * 10000), ".17g") if total else "",
            "tpm": format(raw / length / rate_total * 1e6, ".17g") if rate_total else "",
            "expression_unit": "TPM", "expression_source": "featureCounts_gene_length",
        })
    return result


def import_tpm(rows):
    result = []
    seen = set()
    for row in rows:
        if not row.get("patient_id") or not row.get("cell_id") or not (row.get("gene_name") or row.get("gene_id")):
            raise ValueError("TPM input needs patient_id, cell_id and gene_name or gene_id")
        if "tpm" not in row:
            raise ValueError("TPM input requires a tpm column on the unlogged scale")
        key = (row["patient_id"], row["cell_id"], row.get("gene_id") or row["gene_name"])
        if key in seen:
            raise ValueError(f"Duplicate TPM measurement {key}")
        seen.add(key)
        value = row["tpm"].strip()
        if value:
            tpm = float(value)
            if not math.isfinite(tpm) or tpm < 0:
                raise ValueError(f"Invalid TPM at {key}")
            value = format(tpm, ".17g")
        result.append({**row, "raw_count": "", "gene_length": "",
                       "normalized_expression": "", "tpm": value, "expression_unit": "TPM",
                       "expression_source": "supplied_gene_TPM"})
    return result


def write_matrix(path, rows, value_field):
    cells = sorted({(str(row["patient_id"]), str(row["cell_id"])) for row in rows})
    cell_columns = [f"{patient}|{cell}" for patient, cell in cells]
    values = defaultdict(dict)
    for row in rows:
        key = (str(row.get("gene_id", "")), str(row.get("gene_name", "")))
        column = f"{row['patient_id']}|{row['cell_id']}"
        if column in values[key]:
            raise ValueError(f"Duplicate matrix measurement {key}/{column}")
        values[key][column] = row.get(value_field, "")
    matrix_rows = []
    for gene_id, gene_name in sorted(values):
        matrix = {"gene_id": gene_id, "gene_name": gene_name}
        matrix.update({cell: values[(gene_id, gene_name)].get(cell, "") for cell in cell_columns})
        matrix_rows.append(matrix)
    write_tsv(path, matrix_rows, ["gene_id", "gene_name", *cell_columns])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--counts", nargs="+")
    source.add_argument("--tpm")
    parser.add_argument("--out-long", required=True)
    parser.add_argument("--out-raw-matrix", required=True)
    parser.add_argument("--out-normalized-matrix", required=True)
    parser.add_argument("--out-tpm-matrix")
    args = parser.parse_args()
    if args.tpm:
        normalized = import_tpm(read_tsv(args.tpm))
    else:
        rows = [row for path in args.counts for row in read_tsv(path)]
        normalized = normalize_rows(rows)
    write_tsv(args.out_long, normalized, LONG_FIELDS)
    write_matrix(args.out_raw_matrix, normalized, "raw_count")
    write_matrix(args.out_normalized_matrix, normalized, "normalized_expression")
    if args.out_tpm_matrix:
        write_matrix(args.out_tpm_matrix, normalized, "tpm")


if __name__ == "__main__":
    main()
