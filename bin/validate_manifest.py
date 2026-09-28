#!/usr/bin/env python3
"""Validate and normalize the Smart-seq2 manifest."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path

from scneo_utils import die, read_tsv, write_tsv


ALLOWED_CELL_TYPES = {"tumor", "normal", "tcell"}
REQUIRED = ["patient_id", "cell_id", "cell_type", "fastq_1", "fastq_2"]
OUTPUT_COLUMNS = REQUIRED + ["read_layout"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--out-manifest", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--min-cells-per-type", type=int, default=3)
    parser.add_argument("--allow-single-end", action="store_true")
    args = parser.parse_args()

    rows = read_tsv(args.manifest)
    if not rows:
        die("manifest is empty")

    missing = [col for col in REQUIRED if col not in rows[0]]
    if missing:
        die(f"manifest missing required columns: {','.join(missing)}")

    seen = set()
    errors = []
    normalized = []
    counts = Counter()
    patient_counts = defaultdict(Counter)

    for line_no, row in enumerate(rows, start=2):
        cleaned = {key: (row.get(key) or "").strip() for key in row}
        for col in REQUIRED:
            if col != "fastq_2" or not args.allow_single_end:
                if not cleaned.get(col):
                    errors.append(f"line {line_no}: missing {col}")

        cell_id = cleaned.get("cell_id", "")
        patient_id = cleaned.get("patient_id", "")
        cell_type = cleaned.get("cell_type", "")
        if "." in cell_id:
            errors.append(f"line {line_no}: cell_id must not contain '.' because TRUST4 uses the text before the first dot")
        if cell_id in seen:
            errors.append(f"line {line_no}: duplicated cell_id {cell_id}")
        seen.add(cell_id)
        if cell_type not in ALLOWED_CELL_TYPES:
            errors.append(f"line {line_no}: invalid cell_type {cell_type}")

        for fq_col in ["fastq_1", "fastq_2"]:
            fq = cleaned.get(fq_col, "")
            if not fq and fq_col == "fastq_2" and args.allow_single_end:
                continue
            fq_path = Path(fq)
            if not fq_path.exists():
                errors.append(f"line {line_no}: {fq_col} does not exist: {fq}")
            elif not fq_path.is_file():
                errors.append(f"line {line_no}: {fq_col} is not a file: {fq}")

        counts[cell_type] += 1
        patient_counts[patient_id][cell_type] += 1
        normalized.append(
            {
                "patient_id": patient_id,
                "cell_id": cell_id,
                "cell_type": cell_type,
                "fastq_1": str(Path(cleaned.get("fastq_1", "")).resolve()) if cleaned.get("fastq_1") else "",
                "fastq_2": str(Path(cleaned.get("fastq_2", "")).resolve()) if cleaned.get("fastq_2") else "",
                "read_layout": "paired" if cleaned.get("fastq_2") else "single",
            }
        )

    for patient_id, patient_count in sorted(patient_counts.items()):
        for cell_type in sorted(ALLOWED_CELL_TYPES):
            if patient_count[cell_type] < args.min_cells_per_type:
                errors.append(
                    f"patient {patient_id} has {patient_count[cell_type]} {cell_type} cells; "
                    f"minimum is {args.min_cells_per_type}"
                )

    if errors:
        for error in errors:
            print(error)
        die(f"manifest validation failed with {len(errors)} error(s)")

    write_tsv(args.out_manifest, normalized, OUTPUT_COLUMNS)

    summary_rows = []
    for patient_id in sorted(patient_counts):
        row = {
            "patient_id": patient_id,
            "tumor_cells": patient_counts[patient_id]["tumor"],
            "normal_cells": patient_counts[patient_id]["normal"],
            "tcell_cells": patient_counts[patient_id]["tcell"],
        }
        row["patient_status"] = (
            "ok"
            if all(row[f"{ct}_cells"] >= args.min_cells_per_type for ct in ["tumor", "normal", "tcell"])
            else "insufficient_cells"
        )
        summary_rows.append(row)
    write_tsv(args.summary, summary_rows, ["patient_id", "tumor_cells", "normal_cells", "tcell_cells", "patient_status"])


if __name__ == "__main__":
    main()
