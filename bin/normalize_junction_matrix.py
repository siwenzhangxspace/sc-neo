#!/usr/bin/env python3
"""Normalize AltAnalyze count-matrix sample columns to manifest cell IDs."""

from __future__ import annotations

import argparse
import csv

from scneo_utils import read_manifest


SUFFIXES = (".sorted.bam.bed", ".sorted.bed", ".bam.bed", ".sorted.bam", ".bam", ".bed")


def normalize_sample_name(name: str, expected: set[str]) -> str | None:
    name = name.strip()
    if name in expected:
        return name
    for suffix in SUFFIXES:
        if name.endswith(suffix) and name[: -len(suffix)] in expected:
            return name[: -len(suffix)]
    return None


def normalize_matrix(input_path: str, manifest_path: str, output_path: str) -> None:
    manifest = read_manifest(manifest_path)
    expected_order = [row["cell_id"] for row in manifest if row["cell_type"] == "tumor"]
    expected = set(expected_order)
    with open(input_path, newline="") as source:
        reader = csv.reader(source, delimiter="\t")
        header = next(reader, None)
        if not header or len(header) < 2:
            raise ValueError("AltAnalyze matrix has no sample columns")
        normalized = []
        unknown = []
        for name in header[1:]:
            cell_id = normalize_sample_name(name, expected)
            if cell_id is None:
                unknown.append(name)
            normalized.append(cell_id)
        if unknown:
            raise ValueError("AltAnalyze matrix contains unknown sample columns: " + ", ".join(unknown))
        duplicates = sorted({cell for cell in normalized if normalized.count(cell) > 1})
        if duplicates:
            raise ValueError("AltAnalyze matrix maps multiple columns to: " + ", ".join(duplicates))
        missing = [cell for cell in expected_order if cell not in normalized]
        if missing:
            raise ValueError("AltAnalyze matrix is missing manifest cells: " + ", ".join(missing))
        column_by_cell = {cell: index + 1 for index, cell in enumerate(normalized)}
        with open(output_path, "w", newline="") as destination:
            writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
            writer.writerow([header[0]] + expected_order)
            for line_number, fields in enumerate(reader, start=2):
                if len(fields) != len(header):
                    raise ValueError(f"AltAnalyze matrix line {line_number} has {len(fields)} columns; expected {len(header)}")
                for index, value in enumerate(fields[1:], start=2):
                    try:
                        float(value)
                    except ValueError as exc:
                        raise ValueError(f"AltAnalyze matrix line {line_number}, column {index} is not numeric: {value}") from exc
                writer.writerow([fields[0]] + [fields[column_by_cell[cell]] for cell in expected_order])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    normalize_matrix(args.matrix, args.manifest, args.out)


if __name__ == "__main__":
    main()
