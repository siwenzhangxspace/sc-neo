#!/usr/bin/env python3
"""Convert featureCounts gene-level output to long per-cell raw counts."""

from __future__ import annotations

import argparse
import csv

from scneo_utils import split_list, write_tsv


FIELDS = ["patient_id", "cell_id", "gene_id", "gene_name", "raw_count"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--featurecounts", required=True)
    parser.add_argument("--patient-id", required=True)
    parser.add_argument("--cell-ids", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    cell_ids = split_list(args.cell_ids)
    rows = parse_featurecounts(args.featurecounts, args.patient_id, cell_ids)
    write_tsv(args.out, rows, FIELDS)


def parse_featurecounts(path: str, patient_id: str, cell_ids: list[str]) -> list[dict[str, object]]:
    with open(path) as handle:
        data_lines = (line for line in handle if line.strip() and not line.startswith("#"))
        reader = csv.DictReader(data_lines, delimiter="\t")
        fieldnames = reader.fieldnames or []
        if not fieldnames:
            return []
        if "Geneid" not in fieldnames:
            raise ValueError("featureCounts output is missing Geneid")

        metadata_columns = {"Geneid", "Chr", "Start", "End", "Strand", "Length", "gene_name"}
        count_columns = [column for column in fieldnames if column not in metadata_columns]
        if len(count_columns) != len(cell_ids):
            raise ValueError(
                f"featureCounts sample columns ({len(count_columns)}) do not match cell ids ({len(cell_ids)})"
            )

        rows = []
        for record in reader:
            gene_id = (record.get("Geneid") or "").split(".")[0]
            gene_name = record.get("gene_name") or gene_id
            for cell_id, column in zip(cell_ids, count_columns):
                value = record.get(column, "0")
                try:
                    count = int(float(value or 0))
                except ValueError:
                    count = 0
                rows.append(
                    {
                        "patient_id": patient_id,
                        "cell_id": cell_id,
                        "gene_id": gene_id,
                        "gene_name": gene_name,
                        "raw_count": count,
                    }
                )
        return rows


if __name__ == "__main__":
    main()
