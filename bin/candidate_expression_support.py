#!/usr/bin/env python3
"""Attach candidate-level RNA expression support from supporting tumor cells."""

from __future__ import annotations

import argparse
import math
from collections import defaultdict

from scneo_utils import read_tsv, safe_float, write_tsv
from candidate_utils import IDENTITY_FIELDS, grouped_candidates, supporting_cells


KEY_FIELDS = IDENTITY_FIELDS
EXPRESSION_FIELDS = KEY_FIELDS + [
    "RNA_expression_support",
    "expression_mapping_status",
    "expression_supporting_cell_ids",
    "expression_missing_cell_ids",
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pmhc", required=True)
    parser.add_argument("--expression", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    pmhc = read_tsv(args.pmhc)
    expression = read_tsv(args.expression)
    rows = attach_expression_support(pmhc, expression)
    write_tsv(args.out, rows, EXPRESSION_FIELDS)


def attach_expression_support(
    pmhc_rows: list[dict[str, str]],
    expression_rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    by_name: dict[tuple[str, str, str], float] = {}
    by_id: dict[tuple[str, str, str], float] = {}
    patient_cells: dict[str, set[str]] = defaultdict(set)
    patient_gene_names: dict[str, set[str]] = defaultdict(set)
    patient_gene_ids: dict[str, set[str]] = defaultdict(set)

    for row in expression_rows:
        patient = row.get("patient_id", "")
        cell = row.get("cell_id", "")
        value = safe_float(row.get("normalized_expression"), None)
        gene_name = (row.get("gene_name") or "").strip()
        gene_id = (row.get("gene_id") or "").split(".")[0]
        patient_cells[patient].add(cell)
        if value is None or not math.isfinite(value):
            continue
        if gene_name:
            key = (patient, cell, gene_name.upper())
            if key in by_name and by_name[key] != value:
                raise ValueError(f"Conflicting expression for gene name {key}")
            by_name[key] = value
            patient_gene_names[patient].add(gene_name.upper())
        if gene_id:
            key = (patient, cell, gene_id.upper())
            if key in by_id and by_id[key] != value:
                raise ValueError(f"Conflicting expression for gene ID {key}")
            by_id[key] = value
            patient_gene_ids[patient].add(gene_id.upper())

    out = []
    for group in grouped_candidates(pmhc_rows).values():
        pmhc = group[0]
        patient = pmhc.get("patient_id", "")
        gene = (pmhc.get("gene") or "").strip()
        supporting = supporting_cells(group)
        present_cells = [cell for cell in supporting if cell in patient_cells.get(patient, set())]
        missing_cells = [cell for cell in supporting if cell not in patient_cells.get(patient, set())]

        gene_key = gene.upper()
        use_name = bool(gene_key and gene_key in patient_gene_names.get(patient, set()))
        use_id = bool(gene_key and gene_key in patient_gene_ids.get(patient, set()))
        status = "ok"
        score: object = ""
        scored_cells: list[str] = []

        if not gene:
            status = "gene_missing"
        elif not use_name and not use_id:
            status = "gene_not_mapped"
        elif not supporting:
            status = "supporting_cells_missing"
        elif not present_cells:
            status = "supporting_cells_not_quantified"
        else:
            values = []
            for cell in present_cells:
                if use_name:
                    key = (patient, cell, gene_key)
                    if key in by_name:
                        values.append(by_name[key])
                        scored_cells.append(cell)
                else:
                    key = (patient, cell, gene_key)
                    if key in by_id:
                        values.append(by_id[key])
                        scored_cells.append(cell)
            if values:
                score = f"{max(values):.10g}"
                missing_cells = sorted(set(supporting) - set(scored_cells))
                status = "ok" if not missing_cells else "partial_cells_missing"
            else:
                status = "gene_not_mapped_in_supporting_cells"

        missing_cells = sorted(set(supporting) - set(scored_cells))
        row = {field: pmhc.get(field, "") for field in KEY_FIELDS}
        row.update(
            {
                "gene": gene,
                "RNA_expression_support": score,
                "expression_mapping_status": status,
                "expression_supporting_cell_ids": ",".join(scored_cells),
                "expression_missing_cell_ids": ",".join(missing_cells),
            }
        )
        out.append(row)
    return out


if __name__ == "__main__":
    main()
