#!/usr/bin/env python3
"""Attach mean source-gene TPM across every sampled malignant cell of a patient."""

from __future__ import annotations
import argparse
import math
from collections import defaultdict
from scneo_utils import read_tsv, write_tsv
from candidate_utils import IDENTITY_FIELDS, grouped_candidates

EXPRESSION_FIELDS = IDENTITY_FIELDS + [
    "RNA_expression_support", "mean_malignant_TPM", "expression_unit", "expression_source",
    "malignant_cells", "expression_measured_cells", "expression_detected_cells",
    "expression_mapping_status", "expression_supporting_cell_ids", "expression_missing_cell_ids",
]


def attach_expression_support(pmhc_rows, expression_rows, manifest=None):
    patients = defaultdict(set)
    if manifest is not None:
        for row in manifest:
            if row["cell_type"] == "tumor":
                patients[row["patient_id"]].add(row["cell_id"])
    by_name, by_id = defaultdict(dict), defaultdict(dict)
    sources = defaultdict(set)
    for row in expression_rows:
        patient, cell = row["patient_id"], row["cell_id"]
        if manifest is not None and cell not in patients[patient]:
            raise ValueError(f"Expression cell outside the malignant-cell manifest {patient}/{cell}")
        patients[patient].add(cell)
        if "tpm" not in row:
            raise ValueError("Expression needs unlogged tpm. Regenerate or supply measured TPM.")
        sources[patient].add(row.get("expression_source", "supplied_gene_TPM"))
        value = row["tpm"]
        if value in ("", None):
            continue
        value = float(value)
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"Invalid TPM for {patient}/{cell}")
        identifiers = [(by_name, row.get("gene_name", "").upper()),
                       (by_id, row.get("gene_id", "").split(".")[0].upper())]
        for index, gene in identifiers:
            if gene:
                values = index[(patient, gene)]
                if cell in values and values[cell] != value:
                    raise ValueError(f"Ambiguous expression mapping {patient}/{cell}/{gene}")
                values[cell] = value
    out = []
    for group in grouped_candidates(pmhc_rows).values():
        candidate = group[0]
        patient, gene = candidate["patient_id"], candidate["gene"].upper()
        cells = patients[patient]
        measured = by_name.get((patient, gene), by_id.get((patient, gene), {}))
        missing = cells - measured.keys()
        value = ""
        status = "gene_not_mapped"
        if not cells:
            status = "malignant_cells_missing"
        elif measured:
            status = "incomplete_gene_measurement" if missing else "ok"
            if not missing:
                value = format(math.fsum(measured[cell] for cell in sorted(cells)) / len(cells), ".17g")
        row = {field: candidate.get(field, "") for field in IDENTITY_FIELDS}
        row.update({
            "RNA_expression_support": value, "mean_malignant_TPM": value, "expression_unit": "TPM",
            "expression_source": ";".join(sorted(sources[patient])), "malignant_cells": len(cells),
            "expression_measured_cells": len(measured), "expression_detected_cells": sum(v > 0 for v in measured.values()),
            "expression_mapping_status": status, "expression_supporting_cell_ids": ",".join(sorted(measured)),
            "expression_missing_cell_ids": ",".join(sorted(missing)),
        })
        out.append(row)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pmhc", required=True)
    parser.add_argument("--expression", required=True)
    parser.add_argument("--manifest")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    rows = attach_expression_support(read_tsv(args.pmhc), read_tsv(args.expression),
                                     read_tsv(args.manifest) if args.manifest else None)
    write_tsv(args.out, rows, EXPRESSION_FIELDS)


if __name__ == "__main__":
    main()
