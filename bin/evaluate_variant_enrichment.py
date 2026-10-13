#!/usr/bin/env python3
"""Evaluate tumor enrichment of recurrent variants using cell-level prevalence."""

from __future__ import annotations

import argparse
from collections import defaultdict

from scneo_utils import (
    bh_adjust,
    die,
    fisher_greater,
    infer_cell_id_from_path,
    iter_vcf_records,
    parse_info,
    read_manifest,
    read_tsv,
    variant_key,
    write_tsv,
    write_vcf,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--recurrent-tsv", required=True)
    parser.add_argument("--vcfs", nargs="+", required=True)
    parser.add_argument("--fdr", type=float, default=0.1)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--out-vcf", required=True)
    args = parser.parse_args()

    manifest_rows = read_manifest(args.manifest)
    by_patient_type = defaultdict(list)
    for row in manifest_rows:
        by_patient_type[(row["patient_id"], row["cell_type"])].append(row["cell_id"])

    alt_cells = collect_alt_cells(args.vcfs)
    recurrent = read_tsv(args.recurrent_tsv)

    rows = []
    p_values_by_patient = defaultdict(list)
    for var in recurrent:
        if not var.get("patient_id"):
            die("recurrent-tsv missing required patient_id column or value")
        key = var["variant_key"]
        patient_id = var["patient_id"]
        detected = alt_cells.get(key, set())
        tumor_cells = by_patient_type[(patient_id, "tumor")]
        normal_cells = by_patient_type[(patient_id, "normal")]
        tumor_alt = [cell_id for cell_id in tumor_cells if cell_id in detected]
        normal_alt = [cell_id for cell_id in normal_cells if cell_id in detected]
        tumor_ref = len(tumor_cells) - len(tumor_alt)
        normal_ref = len(normal_cells) - len(normal_alt)
        p_value = fisher_greater(len(tumor_alt), tumor_ref, len(normal_alt), normal_ref)
        p_values_by_patient[patient_id].append((len(rows), p_value))
        row = dict(var)
        row.update(
            {
                "patient_id": patient_id,
                "tumor_alt": len(tumor_alt),
                "tumor_ref": tumor_ref,
                "normal_alt": len(normal_alt),
                "normal_ref": normal_ref,
                # Preserve the schema without claiming unmeasured BAM callability.
                "tumor_callable_cells": "",
                "normal_callable_cells": "",
                "tumor_alt_cell_ids": ",".join(tumor_alt),
                "tumor_callable_cell_ids": "",
                "normal_alt_cell_ids": ",".join(normal_alt),
                "normal_callable_cell_ids": "",
                "not_callable_cell_ids": "",
                "fisher_p_value": p_value,
                "fisher_fdr": "",
                "tumor_enriched_candidate": "false",
                "tumor_specific_candidate": "false",
                "normal_evidence_status": "normal_callability_not_measured",
                "passes_enrichment": "false",
            }
        )
        rows.append(row)

    for patient_values in p_values_by_patient.values():
        q_values = bh_adjust([p_value for _, p_value in patient_values])
        for (row_index, _), q_value in zip(patient_values, q_values):
            row = rows[row_index]
            row["fisher_fdr"] = f"{q_value:.6g}"
            enriched = q_value <= args.fdr
            row["tumor_enriched_candidate"] = "true" if enriched else "false"
            row["passes_enrichment"] = "true" if enriched else "false"

    fields = [
        "variant_key",
        "patient_id",
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "varType",
        "supporting_tumor_cells",
        "supporting_tumor_cell_ids",
        "tumor_alt",
        "tumor_ref",
        "normal_alt",
        "normal_ref",
        "tumor_callable_cells",
        "normal_callable_cells",
        "tumor_alt_cell_ids",
        "tumor_callable_cell_ids",
        "normal_alt_cell_ids",
        "normal_callable_cell_ids",
        "not_callable_cell_ids",
        "fisher_p_value",
        "fisher_fdr",
        "tumor_enriched_candidate",
        "tumor_specific_candidate",
        "normal_evidence_status",
        "passes_enrichment",
    ]
    write_tsv(args.out_tsv, rows, fields)

    passed = defaultdict(list)
    for row in rows:
        if row["passes_enrichment"] == "true":
            passed[row["variant_key"]].append(row)
    write_vcf(args.out_vcf, summarize_vcf_rows(passed))


def summarize_vcf_rows(
    grouped: dict[str, list[dict[str, str]]],
) -> list[dict[str, str]]:
    out = []
    for key in sorted(grouped):
        patient_rows = sorted(grouped[key], key=lambda row: row["patient_id"])
        row = dict(patient_rows[0])
        row["patient_ids"] = ",".join(row["patient_id"] for row in patient_rows)
        row["patient_support_counts"] = ",".join(
            f"{patient_row['patient_id']}:{patient_row.get('supporting_tumor_cells', '')}"
            for patient_row in patient_rows
        )
        row["patient_supporting_tumor_cell_ids"] = "|".join(
            f"{patient_row['patient_id']}:{patient_row.get('supporting_tumor_cell_ids', '')}"
            for patient_row in patient_rows
        )
        out.append(row)
    return out


def collect_alt_cells(vcfs: list[str]) -> dict[str, set[str]]:
    out = defaultdict(set)
    for vcf in vcfs:
        cell_id = infer_cell_id_from_path(vcf)
        for rec in iter_vcf_records(vcf):
            info = parse_info(rec["INFO"])
            vartype = (
                info.get("varType")
                or info.get("TYPE")
                or info.get("SVTYPE")
                or infer_vartype(rec["REF"], rec["ALT"])
            )
            for alt in rec["ALT"].split(","):
                out[
                    variant_key(rec["CHROM"], rec["POS"], rec["REF"], alt, vartype)
                ].add(cell_id)
    return out


def infer_vartype(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV"
    return "INDEL"


if __name__ == "__main__":
    main()
