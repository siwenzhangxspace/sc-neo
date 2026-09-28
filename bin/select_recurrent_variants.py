#!/usr/bin/env python3
"""Merge tumor-cell Scallele VCFs and retain recurrent tumor variants."""

from __future__ import annotations

import argparse
from collections import defaultdict

from scneo_utils import (
    infer_cell_id_from_path,
    iter_vcf_records,
    parse_info,
    read_manifest,
    uniq,
    variant_key,
    write_tsv,
    write_vcf,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--vcfs", nargs="+", required=True)
    parser.add_argument("--min-tumor-cells", type=int, default=3)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--out-vcf", required=True)
    args = parser.parse_args()

    manifest = {row["cell_id"]: row for row in read_manifest(args.manifest)}
    variants = {}
    patient_support = defaultdict(list)

    for vcf in args.vcfs:
        cell_id = infer_cell_id_from_path(vcf)
        meta = manifest.get(cell_id)
        if not meta or meta["cell_type"] != "tumor":
            continue
        for rec in iter_vcf_records(vcf):
            info = parse_info(rec["INFO"])
            vartype = (
                info.get("varType")
                or info.get("TYPE")
                or info.get("SVTYPE")
                or infer_vartype(rec["REF"], rec["ALT"])
            )
            for alt in rec["ALT"].split(","):
                key = variant_key(rec["CHROM"], rec["POS"], rec["REF"], alt, vartype)
                variants.setdefault(
                    key,
                    {
                        "CHROM": rec["CHROM"],
                        "POS": rec["POS"],
                        "ID": rec.get("ID", "."),
                        "REF": rec["REF"],
                        "ALT": alt,
                        "varType": vartype,
                    },
                )
                patient_support[(meta["patient_id"], key)].append(cell_id)

    rows = []
    for patient_id, key in sorted(patient_support):
        cells = uniq(patient_support[(patient_id, key)])
        if len(cells) < args.min_tumor_cells:
            continue
        row = dict(variants[key])
        row.update(
            {
                "patient_id": patient_id,
                "variant_key": key,
                "supporting_tumor_cells": len(cells),
                "supporting_tumor_cell_ids": ",".join(cells),
            }
        )
        rows.append(row)

    fields = [
        "patient_id",
        "variant_key",
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "varType",
        "supporting_tumor_cells",
        "supporting_tumor_cell_ids",
    ]
    write_tsv(args.out_tsv, rows, fields)
    write_vcf(args.out_vcf, summarize_vcf_rows(rows))


def summarize_vcf_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["variant_key"]].append(row)

    out = []
    for key in sorted(grouped):
        patient_rows = sorted(grouped[key], key=lambda row: row["patient_id"])
        row = dict(patient_rows[0])
        patient_ids = [patient_row["patient_id"] for patient_row in patient_rows]
        row["patient_ids"] = ",".join(patient_ids)
        row["patient_support_counts"] = ",".join(
            f"{patient_row['patient_id']}:{patient_row['supporting_tumor_cells']}"
            for patient_row in patient_rows
        )
        row["patient_supporting_tumor_cell_ids"] = "|".join(
            f"{patient_row['patient_id']}:{patient_row['supporting_tumor_cell_ids']}"
            for patient_row in patient_rows
        )
        out.append(row)
    return out


def infer_vartype(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV"
    return "INDEL"


if __name__ == "__main__":
    main()
