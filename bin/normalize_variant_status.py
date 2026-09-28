#!/usr/bin/env python3
"""Normalize legacy mutation classification fields to the final workflow terminology."""

from __future__ import annotations

import argparse

from scneo_utils import read_tsv, write_tsv


STATUS_MAP = {
    "tumor_variant": "somatic_candidate",
    "somatic_candidate": "somatic_candidate",
    "novel_variant": "novel_variant",
    "uncertain": "uncertain",
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maf", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    rows = read_tsv(args.maf)
    for row in rows:
        status = STATUS_MAP.get(row.get("variant_status") or row.get("variant_class", ""), "uncertain")
        row["variant_status"] = status
        row["somatic_candidate"] = "true" if status == "somatic_candidate" else "false"
        row["somatic_mutation_candidate"] = row["somatic_candidate"]
        row["novel_variant"] = "true" if status == "novel_variant" else "false"
        row["uncertain_variant"] = "true" if status == "uncertain" else "false"

    fields = list(rows[0].keys()) if rows else [
        "Hugo_Symbol",
        "Chromosome",
        "Start_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
    ]
    for field in [
        "variant_status",
        "somatic_candidate",
        "somatic_mutation_candidate",
        "novel_variant",
        "uncertain_variant",
    ]:
        if field not in fields:
            fields.append(field)
    write_tsv(args.out, rows, fields)


if __name__ == "__main__":
    main()
