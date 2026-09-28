#!/usr/bin/env python3
"""Merge and filter netMHCpan mutation/splicing predictions."""

from __future__ import annotations

import argparse

from scneo_utils import read_tsv, safe_float, split_list, write_tsv


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mutation", required=True)
    parser.add_argument("--splicing", required=True)
    parser.add_argument("--affinity-max", type=float, default=150.0)
    parser.add_argument("--rank-el-max", type=float, default=2.0)
    parser.add_argument("--bind-levels", default="wb,sb")
    parser.add_argument("--ignore-bind-level", action="store_true")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    allowed = {value.lower() for value in split_list(args.bind_levels)}
    rows = []
    for source_path in [args.mutation, args.splicing]:
        for row in read_tsv(source_path):
            if row.get("netmhcpan_status") != "ok":
                continue
            aff = safe_float(row.get("Aff.nM"), None)
            if aff is None or aff >= args.affinity_max:
                continue
            rank_el = safe_float(row.get("rank_el"), None)
            if rank_el is None or rank_el >= args.rank_el_max:
                continue
            row["bind_level"] = row.get("bind_level", "").lower()
            if not args.ignore_bind_level and row["bind_level"] not in allowed:
                continue
            row["normal_evidence_cell_ids"] = row.get("normal_callable_cell_ids", "")
            rows.append(row)

    fields = [
        "patient_id",
        "source_cell_ids",
        "supporting_tumor_cell_ids",
        "normal_evidence_cell_ids",
        "normal_alt_cell_ids",
        "normal_callable_cell_ids",
        "source_type",
        "peptide",
        "peptide_length",
        "hla_allele",
        "Aff.nM",
        "rank_el",
        "bind_level",
        "variant_or_junction_id",
        "variant_class",
        "gene",
        "transcript_id",
        "junction_annotation",
        "splice_site",
        "translation_frame",
        "n_from_first",
        "phase_supported",
        "phase_evidence",
        "reference_protein_match",
        "gtex_mean",
        "tumor_gtex_delta",
        "supporting_tumor_cells",
        "normal_evidence_status",
        "tumor_specific_candidate",
        "tumor_enriched_candidate",
        "annotation_summary",
        "hla_batch_id",
        "hla_batch_alleles",
    ]
    for row in rows:
        row["supporting_tumor_cells"] = len(split_list(row.get("supporting_tumor_cell_ids", "")))
    write_tsv(args.out, rows, fields)


if __name__ == "__main__":
    main()
