#!/usr/bin/env python3
"""Integrate candidate pMHC, immunogenicity, immunopeptidome, TABR-BERT and RNA evidence."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path

from scneo_utils import read_tsv, split_list, write_tsv


FINAL_FIELDS = [
    "patient_id",
    "patient_status",
    "source_type",
    "gene",
    "transcript_id",
    "variant_or_junction_id",
    "variant_status",
    "variant_class",
    "Variant_Classification",
    "peptide_sequence",
    "peptide_length",
    "HLA_allele",
    "Aff.nM",
    "%Rank_EL",
    "bind_level",
    "IEDB_immunogenicity_score",
    "PCI_DB_enabled",
    "healthy_exact_match",
    "tumor_immune_peptide",
    "best_match_peptide",
    "blast_identity",
    "blast_evalue",
    "max_TABR_BERT_rank",
    "best_tcr_cell_id",
    "best_cdr3_beta_aa",
    "RNA_expression_support",
    "expression_mapping_status",
    "expression_supporting_cell_ids",
    "expression_missing_cell_ids",
    "supporting_tumor_cell_ids",
    "normal_alt_cell_ids",
    "normal_callable_cell_ids",
    "normal_evidence_status",
    "tumor_specific_candidate",
    "tumor_enriched_candidate",
    "annotation_summary",
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pmhc", required=True)
    parser.add_argument("--mutation-maf")
    parser.add_argument("--iedb", required=True)
    parser.add_argument("--immunopeptidome")
    parser.add_argument("--pcidb")  # Backward-compatible legacy input.
    parser.add_argument("--tabr", required=True)
    parser.add_argument("--expression")
    parser.add_argument("--manifest-summary")
    parser.add_argument("--out", required=True)
    parser.add_argument("--by-patient-dir", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    pmhc_rows = read_tsv(args.pmhc)
    mutation_maf = (
        {row.get("variant_key", ""): row for row in read_tsv(args.mutation_maf) if row.get("variant_key")}
        if args.mutation_maf
        else {}
    )
    iedb = index_iedb(read_tsv(args.iedb))
    if args.immunopeptidome:
        immunopeptidome = {
            row.get("candidate_peptide", ""): row
            for row in read_tsv(args.immunopeptidome)
            if row.get("candidate_peptide")
        }
    elif args.pcidb:
        immunopeptidome = legacy_pcidb_evidence(read_tsv(args.pcidb))
    else:
        immunopeptidome = {}
    tabr = {
        tabr_key(row): row
        for row in read_tsv(args.tabr)
    }
    expression = (
        {candidate_key(row): row for row in read_tsv(args.expression)}
        if args.expression
        else {}
    )
    patient_status = (
        {
            row.get("patient_id", ""): row.get("patient_status", "")
            for row in read_tsv(args.manifest_summary)
        }
        if args.manifest_summary
        else {}
    )

    final_rows = []
    seen = set()
    for pmhc in pmhc_rows:
        key = candidate_key(pmhc)
        if key in seen:
            continue
        seen.add(key)
        peptide = pmhc.get("peptide", "")
        mutation_row = mutation_maf.get(pmhc.get("variant_or_junction_id", ""), {}) if pmhc.get("source_type") == "mutation" else {}
        legacy_status = (
            mutation_row.get("variant_status")
            or mutation_row.get("variant_class")
            or pmhc.get("variant_status")
            or pmhc.get("variant_class", "")
        )
        variant_status = "somatic_candidate" if legacy_status == "tumor_variant" else legacy_status
        iedb_row = iedb.get(key) or iedb.get(("", "", "", peptide, "")) or {}
        immune_row = immunopeptidome.get(peptide, {})
        tabr_row = tabr.get(tabr_key(pmhc), {})
        expression_row = expression.get(key, {})

        final_rows.append(
            {
                "patient_id": pmhc.get("patient_id", ""),
                "patient_status": patient_status.get(pmhc.get("patient_id", ""), ""),
                "source_type": pmhc.get("source_type", ""),
                "gene": pmhc.get("gene", ""),
                "transcript_id": pmhc.get("transcript_id", ""),
                "variant_or_junction_id": pmhc.get("variant_or_junction_id", ""),
                "variant_status": pmhc.get("variant_status") or variant_status,
                "variant_class": pmhc.get("variant_class") or mutation_row.get("variant_class", ""),
                "Variant_Classification": pmhc.get("Variant_Classification") or mutation_row.get("Variant_Classification", ""),
                "peptide_sequence": peptide,
                "peptide_length": pmhc.get("peptide_length", len(peptide)),
                "HLA_allele": pmhc.get("hla_allele", ""),
                "Aff.nM": pmhc.get("Aff.nM", ""),
                "%Rank_EL": pmhc.get("rank_el", ""),
                "bind_level": pmhc.get("bind_level", ""),
                "IEDB_immunogenicity_score": iedb_row.get("immunogenicity_score", ""),
                "PCI_DB_enabled": immune_row.get("PCI_DB_enabled", "false"),
                "healthy_exact_match": immune_row.get("healthy_exact_match", "false"),
                "tumor_immune_peptide": immune_row.get("tumor_immune_peptide", ""),
                "best_match_peptide": immune_row.get("best_match_peptide", ""),
                "blast_identity": immune_row.get("blast_identity", ""),
                "blast_evalue": immune_row.get("blast_evalue", ""),
                "max_TABR_BERT_rank": tabr_row.get("max_TABR_BERT_rank", ""),
                "best_tcr_cell_id": tabr_row.get("best_tcr_cell_id", ""),
                "best_cdr3_beta_aa": tabr_row.get("best_cdr3_beta_aa", ""),
                "RNA_expression_support": expression_row.get("RNA_expression_support", ""),
                "expression_mapping_status": expression_row.get("expression_mapping_status", ""),
                "expression_supporting_cell_ids": expression_row.get("expression_supporting_cell_ids", ""),
                "expression_missing_cell_ids": expression_row.get("expression_missing_cell_ids", ""),
                "supporting_tumor_cell_ids": pmhc.get("supporting_tumor_cell_ids", ""),
                "normal_alt_cell_ids": pmhc.get("normal_alt_cell_ids", ""),
                "normal_callable_cell_ids": pmhc.get("normal_callable_cell_ids", ""),
                "normal_evidence_status": pmhc.get("normal_evidence_status", ""),
                "tumor_specific_candidate": pmhc.get("tumor_specific_candidate", ""),
                "tumor_enriched_candidate": pmhc.get("tumor_enriched_candidate", ""),
                "annotation_summary": pmhc.get("annotation_summary", ""),
            }
        )

    write_tsv(args.out, final_rows, FINAL_FIELDS)
    by_patient_dir = Path(args.by_patient_dir)
    by_patient_dir.mkdir(parents=True, exist_ok=True)
    by_patient = defaultdict(list)
    for row in final_rows:
        by_patient[row["patient_id"]].append(row)
    for patient_id, rows in by_patient.items():
        write_tsv(by_patient_dir / f"{patient_id}.tsv", rows, FINAL_FIELDS)

    counts = Counter((row["patient_id"], row["source_type"]) for row in final_rows)
    summary_rows = []
    for patient_id in sorted(by_patient):
        rows = by_patient[patient_id]
        summary_rows.append(
            {
                "patient_id": patient_id,
                "patient_status": patient_status.get(patient_id, ""),
                "integrated_candidates": len(rows),
                "mutation_candidates": counts[(patient_id, "mutation")],
                "splicing_candidates": counts[(patient_id, "splicing")],
                "immunopeptidome_enabled": (
                    "true" if any(row["PCI_DB_enabled"] == "true" for row in rows) else "false"
                ),
            }
        )
    write_tsv(
        args.summary,
        summary_rows,
        [
            "patient_id",
            "patient_status",
            "integrated_candidates",
            "mutation_candidates",
            "splicing_candidates",
            "immunopeptidome_enabled",
        ],
    )


def candidate_key(row: dict[str, str]) -> tuple[str, str, str, str, str]:
    return (
        row.get("patient_id", ""),
        row.get("source_type", ""),
        row.get("variant_or_junction_id", ""),
        row.get("peptide") or row.get("peptide_sequence", ""),
        row.get("hla_allele") or row.get("HLA_allele", ""),
    )


def tabr_key(row: dict[str, str]) -> tuple[str, str, str]:
    return (
        row.get("patient_id", ""),
        row.get("peptide") or row.get("peptide_sequence", ""),
        row.get("hla_allele") or row.get("HLA_allele", ""),
    )


def index_iedb(rows: list[dict[str, str]]) -> dict[tuple[str, str, str, str, str], dict[str, str]]:
    out = {}
    peptide_fallback = {}
    for row in rows:
        if not row.get("peptide") or row.get("immunogenicity_score", "") == "":
            continue
        key = candidate_key(row)
        out.setdefault(key, row)
        peptide_fallback.setdefault(("", "", "", row["peptide"], ""), row)
    out.update(peptide_fallback)
    return out


def legacy_pcidb_evidence(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    """Convert the historical PCI-DB BLAST table to the generic evidence schema."""
    out: dict[str, dict[str, str]] = {}
    for row in rows:
        peptide = row.get("candidate_peptide", "")
        if not peptide:
            continue
        candidate = {
            "candidate_peptide": peptide,
            "PCI_DB_enabled": "true",
            "healthy_exact_match": "false",
            "tumor_immune_peptide": row.get("matched_pcidb_peptide", ""),
            "best_match_peptide": row.get("matched_pcidb_peptide", ""),
            "blast_identity": row.get("percent_identity", ""),
            "blast_evalue": row.get("evalue", ""),
        }
        current = out.get(peptide)
        if current is None:
            out[peptide] = candidate
            continue
        try:
            current_identity = float(current.get("blast_identity") or 0)
        except ValueError:
            current_identity = 0.0
        try:
            candidate_identity = float(candidate.get("blast_identity") or 0)
        except ValueError:
            candidate_identity = 0.0
        if candidate_identity > current_identity:
            out[peptide] = candidate
    return out


if __name__ == "__main__":
    main()
