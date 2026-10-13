#!/usr/bin/env python3
"""Integrate candidate pMHC, immunogenicity, immunopeptidome, TABR-BERT and RNA evidence."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path

from scneo_utils import read_tsv, write_tsv
from candidate_utils import (
    candidate_key, clean, consistent_row, evidence_index,
    grouped_candidates, supporting_cells,
)


FINAL_FIELDS = [
    "patient_id",
    "patient_status",
    "source_type",
    "gene",
    "variant_or_junction_id",
    "var_key",
    "HGVSp_Short",
    "protein_change_status",
    "n_original_rows",
    "n_barcode",
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
    "healthy_blast_match",
    "healthy_match_peptide",
    "healthy_blast_identity",
    "healthy_blast_evalue",
    "healthy_blast_alignment_length",
    "tumor_immune_peptide",
    "best_match_peptide",
    "blast_identity",
    "blast_evalue",
    "max_TABR_BERT_rank",
    "best_tcr_cell_id",
    "best_cdr3_beta_aa",
    "RNA_expression_support",
    "expression_mapping_status",
    "normal_evidence_status",
    "tumor_specific_candidate",
    "tumor_enriched_candidate",
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

    groups = grouped_candidates(read_tsv(args.pmhc))
    mutation_maf = index_maf(read_tsv(args.mutation_maf)) if args.mutation_maf else {}
    iedb = index_iedb(read_tsv(args.iedb))
    if args.immunopeptidome:
        immunopeptidome = evidence_index(
            read_tsv(args.immunopeptidome), lambda row: clean(row.get("candidate_peptide")),
            FINAL_FIELDS[FINAL_FIELDS.index("PCI_DB_enabled"):FINAL_FIELDS.index("max_TABR_BERT_rank")],
        )
    elif args.pcidb:
        immunopeptidome = legacy_pcidb_evidence(read_tsv(args.pcidb))
    else:
        immunopeptidome = {}
    tabr = evidence_index(read_tsv(args.tabr), tabr_key,
                          ["max_TABR_BERT_rank", "best_tcr_cell_id", "best_cdr3_beta_aa"])
    expression = (
        evidence_index(read_tsv(args.expression), candidate_key,
                       ["RNA_expression_support", "expression_mapping_status"])
        if args.expression
        else {}
    )
    patient_status = (
        {key: row["patient_status"] for key, row in evidence_index(
            read_tsv(args.manifest_summary), lambda row: clean(row.get("patient_id")),
            ["patient_status"]).items()}
        if args.manifest_summary
        else {}
    )

    final_rows = []
    for key, group in groups.items():
        enriched = []
        for record in group:
            annotation = match_maf(record, mutation_maf)
            item = dict(record)
            for field in ["variant_class", "variant_status", "Variant_Classification"]:
                if not clean(item.get(field)):
                    item[field] = annotation.get(field, "")
            enriched.append(item)
        pmhc = consistent_row(enriched, [
            "patient_id", "source_type", "gene", "variant_or_junction_id", "var_key", "HGVSp_Short",
            "protein_change_status", "peptide", "peptide_length", "hla_allele", "Aff.nM", "rank_el",
            "bind_level", "variant_status", "variant_class", "Variant_Classification",
            "normal_evidence_status", "tumor_specific_candidate", "tumor_enriched_candidate",
        ], key)
        peptide = pmhc.get("peptide", "")
        legacy_status = (
            pmhc.get("variant_status")
            or pmhc.get("variant_class", "")
        )
        variant_status = "somatic_candidate" if legacy_status == "tumor_variant" else legacy_status
        iedb_row = iedb.get(peptide, {})
        immune_row = immunopeptidome.get(peptide, {})
        tabr_row = tabr.get(tabr_key(pmhc), {})
        if args.expression and key not in expression:
            raise ValueError(f"Missing expression evidence row for key7 {key!r}; regenerate expression support")
        expression_row = expression.get(key, {})

        final_rows.append(
            {
                "patient_id": pmhc.get("patient_id", ""),
                "patient_status": patient_status.get(pmhc.get("patient_id", ""), ""),
                "source_type": pmhc.get("source_type", ""),
                "gene": pmhc.get("gene", ""),
                "variant_or_junction_id": pmhc.get("variant_or_junction_id", ""),
                "var_key": pmhc["var_key"],
                "HGVSp_Short": pmhc["HGVSp_Short"],
                "protein_change_status": pmhc["protein_change_status"],
                "n_original_rows": len(group),
                "n_barcode": len(supporting_cells(group)),
                "variant_status": pmhc.get("variant_status") or variant_status,
                "variant_class": pmhc.get("variant_class", ""),
                "Variant_Classification": pmhc.get("Variant_Classification", ""),
                "peptide_sequence": peptide,
                "peptide_length": pmhc.get("peptide_length", len(peptide)),
                "HLA_allele": pmhc.get("hla_allele", ""),
                "Aff.nM": pmhc.get("Aff.nM", ""),
                "%Rank_EL": pmhc.get("rank_el", ""),
                "bind_level": pmhc.get("bind_level", ""),
                "IEDB_immunogenicity_score": iedb_row.get("immunogenicity_score", ""),
                "PCI_DB_enabled": immune_row.get("PCI_DB_enabled", "false"),
                "healthy_blast_match": immune_row.get("healthy_blast_match", "false"),
                "healthy_match_peptide": immune_row.get("healthy_match_peptide", ""),
                "healthy_blast_identity": immune_row.get("healthy_blast_identity", ""),
                "healthy_blast_evalue": immune_row.get("healthy_blast_evalue", ""),
                "healthy_blast_alignment_length": immune_row.get("healthy_blast_alignment_length", ""),
                "tumor_immune_peptide": immune_row.get("tumor_immune_peptide", ""),
                "best_match_peptide": immune_row.get("best_match_peptide", ""),
                "blast_identity": immune_row.get("blast_identity", ""),
                "blast_evalue": immune_row.get("blast_evalue", ""),
                "max_TABR_BERT_rank": tabr_row.get("max_TABR_BERT_rank", ""),
                "best_tcr_cell_id": tabr_row.get("best_tcr_cell_id", ""),
                "best_cdr3_beta_aa": tabr_row.get("best_cdr3_beta_aa", ""),
                "RNA_expression_support": expression_row.get("RNA_expression_support", ""),
                "expression_mapping_status": expression_row.get("expression_mapping_status", ""),
                "normal_evidence_status": pmhc.get("normal_evidence_status", ""),
                "tumor_specific_candidate": pmhc.get("tumor_specific_candidate", ""),
                "tumor_enriched_candidate": pmhc.get("tumor_enriched_candidate", ""),
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


def index_maf(rows):
    groups = defaultdict(list)
    for row in rows:
        key = tuple(clean(value) for value in (
            row.get("patient_id"), row.get("variant_key"),
            row.get("Hugo_Symbol") or row.get("Gene") or row.get("SYMBOL"), row.get("HGVSp_Short"),
        ))
        groups[key].append(row)
    return groups


def match_maf(pmhc, index):
    if pmhc["source_type"] != "mutation":
        return {}
    key = (pmhc["patient_id"], pmhc["var_key"], pmhc["gene"], pmhc["HGVSp_Short"])
    rows = index.get(key, [])
    tid = clean(pmhc.get("transcript_id"))
    if tid:
        exact = [row for row in rows if maf_transcript(row) == tid]
        if exact:
            rows = exact
        else:
            # Allow the existing unversioned peptide transcript only if unambiguous.
            compatible = [row for row in rows if maf_transcript(row).split(".")[0] == tid.split(".")[0]
                          and ("." not in tid or "." not in maf_transcript(row))]
            if len({maf_transcript(row) for row in compatible}) > 1:
                raise ValueError(f"Ambiguous MAF transcript for {key}: {tid}")
            rows = compatible
    return consistent_row(rows, ["variant_status", "variant_class", "Variant_Classification"], key)


def maf_transcript(row):
    return next((clean(row.get(field)) for field in ["Transcript_ID", "Feature", "transcript_id", "ENST"]
                 if clean(row.get(field))), "")


def tabr_key(row: dict[str, str]) -> tuple[str, str, str]:
    return (
        clean(row.get("patient_id")),
        clean(row.get("peptide") or row.get("peptide_sequence")),
        clean(row.get("hla_allele") or row.get("HLA_allele")),
    )


def index_iedb(rows):
    return evidence_index(rows, lambda row: clean(row.get("peptide")), ["immunogenicity_score"])


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
            "healthy_blast_match": "false",
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
