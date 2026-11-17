#!/usr/bin/env python3
"""Fixed-ideal scNeo TOPSIS and stable cohort-wide top-decile selection."""

from __future__ import annotations

import argparse
import csv
import math

from scneo_utils import _tsv_data_lines, open_text, write_tsv
from candidate_utils import candidate_key

RNA_REFERENCE = 198.738794605166
MODULE_FIELDS = {
    "hla": "topsis_hla_presentation",
    "immunopeptidome": "topsis_immunopeptidome",
    "tabr": "topsis_tabr",
    "rna": "topsis_rna_expression",
    "iedb": "topsis_iedb_immunogenicity",
}
WEIGHTS = {
    "five_modules": [0.25, 0.20, 0.25, 0.15, 0.15],
    "TCR_missing_40_30_15_15": [0.40, 0.30, 0.00, 0.15, 0.15],
    "PCI_missing_35_35_15_15": [0.35, 0.00, 0.35, 0.15, 0.15],
}
FINAL_FIELDS_SUFFIX = list(MODULE_FIELDS.values()) + [
    "active_module_n", "available_module_n", "missing_modules", "score_eligible",
    "evidence_completeness", "topsis_relative_closeness",
    "topsis_completeness_adjusted_score", "TOPSIS_score", "priority_class",
    "weight_policy", "cohort_rank", "top_decile", "RNA_reference_TPM",
]


def number(value, name, minimum=None, maximum=None):
    if value is None or str(value).strip() == "":
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite or an empty missing value")
    if minimum is not None and result < minimum:
        raise ValueError(f"{name} must be >= {minimum}")
    if maximum is not None and result > maximum:
        raise ValueError(f"{name} must be <= {maximum}")
    return result


def flag(value):
    if value in (True, "true", "True", "1", 1):
        return True
    if value in (False, "false", "False", "0", 0, "", None):
        return False
    raise ValueError(f"Invalid boolean evidence {value!r}")


def pci_module(row):
    if not flag(row.get("PCI_DB_enabled")):
        return None
    status = row.get("pci_status", "")
    if not flag(row.get("search_completed")) or status not in {"hit", "no_hit"}:
        raise ValueError("Completed PCI evidence required. Regenerate immunopeptidome evidence.")
    if flag(row.get("healthy_blast_match")):
        raise ValueError("Healthy full-length exact matches must be excluded before scoring")
    support = 0.0
    if status == "hit":
        identity = number(row.get("blast_identity"), "blast_identity", 0, 100)
        coverage = number(row.get("query_coverage"), "query_coverage", 0, 1)
        evalue = number(row.get("blast_evalue"), "blast_evalue", 0)
        if None in (identity, coverage, evalue):
            raise ValueError("PCI hit lacks identity, coverage or E-value. Regenerate evidence.")
        support = identity / 100 * coverage / (1 + evalue)
    if flag(row.get("tumor_full_length_exact")):
        support = 1.0
    return 0.5 + 0.5 * support


def module_values(row, rna_reference):
    if row.get("tabr_bert_status") == "missing_score":
        raise ValueError("TCR prediction is incomplete. Re-run prediction before scoring.")
    affinity = number(row.get("Aff.nM"), "Aff.nM", 0)
    rank = number(row.get("%Rank_EL"), "%Rank_EL", 0)
    presentation = None
    if affinity is not None and rank is not None:
        if affinity == 0:
            raise ValueError("Aff.nM must be positive")
        presentation = (150 / (150 + affinity) + 0.5 / (0.5 + rank)) / 2
    tpm = number(row.get("mean_malignant_TPM"), "mean_malignant_TPM", 0)
    if row.get("RNA_expression_support") not in ("", None) and "mean_malignant_TPM" not in row:
        raise ValueError("Legacy RNA_expression_support has no TPM scale. Regenerate expression support.")
    if tpm is not None and row.get("expression_unit") != "TPM":
        raise ValueError("mean_malignant_TPM requires expression_unit=TPM")
    expression = min(1, math.log1p(tpm) / math.log1p(rna_reference)) if tpm is not None else None
    intrinsic = number(row.get("IEDB_immunogenicity_score"), "IEDB_immunogenicity_score")
    intrinsic = min(1, max(0, (intrinsic + 1) / 2)) if intrinsic is not None else None
    return [presentation, pci_module(row),
            number(row.get("max_TABR_BERT_rank"), "max_TABR_BERT_rank", 0, 1),
            expression, intrinsic]


def rank_cohort(rows, rna_reference=RNA_REFERENCE):
    if not math.isfinite(rna_reference) or rna_reference <= 0:
        raise ValueError("RNA reference must be positive and finite")
    keys = ["|".join(candidate_key(row)) for row in rows]
    if len(set(keys)) != len(keys):
        raise ValueError("Duplicate key7 in TOPSIS input. Integrate candidates first.")
    out = []
    eligible = []
    for key, source in zip(keys, rows):
        values = module_values(source, rna_reference)
        policy = "five_modules"
        if values[2] is None:
            policy = "TCR_missing_40_30_15_15"
        elif values[1] is None:
            policy = "PCI_missing_35_35_15_15"
        weights = WEIGHTS[policy]
        observed = sum(value is not None for value in values)
        valid = observed >= 4 and all(value is not None for value, weight in zip(values, weights) if weight)
        row = dict(source)
        row.update({field: fmt(value) for field, value in zip(MODULE_FIELDS.values(), values)})
        row.update({
            "active_module_n": sum(weight > 0 for weight in weights),
            "available_module_n": observed,
            "missing_modules": ";".join(name for name, value in zip(MODULE_FIELDS, values) if value is None),
            "score_eligible": str(valid).lower(), "weight_policy": policy,
            "evidence_completeness": fmt(observed / len(values)),
            "topsis_relative_closeness": "", "topsis_completeness_adjusted_score": "",
            "TOPSIS_score": "", "priority_class": "Not scored", "cohort_rank": "",
            "top_decile": "false", "RNA_reference_TPM": fmt(rna_reference),
        })
        if valid:
            d0 = math.sqrt(sum((value * weight) ** 2 for value, weight in zip(values, weights) if weight))
            d1 = math.sqrt(sum(((1 - value) * weight) ** 2 for value, weight in zip(values, weights) if weight))
            closeness = d0 / (d0 + d1)
            row.update(TOPSIS_score=fmt(100 * closeness), topsis_relative_closeness=fmt(closeness),
                       priority_class="Other")
            eligible.append((100 * closeness, key, row))
        out.append(row)
    eligible.sort(key=lambda item: (-item[0], item[1]))
    budget = (len(eligible) + 9) // 10
    for rank, (_, _, row) in enumerate(eligible, 1):
        row["cohort_rank"] = rank
        row["top_decile"] = str(rank <= budget).lower()
        row["priority_class"] = "High" if rank <= budget else "Other"
    return out


def fmt(value):
    return "" if value is None else format(value, ".17g")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--rna-reference", type=float, default=RNA_REFERENCE)
    args = parser.parse_args()
    with open_text(args.candidates) as handle:
        reader = csv.DictReader(_tsv_data_lines(handle), delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    try:
        ranked = rank_cohort(rows, args.rna_reference)
    except ValueError as exc:
        parser.error(str(exc))
    fields.extend(field for field in FINAL_FIELDS_SUFFIX if field not in fields)
    write_tsv(args.out, ranked, fields)


if __name__ == "__main__":
    main()
