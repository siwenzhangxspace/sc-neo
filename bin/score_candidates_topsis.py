#!/usr/bin/env python3
"""Rank integrated neoantigen candidates with per-patient equal-weight TOPSIS."""

from __future__ import annotations

import argparse
import math
import statistics
from collections import defaultdict

from scneo_utils import read_tsv, safe_float, write_tsv


MODULE_OUTPUT_FIELDS = [
    "topsis_hla_presentation",
    "topsis_iedb_immunogenicity",
    "topsis_immunopeptidome",
    "topsis_tabr",
    "topsis_rna_expression",
]
FINAL_FIELDS_SUFFIX = MODULE_OUTPUT_FIELDS + ["TOPSIS_score", "priority_class"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    rows = read_tsv(args.candidates)
    by_patient: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_patient[row.get("patient_id", "")].append(row)

    ranked = []
    for patient_id in sorted(by_patient):
        ranked.extend(rank_patient(by_patient[patient_id]))

    fields = list(rows[0].keys()) if rows else []
    for field in FINAL_FIELDS_SUFFIX:
        if field not in fields:
            fields.append(field)
    write_tsv(args.out, ranked, fields)


def rank_patient(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    pci_enabled = any(row.get("PCI_DB_enabled") == "true" for row in rows)

    affinity = [safe_float(row.get("Aff.nM"), None) for row in rows]
    rank_el = [safe_float(row.get("%Rank_EL"), None) for row in rows]
    aff_benefit = minmax(affinity, benefit=False)
    rank_benefit = minmax(rank_el, benefit=False)

    raw_modules: dict[str, list[float | None]] = {
        "hla": [
            mean_present([aff_benefit[index], rank_benefit[index]])
            for index in range(len(rows))
        ],
        "iedb": [safe_float(row.get("IEDB_immunogenicity_score"), None) for row in rows],
        "tabr": [safe_float(row.get("max_TABR_BERT_rank"), None) for row in rows],
        "rna": [safe_float(row.get("RNA_expression_support"), None) for row in rows],
    }
    if pci_enabled:
        identities = [
            safe_float(row.get("blast_identity"), 0.0)
            if row.get("tumor_immune_peptide")
            else 0.0
            for row in rows
        ]
        evalue_strength = [
            evalue_score(row.get("blast_evalue")) if row.get("tumor_immune_peptide") else 0.0
            for row in rows
        ]
        identity_norm = minmax(identities, benefit=True)
        evalue_norm = minmax(evalue_strength, benefit=True)
        raw_modules["immunopeptidome"] = [
            mean_present([identity_norm[index], evalue_norm[index]])
            for index in range(len(rows))
        ]

    order = ["hla", "iedb"]
    if pci_enabled:
        order.append("immunopeptidome")
    order.extend(["tabr", "rna"])

    modules: dict[str, list[float]] = {}
    for name in order:
        imputed = median_impute(raw_modules[name])
        modules[name] = minmax(imputed, benefit=True, constant=0.5)

    weighted: dict[str, list[float]] = {}
    weight = 1.0 / len(order)
    for name in order:
        norm = vector_normalize(modules[name])
        weighted[name] = [value * weight for value in norm]

    positive = {name: max(weighted[name]) if weighted[name] else 0.0 for name in order}
    negative = {name: min(weighted[name]) if weighted[name] else 0.0 for name in order}

    out = []
    for index, source in enumerate(rows):
        d_positive = math.sqrt(
            sum((weighted[name][index] - positive[name]) ** 2 for name in order)
        )
        d_negative = math.sqrt(
            sum((weighted[name][index] - negative[name]) ** 2 for name in order)
        )
        denom = d_positive + d_negative
        score = 50.0 if denom == 0 else d_negative / denom * 100.0
        row: dict[str, object] = dict(source)
        row["topsis_hla_presentation"] = fmt(modules["hla"][index])
        row["topsis_iedb_immunogenicity"] = fmt(modules["iedb"][index])
        row["topsis_immunopeptidome"] = (
            fmt(modules["immunopeptidome"][index]) if pci_enabled else ""
        )
        row["topsis_tabr"] = fmt(modules["tabr"][index])
        row["topsis_rna_expression"] = fmt(modules["rna"][index])
        row["TOPSIS_score"] = f"{score:.6g}"
        row["priority_class"] = priority_class(score)
        out.append(row)
    return out


def evalue_score(value: object) -> float:
    numeric = safe_float(value, None)
    if numeric is None:
        return 0.0
    if numeric <= 0:
        return 300.0
    return max(0.0, -math.log10(numeric))


def mean_present(values: list[float | None]) -> float | None:
    present = [value for value in values if value is not None]
    return sum(present) / len(present) if present else None


def median_impute(values: list[float | None]) -> list[float]:
    present = [float(value) for value in values if value is not None and math.isfinite(value)]
    fill = statistics.median(present) if present else 0.5
    return [
        float(value) if value is not None and math.isfinite(value) else fill
        for value in values
    ]


def minmax(
    values: list[float | None],
    benefit: bool,
    constant: float = 0.5,
) -> list[float | None]:
    present = [float(value) for value in values if value is not None and math.isfinite(value)]
    if not present:
        return [None for _ in values]
    low, high = min(present), max(present)
    if high == low:
        return [constant if value is not None else None for value in values]

    out = []
    for value in values:
        if value is None or not math.isfinite(value):
            out.append(None)
            continue
        scaled = (float(value) - low) / (high - low)
        out.append(scaled if benefit else 1.0 - scaled)
    return out


def vector_normalize(values: list[float]) -> list[float]:
    norm = math.sqrt(sum(value * value for value in values))
    if norm == 0:
        return [0.0 for _ in values]
    return [value / norm for value in values]


def priority_class(score: float) -> str:
    if score >= 70:
        return "High"
    if score >= 50:
        return "Intermediate"
    return "Low"


def fmt(value: float) -> str:
    return f"{value:.10g}"


if __name__ == "__main__":
    main()
