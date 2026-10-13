#!/usr/bin/env python3
"""Rank one cohort of integrated candidates with equal-weight TOPSIS."""

from __future__ import annotations

import argparse
import csv
import math
import statistics

from scneo_utils import _tsv_data_lines, open_text, safe_float, write_tsv
from candidate_utils import candidate_key


MODULE_FIELDS = {
    "hla": "topsis_hla_presentation",
    "iedb": "topsis_iedb_immunogenicity",
    "immunopeptidome": "topsis_immunopeptidome",
    "tabr": "topsis_tabr",
    "rna": "topsis_rna_expression",
}
MODULE_OUTPUT_FIELDS = list(MODULE_FIELDS.values())
FINAL_FIELDS_SUFFIX = MODULE_OUTPUT_FIELDS + [
    "active_module_n", "available_module_n", "missing_modules", "score_eligible",
    "evidence_completeness", "topsis_relative_closeness",
    "topsis_completeness_adjusted_score", "TOPSIS_score", "priority_class",
]
MIN_AVAILABLE_MODULES = 3
EMPTY_MODULE_MEDIAN = 0.5


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    # Preserve the input schema even when it contains no records.
    with open_text(args.candidates) as handle:
        reader = csv.DictReader(_tsv_data_lines(handle), delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    try:
        seen = set()
        for row in rows:
            key = candidate_key(row)
            if key in seen:
                raise ValueError(f"Duplicate key7 in TOPSIS input: {key!r}; integrate candidates first")
            seen.add(key)
        ranked = rank_cohort(rows)
    except ValueError as exc:
        parser.error(str(exc))
    fields.extend(field for field in FINAL_FIELDS_SUFFIX if field not in fields)
    write_tsv(args.out, ranked, fields)


def numeric(value: object, minimum: float | None = None, positive: bool = False) -> float | None:
    result = safe_float(value, None)
    if result is None or not math.isfinite(result):
        return None
    if positive and result <= 0:
        return None
    if minimum is not None and result < minimum:
        return None
    return result


def build_modules(rows: list[dict[str, str]], pci_enabled: bool) -> dict[str, list[float | None]]:
    affinity = [numeric(row.get("Aff.nM"), positive=True) for row in rows]
    aff_score = minmax([math.log10(x) if x is not None else None for x in affinity], False)
    rank_score = minmax([numeric(row.get("%Rank_EL"), 0) for row in rows], False)
    rna = [numeric(row.get("RNA_expression_support"), 0) for row in rows]
    modules = {
        "hla": [mean_present([a, r]) for a, r in zip(aff_score, rank_score)],
        "iedb": minmax([numeric(row.get("IEDB_immunogenicity_score")) for row in rows], True),
        "tabr": minmax([numeric(row.get("max_TABR_BERT_rank"), 0) for row in rows], True),
        "rna": minmax([math.log1p(x) if x is not None else None for x in rna], True),
    }
    if pci_enabled:
        has_hit = [bool((row.get("tumor_immune_peptide") or "").strip()) for row in rows]
        identities = [numeric(row.get("blast_identity"), 0) if hit else None
                      for row, hit in zip(rows, has_hit)]
        identity_max = max((x for x in identities if x is not None), default=100.0)
        identity_scores = [
            x if x is None or identity_max <= 1 else min(1.0, max(0.0, x / 100))
            for x in identities
        ]
        evalues = minmax([evalue_score(row.get("blast_evalue")) if hit else None
                         for row, hit in zip(rows, has_hit)], True)
        modules["immunopeptidome"] = [
            mean_present([identity, evalue]) if hit else 0.0
            for identity, evalue, hit in zip(identity_scores, evalues, has_hit)
        ]
    return modules


def rank_cohort(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    if not rows:
        return []
    states = {row.get("PCI_DB_enabled") for row in rows}
    if len(states) != 1 or not states <= {"true", "false"}:
        raise ValueError("PCI_DB_enabled must be uniformly true or false within one cohort")
    modules = build_modules(rows, states == {"true"})
    order = [name for name in MODULE_FIELDS if name in modules]
    available = [sum(modules[name][i] is not None for name in order) for i in range(len(rows))]
    eligible = [i for i, n in enumerate(available)
                if n >= MIN_AVAILABLE_MODULES and modules["hla"][i] is not None]

    adjusted = []
    closeness = []
    weighted = {
        name: [x / len(order) for x in vector_normalize(median_impute([modules[name][i] for i in eligible]))]
        for name in order
    }
    if eligible:
        positive = {name: max(weighted[name]) for name in order}
        negative = {name: min(weighted[name]) for name in order}
        for j, i in enumerate(eligible):
            dp = math.sqrt(sum((weighted[name][j] - positive[name]) ** 2 for name in order))
            dn = math.sqrt(sum((weighted[name][j] - negative[name]) ** 2 for name in order))
            c = dn / (dp + dn) if dp + dn else 0.5
            closeness.append(c)
            adjusted.append(c * available[i] / len(order) * 100)
    scores = minmax(adjusted, True, constant=0.5)
    scored = {i: (closeness[j], adjusted[j], scores[j] * 100) for j, i in enumerate(eligible)}

    out = []
    for i, source in enumerate(rows):
        row: dict[str, object] = dict(source)
        for name, field in MODULE_FIELDS.items():
            row[field] = fmt(modules[name][i]) if name in modules else ""
        row.update({
            "active_module_n": len(order), "available_module_n": available[i],
            "missing_modules": ";".join(name for name in order if modules[name][i] is None),
            "score_eligible": "true" if i in scored else "false",
            "evidence_completeness": "", "topsis_relative_closeness": "",
            "topsis_completeness_adjusted_score": "", "TOPSIS_score": "",
            "priority_class": "Not scored",
        })
        if i in scored:
            c, fs, score = scored[i]
            row.update({
                "evidence_completeness": fmt(available[i] / len(order)),
                "topsis_relative_closeness": fmt(c),
                "topsis_completeness_adjusted_score": fmt(fs),
                "TOPSIS_score": fmt(score), "priority_class": priority_class(score),
            })
        out.append(row)
    return out


def evalue_score(value: object) -> float | None:
    result = numeric(value, 0)
    return -math.log10(max(result, 1e-300)) if result is not None else None


def mean_present(values: list[float | None]) -> float | None:
    present = [x for x in values if x is not None]
    return sum(present) / len(present) if present else None


def median_impute(values: list[float | None]) -> list[float]:
    present = [x for x in values if x is not None]
    # Explicitly retain the reference R policy for entirely unavailable modules.
    fill = statistics.median(present) if present else EMPTY_MODULE_MEDIAN
    return [x if x is not None else fill for x in values]


def minmax(values: list[float | None], benefit: bool, constant: float | None = None) -> list[float | None]:
    present = [x for x in values if x is not None and math.isfinite(x)]
    if not present:
        return [None for _ in values]
    low, high = min(present), max(present)
    flat = constant if constant is not None else (0.0 if high == 0 else 0.5)
    out = []
    for x in values:
        if x is None or not math.isfinite(x):
            out.append(None)
            continue
        scaled = flat if high == low else (x - low) / (high - low)
        out.append(scaled if benefit else 1 - scaled)
    return out


def vector_normalize(values: list[float]) -> list[float]:
    norm = math.sqrt(sum(x * x for x in values))
    return [x / norm for x in values] if norm else [0.0 for _ in values]


def priority_class(score: float) -> str:
    if score >= 70:
        return "High"
    if score >= 50:
        return "Intermediate"
    return "Low"


def fmt(value: float | None) -> str:
    return "" if value is None else f"{value:.10g}"


if __name__ == "__main__":
    main()
