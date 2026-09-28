#!/usr/bin/env python3
"""Run TABR-BERT on same-patient pMHC/TCR-beta pairs and keep candidate-level max rank."""

from __future__ import annotations

import argparse
import shlex
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

from scneo_utils import read_tsv, safe_float, write_tsv


OUTPUT_FIELDS = [
    "patient_id",
    "peptide",
    "hla_allele",
    "max_TABR_BERT_rank",
    "best_tcr_cell_id",
    "best_cdr3_beta_aa",
    "tabr_bert_status",
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pmhc", required=True)
    parser.add_argument("--tcr", required=True)
    parser.add_argument("--cmd", default="tabr_bert_predict")
    parser.add_argument("--extra-args", default="")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    pmhc_rows = unique_candidates(read_tsv(args.pmhc))
    tcr_rows = [
        row
        for row in read_tsv(args.tcr)
        if row.get("trust4_status") == "ok" and row.get("cdr3_beta_aa")
    ]
    tcr_by_patient: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in tcr_rows:
        tcr_by_patient[row.get("patient_id", "")].append(row)

    combos = []
    for pmhc in pmhc_rows:
        for tcr_row in tcr_by_patient.get(pmhc.get("patient_id", ""), []):
            combos.append(
                {
                    "patient_id": pmhc.get("patient_id", ""),
                    "peptide": pmhc.get("peptide", ""),
                    "hla_allele": pmhc.get("hla_allele", ""),
                    "cdr3_beta_aa": tcr_row.get("cdr3_beta_aa", ""),
                    "tcr_cell_id": tcr_row.get("cell_id", ""),
                }
            )

    scores = run_wrapper(args.cmd, args.extra_args, combos) if combos else {}
    grouped_scores: dict[tuple[str, str, str], list[tuple[float, dict[str, str]]]] = defaultdict(list)
    for combo in combos:
        value = scores.get(score_key(combo))
        numeric = safe_float(value, None)
        if numeric is not None:
            grouped_scores[candidate_key(combo)].append((numeric, combo))

    rows = []
    for pmhc in pmhc_rows:
        key = candidate_key(pmhc)
        patient_tcrs = tcr_by_patient.get(pmhc.get("patient_id", ""), [])
        candidates = grouped_scores.get(key, [])
        if candidates:
            rank, best = max(candidates, key=lambda item: item[0])
            rows.append(
                {
                    "patient_id": key[0],
                    "peptide": key[1],
                    "hla_allele": key[2],
                    "max_TABR_BERT_rank": f"{rank:.10g}",
                    "best_tcr_cell_id": best.get("tcr_cell_id", ""),
                    "best_cdr3_beta_aa": best.get("cdr3_beta_aa", ""),
                    "tabr_bert_status": "ok",
                }
            )
        else:
            rows.append(
                {
                    "patient_id": key[0],
                    "peptide": key[1],
                    "hla_allele": key[2],
                    "max_TABR_BERT_rank": "",
                    "best_tcr_cell_id": "",
                    "best_cdr3_beta_aa": "",
                    "tabr_bert_status": "missing_score" if patient_tcrs else "no_valid_tcr",
                }
            )
    write_tsv(args.out, rows, OUTPUT_FIELDS)


def unique_candidates(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    out = []
    seen = set()
    for row in rows:
        key = candidate_key(row)
        if key not in seen:
            seen.add(key)
            out.append(row)
    return out


def run_wrapper(cmd: str, extra_args: str, combos: list[dict[str, str]]) -> dict[tuple[str, str, str, str], str]:
    with tempfile.TemporaryDirectory() as tmpdir:
        input_path = Path(tmpdir) / "tabr_bert_input.tsv"
        output_path = Path(tmpdir) / "tabr_bert_output.tsv"
        write_tsv(
            input_path,
            combos,
            ["patient_id", "peptide", "hla_allele", "cdr3_beta_aa", "tcr_cell_id"],
        )
        command = shlex.split(cmd) + [str(input_path), str(output_path)]
        if extra_args:
            command.extend(shlex.split(extra_args))
        proc = subprocess.run(
            command,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        if proc.returncode != 0:
            raise SystemExit(f"TABR-BERT wrapper failed: {proc.stderr}")
        rows = read_tsv(output_path)

    scores = {}
    for row in rows:
        value = first(
            row,
            ["max_TABR_BERT_rank", "max_rank", "tabr_bert_rank", "rank", "Rank",
             "tabr_bert_score", "score", "Score"],
        )
        if value == "":
            continue
        scores[score_key(row)] = value
    return scores


def candidate_key(row: dict[str, str]) -> tuple[str, str, str]:
    return (
        row.get("patient_id", ""),
        row.get("peptide", ""),
        row.get("hla_allele", ""),
    )


def score_key(row: dict[str, str]) -> tuple[str, str, str, str]:
    return (
        row.get("patient_id", ""),
        row.get("peptide", ""),
        row.get("hla_allele", ""),
        row.get("cdr3_beta_aa", ""),
    )


def first(row: dict[str, str], columns: list[str]) -> str:
    for column in columns:
        if row.get(column) not in (None, ""):
            return str(row[column])
    return ""


if __name__ == "__main__":
    main()
