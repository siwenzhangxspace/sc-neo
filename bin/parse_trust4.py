#!/usr/bin/env python3
"""Parse TRUST4 Smart-seq report.tsv files and select each cell's dominant beta CDR3."""

from __future__ import annotations

import argparse
import csv

from collections import defaultdict

from scneo_utils import die, open_text, read_manifest, write_tsv


REPORT_COLUMNS = ["#count", "frequency", "CDR3nt", "CDR3aa", "V", "D", "J", "C", "cid", "cid_full_length"]
OUTPUT_COLUMNS = [
    "patient_id",
    "cell_id",
    "cdr3_beta_aa",
    "beta_count",
    "alpha_chain_detected",
    "alpha_beta_tcell",
    "trust4_status",
]
INVALID_CDR3_AA_MARKERS = {"?", "_", "*"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--reports", nargs="*", default=[])
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    manifest_rows = read_manifest(args.manifest)
    manifest = {row["cell_id"]: row for row in manifest_rows}
    tcell_ids = {cell_id for cell_id, row in manifest.items() if row["cell_type"] == "tcell"}
    report_rows_by_cell = defaultdict(list)
    for report in args.reports:
        for row in read_report(report):
            cell_id = cell_id_from_cid(row["cid"], tcell_ids, report)
            report_rows_by_cell[cell_id].append(row)

    rows = []
    for meta in manifest_rows:
        if meta["cell_type"] != "tcell":
            continue
        cell_id = meta["cell_id"]
        if cell_id in report_rows_by_cell:
            parsed = parse_report(report_rows_by_cell[cell_id])
            parsed.update({"patient_id": meta["patient_id"], "cell_id": cell_id})
            rows.append(parsed)
        else:
            rows.append(
                {
                    "patient_id": meta["patient_id"],
                    "cell_id": cell_id,
                    "cdr3_beta_aa": "",
                    "beta_count": "",
                    "alpha_chain_detected": "false",
                    "alpha_beta_tcell": "false",
                    "trust4_status": "missing_report",
                }
            )

    write_tsv(args.out, rows, OUTPUT_COLUMNS)


def read_report(path: str) -> list[dict[str, object]]:
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        missing = [column for column in REPORT_COLUMNS if column not in fieldnames]
        if missing:
            die(f"TRUST4 report {path} is missing required columns: {','.join(missing)}")

        rows = []
        for line_number, row in enumerate(reader, start=2):
            count = row.get("#count", "")
            try:
                row["#count"] = float(count)
            except (TypeError, ValueError):
                die(f"TRUST4 report {path} line {line_number} has invalid #count: {count}")
            rows.append(row)
        return rows


def cell_id_from_cid(cid: object, tcell_ids: set[str], report: str) -> str:
    cid_text = str(cid or "").strip()
    for cell_id in sorted(tcell_ids, key=len, reverse=True):
        if cid_text.startswith(f"{cell_id}_"):
            return cell_id
    die(f"TRUST4 report {report} cid {cid_text!r} does not match a manifest T-cell")


def parse_report(report_rows: list[dict[str, object]]) -> dict[str, str]:
    beta = []
    alpha_detected = False
    for row in report_rows:
        genes = [str(row.get(column) or "").upper() for column in ["V", "D", "J", "C"]]
        if any(gene.startswith("TRA") for gene in genes):
            alpha_detected = True
        cdr3_aa = str(row.get("CDR3aa") or "")
        if (
            any(gene.startswith("TRB") for gene in genes)
            and cdr3_aa
            and not any(marker in cdr3_aa for marker in INVALID_CDR3_AA_MARKERS)
        ):
            beta.append((row["#count"], cdr3_aa))
    if not beta:
        return {
            "cdr3_beta_aa": "",
            "beta_count": "",
            "alpha_chain_detected": "true" if alpha_detected else "false",
            "alpha_beta_tcell": "false",
            "trust4_status": "no_valid_beta_cdr3",
        }
    read_count, cdr3_aa = max(beta, key=lambda candidate: candidate[0])
    return {
        "cdr3_beta_aa": cdr3_aa,
        "beta_count": f"{read_count:.6g}",
        "alpha_chain_detected": "true" if alpha_detected else "false",
        "alpha_beta_tcell": "true" if alpha_detected else "false",
        "trust4_status": "ok",
    }


if __name__ == "__main__":
    main()
