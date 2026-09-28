#!/usr/bin/env python3
"""Append IEDB immunogenicity scores to filtered netMHCpan pMHC rows."""

from __future__ import annotations

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path

from scneo_utils import read_tsv, write_tsv


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pmhc", required=True)
    parser.add_argument("--script", required=True)
    parser.add_argument("--python", default="python3")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    pmhc = read_tsv(args.pmhc)
    peptides = sorted({row["peptide"] for row in pmhc if row.get("peptide")})
    scores = run_script(args.script, args.python, peptides)
    fields = list(pmhc[0]) if pmhc else ["peptide"]
    if "immunogenicity_score" not in fields:
        fields.append("immunogenicity_score")
    rows = annotate_rows(pmhc, scores)
    write_tsv(args.out, rows, fields)


def annotate_rows(rows: list[dict[str, str]], scores: dict[str, str]) -> list[dict[str, str]]:
    """Preserve each netMHCpan row and append its peptide-level score."""
    return [
        {**row, "immunogenicity_score": scores.get(row.get("peptide", ""), "")}
        for row in rows
    ]


def run_script(script: str, python_cmd: str, peptides: list[str]) -> dict[str, str]:
    with tempfile.TemporaryDirectory() as tmpdir:
        input_path = Path(tmpdir) / "peptides.txt"
        input_path.write_text("\n".join(peptides) + "\n")
        proc = subprocess.run(
            [python_cmd, script, str(input_path)],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(f"immunogenicity script failed: {proc.stderr.strip() or 'no stderr output'}")
    return parse_script_output(proc.stdout)


def parse_script_output(output: str) -> dict[str, str]:
    """Parse the downloaded script's metadata plus comma-separated result table."""
    scores = {}
    for line in output.splitlines():
        if not line or line.startswith(("masking:", "masked variables:", "allele:")):
            continue
        columns = next(csv.reader([line]), [])
        if columns and columns[0] == "peptide" and len(columns) >= 3:
            continue
        if len(columns) >= 3 and columns[0]:
            scores[columns[0]] = columns[2]
    return scores


if __name__ == "__main__":
    main()
