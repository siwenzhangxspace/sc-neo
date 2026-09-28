#!/usr/bin/env python3
"""Compare candidate peptides against PCI-DB HCC MHC-I peptides with BLASTP-short."""

from __future__ import annotations

import argparse
import hashlib
import re
import subprocess
import tempfile
from pathlib import Path

from scneo_utils import parse_lengths, read_tsv, write_tsv


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pmhc", required=True)
    parser.add_argument("--pcidb", required=True)
    parser.add_argument("--makeblastdb-cmd", default="makeblastdb")
    parser.add_argument("--blastp-cmd", default="blastp")
    parser.add_argument("--extra-args", default="")
    parser.add_argument("--lengths", default="8,9,10,11")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    pmhc = read_tsv(args.pmhc)
    candidates = sorted({row["peptide"] for row in pmhc if row.get("peptide")})
    peptide_sources = {}
    for row in pmhc:
        peptide_sources.setdefault(row.get("peptide", ""), set()).update([cell for cell in row.get("source_cell_ids", "").split(",") if cell])
    pcidb = load_pcidb(args.pcidb, parse_lengths(args.lengths))
    rows = run_blast(candidates, pcidb, args)
    for row in rows:
        row["source_cell_ids"] = ",".join(sorted(peptide_sources.get(row.get("candidate_peptide", ""), set())))
    write_tsv(
        args.out,
        rows,
        [
            "candidate_peptide",
            "source_cell_ids",
            "matched_pcidb_peptide",
            "pcidb_hla_allele",
            "pcidb_tissue_disease",
            "percent_identity",
            "alignment_length",
            "mismatch",
            "gap_open",
            "evalue",
            "bit_score",
        ],
    )


def load_pcidb(path: str, lengths: list[int]) -> list[dict[str, str]]:
    rows = []
    for row in read_tsv(path):
        peptide = first(row, ["Peptide", "peptide", "Sequence", "sequence"])
        disease = first(row, ["Disease", "disease", "Tissue", "tissue"])
        mhc_class = first(row, ["MHC Class", "MHC_Class", "mhc_class"])
        modified = first(row, ["Modified", "Modification", "modified"])
        if not peptide or len(peptide) not in lengths:
            continue
        if "hepatocellular" not in disease.lower() and "hcc" not in disease.lower():
            continue
        mhc_norm = mhc_class.lower().replace("class", "").replace("-", " ").strip()
        if "ii" in mhc_norm or ("i" not in mhc_norm and "1" not in mhc_norm):
            continue
        if modified and modified.lower() not in {"no", "none", "false", "0", "unmodified"}:
            continue
        rows.append(
            {
                "peptide": peptide,
                "hla": first(row, ["HLA", "HLA allele", "MHC", "Allele"]),
                "disease": disease,
            }
        )
    return rows


def run_blast(candidates: list[str], pcidb: list[dict[str, str]], args: argparse.Namespace) -> list[dict[str, str]]:
    if not candidates or not pcidb:
        return []
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        db_fasta = tmp / "pcidb.fa"
        query_fasta = tmp / "candidates.fa"
        db_fasta.write_text("".join(f">pcidb_{i}|{row['hla']}|{row['disease']}\n{row['peptide']}\n" for i, row in enumerate(pcidb)))
        query_fasta.write_text("".join(f">{sanitize(peptide)}\n{peptide}\n" for peptide in candidates))
        db_prefix = tmp / "pcidb"
        make = subprocess.run([args.makeblastdb_cmd, "-in", str(db_fasta), "-dbtype", "prot", "-out", str(db_prefix)], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        if make.returncode != 0:
            raise SystemExit(f"makeblastdb failed: {make.stderr}")
        cmd = [
            args.blastp_cmd,
            "-task",
            "blastp-short",
            "-query",
            str(query_fasta),
            "-db",
            str(db_prefix),
            "-outfmt",
            "6 qseqid sseqid pident length mismatch gapopen evalue bitscore",
        ]
        if args.extra_args:
            cmd.extend(args.extra_args.split())
        blast = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        if blast.returncode != 0:
            raise SystemExit(f"blastp failed: {blast.stderr}")
        id_to_meta = {f"pcidb_{i}": row for i, row in enumerate(pcidb)}
        rows = []
        for line in blast.stdout.splitlines():
            qseqid, sseqid, pident, length, mismatch, gapopen, evalue, bitscore = line.split("\t")
            subject = sseqid.split("|")[0]
            meta = id_to_meta.get(subject, {})
            peptide = unsanitize(qseqid, candidates)
            rows.append(
                {
                    "candidate_peptide": peptide,
                    "source_cell_ids": "",
                    "matched_pcidb_peptide": meta.get("peptide", ""),
                    "pcidb_hla_allele": meta.get("hla", ""),
                    "pcidb_tissue_disease": meta.get("disease", ""),
                    "percent_identity": pident,
                    "alignment_length": length,
                    "mismatch": mismatch,
                    "gap_open": gapopen,
                    "evalue": evalue,
                    "bit_score": bitscore,
                }
            )
    return rows


def first(row: dict[str, str], columns: list[str]) -> str:
    for col in columns:
        if row.get(col):
            return row[col]
    return ""


def sanitize(peptide: str) -> str:
    return "pep_" + hashlib.md5(peptide.encode()).hexdigest()


def unsanitize(qseqid: str, candidates: list[str]) -> str:
    mapping = {sanitize(peptide): peptide for peptide in candidates}
    return mapping.get(qseqid, qseqid)


if __name__ == "__main__":
    main()
