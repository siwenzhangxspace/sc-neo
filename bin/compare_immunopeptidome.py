#!/usr/bin/env python3
"""Compare mutation/splicing candidates with optional immunopeptidome references."""

from __future__ import annotations

import argparse
import hashlib
import shlex
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

from scneo_utils import parse_lengths, read_tsv, safe_float, write_tsv


EVIDENCE_FIELDS = [
    "candidate_peptide",
    "peptide_length",
    "PCI_DB_enabled",
    "healthy_exact_match",
    "tumor_immune_peptide",
    "best_match_peptide",
    "blast_identity",
    "blast_evalue",
    "blast_bit_score",
    "reference_hla_allele",
    "reference_disease_tissue",
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mutation", required=True)
    parser.add_argument("--splicing", required=True)
    parser.add_argument("--healthy")
    parser.add_argument("--tumor")
    parser.add_argument("--disabled", action="store_true")
    parser.add_argument("--makeblastdb-cmd", default="makeblastdb")
    parser.add_argument("--blastp-cmd", default="blastp")
    parser.add_argument("--extra-args", default="")
    parser.add_argument("--lengths", default="8,9,10,11")
    parser.add_argument("--out-mutation", required=True)
    parser.add_argument("--out-splicing", required=True)
    parser.add_argument("--out-evidence", required=True)
    args = parser.parse_args()

    mutation_rows = read_tsv(args.mutation)
    splicing_rows = read_tsv(args.splicing)
    all_rows = mutation_rows + splicing_rows
    lengths = set(parse_lengths(args.lengths))
    candidates = sorted(
        {
            peptide(row)
            for row in all_rows
            if peptide(row) and len(peptide(row)) in lengths
        },
        key=lambda value: (len(value), value),
    )

    if args.disabled:
        evidence = [
            empty_evidence(candidate, enabled=False)
            for candidate in candidates
        ]
        write_filtered(args.out_mutation, mutation_rows, set())
        write_filtered(args.out_splicing, splicing_rows, set())
        write_tsv(args.out_evidence, evidence, EVIDENCE_FIELDS)
        return

    if not args.healthy or not args.tumor:
        raise SystemExit("--healthy and --tumor are required unless --disabled is set")

    healthy = load_reference(args.healthy, lengths)
    tumor = load_reference(args.tumor, lengths)

    healthy_exact = {
        candidate
        for candidate in candidates
        if candidate in {row["peptide"] for row in healthy.get(len(candidate), [])}
    }
    retained = [candidate for candidate in candidates if candidate not in healthy_exact]
    tumor_hits = blast_best_hits(retained, tumor, args)

    evidence = []
    for candidate in candidates:
        if candidate in healthy_exact:
            row = empty_evidence(candidate, enabled=True)
            row["healthy_exact_match"] = "true"
        else:
            row = empty_evidence(candidate, enabled=True)
            hit = tumor_hits.get(candidate)
            if hit:
                row.update(hit)
        evidence.append(row)

    write_filtered(args.out_mutation, mutation_rows, healthy_exact)
    write_filtered(args.out_splicing, splicing_rows, healthy_exact)
    write_tsv(args.out_evidence, evidence, EVIDENCE_FIELDS)


def peptide(row: dict[str, str]) -> str:
    return (row.get("peptide") or row.get("peptide_seq") or "").strip()


def write_filtered(path: str, rows: list[dict[str, str]], excluded: set[str]) -> None:
    fieldnames = list(rows[0].keys()) if rows else [
        "patient_id",
        "source_type",
        "peptide",
        "peptide_length",
    ]
    write_tsv(path, [row for row in rows if peptide(row) not in excluded], fieldnames)


def empty_evidence(candidate: str, enabled: bool) -> dict[str, object]:
    return {
        "candidate_peptide": candidate,
        "peptide_length": len(candidate),
        "PCI_DB_enabled": "true" if enabled else "false",
        "healthy_exact_match": "false",
        "tumor_immune_peptide": "",
        "best_match_peptide": "",
        "blast_identity": "",
        "blast_evalue": "",
        "blast_bit_score": "",
        "reference_hla_allele": "",
        "reference_disease_tissue": "",
    }


def load_reference(path: str, lengths: set[int]) -> dict[int, list[dict[str, str]]]:
    grouped: dict[int, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(path):
        seq = first(row, ["Peptide", "peptide", "Sequence", "sequence", "peptide_sequence"])
        seq = seq.strip().upper()
        if not seq or len(seq) not in lengths:
            continue
        grouped[len(seq)].append(
            {
                "peptide": seq,
                "hla": first(row, ["HLA", "HLA allele", "MHC", "Allele", "hla_allele"]),
                "disease": first(row, ["Disease", "disease", "Tissue", "tissue", "Source", "source"]),
            }
        )
    return grouped


def blast_best_hits(
    candidates: list[str],
    references: dict[int, list[dict[str, str]]],
    args: argparse.Namespace,
) -> dict[str, dict[str, str]]:
    by_length: dict[int, list[str]] = defaultdict(list)
    for candidate in candidates:
        by_length[len(candidate)].append(candidate)

    best: dict[str, dict[str, str]] = {}
    for length, queries in sorted(by_length.items()):
        refs = references.get(length, [])
        if not refs:
            continue
        for row in run_blast_for_length(queries, refs, args):
            candidate = row["candidate_peptide"]
            current = best.get(candidate)
            if current is None or better_hit(row, current):
                best[candidate] = row
    return best


def run_blast_for_length(
    candidates: list[str],
    refs: list[dict[str, str]],
    args: argparse.Namespace,
) -> list[dict[str, str]]:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        db_fasta = tmp / "tumor.fa"
        query_fasta = tmp / "candidates.fa"
        db_fasta.write_text(
            "".join(f">ref_{idx}\n{row['peptide']}\n" for idx, row in enumerate(refs))
        )
        query_fasta.write_text(
            "".join(f">{query_id(seq)}\n{seq}\n" for seq in candidates)
        )
        db_prefix = tmp / "tumor"
        make = subprocess.run(
            [args.makeblastdb_cmd, "-in", str(db_fasta), "-dbtype", "prot", "-out", str(db_prefix)],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
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
            cmd.extend(shlex.split(args.extra_args))
        blast = subprocess.run(
            cmd,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        if blast.returncode != 0:
            raise SystemExit(f"blastp failed: {blast.stderr}")

        qmap = {query_id(seq): seq for seq in candidates}
        rows = []
        for line in blast.stdout.splitlines():
            fields = line.split("\t")
            if len(fields) != 8:
                continue
            qseqid, sseqid, identity, align_len, _mismatch, _gapopen, evalue, bitscore = fields
            candidate = qmap.get(qseqid, "")
            try:
                ref_idx = int(sseqid.removeprefix("ref_"))
            except ValueError:
                continue
            if not candidate or ref_idx < 0 or ref_idx >= len(refs):
                continue
            # Reference groups are already length-matched, but guard the parsed alignment too.
            ref = refs[ref_idx]
            if len(candidate) != len(ref["peptide"]):
                continue
            rows.append(
                {
                    "candidate_peptide": candidate,
                    "peptide_length": str(len(candidate)),
                    "PCI_DB_enabled": "true",
                    "healthy_exact_match": "false",
                    "tumor_immune_peptide": ref["peptide"],
                    "best_match_peptide": ref["peptide"],
                    "blast_identity": identity,
                    "blast_evalue": evalue,
                    "blast_bit_score": bitscore,
                    "reference_hla_allele": ref.get("hla", ""),
                    "reference_disease_tissue": ref.get("disease", ""),
                }
            )
        return rows


def query_id(sequence: str) -> str:
    return "pep_" + hashlib.md5(sequence.encode()).hexdigest()


def better_hit(candidate: dict[str, str], current: dict[str, str]) -> bool:
    cand_identity_value = safe_float(candidate.get("blast_identity"), -1.0)
    curr_identity_value = safe_float(current.get("blast_identity"), -1.0)
    cand_identity = -1.0 if cand_identity_value is None else cand_identity_value
    curr_identity = -1.0 if curr_identity_value is None else curr_identity_value
    if cand_identity != curr_identity:
        return cand_identity > curr_identity
    cand_evalue = safe_float(candidate.get("blast_evalue"), float("inf"))
    curr_evalue = safe_float(current.get("blast_evalue"), float("inf"))
    if cand_evalue != curr_evalue:
        return cand_evalue < curr_evalue
    cand_bits_value = safe_float(candidate.get("blast_bit_score"), -1.0)
    curr_bits_value = safe_float(current.get("blast_bit_score"), -1.0)
    cand_bits = -1.0 if cand_bits_value is None else cand_bits_value
    curr_bits = -1.0 if curr_bits_value is None else curr_bits_value
    return cand_bits > curr_bits


def first(row: dict[str, str], columns: list[str]) -> str:
    for column in columns:
        if row.get(column):
            return str(row[column])
    return ""


if __name__ == "__main__":
    main()
