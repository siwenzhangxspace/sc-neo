#!/usr/bin/env python3
"""Filter healthy exact peptide matches and collect tumor immunopeptidome evidence."""
from __future__ import annotations
import argparse
import hashlib
import math
import shlex
import subprocess
import tempfile
from pathlib import Path
from scneo_utils import open_text, parse_lengths, read_tsv, write_tsv

EVIDENCE_FIELDS = [
    "candidate_peptide", "peptide_length", "PCI_DB_enabled", "healthy_blast_match",
    "healthy_match_peptide", "healthy_blast_identity", "healthy_blast_evalue",
    "healthy_blast_alignment_length", "tumor_immune_peptide", "best_match_peptide",
    "blast_identity", "blast_evalue", "blast_bit_score", "blast_alignment_length",
    "qstart", "qend", "qlen", "sstart", "send", "slen", "query_coverage",
    "tumor_full_length_exact", "pci_status", "search_completed",
    "reference_id", "reference_hla_allele", "reference_disease_tissue",
]


def load_reference(path):
    with open_text(path) as handle:
        first_line = handle.readline()
        if first_line.startswith(">"):
            rows = []
            identifier, sequence = first_line[1:].strip().split()[0], []
            for line in handle:
                if line.startswith(">"):
                    rows.append({"peptide": "".join(sequence).upper(), "reference_id": identifier})
                    identifier, sequence = line[1:].strip().split()[0], []
                else:
                    sequence.append(line.strip())
            rows.append({"peptide": "".join(sequence).upper(), "reference_id": identifier})
            return [row for row in rows if row["peptide"]]
    rows = []
    for i, row in enumerate(read_tsv(path)):
        sequence = first(row, ["Peptide", "peptide", "Sequence", "sequence", "peptide_sequence"]).upper()
        if sequence:
            rows.append({
                "peptide": sequence,
                "reference_id": first(row, ["reference_id", "id", "ID"]) or f"ref_{i:012d}",
                "hla": first(row, ["HLA", "HLA allele", "MHC", "Allele", "hla_allele"]),
                "disease": first(row, ["Disease", "disease", "Tissue", "tissue", "Source", "source"]),
            })
    return rows


def blast_best_hits(candidates, refs, args):
    if not candidates:
        return {}
    if not refs:
        raise ValueError("Tumor immunopeptidome reference is empty")
    best = {}
    with tempfile.TemporaryDirectory() as directory:
        tmp = Path(directory)
        (tmp / "reference.fa").write_text("".join(f">ref_{i}\n{row['peptide']}\n" for i, row in enumerate(refs)))
        (tmp / "candidates.fa").write_text("".join(f">{query_id(seq)}\n{seq}\n" for seq in candidates))
        database = subprocess.run([args.makeblastdb_cmd, "-in", "reference.fa", "-dbtype", "prot", "-out", "reference"],
                                  cwd=tmp, capture_output=True, text=True)
        if database.returncode:
            raise RuntimeError("makeblastdb failed: " + (database.stderr or database.stdout).strip())
        command = [args.blastp_cmd, "-task", "blastp-short", "-query", "candidates.fa",
                   "-db", "reference", "-outfmt",
                   "6 qseqid sseqid pident length mismatch gapopen evalue bitscore qstart qend sstart send qlen slen"]
        command += shlex.split(args.extra_args)
        result = subprocess.run(command, cwd=tmp, capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError("blastp failed: " + (result.stderr or result.stdout).strip())
    queries = {query_id(seq): seq for seq in candidates}
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) != 14:
            raise ValueError("Expected 14 BLAST columns including query coordinates")
        query, subject, identity, length, mismatch, gaps, evalue, bits, qs, qe, ss, se, qlen, slen = fields
        if query not in queries or not subject.startswith("ref_"):
            raise ValueError("BLAST output contains an unknown identifier")
        index = int(subject[4:])
        if not 0 <= index < len(refs):
            raise ValueError("BLAST reference index out of range")
        ref, sequence = refs[index], queries[query]
        numbers = [float(v) for v in [identity, evalue, bits]]
        if not all(math.isfinite(v) for v in numbers) or not 0 <= numbers[0] <= 100 or min(numbers[1:]) < 0:
            raise ValueError("Invalid BLAST score")
        if not 1 <= int(qs) <= int(qe) <= int(qlen) == len(sequence):
            raise ValueError("Invalid BLAST query coordinates")
        if not 1 <= int(ss) <= int(se) <= int(slen) == len(ref["peptide"]):
            raise ValueError("Invalid BLAST reference coordinates")
        row = {
            "candidate_peptide": sequence, "tumor_immune_peptide": ref["peptide"],
            "best_match_peptide": ref["peptide"], "blast_identity": identity, "blast_evalue": evalue,
            "blast_bit_score": bits, "blast_alignment_length": length,
            "qstart": qs, "qend": qe, "qlen": qlen, "sstart": ss, "send": se, "slen": slen,
            "query_coverage": format((int(qe) - int(qs) + 1) / len(sequence), ".17g"),
            "reference_id": ref["reference_id"], "reference_hla_allele": ref.get("hla", ""),
            "reference_disease_tissue": ref.get("disease", ""),
        }
        if sequence not in best or better_hit(row, best[sequence]):
            best[sequence] = row
    return best


def better_hit(candidate, current):
    def order(row):
        return (float(row["blast_evalue"]), -float(row["blast_bit_score"]),
                -float(row["blast_identity"]), -int(row["blast_alignment_length"]),
                row["reference_id"])
    return order(candidate) < order(current)


def query_id(sequence):
    return "pep_" + hashlib.md5(sequence.encode()).hexdigest()


def first(row, columns):
    return next((str(row[column]).strip() for column in columns if row.get(column)), "")


def peptide(row):
    return first(row, ["peptide", "peptide_seq", "peptide_sequence"])


def write_filtered(path, rows, excluded):
    fields = list(rows[0]) if rows else ["patient_id", "source_type", "peptide", "peptide_length"]
    write_tsv(path, [row for row in rows if peptide(row) not in excluded], fields)


def main():
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
    mutations, splicing = read_tsv(args.mutation), read_tsv(args.splicing)
    lengths = set(parse_lengths(args.lengths))
    candidates = sorted({peptide(row) for row in mutations + splicing if len(peptide(row)) in lengths})
    healthy, tumor, hits, excluded = set(), set(), {}, set()
    if not args.disabled:
        if not args.healthy or not args.tumor:
            parser.error("--healthy and --tumor are required unless --disabled")
        healthy_refs, tumor_refs = load_reference(args.healthy), load_reference(args.tumor)
        if not healthy_refs or not tumor_refs:
            parser.error("Enabled immunopeptidome references must contain peptide sequences")
        healthy = {row["peptide"] for row in healthy_refs}
        tumor = {row["peptide"] for row in tumor_refs}
        excluded = set(candidates) & healthy
        hits = blast_best_hits([seq for seq in candidates if seq not in excluded], tumor_refs, args)
    evidence = []
    for sequence in candidates:
        row = {field: "" for field in EVIDENCE_FIELDS}
        row.update(candidate_peptide=sequence, peptide_length=len(sequence),
                   PCI_DB_enabled=str(not args.disabled).lower(),
                   healthy_blast_match=str(sequence in excluded).lower(),
                   tumor_full_length_exact=str(sequence in tumor).lower(),
                   search_completed=str(not args.disabled and sequence not in excluded).lower(),
                   pci_status="not_measured" if args.disabled else "no_hit")
        if sequence in excluded:
            row.update(healthy_match_peptide=sequence, healthy_blast_identity=100,
                       healthy_blast_alignment_length=len(sequence), pci_status="excluded_healthy")
        if sequence in hits:
            row.update(hits[sequence], pci_status="hit")
        evidence.append(row)
    write_filtered(args.out_mutation, mutations, excluded)
    write_filtered(args.out_splicing, splicing, excluded)
    write_tsv(args.out_evidence, evidence, EVIDENCE_FIELDS)


if __name__ == "__main__":
    main()
