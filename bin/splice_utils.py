#!/usr/bin/env python3
"""Shared AltAnalyze junction helpers."""

from __future__ import annotations

import re
import csv
from collections import defaultdict
from itertools import chain


CODON = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def parse_uid(uid: str) -> tuple[str, str, str] | None:
    parts = str(uid).split(":")
    if len(parts) < 2:
        return None
    ensg = parts[0].split(".")[0]
    exon_part = ":".join(parts[1:])
    if "-" not in exon_part:
        return None
    subexon1, subexon2 = exon_part.rsplit("-", 1)
    if not ensg or not subexon1 or not subexon2:
        return None
    return ensg, subexon1, subexon2


def split_subexon(subexon: str, default_gene: str) -> tuple[str, str, str | None]:
    gene = default_gene
    actual = subexon
    if ":" in actual:
        gene, actual = actual.split(":", 1)
        gene = gene.split(".")[0]
    suffix = None
    if "_" in actual:
        actual, suffix = actual.split("_", 1)
    return gene, actual, suffix


def clean_subexon_for_lookup(subexon: str) -> str:
    _gene, actual, _suffix = split_subexon(subexon, "")
    return actual


def load_exon_coords(path: str) -> dict[str, dict[str, tuple[str, str, int, int]]]:
    coords: dict[str, dict[str, tuple[str, str, int, int]]] = defaultdict(dict)
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        positions = {name: idx for idx, name in enumerate(header)}
        required = ["EnsGID", "SubExonID", "Chromosome", "Strand", "Start", "End"]
        if all(name in positions for name in required):
            column_indexes = [positions[name] for name in required]
        elif len(header) >= 6:
            # Native AltAnalyze/SNAF Hs_Ensembl_exon.txt layout: the first six
            # columns are gene, subexon, chromosome, strand, start and end.
            column_indexes = list(range(6))
        else:
            raise ValueError("exon table must be a six-column normalized table or AltAnalyze Hs_Ensembl_exon.txt")
        for line_number, line in enumerate(handle, start=2):
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(column_indexes):
                raise ValueError(f"exon table line {line_number} has too few columns")
            gene, subexon, chrom, strand, start, end = (fields[index] for index in column_indexes)
            gene = gene.split(".")[0]
            coords[gene][subexon] = (
                chrom,
                strand,
                int(start),
                int(end),
            )
    return dict(coords)


def load_documented_junctions(path: str) -> set[str]:
    """Read headered or SNAF mRNA-ExonIDs transcript paths."""
    documented = set()
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        first_row = next(reader, None)
        if not first_row:
            raise ValueError("transcript table is empty")
        header_names = {value.lower() for value in first_row}
        has_header = bool(header_names & {"ensgid", "gene", "gene_id"})
        if has_header:
            gene_index = next(
                index
                for index, value in enumerate(first_row)
                if value.lower() in {"ensgid", "gene", "gene_id"}
            )
            rows = reader
        else:
            if len(first_row) < 4:
                raise ValueError("headerless transcript table must have at least four columns")
            gene_index = 0
            rows = chain([first_row], reader)
        for fields in rows:
            if len(fields) <= gene_index:
                continue
            gene = fields[gene_index].split(".")[0]
            exon_path = next((value for value in fields if value and "|" in value), "")
            exons = [item for item in exon_path.split("|") if item]
            for left, right in zip(exons, exons[1:]):
                documented.add(f"{gene}:{left}-{right}")
    return documented


def translate_to_stop(sequence: str) -> str:
    sequence = re.sub(r"[^ACGT]", "N", sequence.upper())
    amino_acids = []
    for index in range(0, len(sequence) - 2, 3):
        amino_acid = CODON.get(sequence[index : index + 3], "X")
        if amino_acid == "*":
            break
        amino_acids.append(amino_acid)
    return "".join(amino_acids)


def translate_junction_to_peptides(
    seq1: str,
    seq2: str,
    lengths: tuple[int, ...] | list[int] = (8, 9, 10, 11),
) -> list[dict[str, int | str]]:
    peptides: list[dict[str, int | str]] = []
    for frame in range(3):
        first = seq1[frame:]
        if len(first) < 3:
            continue
        extra = len(first) % 3
        complete_codons = len(first) // 3
        aa_first = translate_to_stop(first)
        if not aa_first or len(aa_first) != complete_codons:
            continue
        if extra == 0:
            second = seq2
        else:
            second = first[-extra:] + seq2
        aa_second = translate_to_stop(second)
        if not aa_second:
            continue
        for length in lengths:
            for from_second in range(min(length, len(aa_second)), 0, -1):
                from_first = length - from_second
                if from_first == 0 and extra == 0:
                    continue
                if from_first > len(aa_first):
                    continue
                peptide = (aa_first[-from_first:] if from_first else "") + aa_second[:from_second]
                if len(peptide) == length:
                    peptides.append(
                        {
                            "peptide": peptide,
                            "peptide_length": length,
                            "n_from_first": from_first,
                            "translation_frame": frame,
                        }
                    )
    unique = []
    seen = set()
    for item in peptides:
        key = (item["peptide"], item["peptide_length"], item["translation_frame"])
        if key not in seen:
            seen.add(key)
            unique.append(item)
    return unique
