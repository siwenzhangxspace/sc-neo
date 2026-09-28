#!/usr/bin/env python3

from __future__ import annotations
import ast
import csv
import re
from bisect import bisect, bisect_left

def load_transcripts(gtf_path: str):
    result = {}
    with open(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            attrs = dict(re.findall(r'(\S+)\s+"([^"]*)"', fields[8]))
            gene = attrs.get("gene_id", "").split(".")[0]
            transcript = attrs.get("transcript_id", "").split(".")[0]
            if gene and transcript:
                result.setdefault(gene, {}).setdefault(transcript, []).append((int(fields[3]), int(fields[4])))
    for transcripts in result.values():
        for exons in transcripts.values():
            exons.sort()
    return result

def load_start_codons(path: str):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "start_codon" not in reader.fieldnames:
            raise ValueError("start-codon table must contain a start_codon column")
        id_column = reader.fieldnames[0]
        result = {}
        for row in reader:
            gene = row[id_column].split(".")[0]
            try:
                values = ast.literal_eval(row["start_codon"])
            except (SyntaxError, ValueError) as exc:
                raise ValueError(f"invalid start_codon value for {gene}") from exc
            result[gene] = [int(value) for value in values]
        return result

def infer_support_phases(gene, boundary, strand, first_length, transcripts, start_codons):
    support = {}
    for start in start_codons.get(gene, []):
        for phase, pssc, transcript, support_strand in get_support_phase(
            transcripts.get(gene, {}), boundary, start, strand, first_length
        ):
            support.setdefault(phase, []).append((pssc, transcript, support_strand))
    return support

def get_support_phase(all_trans, boundary, pssc, strand, first_length):
    boundary, pssc, first_length = int(boundary), int(pssc), int(first_length)
    supports = []
    for transcript, original in all_trans.items():
        exons = sorted(original)
        flat = [coordinate for exon in exons for coordinate in exon]
        if strand == "+":
            start_pos = bisect(flat, pssc)
            first_base = boundary - first_length + 1
            junction_pos = bisect(flat, first_base)
            if start_pos % 2 != 1 or junction_pos % 2 != 1 or junction_pos < start_pos:
                continue
            n_bases = _path_length(exons, start_pos, junction_pos, pssc, first_base, "+")
        else:
            start_pos = bisect_left(flat, pssc)
            first_base = boundary + first_length - 1
            junction_pos = bisect_left(flat, first_base)
            if start_pos % 2 != 1 or junction_pos % 2 != 1 or junction_pos > start_pos:
                continue
            n_bases = _path_length(exons, start_pos, junction_pos, pssc, first_base, "-")
        phase = {1: 0, 2: 2, 0: 1}[n_bases % 3]
        supports.append((phase, pssc, transcript, strand))
    return supports

def _path_length(exons, start_pos, junction_pos, pssc, first_base, strand):
    if start_pos == junction_pos:
        return first_base - pssc + 1 if strand == "+" else pssc - first_base + 1
    start_index, end_index = (start_pos - 1) // 2, (junction_pos - 1) // 2
    low, high = sorted((start_index, end_index))
    total = 0
    for index in range(low, high + 1):
        exon_start, exon_end = exons[index]
        if strand == "+":
            left = pssc if index == start_index else exon_start
            right = first_base if index == end_index else exon_end
        else:
            left = first_base if index == end_index else exon_start
            right = pssc if index == start_index else exon_end
        total += right - left + 1
    return total
