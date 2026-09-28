#!/usr/bin/env python3
"""Generate strict junction-spanning MHC-I peptides from AltAnalyze UIDs."""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter

from scneo_utils import parse_lengths, read_tsv, write_tsv
from splice_utils import load_exon_coords, parse_uid, split_subexon, translate_junction_to_peptides
from splice_phase import infer_support_phases, load_start_codons, load_transcripts


SPLICING_VARIANT_CLASS = "tumor_as"


FIELDS = [
    "patient_id", "source_type", "junction_id", "variant_or_junction_id", "gene",
    "transcript_id", "junction_annotation", "splice_site", "translation_frame",
    "n_from_first", "phase_supported", "phase_evidence", "peptide", "peptide_length", "reference_protein_match",
    "gtex_mean", "tumor_gtex_delta", "junction_supporting_tumor_cell_ids",
    "junction_supporting_normal_cell_ids", "supporting_tumor_cell_ids", "source_cell_ids",
    "normal_alt_cell_ids", "normal_callable_cell_ids", "normal_evidence_status",
    "tumor_specific_candidate", "tumor_enriched_candidate", "variant_class", "annotation_summary",
]
REJECT_FIELDS = ["patient_id", "junction_id", "reason"]


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))[::-1]


def query_genome(fasta, chrom: str, start: int, end: int, strand: str) -> str | None:
    if start < 1 or end < start:
        return None
    try:
        sequence = fasta.fetch(chrom, start - 1, end).upper()
    except (KeyError, ValueError, OSError):
        return None
    return reverse_complement(sequence) if strand == "-" else sequence


def gene_chrom_strand(gene: str, coords) -> tuple[str, str] | None:
    entries = coords.get(gene, {})
    if not entries:
        return None
    chrom, strand, _start, _end = next(iter(entries.values()))
    return chrom, strand


def unknown_coordinate(actual: str, suffix: str | None) -> int | None:
    value = suffix
    if value is None and actual.startswith("U"):
        value = actual[1:]
    try:
        return int(value) if value else None
    except ValueError:
        return None


def get_subexon_sequence(subexon: str, default_gene: str, coords, fasta, is_first: bool, consecutive: bool):
    gene, actual, suffix = split_subexon(subexon, default_gene)
    if actual.startswith("U"):
        chrom_strand = gene_chrom_strand(gene, coords)
        coordinate = unknown_coordinate(actual, suffix)
        if chrom_strand is None or coordinate is None:
            return None, "missing_exon"
        chrom, strand = chrom_strand
        flank = 100
        if is_first:
            start, end = ((coordinate - flank + 1, coordinate) if strand == "+" else (coordinate, coordinate + flank - 1))
        else:
            start, end = ((coordinate, coordinate + flank - 1) if strand == "+" else (coordinate - flank + 1, coordinate))
        sequence = query_genome(fasta, chrom, start, end, strand)
        return (sequence, "ok") if sequence is not None else (None, "fasta_fail")
    if actual not in coords.get(gene, {}):
        return None, "missing_exon"
    chrom, strand, start, end = coords[gene][actual]
    if suffix is not None:
        try:
            coordinate = int(suffix)
        except ValueError:
            return None, "missing_exon"
        if is_first:
            start, end = ((start, coordinate) if strand == "+" else (coordinate, end))
        else:
            start, end = ((coordinate, end) if strand == "+" else (start, coordinate))
    elif not is_first and consecutive:
        if strand == "+":
            start += 1
        else:
            end -= 1
    sequence = query_genome(fasta, chrom, start, end, strand)
    return (sequence, "ok") if sequence is not None else (None, "fasta_fail")


def is_consecutive(left: str, right: str) -> bool:
    import re

    left_match = re.match(r"^E(\d+)\.(\d+)$", left)
    right_match = re.match(r"^E(\d+)\.(\d+)$", right)
    return bool(
        left_match
        and right_match
        and left_match.group(1) == right_match.group(1)
        and int(left_match.group(2)) + 1 == int(right_match.group(2))
    )


def retrieve_junction_sequence(uid: str, coords, fasta):
    parsed = parse_uid(uid)
    if parsed is None:
        return None, None, "invalid_uid"
    gene, left, right = parsed
    consecutive = is_consecutive(left, right)
    seq1, status1 = get_subexon_sequence(left, gene, coords, fasta, True, consecutive)
    if seq1 is None:
        return None, None, status1
    seq2, status2 = get_subexon_sequence(right, gene, coords, fasta, False, consecutive)
    if seq2 is None:
        return None, None, status2
    return seq1, seq2, "ok"


def subexon_boundary(subexon: str, default_gene: str, coords, is_first: bool):
    gene, actual, suffix = split_subexon(subexon, default_gene)
    coordinate = unknown_coordinate(actual, suffix) if actual.startswith("U") else None
    if actual.startswith("U"):
        chrom_strand = gene_chrom_strand(gene, coords)
        if chrom_strand is None or coordinate is None:
            return "?", "?", "?"
        return chrom_strand[0], chrom_strand[1], str(coordinate)
    if actual not in coords.get(gene, {}):
        return "?", "?", "?"
    chrom, strand, start, end = coords[gene][actual]
    if suffix is not None:
        return chrom, strand, suffix
    coordinate = end if (is_first and strand == "+") else start
    if not is_first:
        coordinate = start if strand == "+" else end
    return chrom, strand, str(coordinate)


def uid_to_coordinate(uid: str, coords) -> str:
    parsed = parse_uid(uid)
    if parsed is None:
        return "unknown"
    gene, left, right = parsed
    chrom, strand, left_coordinate = subexon_boundary(left, gene, coords, True)
    _chrom2, _strand2, right_coordinate = subexon_boundary(right, gene, coords, False)
    if strand == "-":
        return f"{chrom}:{right_coordinate}-{left_coordinate}(-)"
    return f"{chrom}:{left_coordinate}-{right_coordinate}(+)"


def load_gene_names(gtf_path: str) -> dict[str, str]:
    names = {}
    with open(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = dict(re.findall(r'(\S+)\s+"([^"]*)"', fields[8]))
            gene = attrs.get("gene_id", "").split(".")[0]
            if gene:
                names[gene] = attrs.get("gene_name", gene)
    return names


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--junctions", required=True)
    parser.add_argument("--exon-table", required=True)
    parser.add_argument("--genome-fasta", required=True)
    parser.add_argument("--reference-gtf", required=True)
    parser.add_argument("--start-codon-table", required=True)
    parser.add_argument("--strict-phase", action="store_true")
    parser.add_argument("--lengths", default="8,9,10,11")
    parser.add_argument("--out", required=True)
    parser.add_argument("--rejected-out", required=True)
    args = parser.parse_args()
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError("splicing peptide generation requires pysam") from exc
    coords = load_exon_coords(args.exon_table)
    gene_names = load_gene_names(args.reference_gtf)
    transcripts = load_transcripts(args.reference_gtf)
    start_codons = load_start_codons(args.start_codon_table)
    fasta = pysam.FastaFile(args.genome_fasta)
    failures = Counter()
    output = []
    rejected = []
    seen = set()
    try:
        for junction in read_tsv(args.junctions):
            uid = junction.get("uid") or junction.get("junction_id", "")
            parsed = parse_uid(uid)
            if parsed is None:
                failures["invalid_uid"] += 1
                rejected.append({"patient_id": junction.get("patient_id", ""), "junction_id": uid, "reason": "invalid_uid"})
                continue
            seq1, seq2, status = retrieve_junction_sequence(uid, coords, fasta)
            if status != "ok" or seq1 is None or seq2 is None:
                failures[status] += 1
                rejected.append({"patient_id": junction.get("patient_id", ""), "junction_id": uid, "reason": status})
                continue
            peptides = translate_junction_to_peptides(seq1, seq2, parse_lengths(args.lengths))
            _chrom, strand, boundary = subexon_boundary(parsed[1], parsed[0], coords, True)
            try:
                support_phases = infer_support_phases(
                    parsed[0], int(boundary), strand, len(seq1), transcripts, start_codons
                )
            except (KeyError, TypeError, ValueError):
                support_phases = {}
            if args.strict_phase:
                peptides = [item for item in peptides if item["translation_frame"] in support_phases]
            if not peptides:
                reason = "no_supported_phase" if args.strict_phase else "no_translation"
                failures[reason] += 1
                rejected.append({"patient_id": junction.get("patient_id", ""), "junction_id": uid, "reason": reason})
                continue
            gene_id = parsed[0]
            splice_site = uid_to_coordinate(uid, coords)
            for peptide in peptides:
                evidences = support_phases.get(peptide["translation_frame"], [])
                transcript_ids = sorted({item[1] for item in evidences})
                key = (
                    junction.get("patient_id", ""), uid, peptide["peptide"], peptide["translation_frame"]
                )
                if key in seen:
                    continue
                seen.add(key)
                tumor_cells = junction.get("junction_supporting_tumor_cell_ids", "")
                normal_cells = junction.get("junction_supporting_normal_cell_ids", "")
                output.append(
                    {
                        "patient_id": junction.get("patient_id", ""),
                        "source_type": "splicing",
                        "junction_id": uid,
                        "variant_or_junction_id": uid,
                        "gene": gene_names.get(gene_id, gene_id),
                        "transcript_id": ",".join(transcript_ids),
                        "junction_annotation": "altanalyze_uid",
                        "splice_site": splice_site,
                        "translation_frame": peptide["translation_frame"],
                        "n_from_first": peptide["n_from_first"],
                        "phase_supported": "true" if evidences else "false",
                        "phase_evidence": ";".join(f"{start}:{transcript}:{strand_}" for start, transcript, strand_ in evidences),
                        "peptide": peptide["peptide"],
                        "peptide_length": peptide["peptide_length"],
                        "reference_protein_match": "not_evaluated",
                        "gtex_mean": junction.get("gtex_mean", ""),
                        "tumor_gtex_delta": junction.get("tumor_gtex_delta", ""),
                        "junction_supporting_tumor_cell_ids": tumor_cells,
                        "junction_supporting_normal_cell_ids": normal_cells,
                        "supporting_tumor_cell_ids": tumor_cells,
                        "source_cell_ids": tumor_cells,
                        "normal_alt_cell_ids": normal_cells,
                        "normal_callable_cell_ids": junction.get("normal_callable_cell_ids", ""),
                        "normal_evidence_status": junction.get("normal_evidence_status", ""),
                        "tumor_specific_candidate": junction.get("tumor_specific_junction", ""),
                        "tumor_enriched_candidate": junction.get("tumor_enriched_junction", ""),
                        "variant_class": SPLICING_VARIANT_CLASS,
                        "annotation_summary": (
                            f"junction={uid};splice_site={splice_site};filter=external_gtex"
                        ),
                    }
                )
    finally:
        fasta.close()
    print(f"Generated {len(output)} peptide rows; skipped={dict(failures)}", file=sys.stderr)
    write_tsv(args.out, output, FIELDS)
    write_tsv(args.rejected_out, rejected, REJECT_FIELDS)


if __name__ == "__main__":
    main()
