#!/usr/bin/env python3
"""Generate 8-11mer mutation peptides covering the annotated protein change."""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from contextlib import ExitStack
from pathlib import Path

from reconstruct_frameshift import (RejectedVariant, base_id, is_frameshift, load_models,
                                   reconstruct, transcript_id)

from scneo_utils import parse_lengths, read_tsv, write_tsv


PEPTIDE_LENGTHS = [8, 9, 10, 11]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maf", required=True)
    parser.add_argument("--protein-fasta", help="Protein reference required for non-frameshift peptide generation")
    parser.add_argument("--lengths", default="8,9,10,11")
    parser.add_argument("--out", required=True)
    parser.add_argument("--genome-fasta")
    parser.add_argument("--reference-gtf")
    parser.add_argument("--rejected-out", required=True)
    args = parser.parse_args()
    if bool(args.genome_fasta) != bool(args.reference_gtf):
        parser.error("--genome-fasta and --reference-gtf must be supplied together")
    with ExitStack() as stack:
        generate(args, stack)


def generate(args, stack):

    proteins = None
    lengths = parse_lengths(args.lengths)
    rows = []
    rejected = []
    successes = 0
    maf_rows = read_tsv(args.maf)
    fasta = None
    models = {}
    if args.genome_fasta:
        import pysam
        if not Path(args.genome_fasta + ".fai").is_file():
            raise FileNotFoundError(args.genome_fasta + ".fai")
        fasta = stack.enter_context(pysam.FastaFile(args.genome_fasta))
        models = load_models(args.reference_gtf, {base_id(transcript_id(row)) for row in maf_rows if is_frameshift(row)})
    for maf in maf_rows:
        consequence = maf.get("Consequence") or maf.get("Variant_Classification") or ""
        fs = is_frameshift(maf)
        if not fs and not coding_consequence(consequence):
            continue
        parsed = parse_hgvsp(maf.get("HGVSp_Short") or maf.get("HGVSp", ""))
        if not fs and (parsed is None or parsed.get("parse_error")):
            continue
        gene = str(maf.get("Hugo_Symbol", "")).strip()
        termination, resolved_tid = "", transcript_id(maf)
        if fs:
            try:
                if fasta is None:
                    raise RejectedVariant("missing_nucleotide_reference", "Supply --genome-fasta and --reference-gtf")
                mutant, start, end, termination, resolved_tid = reconstruct(maf, models, fasta)
            except RejectedVariant as exc:
                rejected.append({"patient_id": maf.get("patient_id", ""), "variant_key": maf.get("variant_key", ""),
                                 "transcript_id": transcript_id(maf), "HGVSc": maf.get("HGVSc", ""),
                                 "HGVSp": maf.get("HGVSp_Short") or maf.get("HGVSp", ""),
                                 "reason": exc.reason, "detail": str(exc)})
                continue
            successes += 1
            mutation_type = "FS"
        else:
            if proteins is None:
                if not args.protein_fasta:
                    raise ValueError("--protein-fasta is required for non-frameshift peptide generation")
                proteins = load_protein_sequences(args.protein_fasta)
            wild_type = proteins.get(gene, "") or find_protein_sequence(proteins, transcript(maf), maf)
            built = build_mutant_sequence(wild_type, parsed)
            if built is None:
                continue
            mutant, start, end = built
            mutation_type = parsed["type"]
        peptides = sliding_window_peptides(mutant, start, end, lengths) if start <= end else []
        for peptide, length in peptides:
            # Preserve annotation and cell evidence until patient-level key7 aggregation.
            output = make_output_row(maf, consequence, peptide, length, start, mutation_type)
            output["frameshift_translation_end"] = termination
            if fs:
                output["transcript_id"] = resolved_tid
            rows.append(output)

    write_tsv(args.rejected_out, rejected, ["patient_id", "variant_key", "transcript_id", "HGVSc", "HGVSp", "reason", "detail"])
    print(f"Frameshift reconstructed={successes} rejected={len(rejected)} reasons={dict(Counter(row['reason'] for row in rejected))}", file=sys.stderr)

    fields = [
        "patient_id",
        "source_type",
        "variant_key",
        "variant_or_junction_id",
        "HGVSp_Short",
        "Variant_Classification",
        "variant_status",
        "gene",
        "transcript_id",
        "consequence",
        "peptide",
        "peptide_seq",
        "peptide_length",
        "mutation_aa_position",
        "mutation_mode",
        "mutation_type",
        "variant_class",
        "supporting_tumor_cell_ids",
        "source_cell_ids",
        "normal_alt_cell_ids",
        "normal_callable_cell_ids",
        "normal_evidence_status",
        "tumor_specific_candidate",
        "tumor_enriched_candidate",
        "annotation_summary",
        "frameshift_translation_end",
    ]
    write_tsv(args.out, rows, fields)


def load_protein_sequences(path: str) -> dict[str, str]:
    """Load either the reference TSV (gene/alias/sequence) or FASTA format."""
    with open(path) as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    if not lines:
        return {}
    if lines[0].startswith(">"):
        return load_fasta_lines(lines)

    sequences = {}
    for line in lines:
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        gene, alias, sequence = (part.strip() for part in parts[:3])
        if gene:
            sequences.setdefault(gene, sequence)
        if alias:
            sequences.setdefault(alias, sequence)
    return sequences


def load_fasta(path: str) -> dict[str, str]:
    """Backward-compatible name for loading the configured protein reference."""
    return load_protein_sequences(path)


def load_fasta_lines(lines: list[str]) -> dict[str, str]:
    sequences = {}
    name = None
    chunks = []
    for line in lines:
        if line.startswith(">"):
            if name:
                store_sequence(sequences, name, "".join(chunks))
            name = line[1:]
            chunks = []
        else:
            chunks.append(line.strip())
    if name:
        store_sequence(sequences, name, "".join(chunks))
    return sequences


def store_sequence(sequences: dict[str, str], header: str, sequence: str) -> None:
    primary = header.split()[0]
    for token in [primary] + re.split(r"[|\s]", header):
        token = token.strip()
        if token:
            sequences.setdefault(token, sequence)
            sequences.setdefault(token.split(".")[0], sequence)


def parse_hgvsp(value: str) -> dict[str, object] | None:
    """Parse the HGVSp forms used by the supplied peptide script."""
    text = str(value or "").strip()
    if not text or text == ".":
        return None
    result: dict[str, object] = {"raw": text}
    if "fs" in text:
        result["type"] = "FS"
    elif "delins" in text:
        result["type"] = "DELINS"
    elif "dup" in text:
        result["type"] = "DUP"
    elif "ins" in text:
        result["type"] = "INS"
    elif "del" in text:
        result["type"] = "DEL"
    else:
        result["type"] = "SNV"

    try:
        body = text[2:] if text.startswith("p.") else text
        mutation_type = result["type"]
        if mutation_type == "SNV":
            match = re.fullmatch(r"([A-Z])([0-9]+)([A-Z*])", body)
            if not match:
                raise ValueError("invalid SNV HGVSp")
            result.update(ref1=match.group(1), num1=int(match.group(2)), alt=match.group(3))
        elif mutation_type == "FS":
            prefix, suffix = body.split("fs", 1)
            match = re.fullmatch(r"([A-Z])([0-9]+)([A-Z*])", prefix)
            if not match:
                raise ValueError("invalid frameshift HGVSp")
            result.update(ref1=match.group(1), num1=int(match.group(2)), alt=match.group(3))
            result["count"] = int(suffix[1:]) if suffix.startswith("*") and suffix[1:] != "?" else None
        elif mutation_type in {"DEL", "DUP"}:
            marker = "del" if mutation_type == "DEL" else "dup"
            segment = body.split(marker, 1)[0]
            first, second = parse_segment(segment)
            result.update(ref1=first[0], num1=first[1], num2=second[1] if second else None)
            if second:
                result["ref2"] = second[0]
        elif mutation_type in {"INS", "DELINS"}:
            marker = "ins" if mutation_type == "INS" else "delins"
            segment, inserted = body.split(marker, 1)
            first, second = parse_segment(segment)
            result.update(ref1=first[0], num1=first[1], ins=inserted)
            if second:
                result.update(ref2=second[0], num2=second[1])
            else:
                result["num2"] = None
    except (ValueError, IndexError, KeyError):
        result["parse_error"] = True
    return result


def parse_segment(segment: str) -> tuple[tuple[str, int], tuple[str, int] | None]:
    parts = segment.split("_")
    first_match = re.fullmatch(r"([A-Z])([0-9]+)", parts[0])
    if not first_match:
        raise ValueError("invalid HGVSp segment")
    first = (first_match.group(1), int(first_match.group(2)))
    second = None
    if len(parts) > 1:
        second_match = re.fullmatch(r"([A-Z])([0-9]+)", parts[1])
        if not second_match:
            raise ValueError("invalid HGVSp range")
        second = (second_match.group(1), int(second_match.group(2)))
    return first, second


def build_mutant_sequence(
    wild_type: str, parsed: dict[str, object] | None
) -> tuple[str, int, int] | None:
    """Apply one parsed protein change and return 1-indexed affected bounds."""
    if not wild_type or parsed is None:
        return None
    mutation_type = parsed["type"]
    length = len(wild_type)
    try:
        n1 = int(parsed["num1"])
        if n1 > length:
            return None
        if mutation_type == "SNV":
            mutant = wild_type[: n1 - 1] + str(parsed["alt"]) + wild_type[n1:]
            return mutant, n1, n1
        n2 = int(parsed.get("num2") or n1)
        if mutation_type == "DEL":
            return wild_type[: n1 - 1] + wild_type[n2:], n1, n1
        if mutation_type == "INS":
            inserted = str(parsed.get("ins", ""))
            return wild_type[:n1] + inserted + wild_type[n1:], n1 + 1, n1 + len(inserted)
        if mutation_type == "DELINS":
            inserted = str(parsed.get("ins", ""))
            mutant = wild_type[: n1 - 1] + inserted + wild_type[n2:]
            return mutant, n1, n1 + len(inserted) - 1
        if mutation_type == "DUP":
            duplicate = wild_type[n1 - 1 : n2]
            return wild_type[:n2] + duplicate + wild_type[n2:], n2 + 1, n2 + len(duplicate)
    except (IndexError, KeyError, TypeError, ValueError):
        return None
    return None


def sliding_window_peptides(
    mutant: str,
    start: int,
    end: int,
    lengths: list[int] = PEPTIDE_LENGTHS,
    stop_pos: int | None = None,
) -> list[tuple[str, int]]:
    """Generate unique windows containing the affected 1-indexed region."""
    peptides = []
    for length in lengths:
        first = max(1, start - length + 1)
        last = min(len(mutant) - length + 1, end)
        for position in range(first, last + 1):
            peptide = mutant[position - 1 : position - 1 + length]
            if len(peptide) != length or "*" in peptide:
                continue
            if stop_pos is not None and position <= stop_pos <= position + length - 1:
                continue
            if start <= position + length - 1 and position <= end:
                peptides.append((peptide, length))
    unique = []
    seen = set()
    for peptide in peptides:
        if peptide[0] not in seen:
            seen.add(peptide[0])
            unique.append(peptide)
    return unique


def find_protein_sequence(proteins: dict[str, str], transcript_id: str, row: dict[str, str]) -> str:
    for key in [transcript_id, row.get("ENSP", ""), row.get("Protein_ID", "")]:
        key = (key or "").split(".")[0]
        if key and key in proteins:
            return proteins[key]
    for column in ["Mutant_Protein_Sequence", "Protein_Sequence", "protein_sequence"]:
        if row.get(column):
            return row[column]
    return ""


def transcript(row: dict[str, str]) -> str:
    for column in ["Transcript_ID", "Feature", "transcript_id", "ENST"]:
        if row.get(column):
            return row[column].split(".")[0]
    return ""


def coding_consequence(consequence: str) -> bool:
    terms = set(re.split(r"[,;& ]+", consequence))
    return bool(
        terms
        & {
            "missense_variant",
            "inframe_insertion",
            "inframe_deletion",
            "frameshift_variant",
            "stop_lost",
            "protein_altering_variant",
            "Missense_Mutation",
            "In_Frame_Ins",
            "In_Frame_Del",
            "Frame_Shift_Ins",
            "Frame_Shift_Del",
        }
    )


def make_output_row(
    maf: dict[str, str], consequence: str, peptide: str, length: int, position: int, mutation_type: str
) -> dict[str, object]:
    variant_key = maf.get("variant_key", "")
    return {
        "patient_id": maf.get("patient_id", ""),
        "source_type": "mutation",
        "variant_key": variant_key,
        "variant_or_junction_id": variant_key,
        "HGVSp_Short": maf.get("HGVSp_Short", ""),
        "Variant_Classification": maf.get("Variant_Classification", ""),
        "variant_status": maf.get("variant_status", ""),
        "gene": maf.get("Hugo_Symbol") or maf.get("Gene") or maf.get("SYMBOL", ""),
        "transcript_id": transcript(maf),
        "consequence": consequence,
        "peptide": peptide,
        "peptide_seq": peptide,
        "peptide_length": length,
        "mutation_aa_position": position,
        "mutation_mode": mutation_type.lower(),
        "mutation_type": mutation_type,
        "variant_class": maf.get("variant_class", ""),
        "supporting_tumor_cell_ids": maf.get("supporting_tumor_cell_ids", ""),
        "source_cell_ids": maf.get("supporting_tumor_cell_ids", ""),
        "normal_alt_cell_ids": maf.get("normal_alt_cell_ids", ""),
        "normal_callable_cell_ids": maf.get("normal_callable_cell_ids", ""),
        "normal_evidence_status": maf.get("normal_evidence_status", ""),
        "tumor_specific_candidate": maf.get("tumor_specific_candidate", ""),
        "tumor_enriched_candidate": maf.get("tumor_enriched_candidate", ""),
        "annotation_summary": summarize_annotation(maf),
    }


def summarize_annotation(row: dict[str, str]) -> str:
    parts = []
    for column in ["Hugo_Symbol", "Consequence", "HGVSp_Short", "HGVSc", "known_cancer_mutation", "novel_variant"]:
        if row.get(column):
            parts.append(f"{column}={row[column]}")
    return ";".join(parts)


if __name__ == "__main__":
    main()
