#!/usr/bin/env python3
"""Apply population/germline filters to MAF output and preserve trace fields."""

from __future__ import annotations

import argparse
from collections import defaultdict

from scneo_utils import read_tsv, safe_float, safe_int, variant_key, write_tsv


AF_COLUMNS = [
    "gnomADg_AF",
    "gnomADe_AF",
    "AF",
    "AFR_AF",
    "AMR_AF",
    "EAS_AF",
    "EUR_AF",
    "SAS_AF",
    "TOPMED",
]

CODING_VARIANT_CLASSIFICATIONS = frozenset(
    {
        "Missense_Mutation",
        "Frame_Shift_Del",
        "Frame_Shift_Ins",
        "In_Frame_Del",
        "In_Frame_Ins",
    }
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maf", required=True)
    parser.add_argument("--trace-tsv", required=True)
    parser.add_argument("--population-af-max", type=float, default=0.001)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    maf_rows = read_tsv(args.maf)
    trace_rows = [row for row in read_tsv(args.trace_tsv) if row.get("passes_enrichment") == "true"]
    trace_by_locus = defaultdict(list)
    for row in trace_rows:
        trace_by_locus[locus_key(row.get("CHROM", ""), row.get("POS", ""), row.get("REF", ""), row.get("ALT", ""))].append(row)

    out_rows = []
    for maf in maf_rows:
        chrom = maf.get("Chromosome") or maf.get("CHROM") or maf.get("#CHROM")
        pos = maf.get("Start_Position") or maf.get("POS")
        ref = maf.get("Reference_Allele") or maf.get("REF")
        alt = maf.get("Tumor_Seq_Allele2") or maf.get("Tumor_Seq_Allele1") or maf.get("ALT")
        traces = trace_by_locus.get(locus_key(chrom, pos, ref, alt), [])
        if not traces:
            continue
        if is_pre_filtered(maf, args.population_af_max):
            continue
        variant_class = classify_variant(maf)
        if variant_class is None:
            continue
        if not is_coding_variant(maf):
            continue
        for trace in traces:
            row = dict(maf)
            cosmic_support = cosmic_count(maf)
            row.update(
                {
                    "patient_id": trace.get("patient_id", ""),
                    "variant_key": trace.get("variant_key", variant_key(chrom, pos, ref, alt, trace.get("varType", ""))),
                    "varType": trace.get("varType", ""),
                    "supporting_tumor_cell_ids": trace.get("supporting_tumor_cell_ids", ""),
                    "tumor_callable_cell_ids": trace.get("tumor_callable_cell_ids", ""),
                    "normal_alt_cell_ids": trace.get("normal_alt_cell_ids", ""),
                    "normal_callable_cell_ids": trace.get("normal_callable_cell_ids", ""),
                    "normal_evidence_status": trace.get("normal_evidence_status", ""),
                    "tumor_specific_candidate": trace.get("tumor_specific_candidate", ""),
                    "tumor_enriched_candidate": trace.get("tumor_enriched_candidate", ""),
                    "somatic_candidate": "true",
                    "somatic_mutation_candidate": "true",
                    "known_cancer_mutation": "true" if cosmic_support >= 5 else "false",
                    "high_confidence_somatic_candidate": "true" if cosmic_support >= 5 else "false",
                    "novel_variant": "true" if variant_class == "novel_variant" else "false",
                    "variant_class": variant_class,
                }
            )
            out_rows.append(row)

    fields = list(maf_rows[0].keys()) if maf_rows else ["Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"]
    for field in [
        "patient_id",
        "variant_key",
        "varType",
        "supporting_tumor_cell_ids",
        "tumor_callable_cell_ids",
        "normal_alt_cell_ids",
        "normal_callable_cell_ids",
        "normal_evidence_status",
        "tumor_specific_candidate",
        "tumor_enriched_candidate",
        "somatic_candidate",
        "somatic_mutation_candidate",
        "known_cancer_mutation",
        "high_confidence_somatic_candidate",
        "novel_variant",
        "variant_class",
    ]:
        if field not in fields:
            fields.append(field)
    write_tsv(args.out, out_rows, fields)


def locus_key(chrom: str, pos: str, ref: str, alt: str) -> str:
    return "|".join([str(chrom).replace("chr", ""), str(pos), str(ref), str(alt)])


def exceeds_population_af(row: dict[str, str], max_af: float) -> bool:
    for col in AF_COLUMNS:
        for value in str(row.get(col, "")).replace("&", ",").split(","):
            af = safe_float(value, None)
            if af is not None and af > max_af:
                return True
    return False


def is_pre_filtered(row: dict[str, str], population_af_max: float) -> bool:
    """Apply the reference script's two unconditional exclusion rules."""
    return exceeds_population_af(row, population_af_max) or safe_float(row.get("dbSNP_SAO"), None) == 1


def is_coding_variant(row: dict[str, str]) -> bool:
    """Keep only the explicit protein-altering MAF classifications requested."""
    return row.get("Variant_Classification", "") in CODING_VARIANT_CLASSIFICATIONS


def cosmic_count(row: dict[str, str]) -> int:
    targeted = safe_int(row.get("COSMIC_TARGETED_SCREEN_SAMPLE_COUNT"), 0)
    genome = safe_int(row.get("COSMIC_GENOME_SCREEN_SAMPLE_COUNT"), 0)
    if "COSMIC_TARGETED_SCREEN_SAMPLE_COUNT" in row or "COSMIC_GENOME_SCREEN_SAMPLE_COUNT" in row:
        return targeted + genome
    for col in ["COSMIC_CNT", "COSMIC_COUNT", "CNT"]:
        if col in row:
            return safe_int(row.get(col), 0)
    return 0


def is_valid(value: object) -> bool:
    """Return whether an annotation value is present and meaningful."""
    if value is None:
        return False
    return str(value).strip().lower() not in {"", ".", "nan", "none"}


def classify_variant(row: dict[str, str]) -> str | None:
    """Apply the final variant classification decision tree."""
    cosmic_support = cosmic_count(row)
    sao = safe_float(row.get("dbSNP_SAO"), None)
    in_cosmic = is_valid(row.get("COSMIC"))
    in_dbsnp = is_valid(row.get("dbSNP_RS"))

    if cosmic_support >= 5:
        return "tumor_variant"
    if sao == 2:
        return "tumor_variant"
    if sao == 3:
        return None
    if not in_cosmic and in_dbsnp and sao is None:
        return "uncertain"
    if in_dbsnp and sao in (None, 0):
        return "uncertain"
    if not in_cosmic and not in_dbsnp:
        return "novel_variant"
    return "uncertain"


if __name__ == "__main__":
    main()
