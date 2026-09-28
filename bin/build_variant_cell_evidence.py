#!/usr/bin/env python3
"""Measure per-cell BAM evidence for recurrent small variants."""
from __future__ import annotations
import argparse
from scneo_utils import infer_cell_id_from_path, read_manifest, read_tsv, variant_key, write_tsv

FIELDS = ["patient_id", "cell_id", "cell_type", "variant_key", "CHROM", "POS", "REF", "ALT", "depth", "ref_reads", "alt_reads", "vaf", "callability_status"]

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True); parser.add_argument("--variants", required=True)
    parser.add_argument("--bams", nargs="+", required=True); parser.add_argument("--out", required=True)
    parser.add_argument("--min-depth", type=int, default=8); parser.add_argument("--min-alt-reads", type=int, default=3)
    parser.add_argument("--min-vaf", type=float, default=.05); parser.add_argument("--min-base-quality", type=int, default=20)
    parser.add_argument("--min-mapping-quality", type=int, default=20)
    args = parser.parse_args()
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError("BAM callability requires pysam") from exc
    manifest = {row["cell_id"]: row for row in read_manifest(args.manifest)}
    bam_by_cell = {infer_cell_id_from_path(path): path for path in args.bams if str(path).endswith(".bam")}
    rows = []
    for variant in read_tsv(args.variants):
        chrom, pos, ref, alt = variant["CHROM"], int(variant["POS"]), variant["REF"], variant["ALT"]
        key = variant.get("variant_key") or variant_key(chrom, pos, ref, alt, variant.get("varType", ""))
        for cell_id, meta in manifest.items():
            if meta["patient_id"] != variant["patient_id"] or meta["cell_type"] not in {"tumor", "normal"}:
                continue
            if cell_id not in bam_by_cell:
                raise ValueError(f"missing BAM for manifest cell {cell_id}")
            with pysam.AlignmentFile(bam_by_cell[cell_id]) as bam:
                evidence = count_evidence(bam, chrom, pos, ref, alt, args)
            rows.append({"patient_id": variant["patient_id"], "cell_id": cell_id, "cell_type": meta["cell_type"], "variant_key": key, "CHROM": chrom, "POS": pos, "REF": ref, "ALT": alt, **evidence})
    write_tsv(args.out, rows, FIELDS)

def count_evidence(bam, chrom, pos, ref, alt, args):
    depth = ref_reads = alt_reads = 0
    try:
        pileup_iterator = bam.pileup(
            chrom,
            pos - 1,
            pos,
            truncate=True,
            min_base_quality=args.min_base_quality,
            min_mapping_quality=args.min_mapping_quality,
            stepper="all",
        )
        column = next(pileup_iterator, None)
    except ValueError:
        column = None
    if column is not None:
        for read in column.pileups:
            if read.is_refskip or read.is_del or read.query_position is None or read.alignment.mapping_quality < args.min_mapping_quality:
                continue
            if read.alignment.query_qualities[read.query_position] < args.min_base_quality:
                continue
            depth += 1
            base = read.alignment.query_sequence[read.query_position].upper()
            is_alt = (len(ref) == len(alt) == 1 and base == alt.upper()) or (len(alt) > len(ref) and read.indel == len(alt) - len(ref)) or (len(ref) > len(alt) and read.indel == -(len(ref) - len(alt)))
            if is_alt: alt_reads += 1
            else: ref_reads += 1
    vaf = alt_reads / depth if depth else 0.0
    status = "not_callable" if depth < args.min_depth else "alt" if alt_reads >= args.min_alt_reads and vaf >= args.min_vaf else "ref"
    return {"depth": depth, "ref_reads": ref_reads, "alt_reads": alt_reads, "vaf": f"{vaf:.6g}", "callability_status": status}

if __name__ == "__main__": main()
