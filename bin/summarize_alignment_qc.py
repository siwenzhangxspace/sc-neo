#!/usr/bin/env python3
"""Create one cell-level QC row from STAR and optional fastp output."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from scneo_utils import safe_int, write_tsv


def parse_star_log(path: str) -> dict:
    metrics = {}
    if not Path(path).exists():
        return metrics
    for line in Path(path).read_text().splitlines():
        if "|" not in line:
            continue
        key, value = [part.strip() for part in line.split("|", 1)]
        metrics[key] = value
    return metrics


def parse_fastp(path: str | None) -> dict:
    if not path or not Path(path).exists() or Path(path).stat().st_size == 0:
        return {}
    try:
        data = json.loads(Path(path).read_text())
    except json.JSONDecodeError:
        return {}
    before = data.get("summary", {}).get("before_filtering", {})
    after = data.get("summary", {}).get("after_filtering", {})
    before_reads = safe_int(before.get("total_reads"), 0)
    after_reads = safe_int(after.get("total_reads"), 0)
    return {
        "raw_reads": before_reads,
        "trimmed_reads": after_reads,
        "trimming_retention": f"{after_reads / before_reads:.6g}" if before_reads else "",
        "gc_content": after.get("gc_content", before.get("gc_content", "")),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patient-id", required=True)
    parser.add_argument("--cell-id", required=True)
    parser.add_argument("--cell-type", required=True)
    parser.add_argument("--star-log", required=True)
    parser.add_argument("--fastp-json")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    star = parse_star_log(args.star_log)
    fastp = parse_fastp(args.fastp_json)
    uniquely = star.get("Uniquely mapped reads %", "")
    multimapped = star.get("% of reads mapped to multiple loci", "")
    input_reads = star.get("Number of input reads", "")
    row = {
        "patient_id": args.patient_id,
        "cell_id": args.cell_id,
        "cell_type": args.cell_type,
        "total_reads": input_reads or fastp.get("raw_reads", ""),
        "raw_reads": fastp.get("raw_reads", ""),
        "trimmed_reads": fastp.get("trimmed_reads", ""),
        "trimming_retention": fastp.get("trimming_retention", ""),
        "gc_content": fastp.get("gc_content", ""),
        "uniquely_mapped_reads_pct": uniquely,
        "multi_mapped_reads_pct": multimapped,
        "exon_reads_pct": "",
        "intron_reads_pct": "",
        "intergenic_reads_pct": "",
        "detected_genes": "",
        "basic_qc_pass": "true" if uniquely else "unknown",
    }
    fields = [
        "patient_id",
        "cell_id",
        "cell_type",
        "total_reads",
        "raw_reads",
        "trimmed_reads",
        "trimming_retention",
        "gc_content",
        "uniquely_mapped_reads_pct",
        "multi_mapped_reads_pct",
        "exon_reads_pct",
        "intron_reads_pct",
        "intergenic_reads_pct",
        "detected_genes",
        "basic_qc_pass",
    ]
    write_tsv(args.out, [row], fields)


if __name__ == "__main__":
    main()
