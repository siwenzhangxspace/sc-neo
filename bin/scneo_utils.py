#!/usr/bin/env python3
"""Shared utilities for the Smart-seq2 neoantigen workflow."""

from __future__ import annotations

import csv
import gzip
import math
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence


LIST_SEP = ","


def _tsv_data_lines(handle: Iterable[str]) -> Iterator[str]:
    skipping_header_comments = True
    for line in handle:
        if skipping_header_comments and line.startswith("#"):
            continue
        skipping_header_comments = False
        yield line


def die(message: str, exit_code: int = 2) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(exit_code)


def open_text(path: os.PathLike | str, mode: str = "rt"):
    path = str(path)
    if path == "-":
        return sys.stdin if "r" in mode else sys.stdout
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, newline="" if "t" in mode else None)


def read_tsv(path: os.PathLike | str) -> List[Dict[str, str]]:
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(_tsv_data_lines(handle), delimiter="\t")
        if reader.fieldnames is None:
            return []
        return [dict(row) for row in reader]


def iter_tsv(path: os.PathLike | str) -> Iterator[Dict[str, str]]:
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(_tsv_data_lines(handle), delimiter="\t")
        if reader.fieldnames is None:
            return
        for row in reader:
            yield dict(row)


def write_tsv(path: os.PathLike | str, rows: Iterable[Mapping[str, object]], fieldnames: Sequence[str]) -> None:
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(fieldnames), lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: stringify(row.get(field, "")) for field in fieldnames})


def stringify(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, (list, tuple, set)):
        return LIST_SEP.join(str(v) for v in value if str(v) != "")
    return str(value)


def split_list(value: str | None) -> List[str]:
    if not value:
        return []
    parts = re.split(r"[,;|]", value)
    return [p.strip() for p in parts if p.strip()]


def uniq(values: Iterable[str]) -> List[str]:
    seen = set()
    out = []
    for value in values:
        value = str(value)
        if value and value not in seen:
            seen.add(value)
            out.append(value)
    return out


def safe_float(value: object, default: Optional[float] = None) -> Optional[float]:
    try:
        if value is None or str(value).strip() == "":
            return default
        return float(str(value).replace(",", ""))
    except ValueError:
        return default


def safe_int(value: object, default: int = 0) -> int:
    try:
        if value is None or str(value).strip() == "":
            return default
        return int(float(str(value).replace(",", "")))
    except ValueError:
        return default


def read_manifest(path: os.PathLike | str) -> List[Dict[str, str]]:
    required = ["patient_id", "cell_id", "cell_type", "fastq_1", "fastq_2"]
    rows = read_tsv(path)
    if not rows:
        die(f"manifest is empty: {path}")
    missing = [field for field in required if field not in rows[0]]
    if missing:
        die(f"manifest is missing required columns: {','.join(missing)}")
    return rows


def group_by(rows: Iterable[Mapping[str, str]], key: str) -> Dict[str, List[Mapping[str, str]]]:
    grouped: Dict[str, List[Mapping[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[str(row.get(key, ""))].append(row)
    return dict(grouped)


def parse_info(info: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not info or info == ".":
        return out
    for item in info.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = value
        else:
            out[item] = "true"
    return out


def format_info(info: Mapping[str, object]) -> str:
    parts = []
    for key, value in info.items():
        if value is None or value == "":
            continue
        if value is True:
            parts.append(str(key))
        else:
            parts.append(f"{key}={value}")
    return ";".join(parts) if parts else "."


def variant_key(chrom: str, pos: str | int, ref: str, alt: str, vartype: str = "") -> str:
    return "|".join([str(chrom), str(pos), str(ref), str(alt), str(vartype or "")])


def iter_vcf_records(path: os.PathLike | str) -> Iterator[Dict[str, str]]:
    samples: List[str] = []
    with open_text(path, "rt") as handle:
        for line in handle:
            if not line or line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").lstrip("#").split("\t")
                samples = header[9:]
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            record = {
                "CHROM": fields[0],
                "POS": fields[1],
                "ID": fields[2] if len(fields) > 2 else ".",
                "REF": fields[3],
                "ALT": fields[4],
                "QUAL": fields[5] if len(fields) > 5 else ".",
                "FILTER": fields[6] if len(fields) > 6 else ".",
                "INFO": fields[7],
                "FORMAT": fields[8] if len(fields) > 8 else "",
            }
            for sample, sample_value in zip(samples, fields[9:]):
                record[f"SAMPLE:{sample}"] = sample_value
            yield record


def write_vcf(path: os.PathLike | str, records: Iterable[Mapping[str, str]], source: str = "scneo-nextflow") -> None:
    with open_text(path, "wt") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write(f"##source={source}\n")
        handle.write("##INFO=<ID=varType,Number=1,Type=String,Description=\"Variant type\">\n")
        handle.write("##INFO=<ID=patient_ids,Number=.,Type=String,Description=\"Supporting patient ids\">\n")
        handle.write("##INFO=<ID=supporting_tumor_cell_ids,Number=.,Type=String,Description=\"Tumor cells supporting this variant\">\n")
        handle.write("##INFO=<ID=patient_support_counts,Number=.,Type=String,Description=\"Per-patient tumor-cell support counts\">\n")
        handle.write("##INFO=<ID=patient_supporting_tumor_cell_ids,Number=.,Type=String,Description=\"Per-patient tumor cells supporting this variant\">\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for row in records:
            info = row.get("INFO")
            if not info:
                info = format_info(
                    {
                        "varType": row.get("varType", ""),
                        "patient_ids": row.get("patient_ids", ""),
                        "supporting_tumor_cell_ids": row.get("supporting_tumor_cell_ids", ""),
                        "patient_support_counts": row.get("patient_support_counts", ""),
                        "patient_supporting_tumor_cell_ids": row.get("patient_supporting_tumor_cell_ids", ""),
                    }
                )
            handle.write(
                "\t".join(
                    [
                        str(row.get("CHROM", "")),
                        str(row.get("POS", "")),
                        str(row.get("ID", ".")) or ".",
                        str(row.get("REF", "")),
                        str(row.get("ALT", "")),
                        str(row.get("QUAL", ".")) or ".",
                        str(row.get("FILTER", "PASS")) or "PASS",
                        str(info),
                    ]
                )
                + "\n"
            )


def fisher_greater(a: int, b: int, c: int, d: int) -> float:
    """One-sided Fisher exact test for enrichment in row 1."""
    a, b, c, d = [max(0, int(x)) for x in (a, b, c, d)]
    n = a + b + c + d
    row1 = a + b
    col1 = a + c
    if n == 0 or row1 == 0 or col1 == 0:
        return 1.0
    min_x = max(0, row1 - (n - col1))
    max_x = min(row1, col1)

    def hypergeom(x: int) -> float:
        return math.comb(col1, x) * math.comb(n - col1, row1 - x) / math.comb(n, row1)

    return min(1.0, sum(hypergeom(x) for x in range(max(a, min_x), max_x + 1)))


def bh_adjust(p_values: Sequence[float]) -> List[float]:
    n = len(p_values)
    if n == 0:
        return []
    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    adjusted = [1.0] * n
    prev = 1.0
    for rank, (idx, p_value) in reversed(list(enumerate(indexed, start=1))):
        q_value = min(prev, float(p_value) * n / rank)
        adjusted[idx] = min(1.0, q_value)
        prev = q_value
    return adjusted


def run_command(cmd: Sequence[str], cwd: str | None = None, check: bool = True) -> subprocess.CompletedProcess:
    print("+ " + " ".join(str(x) for x in cmd), file=sys.stderr)
    return subprocess.run(cmd, cwd=cwd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=check)


def infer_cell_id_from_path(path: os.PathLike | str) -> str:
    name = Path(path).name
    for suffix in [
        ".sorted.bam",
        ".filtered.vcf.gz",
        ".filtered.vcf",
        ".vcf.gz",
        ".vcf",
        ".SJ.out.tab",
        ".report.tsv",
        ".tsv",
        ".txt",
    ]:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(path).stem


def require_columns(rows: Sequence[Mapping[str, str]], columns: Sequence[str], label: str) -> None:
    if not rows:
        die(f"{label} has no rows")
    missing = [col for col in columns if col not in rows[0]]
    if missing:
        die(f"{label} missing columns: {','.join(missing)}")


def parse_lengths(value: str) -> List[int]:
    lengths = [safe_int(item, 0) for item in split_list(value)]
    lengths = [length for length in lengths if length > 0]
    return lengths or [8, 9, 10, 11]
