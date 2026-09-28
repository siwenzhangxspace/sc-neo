#!/usr/bin/env python3
"""Merge per-cell FASTQ files into one patient pool, keeping the first N reads of each file.

HLA-I typing only needs reads that map to the HLA reference; the read order within a
FASTQ file is unrelated to HLA content, so taking the leading reads of every file is a
safe and reproducible downsample strategy (no randomness, no seed). Each input file
contributes at most ceil(target_reads / n_files) records, so every cell is represented.

Paired-end mode (`--r2` given): R1 and R2 files are truncated with the same per-file
limit and concatenated in the same file order, so read pairs stay aligned between the
two output pools. Gzipped inputs are decompressed on the fly and truncated reads are
dropped with a warning on stderr.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path
from typing import IO, Iterator, List, Optional, TextIO

from scneo_utils import die


def open_maybe_gzip(path: Path) -> IO[str]:
    """Open a FASTQ file transparently, decompressing .gz inputs."""
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def read_records(handle: TextIO, limit: int) -> Iterator[Optional[List[str]]]:
    """Yield up to `limit` FASTQ records (4 lines each) from an open text handle.

    A trailing truncated record (file length not a multiple of 4) is reported by
    yielding None and is not emitted.
    """
    taken = 0
    while taken < limit:
        lines = [handle.readline() for _ in range(4)]
        if lines[0] == "":
            return
        if any(line == "" for line in lines[1:]):
            yield None
            return
        yield lines
        taken += 1


def _truncate_one(path: Path, limit: int, out: TextIO, truncated: List[str]) -> int:
    """Write up to `limit` records of one file to `out`; return records written."""
    written = 0
    with open_maybe_gzip(path) as handle:
        for record in read_records(handle, limit):
            if record is None:
                truncated.append(path.name)
                break
            out.write("".join(record))
            written += 1
    print(f"{path.name}: {written} reads", file=sys.stderr)
    return written


def merge_prefix(
    r1_files: List[Path],
    target_reads: int,
    out_r1: Path,
    out_r2: Optional[Path] = None,
    r2_files: Optional[List[Path]] = None,
) -> None:
    """Merge leading reads of every file into plain FASTQ pool(s)."""
    if not r1_files:
        die("no input files given")
    if target_reads <= 0:
        die(f"--target-reads must be positive, got {target_reads}")
    paired = out_r2 is not None
    if paired:
        if r2_files is None or len(r2_files) != len(r1_files):
            die(f"paired mode needs one --r2 file per --r1 file (got {len(r1_files)} R1, {len(r2_files or [])} R2)")
    n_files = len(r1_files)
    per_file = (target_reads + n_files - 1) // n_files
    truncated: List[str] = []
    total_written = 0
    with open(out_r1, "w", encoding="utf-8") as out1:
        for path in r1_files:
            total_written += _truncate_one(path, per_file, out1, truncated)
        if paired:
            with open(out_r2, "w", encoding="utf-8") as out2:
                for path in r2_files:
                    total_written += _truncate_one(path, per_file, out2, truncated)
    print(f"total: {total_written} reads ({per_file} per file target)", file=sys.stderr)
    if truncated:
        print(
            f"warning: dropped truncated final records in: {', '.join(truncated)}",
            file=sys.stderr,
        )


def _flatten(groups: List[List[str]]) -> List[Path]:
    return [Path(p) for group in groups for p in group]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r1", action="append", nargs="+", required=True, help="R1 FASTQ files (plain or .gz); may be repeated")
    parser.add_argument("--r2", action="append", nargs="+", default=[], help="R2 FASTQ files for paired-end mode; one per --r1 file")
    parser.add_argument("--target-reads", type=int, required=True, help="target total reads in the pool")
    parser.add_argument("--out-r1", required=True, help="output R1 plain FASTQ pool")
    parser.add_argument("--out-r2", help="output R2 plain FASTQ pool (paired-end mode)")
    args = parser.parse_args()

    r1_files = _flatten(args.r1)
    r2_files = _flatten(args.r2)
    missing = [str(f) for f in r1_files + r2_files if not f.is_file()]
    if missing:
        die(f"input files not found: {', '.join(missing)}")
    merge_prefix(
        r1_files,
        args.target_reads,
        Path(args.out_r1),
        out_r2=Path(args.out_r2) if args.out_r2 else None,
        r2_files=r2_files if args.out_r2 else None,
    )


if __name__ == "__main__":
    main()
