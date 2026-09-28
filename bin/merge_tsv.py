#!/usr/bin/env python3
"""Merge TSV files with matching or partially overlapping headers."""

from __future__ import annotations

import argparse
import csv

from scneo_utils import open_text, read_tsv, write_tsv


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True)
    parser.add_argument("inputs", nargs="+")
    args = parser.parse_args()

    rows = []
    fields = []
    for path in args.inputs:
        part = read_tsv(path)
        if not part:
            with open_text(path, "rt") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if reader.fieldnames:
                    for field in reader.fieldnames:
                        if field not in fields:
                            fields.append(field)
            continue
        for field in part[0].keys():
            if field not in fields:
                fields.append(field)
        rows.extend(part)
    if not fields:
        fields = ["status"]
        rows = [{"status": "no_input_rows"}]
    write_tsv(args.out, rows, fields)


if __name__ == "__main__":
    main()
