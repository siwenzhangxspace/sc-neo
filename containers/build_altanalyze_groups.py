#!/usr/bin/env python

from __future__ import print_function

import argparse
import os


JUNCTION_SUFFIX = "__junction.bed"


def load_junction_beds(bed_dir):
    filenames = []
    for filename in sorted(os.listdir(bed_dir)):
        if not filename.endswith(JUNCTION_SUFFIX):
            continue
        sample_name = filename[: -len(JUNCTION_SUFFIX)]
        if not sample_name or sample_name.endswith(".bam"):
            raise ValueError("unexpected AltAnalyze junction BED name: " + filename)
        filenames.append(filename)
    if not filenames:
        raise ValueError("no AltAnalyze junction BED files found in: " + bed_dir)
    if len(filenames) < 2:
        raise ValueError("AltAnalyze requires at least two junction BED files")
    return filenames


def build_groups(bed_dir, group_path, comparison_path):
    beds = load_junction_beds(bed_dir)
    group_rows = []
    with open(group_path, "w") as group_handle:
        for index, filename in enumerate(beds):
            group_number, group_label = ("1", "exp") if index % 2 == 0 else ("2", "ctl")
            stream = filename[: -len(JUNCTION_SUFFIX)] + ".bed"
            group_handle.write("{0}\t{1}\t{2}\n".format(stream, group_number, group_label))
            group_rows.append((stream, group_number, group_label))
    with open(comparison_path, "w") as comparison_handle:
        comparison_handle.write("1\t2\n")
    return group_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bed-dir", required=True)
    parser.add_argument("--groups", required=True)
    parser.add_argument("--comparisons", required=True)
    args = parser.parse_args()
    build_groups(args.bed_dir, args.groups, args.comparisons)


if __name__ == "__main__":
    main()
