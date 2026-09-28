#!/usr/bin/env python
"""Prepare a task-local AltAnalyze application with writable runtime data."""

from __future__ import print_function

import argparse
import os
import shutil
import sys


def require_directory(path, label):
    if not os.path.isdir(path):
        raise ValueError("{0} directory not found: {1}".format(label, path))


def link_children(source, destination, excluded_names):
    os.mkdir(destination)
    for name in sorted(os.listdir(source)):
        if name in excluded_names:
            continue
        os.symlink(os.path.join(source, name), os.path.join(destination, name))


def mutable_rnaseq_entries(species):
    return {
        "chr",
        species + "_Ensembl_exons.txt",
        species + "_Ensembl_exons.bed",
        species + "_Ensembl_junctions.txt",
        species + "_Ensembl_junctions.bed",
        species + "_junction_comps_updated.txt",
    }


def prepare_runtime(source_app, runtime_app, version, species):
    source_app = os.path.abspath(source_app)
    runtime_app = os.path.abspath(runtime_app)
    source_config = os.path.join(source_app, "Config")
    source_database = os.path.join(source_app, "AltDatabase")
    source_version = os.path.join(source_database, version)
    source_species = os.path.join(source_version, species)
    source_rnaseq = os.path.join(source_species, "RNASeq")

    require_directory(source_app, "AltAnalyze application")
    require_directory(source_config, "AltAnalyze Config")
    require_directory(source_database, "AltAnalyze database")
    require_directory(source_version, "AltAnalyze database version")
    require_directory(source_species, "AltAnalyze species")
    require_directory(source_rnaseq, "AltAnalyze RNASeq")
    if os.path.lexists(runtime_app):
        raise ValueError("AltAnalyze runtime path already exists: " + runtime_app)

    link_children(source_app, runtime_app, {"AltDatabase", "Config"})
    shutil.copytree(source_config, os.path.join(runtime_app, "Config"), symlinks=True)

    runtime_database = os.path.join(runtime_app, "AltDatabase")
    link_children(source_database, runtime_database, {version})

    runtime_version = os.path.join(runtime_database, version)
    link_children(source_version, runtime_version, {species})

    runtime_species = os.path.join(runtime_version, species)
    link_children(source_species, runtime_species, {"RNASeq"})

    runtime_rnaseq = os.path.join(runtime_species, "RNASeq")
    link_children(source_rnaseq, runtime_rnaseq, mutable_rnaseq_entries(species))
    os.mkdir(os.path.join(runtime_rnaseq, "chr"))
    return runtime_app


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-app", required=True)
    parser.add_argument("--runtime-app", required=True)
    parser.add_argument("--version", required=True)
    parser.add_argument("--species", required=True)
    args = parser.parse_args()
    try:
        prepare_runtime(args.source_app, args.runtime_app, args.version, args.species)
    except (OSError, ValueError) as exc:
        raise SystemExit("ERROR: " + str(exc))
    print("Prepared task-local AltAnalyze runtime", file=sys.stderr)


if __name__ == "__main__":
    main()
