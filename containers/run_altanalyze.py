#!/usr/bin/env python
"""Run AltAnalyze with a bounded multiprocessing CPU count."""

from __future__ import print_function

import multiprocessing
import os
import runpy
import sys


def parse_core_limit(value):
    try:
        cores = int(value)
    except ValueError:
        raise ValueError("AltAnalyze core limit must be a positive integer: " + value)
    if cores < 1:
        raise ValueError("AltAnalyze core limit must be a positive integer: " + value)
    return cores


def run_altanalyze(core_limit, script, arguments):
    multiprocessing.cpu_count = lambda: core_limit
    sys.argv = [script] + list(arguments)
    sys.path.insert(0, os.path.dirname(os.path.abspath(script)))
    runpy.run_path(script, run_name="__main__")


def main():
    if len(sys.argv) < 3:
        raise SystemExit("Usage: run-altanalyze.py CORES ALTANALYZE_SCRIPT [ARGS...]")
    try:
        core_limit = parse_core_limit(sys.argv[1])
    except ValueError as exc:
        raise SystemExit("ERROR: " + str(exc))
    script = sys.argv[2]
    if not os.path.isfile(script):
        raise SystemExit("ERROR: AltAnalyze script not found: " + script)
    run_altanalyze(core_limit, script, sys.argv[3:])


if __name__ == "__main__":
    main()
