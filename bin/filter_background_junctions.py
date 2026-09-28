#!/usr/bin/env python3
"""Filter AltAnalyze junctions against GTEx and expand candidates by patient."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict

from scneo_utils import read_manifest, write_tsv
from splice_utils import load_documented_junctions


FIELDS = [
    "patient_id", "junction_id", "uid", "gtex_mean", "tumor_max",
    "patient_tumor_max", "tumor_gtex_delta", "junction_supporting_tumor_cell_ids",
    "junction_supporting_normal_cell_ids", "normal_callable_cell_ids",
    "normal_evidence_status", "tumor_enriched_junction", "tumor_specific_junction",
    "control_means", "control_deltas", "in_ensembl_transcript",
]


def read_count_matrix(path: str) -> tuple[list[str], list[tuple[str, dict[str, float]]]]:
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header or len(header) < 2:
            raise ValueError("count matrix has no sample columns")
        samples = header[1:]
        rows = []
        for line_number, fields in enumerate(reader, start=2):
            if len(fields) != len(header):
                raise ValueError(f"count matrix line {line_number} has the wrong number of columns")
            try:
                counts = {sample: float(value) for sample, value in zip(samples, fields[1:])}
            except ValueError as exc:
                raise ValueError(f"count matrix line {line_number} contains a non-numeric value") from exc
            rows.append((fields[0], counts))
    return samples, rows


def load_gtex_means(path: str) -> dict[str, float]:
    try:
        import anndata as ad
        import numpy as np
    except ImportError as exc:
        raise RuntimeError("reading GTEx H5AD requires anndata and numpy") from exc
    adata = ad.read_h5ad(path)
    means = np.asarray(adata.X.mean(axis=1)).reshape(-1)
    out = {}
    for uid, mean in zip(adata.obs_names.astype(str), means):
        out.setdefault(uid, float(mean))
    return out


def filter_rows(
    matrix_rows: list[tuple[str, dict[str, float]]],
    manifest: list[dict[str, str]],
    gtex_means: dict[str, float],
    n_max: float,
    t_min: float,
    min_reads: float,
    additional_controls: dict[str, dict[str, float]] | None = None,
    documented_junctions: set[str] | None = None,
    exclude_documented: bool = False,
) -> list[dict[str, object]]:
    controls = {"gtex": gtex_means, **(additional_controls or {})}
    by_patient_type: dict[tuple[str, str], list[str]] = defaultdict(list)
    tumor_cells = []
    for row in manifest:
        if row["cell_type"] == "tumor":
            by_patient_type[(row["patient_id"], row["cell_type"])].append(row["cell_id"])
        if row["cell_type"] == "tumor":
            tumor_cells.append(row["cell_id"])
    if not tumor_cells:
        raise ValueError("manifest has no tumor cells")
    patients = sorted({row["patient_id"] for row in manifest})
    output = []
    for uid, counts in matrix_rows:
        in_ensembl = uid in (documented_junctions or set())
        if exclude_documented and in_ensembl:
            continue
        missing = [cell for cell in tumor_cells if cell not in counts]
        if missing:
            raise ValueError("count matrix is missing tumor cells: " + ", ".join(missing))
        control_means = {name: means.get(uid, 0.0) for name, means in controls.items()}
        gtex_mean = control_means["gtex"]
        tumor_max = max(counts[cell] for cell in tumor_cells)
        for patient in patients:
            patient_tumor = by_patient_type[(patient, "tumor")]
            patient_tumor_max = max(counts.get(cell, 0.0) for cell in patient_tumor)
            control_deltas = {name: patient_tumor_max - mean for name, mean in control_means.items()}
            if any(mean >= n_max or control_deltas[name] <= t_min for name, mean in control_means.items()):
                continue
            supporting_tumor = [
                cell for cell in patient_tumor
                if counts.get(cell, 0.0) >= min_reads
                and all(counts.get(cell, 0.0) - mean > t_min for mean in control_means.values())
            ]
            if not supporting_tumor:
                continue
            output.append(
                {
                    "patient_id": patient,
                    "junction_id": uid,
                    "uid": uid,
                    "gtex_mean": f"{gtex_mean:.12g}",
                    "tumor_max": f"{tumor_max:.12g}",
                    "patient_tumor_max": f"{patient_tumor_max:.12g}",
                    "tumor_gtex_delta": f"{control_deltas['gtex']:.12g}",
                    "junction_supporting_tumor_cell_ids": ",".join(supporting_tumor),
                    "junction_supporting_normal_cell_ids": "",
                    "normal_callable_cell_ids": "",
                    "normal_evidence_status": "external_gtex_filter",
                    "tumor_enriched_junction": "true",
                    "tumor_specific_junction": "false",
                    "control_means": ";".join(f"{name}={mean:.12g}" for name, mean in control_means.items()),
                    "control_deltas": ";".join(f"{name}={delta:.12g}" for name, delta in control_deltas.items()),
                    "in_ensembl_transcript": "true" if in_ensembl else "false",
                }
            )
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--gtex", required=True)
    parser.add_argument("--additional-control", action="append", default=[], metavar="NAME=H5AD")
    parser.add_argument("--transcript-table", required=True)
    parser.add_argument("--exclude-documented-junctions", action="store_true")
    parser.add_argument("--n-max", type=float, default=3.0)
    parser.add_argument("--t-min", type=float, default=20.0)
    parser.add_argument("--min-reads", type=float, default=3.0)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    samples, matrix_rows = read_count_matrix(args.matrix)
    manifest = read_manifest(args.manifest)
    expected = [row["cell_id"] for row in manifest if row["cell_type"] == "tumor"]
    if samples != expected:
        raise ValueError("count matrix sample columns must exactly match manifest tumor/normal cell order")
    controls = {}
    for item in args.additional_control:
        name, separator, path = item.partition("=")
        if not separator or not name or not path:
            raise ValueError("additional-control must use NAME=H5AD format")
        if name == "gtex" or name in controls:
            raise ValueError(f"duplicate control name: {name}")
        controls[name] = load_gtex_means(path)
    rows = filter_rows(matrix_rows, manifest, load_gtex_means(args.gtex), args.n_max, args.t_min, args.min_reads, controls, load_documented_junctions(args.transcript_table), args.exclude_documented_junctions)
    print(f"GTEx filter retained {len(rows)} patient-junction rows", file=sys.stderr)
    write_tsv(args.out, rows, FIELDS)


if __name__ == "__main__":
    main()
