#!/usr/bin/env python3
"""Run netMHCpan per patient with HLA batches bounded by -a string length."""

from __future__ import annotations

import argparse
import shlex
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

from scneo_utils import read_tsv, safe_float, split_list, uniq, write_tsv


ALLELE_COLUMNS = ("MHC", "Allele", "HLA")
AFFINITY_COLUMNS = ("Aff(nM)", "Aff.nM")
RANK_EL_COLUMNS = ("%Rank_EL", "%Rank.EL", "Rank_EL", "Rnk_EL")
BIND_COLUMNS = ("BindLevel", "Bind")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peptides", required=True)
    parser.add_argument("--hla", required=True)
    parser.add_argument("--source-type", required=True, choices=["mutation", "splicing"])
    parser.add_argument("--cmd", default="netMHCpan")
    parser.add_argument("--extra-args", default="")
    parser.add_argument("--max-a-chars", type=int, default=900)
    parser.add_argument("--out-predictions", required=True)
    parser.add_argument("--out-batches", required=True)
    args = parser.parse_args()

    peptide_rows = [row for row in read_tsv(args.peptides) if row.get("peptide") or row.get("peptide_seq")]
    hla_rows = {row["patient_id"]: row for row in read_tsv(args.hla)}
    by_patient_length = defaultdict(list)
    for row in peptide_rows:
        row["peptide"] = row.get("peptide") or row.get("peptide_seq", "")
        length = peptide_length(row)
        if length is not None:
            by_patient_length[(row.get("patient_id", ""), length)].append(row)

    predictions = []
    batches = []
    for patient_id, length in sorted(by_patient_length):
        patient_peptides = by_patient_length[(patient_id, length)]
        hla_row = hla_rows.get(patient_id, {})
        alleles = uniq(format_allele(a) for a in split_list(hla_row.get("hla_alleles", "")) if format_allele(a))
        if not alleles:
            for pep in patient_peptides:
                predictions.append(status_row(pep, "", args.source_type, "missing_hla"))
            batches.append(
                {
                    "patient_id": patient_id,
                    "source_type": args.source_type,
                    "hla_batch_id": "",
                    "hla_batch_alleles": "",
                    "batch_status": "missing_hla",
                }
            )
            continue
        batch_alleles = make_batches(alleles, args.max_a_chars)
        unique_peptides = sorted({row["peptide"] for row in patient_peptides})
        for idx, batch in enumerate(batch_alleles, start=1):
            batch_id = f"{patient_id}_{args.source_type}_L{length}_hla_batch_{idx}"
            run_batch_recursive(
                patient_id,
                batch_id,
                batch,
                unique_peptides,
                patient_peptides,
                args,
                predictions,
                batches,
                peptide_length=length,
            )

    pred_fields = [
        "patient_id",
        "source_type",
        "peptide",
        "peptide_length",
        "hla_allele",
        "Aff.nM",
        "rank_el",
        "bind_level",
        "netmhcpan_status",
        "hla_batch_id",
        "hla_batch_alleles",
        "source_cell_ids",
        "supporting_tumor_cell_ids",
        "normal_alt_cell_ids",
        "normal_callable_cell_ids",
        "normal_evidence_status",
        "tumor_specific_candidate",
        "tumor_enriched_candidate",
        "variant_or_junction_id",
        "variant_key",
        "junction_id",
        "HGVSp_Short",
        "Variant_Classification",
        "variant_status",
        "variant_class",
        "gene",
        "transcript_id",
        "junction_annotation",
        "splice_site",
        "translation_frame",
        "n_from_first",
        "phase_supported",
        "phase_evidence",
        "reference_protein_match",
        "gtex_mean",
        "tumor_gtex_delta",
        "annotation_summary",
    ]
    batch_fields = ["patient_id", "source_type", "hla_batch_id", "hla_batch_alleles", "batch_status", "stderr"]
    write_tsv(args.out_predictions, predictions, pred_fields)
    write_tsv(args.out_batches, batches, batch_fields)


def run_batch_recursive(
    patient_id: str,
    batch_id: str,
    alleles: list[str],
    unique_peptides: list[str],
    peptide_rows: list[dict[str, str]],
    args: argparse.Namespace,
    predictions: list[dict[str, str]],
    batches: list[dict[str, str]],
    peptide_length: int | None = None,
) -> None:
    result, stderr = run_netmhcpan(args.cmd, args.extra_args, unique_peptides, alleles, peptide_length)
    if result is None:
        if not is_allele_batch_too_large(stderr):
            raise RuntimeError(
                f"netMHCpan failed for {patient_id} {batch_id}: {stderr.strip() or 'no stderr output'}"
            )
        if len(alleles) == 1:
            allele = alleles[0]
            for pep in peptide_rows:
                predictions.append(status_row(pep, allele, args.source_type, "netmhcpan_failed_allele", batch_id, alleles))
            batches.append(batch_row(patient_id, args.source_type, batch_id, alleles, "failed_single_allele", stderr))
            return
        mid = max(1, len(alleles) // 2)
        run_batch_recursive(patient_id, batch_id + "a", alleles[:mid], unique_peptides, peptide_rows, args, predictions, batches, peptide_length)
        run_batch_recursive(patient_id, batch_id + "b", alleles[mid:], unique_peptides, peptide_rows, args, predictions, batches, peptide_length)
        return

    parsed = parse_netmhcpan_output(result, set(unique_peptides), set(alleles))
    if unique_peptides and not parsed:
        raise RuntimeError(
            f"netMHCpan stdout format was not recognized for {patient_id} {batch_id}; "
            "the tool exited successfully and produced non-empty stdout, but no prediction rows "
            "matched the expected columns; "
            f"{summarize_netmhcpan_output(result)}; "
            "required columns: Peptide, one of MHC/Allele/HLA, and one of Aff(nM)/Aff.nM"
        )
    append_predictions(parsed, peptide_rows, args.source_type, batch_id, alleles, predictions)
    batches.append(batch_row(patient_id, args.source_type, batch_id, alleles, "ok", stderr))


def run_netmhcpan(
    cmd: str,
    extra_args: str,
    peptides: list[str],
    alleles: list[str],
    peptide_length: int | None = None,
) -> tuple[str | None, str]:
    with tempfile.TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir)
        pep_path = workdir / "peptides.txt"
        pep_path.write_text("\n".join(peptides) + "\n")
        command = [
            cmd,
            "-BA",
            "-p",
            str(pep_path),
            "-a",
            ",".join(alleles),
            "-s",
        ]
        if peptide_length is not None:
            command.extend(["-l", str(peptide_length)])
        if extra_args:
            command.extend(shlex.split(extra_args))
        try:
            proc = subprocess.run(
                command,
                cwd=workdir,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
            )
        except FileNotFoundError as exc:
            return None, str(exc)
        if proc.returncode != 0 or "PTYPE_LINE parameter option -a too long" in proc.stderr:
            return None, proc.stderr
        if not proc.stdout.strip():
            diagnostic = proc.stderr.strip() or "netMHCpan produced empty stdout"
            return None, diagnostic
        return proc.stdout, proc.stderr


def is_allele_batch_too_large(stderr: str) -> bool:
    return "PTYPE_LINE parameter option -a too long" in stderr


def parse_netmhcpan_output(text: str, peptides: set[str], alleles: set[str]) -> list[dict[str, str]]:
    header = None
    rows = []
    normalized_alleles = {format_allele(allele) or allele for allele in alleles}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "-", "Number")):
            continue
        fields = line.split()
        if "Peptide" in fields and any(column in fields for column in ALLELE_COLUMNS):
            header = fields
            continue
        if header is None:
            continue
        fields = align_netmhcpan_fields(header, fields)
        if fields is None:
            continue
        values = dict(zip(header, fields))
        peptide = values.get("Peptide", "")
        allele = first_column_value(values, ALLELE_COLUMNS)
        if not peptide or not allele:
            continue
        normalized_allele = format_allele(allele) or allele
        if peptide not in peptides or normalized_allele not in normalized_alleles:
            continue
        bind = first_column_value(values, BIND_COLUMNS).replace("<=", "").strip().lower()
        aff = first_column_value(values, AFFINITY_COLUMNS)
        if safe_float(aff, None) is None:
            continue
        rank_el = first_column_value(values, RANK_EL_COLUMNS)
        rows.append(
            {
                "peptide": peptide,
                "hla_allele": normalized_allele,
                "Aff.nM": aff,
                "rank_el": rank_el,
                "bind_level": bind,
                "netmhcpan_status": "ok",
            }
        )
    return rows


def align_netmhcpan_fields(header: list[str], fields: list[str]) -> list[str] | None:
    if header[-1] not in BIND_COLUMNS:
        return fields if len(fields) == len(header) else None
    if len(fields) == len(header) - 1:
        return fields + [""]
    if len(fields) < len(header):
        return None
    # netMHCpan may print thresholds such as "<= SB". Since the table is
    # whitespace-delimited, keep all trailing bind-level tokens together.
    return fields[: len(header) - 1] + [" ".join(fields[len(header) - 1 :])]


def first_column_value(values: dict[str, str], columns: tuple[str, ...]) -> str:
    return next((values[column] for column in columns if values.get(column)), "")


def summarize_netmhcpan_output(text: str) -> str:
    nonempty_lines = [line.strip() for line in text.splitlines() if line.strip()]
    candidate_header = None
    candidate_index = None
    for index, line in enumerate(nonempty_lines):
        if line.startswith(("#", "-", "Number")):
            continue
        fields = line.split()
        if "Peptide" in fields:
            candidate_header = fields
            candidate_index = index
            break
    if candidate_header is None or candidate_index is None:
        return f"nonempty_lines={len(nonempty_lines)}, candidate_header=<not found>"

    following_field_counts = set()
    for line in nonempty_lines[candidate_index + 1 :]:
        if line.startswith(("#", "-", "Number")):
            continue
        fields = line.split()
        if "Peptide" not in fields:
            following_field_counts.add(len(fields))
    columns = ",".join(candidate_header)
    counts = ",".join(str(count) for count in sorted(following_field_counts)) or "<none>"
    return (
        f"nonempty_lines={len(nonempty_lines)}, candidate_header_columns={columns}, "
        f"following_field_counts={counts}"
    )


def append_predictions(
    parsed: list[dict[str, str]],
    peptide_rows: list[dict[str, str]],
    source_type: str,
    batch_id: str,
    alleles: list[str],
    predictions: list[dict[str, str]],
) -> None:
    by_peptide = defaultdict(list)
    for row in peptide_rows:
        by_peptide[row["peptide"]].append(row)
    for pred in parsed:
        for meta in by_peptide[pred["peptide"]]:
            merged = base_prediction(meta, pred.get("hla_allele", ""), source_type, pred.get("netmhcpan_status", "ok"), batch_id, alleles)
            merged.update(
                {
                    "Aff.nM": pred.get("Aff.nM", ""),
                    "rank_el": pred.get("rank_el", ""),
                    "bind_level": pred.get("bind_level", ""),
                }
            )
            predictions.append(merged)


def base_prediction(meta: dict[str, str], allele: str, source_type: str, status: str, batch_id: str = "", alleles: list[str] | None = None) -> dict[str, str]:
    return {
        "patient_id": meta.get("patient_id", ""),
        "source_type": source_type,
        "peptide": meta.get("peptide", ""),
        "peptide_length": meta.get("peptide_length", len(meta.get("peptide", ""))),
        "hla_allele": allele,
        "Aff.nM": "",
        "rank_el": "",
        "bind_level": "",
        "netmhcpan_status": status,
        "hla_batch_id": batch_id,
        "hla_batch_alleles": ",".join(alleles or []),
        "source_cell_ids": meta.get("source_cell_ids", ""),
        "supporting_tumor_cell_ids": meta.get("supporting_tumor_cell_ids", ""),
        "normal_alt_cell_ids": meta.get("normal_alt_cell_ids", ""),
        "normal_callable_cell_ids": meta.get("normal_callable_cell_ids", ""),
        "normal_evidence_status": meta.get("normal_evidence_status", ""),
        "tumor_specific_candidate": meta.get("tumor_specific_candidate", ""),
        "tumor_enriched_candidate": meta.get("tumor_enriched_candidate", ""),
        "variant_or_junction_id": meta.get("variant_or_junction_id") or meta.get("variant_key") or meta.get("junction_id", ""),
        "variant_key": meta.get("variant_key", ""),
        "junction_id": meta.get("junction_id", ""),
        "HGVSp_Short": meta.get("HGVSp_Short", ""),
        "Variant_Classification": meta.get("Variant_Classification", ""),
        "variant_status": meta.get("variant_status", ""),
        "variant_class": meta.get("variant_class", ""),
        "gene": meta.get("gene", ""),
        "transcript_id": meta.get("transcript_id", ""),
        "junction_annotation": meta.get("junction_annotation", ""),
        "splice_site": meta.get("splice_site", ""),
        "translation_frame": meta.get("translation_frame", ""),
        "n_from_first": meta.get("n_from_first", ""),
        "phase_supported": meta.get("phase_supported", ""),
        "phase_evidence": meta.get("phase_evidence", ""),
        "reference_protein_match": meta.get("reference_protein_match", ""),
        "gtex_mean": meta.get("gtex_mean", ""),
        "tumor_gtex_delta": meta.get("tumor_gtex_delta", ""),
        "annotation_summary": meta.get("annotation_summary", ""),
    }


def status_row(meta: dict[str, str], allele: str, source_type: str, status: str, batch_id: str = "", alleles: list[str] | None = None) -> dict[str, str]:
    return base_prediction(meta, allele, source_type, status, batch_id, alleles)


def batch_row(patient_id: str, source_type: str, batch_id: str, alleles: list[str], status: str, stderr: str) -> dict[str, str]:
    return {
        "patient_id": patient_id,
        "source_type": source_type,
        "hla_batch_id": batch_id,
        "hla_batch_alleles": ",".join(alleles),
        "batch_status": status,
        "stderr": stderr.replace("\n", " ")[:1000],
    }


def make_batches(alleles: list[str], max_chars: int) -> list[list[str]]:
    batches = []
    current = []
    for allele in alleles:
        proposed = current + [allele]
        if current and len(",".join(proposed)) > max_chars:
            batches.append(current)
            current = [allele]
        else:
            current = proposed
    if current:
        batches.append(current)
    return batches


def peptide_length(row: dict[str, str]) -> int | None:
    peptide = row.get("peptide") or row.get("peptide_seq", "")
    if not peptide:
        return None
    try:
        length = int(row.get("peptide_length", "") or len(peptide))
    except ValueError:
        length = len(peptide)
    return length if length > 0 else None


def format_allele(value: str) -> str:
    value = (value or "").strip()
    if not value:
        return ""
    value = value.replace("HLA-", "")
    if "*" in value:
        gene, rest = value.split("*", 1)
        return f"HLA-{gene}{rest}"
    if value.startswith(("A", "B", "C")):
        return "HLA-" + value
    if value.startswith("HLA"):
        return value
    return ""


if __name__ == "__main__":
    main()
