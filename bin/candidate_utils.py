"""Patient-level candidate identity and explicit evidence conflict checks."""

from collections import defaultdict
from decimal import Decimal, InvalidOperation

from scneo_utils import split_list


IDENTITY_FIELDS = [
    "patient_id", "source_type", "gene", "peptide", "peptide_length", "hla_allele",
    "variant_or_junction_id", "var_key", "HGVSp_Short", "protein_change_status",
]
NUMERIC_FIELDS = frozenset({
    "Aff.nM", "rank_el", "%Rank_EL", "immunogenicity_score", "IEDB_immunogenicity_score",
    "healthy_blast_identity", "healthy_blast_evalue", "healthy_blast_alignment_length",
    "blast_identity", "blast_evalue", "max_TABR_BERT_rank", "RNA_expression_support",
    "peptide_length", "n_original_rows", "n_barcode",
})


def clean(value):
    return "" if value is None else str(value).strip()


def alias(row, *fields):
    values = {clean(row.get(field)) for field in fields if clean(row.get(field))}
    if len(values) > 1:
        raise ValueError(f"Conflicting identity aliases {fields}: {sorted(values)}")
    return next(iter(values), "")


def candidate_identity(row):
    source = clean(row.get("source_type"))
    peptide = alias(row, "peptide", "peptide_sequence", "peptide_seq")
    length = clean(row.get("peptide_length"))
    if length and (not length.isdigit() or int(length) != len(peptide)):
        raise ValueError(f"peptide_length disagrees with sequence: {length!r}, {peptide!r}")
    if source == "mutation":
        event = alias(row, "var_key", "variant_key", "variant_or_junction_id")
        change = clean(row.get("HGVSp_Short"))
        event_id = event
        status = "present" if change else "missing"
    elif source == "splicing":
        event = alias(row, "var_key", "splice_site")
        change = alias(row, "junction_id", "HGVSp_Short", "variant_or_junction_id")
        event_id = change
        status = "junction_id"
    else:
        raise ValueError(f"Unknown candidate source_type: {source!r}")
    identity = dict(patient_id=clean(row.get("patient_id")), source_type=source,
                    gene=clean(row.get("gene")), peptide=peptide, peptide_length=str(len(peptide)),
                    hla_allele=alias(row, "hla_allele", "HLA_allele"), var_key=event,
                    HGVSp_Short=change, variant_or_junction_id=event_id,
                    protein_change_status=status)
    required = ["patient_id", "gene", "peptide", "hla_allele", "var_key"]
    if source == "splicing":
        required.append("HGVSp_Short")
    missing = [field for field in required if not identity[field]]
    if missing:
        raise ValueError(f"Missing candidate identity fields {missing}: patient={identity['patient_id']!r}")
    return identity


def candidate_key(row):
    identity = candidate_identity(row)
    return tuple(identity[field] for field in (
        "patient_id", "gene", "hla_allele", "peptide_length", "peptide", "var_key", "HGVSp_Short"
    ))


def grouped_candidates(rows):
    groups = defaultdict(list)
    for row in rows:
        identity = candidate_identity(row)
        groups[candidate_key(identity)].append({**row, **identity})
    return dict(sorted(groups.items()))


def supporting_cells(rows):
    return sorted({cell for row in rows for cell in split_list(row.get("supporting_tumor_cell_ids", ""))})


def normalized_value(field, value):
    value = clean(value)
    if not value or field not in NUMERIC_FIELDS:
        return value
    try:
        number = Decimal(value)
    except InvalidOperation as exc:
        raise ValueError(f"Invalid numeric evidence {field}={value!r}") from exc
    if not number.is_finite():
        raise ValueError(f"Non-finite numeric evidence {field}={value!r}")
    return format(number.normalize(), "f") if number else "0"


def consistent_row(rows, fields, context):
    result = {}
    for field in fields:
        values = {normalized_value(field, row.get(field)) for row in rows}
        if len(values) > 1:
            raise ValueError(f"Conflicting evidence for {context!r}: {field}={sorted(values)!r}")
        result[field] = next(iter(values), "")
    return result


def evidence_index(rows, key_function, fields):
    groups = defaultdict(list)
    for row in rows:
        groups[key_function(row)].append(row)
    return {key: consistent_row(group, fields, key) for key, group in groups.items()}
