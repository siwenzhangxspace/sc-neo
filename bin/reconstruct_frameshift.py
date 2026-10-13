"""Transcript-specific nucleotide reconstruction for small frameshift variants."""

from __future__ import annotations

import re

from splice_utils import CODON


class RejectedVariant(ValueError):
    def __init__(self, reason: str, detail: str):
        super().__init__(detail)
        self.reason = reason


def reject(reason, detail):
    raise RejectedVariant(reason, detail)


def base_id(identifier):
    return re.sub(r"\.\d+$", "", identifier)


def resolve(identifier, keys):
    if identifier in keys:
        return identifier
    matches = [key for key in keys if base_id(key) == base_id(identifier)
               and (key == base_id(key) or identifier == base_id(identifier))]
    if len(matches) != 1:
        reject("reference_id_unresolved", f"{identifier}: {len(matches)} compatible matches")
    return matches[0]


def transcript_id(row):
    return next((row[k] for k in ("Transcript_ID", "Feature", "transcript_id", "ENST")
                 if row.get(k) not in (None, "", ".")), "")


def is_frameshift(row):
    terms = re.split(r"[,;& ]+", row.get("Consequence", "") + " " + row.get("Variant_Classification", ""))
    return bool(set(terms) & {"frameshift_variant", "Frame_Shift_Ins", "Frame_Shift_Del"}) or "fs" in (row.get("HGVSp_Short", "") + row.get("HGVSp", ""))


def load_models(path, wanted):
    models = {}
    feature_count = 0
    with open(path) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"GTF line {number}: expected nine fields")
            feature_count += 1
            chrom, _, feature, start, end, _, strand, phase, attributes = fields
            start, end = int(start), int(end)
            if start < 1 or end < start:
                raise ValueError(f"GTF line {number}: invalid interval")
            attrs = dict(re.findall(r'(\S+)\s+"([^"]*)"', attributes))
            tid = attrs.get("transcript_id", "")
            if base_id(tid) not in wanted:
                continue
            if attrs.get("transcript_version") and tid == base_id(tid):
                tid += "." + attrs["transcript_version"]
            model = models.setdefault(tid, {"chrom": chrom, "strand": strand, "exon": [], "CDS": [], "partial": False})
            if (chrom, strand) != (model["chrom"], model["strand"]):
                model["partial"] = True
            if feature in ("exon", "CDS"):
                model[feature].append((start, end, phase))
            if any(tag in attributes for tag in ("cds_start_NF", "cds_end_NF", "selenocysteine")) or feature == "Selenocysteine":
                model["partial"] = True
    if not feature_count:
        raise ValueError("GTF contains no features")
    return models


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def translate(sequence):
    protein = []
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i + 3]
        if codon not in CODON:
            reject("unknown_codon", f"Ambiguous codon at nucleotide {i + 1}")
        aa = CODON[codon]
        if aa == "*":
            return "".join(protein), "stop_codon"
        protein.append(aa)
    return "".join(protein), "transcript_end"


def reconstruct(row, models, fasta):
    tid = resolve(transcript_id(row), models)
    model = models[tid]
    chrom, strand = model["chrom"], model["strand"]
    if chrom.removeprefix("chr") in {"M", "MT"} or strand not in {"+", "-"} or model["partial"]:
        reject("unsupported_transcript", tid)
    if not model["exon"] or not model["CDS"]:
        reject("incomplete_cds", tid)
    if chrom not in fasta.references:
        reject("missing_contig", chrom)
    positions, chunks = [], []
    exons = sorted(model["exon"], reverse=strand == "-")
    exon_index = {}
    for index, (start, end, _) in enumerate(exons):
        seq = fasta.fetch(chrom, start - 1, end).upper()
        if len(seq) != end - start + 1:
            reject("reference_interval_out_of_bounds", tid)
        coords = list(range(start, end + 1))
        if strand == "-":
            coords.reverse()
            seq = reverse_complement(seq)
        positions.extend(coords)
        chunks.append(seq)
        exon_index.update({p: index for p in coords})
    if len(set(positions)) != len(positions):
        reject("overlapping_exons", tid)
    offsets = {p: i for i, p in enumerate(positions)}
    sequence = "".join(chunks)
    cds_offsets = []
    for start, end, phase in sorted(model["CDS"], reverse=strand == "-"):
        if phase not in {"0", "1", "2"} or int(phase) != (-len(cds_offsets)) % 3:
            reject("invalid_cds_phase", tid)
        coords = range(start, end + 1) if strand == "+" else range(end, start - 1, -1)
        if any(p not in offsets for p in coords):
            reject("cds_outside_exons", tid)
        cds_offsets.extend(offsets[p] for p in coords)
    if not cds_offsets or cds_offsets != list(range(cds_offsets[0], cds_offsets[-1] + 1)) or len(cds_offsets) % 3:
        reject("incomplete_cds", tid)
    begin, finish = cds_offsets[0], cds_offsets[-1] + 1
    coding = sequence[begin:finish]
    if not coding.startswith("ATG"):
        reject("incomplete_cds", "CDS does not start with ATG")
    wt, _ = translate(coding)
    if len(wt) * 3 not in {len(coding), len(coding) - 3}:
        reject("internal_reference_stop", tid)
    # A complete annotated CDS must end at a real stop, included in CDS or just after it.
    stop_offset = begin + len(wt) * 3
    if sequence[stop_offset:stop_offset + 3] not in {"TAA", "TAG", "TGA"}:
        reject("incomplete_cds", "Missing reference stop codon")
    parts = row.get("variant_key", "").split("|")
    if len(parts) != 5 or not parts[1].isdigit() or int(parts[1]) < 1 or not all(re.fullmatch("[ACGT]+", v) for v in parts[2:4]):
        reject("invalid_variant_key", row.get("variant_key", ""))
    vchrom, pos, ref, alt, _ = parts
    pos = int(pos)
    if vchrom.removeprefix("chr") != chrom.removeprefix("chr"):
        reject("chromosome_mismatch", tid)
    if fasta.fetch(chrom, pos - 1, pos - 1 + len(ref)).upper() != ref:
        reject("reference_allele_mismatch", row["variant_key"])
    while ref and alt and ref[0] == alt[0]:
        pos, ref, alt = pos + 1, ref[1:], alt[1:]
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    if (len(alt) - len(ref)) % 3 == 0:
        reject("not_frameshift", row["variant_key"])
    coords = list(range(pos, pos + len(ref))) if ref else [pos - 1, pos]
    if any(p not in offsets for p in coords) or len({exon_index[p] for p in coords}) != 1:
        reject("unsupported_splice_boundary", row["variant_key"])
    if ref:
        left = min(offsets[p] for p in coords)
        right = left + len(ref)
    else:
        left = right = max(offsets[p] for p in coords)
    if left < begin + 3 or left >= stop_offset or right > stop_offset:
        reject("unsupported_cds_boundary", row["variant_key"])
    if strand == "-":
        alt = reverse_complement(alt)
    mutant, termination = translate((sequence[:left] + alt + sequence[right:])[begin:])
    first = 0
    while first < min(len(wt), len(mutant)) and wt[first] == mutant[first]:
        first += 1
    return mutant, first + 1, len(mutant), termination, tid
