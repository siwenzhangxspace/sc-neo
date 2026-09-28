#!/usr/bin/env bash
set -euo pipefail

if [[ ${1:-} != "identify" || $# -ne 3 ]]; then
    echo "Usage: altanalyze-nextflow identify BAM_DIR CORES" >&2
    exit 2
fi

root=$PWD
bam_folder="${root}/$2"
cores=$3
app=/usr/src/app
source_app="${app}/altanalyze"
runtime_app="${root}/altanalyze_runtime"
reference_bed="${app}/altanalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs.bed"
reference_root="${root}/altanalyze_reference"

if [[ ! -d ${bam_folder} ]]; then
    echo "BAM directory not found: ${bam_folder}" >&2
    exit 2
fi
if [[ ! ${cores} =~ ^[1-9][0-9]*$ ]]; then
    echo "AltAnalyze cores must be a positive integer: ${cores}" >&2
    exit 2
fi

mkdir tmp
export TMPDIR="${root}/tmp"
python /usr/local/bin/prepare-altanalyze-runtime.py \
    --source-app "${source_app}" \
    --runtime-app "${runtime_app}" \
    --version EnsMart91 \
    --species Hs

run_bam_to_bed() {
    local bam=$1
    local sample_reference_dir="${reference_root}/${bam%.bam}"
    mkdir -p "${sample_reference_dir}"
    ln -sfn "${reference_bed}" "${sample_reference_dir}/Hs.bed"
    python "${app}/altanalyze/import_scripts/BAMtoJunctionBED.py" \
        --i "${bam_folder}/${bam}" \
        --species Hs \
        --r "${app}/altanalyze/AltDatabase/EnsMart91/ensembl/Hs/Hs_Ensembl_exon.txt"
    python "${app}/altanalyze/import_scripts/BAMtoExonBED.py" \
        --i "${bam_folder}/${bam}" \
        --r "${sample_reference_dir}/Hs.bed" \
        --s Hs
}
export -f run_bam_to_bed
export app bam_folder reference_bed reference_root

find -L "${bam_folder}" -maxdepth 1 -type f -name '*.bam' -printf '%f\n' | sort > samples.txt
if [[ $(wc -l < samples.txt) -lt 2 ]]; then
    echo "AltAnalyze requires at least two BAM files" >&2
    exit 2
fi
parallel -P "${cores}" run_bam_to_bed {} < samples.txt

mkdir bed altanalyze_output altanalyze_output/ExpressionInput
find "${bam_folder}" -maxdepth 1 -type f -name '*.bed' -exec mv {} bed/ \;

group_file="${root}/altanalyze_output/ExpressionInput/groups.original.txt"
comparison_file="${root}/altanalyze_output/ExpressionInput/comps.original.txt"
python /usr/local/bin/build-altanalyze-groups.py \
    --bed-dir "${root}/bed" \
    --groups "${group_file}" \
    --comparisons "${comparison_file}"

python /usr/local/bin/run-altanalyze.py "${cores}" "${runtime_app}/AltAnalyze.py" \
    --species Hs \
    --platform RNASeq \
    --version EnsMart91 \
    --bedDir "${root}/bed" \
    --output "${root}/altanalyze_output" \
    --groupdir "${group_file}" \
    --compdir "${comparison_file}" \
    --expname original \
    --runGOElite no

python /usr/local/bin/prune-altanalyze.py "${root}/altanalyze_output"
