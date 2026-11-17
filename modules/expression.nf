process EXPRESSION_QUANTIFY {
    tag "${patient_id}:${cell_id}"
    label "cpu_medium"

    input:
    tuple val(patient_id), val(cell_id), path(bam), path(bai)
    path ref_gtf

    output:
    path "${patient_id}.${cell_id}.gene_counts.tsv", emit: counts
    path "${patient_id}.${cell_id}.featurecounts.txt", emit: raw

    script:
    """
    paired_args=""
    if [ "\$(${params.tool_samtools_cmd} view -c -f 1 ${bam})" -gt 0 ]; then
        paired_args="-p --countReadPairs"
    fi

    ${params.tool_featurecounts_cmd} \
        -T ${task.cpus} \
        -a ${ref_gtf} \
        -o ${patient_id}.${cell_id}.featurecounts.txt \
        -t exon \
        -g gene_id \
        --extraAttributes gene_name \
        \$paired_args \
        ${params.tool_featurecounts_extra_args} \
        ${bam}

    ${params.tool_python_cmd} ${workflow.projectDir}/bin/parse_featurecounts.py \
        --featurecounts ${patient_id}.${cell_id}.featurecounts.txt \
        --patient-id ${patient_id} \
        --cell-ids ${cell_id} \
        --out ${patient_id}.${cell_id}.gene_counts.tsv
    """
}

process EXPRESSION_NORMALIZE {
    tag "tumor_expression"
    label "utility_light"

    input:
    path counts

    output:
    path "tumor_gene_counts.tsv", emit: raw_matrix
    path "tumor_expression_normalized.tsv", emit: normalized
    path "tumor_expression_matrix.tsv", emit: normalized_matrix
    path "tumor_gene_tpm_matrix.tsv", emit: tpm_matrix

    script:
    def countArgs = counts.collect { file -> file }.join(' ')
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/normalize_expression.py \
        --counts ${countArgs} \
        --out-long tumor_expression_normalized.tsv \
        --out-raw-matrix tumor_gene_counts.tsv \
        --out-normalized-matrix tumor_expression_matrix.tsv \
        --out-tpm-matrix tumor_gene_tpm_matrix.tsv
    """
}

process CANDIDATE_EXPRESSION_SUPPORT {
    tag "candidate_expression"
    label "utility_light"

    input:
    path pmhc
    path expression
    path manifest

    output:
    path "candidate_expression_support.tsv", emit: support

    script:
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/candidate_expression_support.py \
        --pmhc ${pmhc} \
        --expression ${expression} \
        --manifest ${manifest} \
        --out candidate_expression_support.tsv
    """
}

process EXPRESSION_IMPORT_TPM {
    tag "supplied_gene_tpm"
    label "utility_light"

    input:
    path gene_tpm

    output:
    path "tumor_gene_counts.tsv", emit: raw_matrix
    path "tumor_expression_normalized.tsv", emit: normalized
    path "tumor_expression_matrix.tsv", emit: normalized_matrix
    path "tumor_gene_tpm_matrix.tsv", emit: tpm_matrix

    script:
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/normalize_expression.py \
        --tpm ${gene_tpm} \
        --out-long tumor_expression_normalized.tsv \
        --out-raw-matrix tumor_gene_counts.tsv \
        --out-normalized-matrix tumor_expression_matrix.tsv \
        --out-tpm-matrix tumor_gene_tpm_matrix.tsv
    """
}
