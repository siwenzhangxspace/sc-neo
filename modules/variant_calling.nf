process VARIANT_CALL {
    tag "${cell_id}"
    label "cpu_medium"

    input:
    tuple val(patient_id), val(cell_id), val(cell_type), path(bam), path(bai), path(sj), path(star_log), path(qc)

    output:
    tuple val(patient_id), val(cell_id), val(cell_type), path("${cell_id}.vcf"), emit: vcf

    script:
    """
    ${params.tool_variant_caller_cmd} \\
        -b ${bam} \\
        -g ${params.ref_variant_caller_fasta ?: params.ref_genome_fasta} \\
        -o ${cell_id} \\
        ${params.tool_variant_caller_extra_args}
    """
}

process VARIANT_FILTER_QUALITY {
    tag "${cell_id}"
    label "cpu_threaded"

    input:
    tuple val(patient_id), val(cell_id), val(cell_type), path(vcf)

    output:
    tuple val(patient_id), val(cell_id), val(cell_type), path("${cell_id}.filtered.vcf"), emit: vcf

    script:
    """
    ${params.tool_bcftools_cmd} filter -i 'QUAL>=10 && FILTER=="PASS"' ${vcf} --threads ${task.cpus} -o ${cell_id}.filtered.vcf
    """
}

process VARIANT_SELECT_RECURRENT {
    tag "recurrent_variants"
    label "utility_light"

    input:
    path input_manifest
    path vcfs

    output:
    path "recurrent_tumor_variants.tsv", emit: tsv
    path "recurrent_tumor_variants.vcf", emit: vcf

    script:
    def inputVcfs = vcfs instanceof List ? vcfs.join(" ") : vcfs
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/select_recurrent_variants.py \\
        --manifest ${input_manifest} \\
        --vcfs ${inputVcfs} \\
        --min-tumor-cells ${params.threshold_min_tumor_variant_cells} \\
        --out-tsv recurrent_tumor_variants.tsv \\
        --out-vcf recurrent_tumor_variants.vcf
    """
}

process VARIANT_TEST_ENRICHMENT {
    tag "variant_enrichment"
    label "cpu_medium"

    input:
    path input_manifest
    path recurrent_tsv
    path vcfs
    path evidence

    output:
    path "tumor_enrichment.tsv", emit: tsv
    path "tumor_enriched_variants.vcf", emit: vcf

    script:
    def inputVcfs = vcfs instanceof List ? vcfs.join(" ") : vcfs
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/test_variant_enrichment.py \\
        --manifest ${input_manifest} \\
        --recurrent-tsv ${recurrent_tsv} \\
        --vcfs ${inputVcfs} \\
        --evidence ${evidence} \\
        --fdr ${params.threshold_variant_enrichment_fdr} \\
        --out-tsv tumor_enrichment.tsv \\
        --out-vcf tumor_enriched_variants.vcf
    """
}

process VARIANT_BUILD_CELL_EVIDENCE {
    tag "variant_evidence"
    label "cpu_medium"

    input:
    path input_manifest
    path variants
    path bams

    output:
    path "variant_cell_evidence.tsv", emit: evidence

    script:
    def inputBams = bams instanceof List ? bams.join(" ") : bams
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/build_variant_cell_evidence.py \\
        --manifest ${input_manifest} \\
        --variants ${variants} \\
        --bams ${inputBams} \\
        --min-depth ${params.threshold_variant_min_depth} \\
        --min-alt-reads ${params.threshold_variant_min_alt_reads} \\
        --min-vaf ${params.threshold_variant_min_vaf} \\
        --min-base-quality ${params.threshold_variant_min_base_quality} \\
        --min-mapping-quality ${params.threshold_variant_min_mapping_quality} \\
        --out variant_cell_evidence.tsv
    """
}

