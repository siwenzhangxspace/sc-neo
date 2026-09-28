process SPLICE_FILTER_MANIFEST {
    tag "splice_tumor_manifest"
    label "utility_light"

    input:
    path input_manifest

    output:
    path "splice_tumor_manifest.tsv", emit: manifest

    script:
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/filter_manifest_cells.py \\
        --manifest ${input_manifest} --cell-type tumor --out splice_tumor_manifest.tsv
    """
}

process SPLICE_BUILD_MATRIX {
    tag "altanalyze_junction_matrix"
    label "mem_high_threaded"

    input:
    path aligned_files

    output:
    path "counts.original.pruned.raw.txt", emit: raw_matrix
    path "groups.original.txt", emit: groups
    path "comps.original.txt", emit: comparisons

    script:
    """
    mkdir -p bam
    for input_bam in *.sorted.bam; do
        cell_id=\${input_bam%.sorted.bam}
        if [[ "\${cell_id}" == *"__"* || "\${cell_id}" =~ [[:space:]] ]]; then
            echo "AltAnalyze cell_id cannot contain whitespace or '__': \${cell_id}" >&2
            exit 2
        fi
        ln -s "\${PWD}/\${input_bam}" "bam/\${cell_id}.bam"
        ln -s "\${PWD}/\${input_bam}.bai" "bam/\${cell_id}.bam.bai"
    done

    ${params.tool_altanalyze_cmd} identify bam ${task.cpus}

    cp altanalyze_output/ExpressionInput/counts.original.pruned.txt \\
        counts.original.pruned.raw.txt
    cp altanalyze_output/ExpressionInput/groups.original.txt groups.original.txt
    cp altanalyze_output/ExpressionInput/comps.original.txt comps.original.txt
    """
}
process SPLICE_NORMALIZE_MATRIX {
    tag "normalize_altanalyze_matrix"
    label "io_light"

    input:
    path input_manifest
    path raw_matrix

    output:
    path "counts.original.pruned.txt", emit: matrix

    script:
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/normalize_junction_matrix.py \\
        --matrix ${raw_matrix} \\
        --manifest ${input_manifest} \\
        --out counts.original.pruned.txt
    """
}

process SPLICE_FILTER_BACKGROUND {
    tag "gtex_junction_filter"
    label "mem_high_threaded"

    input:
    path input_manifest
    path matrix
    path gtex
    path transcript_table
    path tcga_control

    output:
    path "filtered_junctions.tsv", emit: junctions

    script:
    def additionalControl = tcga_control ? "--additional-control tcga=${tcga_control}" : ""
    def excludeDocumented = params.enable_splice_exclude_documented_junctions ? "--exclude-documented-junctions" : ""
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/filter_background_junctions.py \\
        --matrix ${matrix} \\
        --manifest ${input_manifest} \\
        --gtex ${gtex} \\
        ${additionalControl} \\
        --transcript-table ${transcript_table} \\
        ${excludeDocumented} \\
        --n-max ${params.threshold_splice_gtex_n_max} \\
        --t-min ${params.threshold_splice_gtex_t_min} \\
        --min-reads ${params.threshold_splice_min_reads} \\
        --out filtered_junctions.tsv
    """
}
