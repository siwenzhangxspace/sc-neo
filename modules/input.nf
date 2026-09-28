process INPUT_VALIDATE_MANIFEST {
    tag "input_manifest"
    label "utility_light"

    input:
    path input_manifest

    output:
    path "manifest.validated.tsv", emit: manifest
    path "manifest.summary.tsv", emit: summary

    script:
    def singleEnd = params.allow_single_end_reads ? "--allow-single-end" : ""
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/validate_manifest.py \\
        --manifest ${input_manifest} \\
        --out-manifest manifest.validated.tsv \\
        --summary manifest.summary.tsv \\
        --min-cells-per-type ${params.threshold_min_cells_per_type} \\
        ${singleEnd}
    """
}

