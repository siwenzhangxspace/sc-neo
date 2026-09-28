process HLA_CALL {
    tag "${patient_id}"
    label "cpu_threaded"
    errorStrategy { task.attempt < 2 ? 'retry' : 'terminate' }
    maxRetries 1

    input:
    tuple val(patient_id), val(cell_ids), val(read_layouts), path(read_files)

    output:
    tuple val(patient_id), path("${patient_id}.hla.tsv"), emit: hla
    tuple val(patient_id), path("optitype_out/*_result.tsv"), path("optitype_out/*_coverage_plot.pdf"), emit: raw

    script:
    def layouts = (read_layouts as Set).toList()
    if (layouts.size() > 1) {
        error "Patient ${patient_id}: mixed read layouts ${read_layouts}; one patient must share a single layout"
    }
    def layout = layouts[0]
    def n_cells = cell_ids.size()
    def expected_files = layout == "paired" ? n_cells * 2 : n_cells
    if (read_files.size() != expected_files) {
        error "Patient ${patient_id}: expected ${expected_files} FASTQ files for ${n_cells} cells (${layout} layout) but got ${read_files.size()}"
    }
    def merge_args
    def optitype_inputs
    if (layout == "paired") {
        def r1 = []
        def r2 = []
        read_files.eachWithIndex { f, i -> (i % 2 == 0 ? r1 : r2) << f }
        merge_args = "--r1 ${r1.join(" --r1 ")} \\\n        --r2 ${r2.join(" --r2 ")} \\\n        --out-r1 pool.R1.fq \\\n        --out-r2 pool.R2.fq"
        optitype_inputs = "-i pool.R1.fq \\\n        -i pool.R2.fq"
    } else {
        merge_args = "--r1 ${read_files.join(" --r1 ")} \\\n        --out-r1 pool.fq"
        optitype_inputs = "-i pool.fq"
    }
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/merge_fastq_prefix.py \\
        ${merge_args} \\
        --target-reads ${params.threshold_hla_max_reads}
    mkdir -p optitype_out
    ${params.tool_optitype_cmd} \\
        ${optitype_inputs} \\
        --rna \\
        -o optitype_out \\
        -p ${patient_id} \\
        --solver glpk \\
        --threads ${task.cpus} \\
        ${params.tool_optitype_extra_args}
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/parse_optitype.py \\
        --patient-id ${patient_id} \\
        --n-cells ${n_cells} \\
        --input-dir optitype_out \\
        --out ${patient_id}.hla.tsv
    """
}

process HLA_SUMMARIZE {
    tag "hla"
    label "utility_light"

    input:
    path input_manifest
    path per_patient_hla

    output:
    path "patient_hla_alleles.tsv", emit: patient

    script:
    def inputs = per_patient_hla instanceof List ? per_patient_hla.join(" ") : per_patient_hla
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/summarize_hla.py \\
        --manifest ${input_manifest} \\
        --per-patient ${inputs} \\
        --out-patient patient_hla_alleles.tsv
    """
}
