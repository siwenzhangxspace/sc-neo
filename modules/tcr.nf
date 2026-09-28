process TCR_ASSEMBLE {
    tag "trust4_batch"
    label "mem_extreme_batch"

    input:
    path fastqs
    path ref_trust4_sequence, stageAs: 'trust4_sequence.fa'
    path ref_trust4_imgt, stageAs: 'trust4_imgt.fa'

    output:
    path "TRUST_report.tsv", emit: report
    tuple path("TRUST_report.tsv"), path("TRUST_annot.fa"), emit: raw

    script:
    def staged_fastqs = fastqs instanceof List ? fastqs.flatten() : [fastqs]
    def r1_files = staged_fastqs.findAll { f -> f.name.endsWith(".R1.clean.fastq.gz") }.sort { a, b -> a.name <=> b.name }
    def r2_files = staged_fastqs.findAll { f -> f.name.endsWith(".R2.clean.fastq.gz") }.sort { a, b -> a.name <=> b.name }
    if (r1_files.isEmpty()) {
        throw new IllegalArgumentException("TRUST4 received no R1 FASTQ files")
    }
    if (!r2_files.isEmpty() && r1_files.size() != r2_files.size()) {
        throw new IllegalArgumentException("TRUST4 Smart-seq batch requires all T cells to have the same read layout")
     }
    def decompress_fastq = params.trust4_decompress_fastq instanceof Boolean ? params.trust4_decompress_fastq : params.trust4_decompress_fastq.toString().toBoolean()
    def decompress_cmd = params.trust4_decompress_cmd
    def r1_input_files = r1_files.collect { f -> decompress_fastq ? f.name.replaceFirst(/\.fastq\.gz$/, '.fastq') : f.name }
    def r2_input_files = r2_files.collect { f -> decompress_fastq ? f.name.replaceFirst(/\.fastq\.gz$/, '.fastq') : f.name }
    def r1_decompress = decompress_fastq ? r1_files.collect { f ->
        def out = f.name.replaceFirst(/\.fastq\.gz$/, '.fastq')
        "${decompress_cmd} '${f.name}' > '${out}'"
    }.join("\n        ") : ""
    def r2_decompress = decompress_fastq ? r2_files.collect { f ->
        def out = f.name.replaceFirst(/\.fastq\.gz$/, '.fastq')
        "${decompress_cmd} '${f.name}' > '${out}'"
    }.join("\n        ") : ""
    def decompress_script = [r1_decompress, r2_decompress].findAll { s -> s }.join("\n        ")
    def r1_args = r1_input_files.collect { f -> "'${f}'" }.join(" ")
    def r2_args = r2_input_files.collect { f -> "'${f}'" }.join(" ")
    if (r2_files.isEmpty()) {
        """
        ${decompress_script}
        printf '%s\\n' ${r1_args} > read_list.txt

        ${params.tool_tcr_assembler_cmd} \\
            -1 read_list.txt \\
            -t ${task.cpus} \\
            -f '${ref_trust4_sequence}' \\
            --ref '${ref_trust4_imgt}' \\
            -o TRUST \\
            ${params.tool_tcr_assembler_extra_args}
        """
    } else {
        """
        ${decompress_script}
        printf '%s\\n' ${r1_args} > read1_list.txt
        printf '%s\\n' ${r2_args} > read2_list.txt

        ${params.tool_tcr_assembler_cmd} \\
            -1 read1_list.txt \\
            -2 read2_list.txt \\
            -t ${task.cpus} \\
            -f '${ref_trust4_sequence}' \\
            --ref '${ref_trust4_imgt}' \\
            -o TRUST \\
            ${params.tool_tcr_assembler_extra_args}
        """
    }
}

process TCR_PARSE_CDR3 {
    tag "trust4_parse"
    label "utility_light"

    input:
    path input_manifest
    path reports

    output:
    path "tcell_beta_cdr3.tsv", emit: tcr

    script:
    def inputs = reports instanceof List ? reports.join(" ") : reports
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/parse_trust4.py \\
        --manifest ${input_manifest} \\
        --reports ${inputs} \\
        --out tcell_beta_cdr3.tsv
    """
}
