process READS_PREPARE {
    tag "${cell_id}"
    label "io_light"

    input:
    tuple val(patient_id), val(cell_id), val(cell_type), val(read_layout), path(read_files)

    output:
    tuple val(patient_id), val(cell_id), val(cell_type), val(read_layout), path("${cell_id}.R*.clean.fastq.gz"), path("${cell_id}.fastp.json"), emit: reads

    script:
    def staged_reads = read_files instanceof List ? read_files : [read_files]
    def fastq_1 = staged_reads[0]
    def fastq_2 = read_layout == "paired" ? staged_reads[1] : null
    if (params.enable_read_trimming && read_layout == "single") {
        """
        ${params.tool_fastqc_cmd} ${fastq_1}
        ${params.tool_fastp_cmd} \\
            -i ${fastq_1} \\
            -o ${cell_id}.R1.clean.fastq.gz \\
            -j ${cell_id}.fastp.json \\
            -h ${cell_id}.fastp.html \\
            ${params.tool_fastp_extra_args}
        """
    } else if (params.enable_read_trimming) {
        """
        ${params.tool_fastqc_cmd} ${fastq_1} ${fastq_2}
        ${params.tool_fastp_cmd} \\
            -i ${fastq_1} \\
            -I ${fastq_2} \\
            -o ${cell_id}.R1.clean.fastq.gz \\
            -O ${cell_id}.R2.clean.fastq.gz \\
            -j ${cell_id}.fastp.json \\
            -h ${cell_id}.fastp.html \\
            ${params.tool_fastp_extra_args}
        """
    } else if (read_layout == "single") {
        """
        ${params.tool_fastqc_cmd} ${fastq_1}
        if [[ "${fastq_1}" == *.gz ]]; then ln -sf ${fastq_1} ${cell_id}.R1.clean.fastq.gz; else gzip -c ${fastq_1} > ${cell_id}.R1.clean.fastq.gz; fi
        printf '{"summary":{"before_filtering":{},"after_filtering":{}}}\\n' > ${cell_id}.fastp.json
        """
    } else {
        """
        ${params.tool_fastqc_cmd} ${fastq_1} ${fastq_2}
        if [[ "${fastq_1}" == *.gz ]]; then ln -sf ${fastq_1} ${cell_id}.R1.clean.fastq.gz; else gzip -c ${fastq_1} > ${cell_id}.R1.clean.fastq.gz; fi
        if [[ "${fastq_2}" == *.gz ]]; then ln -sf ${fastq_2} ${cell_id}.R2.clean.fastq.gz; else gzip -c ${fastq_2} > ${cell_id}.R2.clean.fastq.gz; fi
        printf '{"summary":{"before_filtering":{},"after_filtering":{}}}\\n' > ${cell_id}.fastp.json
        """
    }
}

process ALIGN_READS {
    tag "${cell_id}"
    label "mem_high_threaded"

    input:
    tuple val(patient_id), val(cell_id), val(cell_type), val(read_layout), path(read_files), path(fastp_json)

    output:
    tuple val(patient_id), val(cell_id), val(cell_type), path("${cell_id}.sorted.bam"), path("${cell_id}.sorted.bam.bai"), path("${cell_id}.SJ.out.tab"), path("${cell_id}.star.Log.final.out"), path("${cell_id}.cell_qc.tsv"), emit: aligned

    script:
    def staged_reads = read_files instanceof List ? read_files : [read_files]
    def fastq_1 = staged_reads[0]
    def fastq_2 = read_layout == "paired" ? staged_reads[1] : null
    def readFilesCommand = staged_reads.every {it -> it.name.endsWith('.gz') } ? '--readFilesCommand zcat' : ''
    if (read_layout == "single") {
        """
        ${params.tool_star_cmd} \\
            --genomeDir ${params.ref_star_index} \\
            --readFilesIn ${fastq_1} \\
            ${readFilesCommand} \\
            --runThreadN ${task.cpus} \\
            --outSAMtype BAM SortedByCoordinate \\
            --outFileNamePrefix ${cell_id}. \\
            ${params.tool_star_extra_args}
        mv ${cell_id}.Aligned.sortedByCoord.out.bam ${cell_id}.sorted.bam
        ${params.tool_samtools_cmd} index ${cell_id}.sorted.bam
        cp ${cell_id}.Log.final.out ${cell_id}.star.Log.final.out
        ${params.tool_python_cmd} ${workflow.projectDir}/bin/summarize_alignment_qc.py \\
            --patient-id ${patient_id} \\
            --cell-id ${cell_id} \\
            --cell-type ${cell_type} \\
            --star-log ${cell_id}.star.Log.final.out \\
            --fastp-json ${fastp_json} \\
            --out ${cell_id}.cell_qc.tsv
        """
    } else {
        """
        ${params.tool_star_cmd} \\
            --genomeDir ${params.ref_star_index} \\
            --readFilesIn ${fastq_1} ${fastq_2} \\
            ${readFilesCommand} \\
            --runThreadN ${task.cpus} \\
            --outSAMtype BAM SortedByCoordinate \\
            --outFileNamePrefix ${cell_id}. \\
            ${params.tool_star_extra_args}
        mv ${cell_id}.Aligned.sortedByCoord.out.bam ${cell_id}.sorted.bam
        ${params.tool_samtools_cmd} index ${cell_id}.sorted.bam
        cp ${cell_id}.Log.final.out ${cell_id}.star.Log.final.out
        ${params.tool_python_cmd} ${workflow.projectDir}/bin/summarize_alignment_qc.py \\
            --patient-id ${patient_id} \\
            --cell-id ${cell_id} \\
            --cell-type ${cell_type} \\
            --star-log ${cell_id}.star.Log.final.out \\
            --fastp-json ${fastp_json} \\
            --out ${cell_id}.cell_qc.tsv
        """
    }
}
