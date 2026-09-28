process VARIANT_ANNOTATE {
    tag "vep_annotate"
    label "annotation_heavy"

    input:
    path enriched_vcf
    path ref_genome_fasta
    path genome_fasta_fai
    path ref_vep_cache_dir
    path ref_dbsnp_vcf
    path ref_dbsnp_vcf_index
    path ref_cosmic_vcf
    path ref_cosmic_vcf_index

    output:
    path "annotated.vcf", emit: vcf

    script:
    """
    if grep -v '^#' ${enriched_vcf} | grep -q .; then
        ${params.tool_variant_annotation_cmd} \
            --input_file ${enriched_vcf} \
            --output_file annotated.vcf \
            --format vcf \
            --vcf \
            --offline \
            --cache \
            --dir_cache ${ref_vep_cache_dir} \
            --species ${params.ref_vep_species} \
            --assembly ${params.ref_vep_assembly} \
            --cache_version ${params.ref_vep_cache_version} \
            --fasta ${ref_genome_fasta} \
            --everything \
            --force_overwrite \
            --no_stats \
            --fork ${task.cpus} \
            --custom file=${ref_dbsnp_vcf},short_name=dbSNP,format=vcf,type=exact,coords=0,fields=RS%SAO%CAF%TOPMED \
            --custom file=${ref_cosmic_vcf},short_name=COSMIC,format=vcf,type=exact,coords=0,fields=GENE%AA%SO_TERM%TARGETED_SCREEN_SAMPLE_COUNT%GENOME_SCREEN_SAMPLE_COUNT%TIER \
            ${params.tool_variant_annotation_extra_args}
    else
        cp ${enriched_vcf} annotated.vcf
    fi
    """
}

process VARIANT_CONVERT_MAF {
    tag "mafsmith_vcf2maf"
    label "io_light"

    input:
    path annotated_vcf

    output:
    path "raw.maf", emit: maf

    script:
    """
    if grep -v '^#' ${annotated_vcf} | grep -q .; then
        ${params.tool_maf_converter_cmd} vcf2maf \
            -i ${annotated_vcf} \
            -o raw.maf \
            --skip-annotation \
            --tumor-id "sample_tumor" \
            --retain-ann "dbSNP,dbSNP_RS,dbSNP_SAO,dbSNP_CAF,dbSNP_TOPMED,COSMIC,COSMIC_GENE,COSMIC_AA,COSMIC_SO_TERM,COSMIC_TARGETED_SCREEN_SAMPLE_COUNT,COSMIC_GENOME_SCREEN_SAMPLE_COUNT,COSMIC_TIER,gnomADg_AF,gnomADe_AF,AF,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF" \
            ${params.tool_maf_converter_extra_args}
    else
        printf "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tVariant_Classification\tConsequence\tTranscript_ID\tFeature\tProtein_position\tAmino_acids\tHGVSp_Short\tHGVSc\n" > raw.maf
    fi
    """
}

process VARIANT_FILTER_SOMATIC {
    tag "filter_maf"
    label "io_light"

    input:
    path raw_maf
    path trace_tsv

    output:
    path "somatic_mutation_candidates.maf", emit: maf

    script:
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/filter_maf.py \
        --maf ${raw_maf} \
        --trace-tsv ${trace_tsv} \
        --population-af-max ${params.threshold_population_af_max} \
        --out somatic_mutation_candidates.pre_status.maf

    ${params.tool_python_cmd} ${workflow.projectDir}/bin/normalize_variant_status.py \
        --maf somatic_mutation_candidates.pre_status.maf \
        --out somatic_mutation_candidates.maf
    """
}
