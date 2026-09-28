process PREDICT_IMMUNOGENICITY {
    tag "iedb"
    label "prediction_medium"

    input:
    path pmhc
    path immunogenicity_script

    output:
    path "netmhcpan_immunogenicity.tsv", emit: netmhcpan

    script:
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/predict_immunogenicity.py \
        --pmhc ${pmhc} \
        --script ${immunogenicity_script} \
        --python ${params.tool_python_cmd} \
        --out netmhcpan_immunogenicity.tsv
    """
}


process PREDICT_MATCH_KNOWN_EPITOPES {
    tag "pcidb"
    label "annotation_heavy"

    input:
    path pmhc
    path pcidb

    output:
    path "pcidb_blast.tsv", emit: pcidb

    script:
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/pcidb_blast.py \
        --pmhc ${pmhc} \
        --pcidb ${pcidb} \
        --makeblastdb-cmd ${params.tool_makeblastdb_cmd} \
        --blastp-cmd ${params.tool_blastp_cmd} \
        --extra-args "${params.tool_blastp_extra_args}" \
        --lengths ${params.peptide_class_i_lengths} \
        --out pcidb_blast.tsv
    """
}

process PREDICT_TCR_RECOGNITION {
    tag "tabr_bert"
    label "annotation_heavy"

    input:
    path pmhc
    path tcr

    output:
    path "tabr_bert_predictions.tsv", emit: tabr

    script:
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/predict_tcr_recognition.py \
        --pmhc ${pmhc} \
        --tcr ${tcr} \
        --cmd "${params.tool_tcr_recognition_cmd}" \
        --extra-args "${params.tool_tcr_recognition_extra_args}" \
        --out tabr_bert_predictions.tsv
    """
}
