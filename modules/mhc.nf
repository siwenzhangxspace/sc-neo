process MHC_PREDICT_MUTATION {
    tag "netmhcpan_mutation"
    label "prediction_medium"

    input:
    path peptides
    path hla

    output:
    path "mutation_netmhcpan.tsv", emit: predictions
    path "mutation_netmhcpan_hla_batches.tsv", emit: batches

    script:
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/predict_mhc_batches.py \
        --peptides ${peptides} \
        --hla ${hla} \
        --source-type mutation \
        --cmd ${params.tool_mhc_predictor_cmd} \
        --extra-args "${params.tool_mhc_predictor_extra_args}" \
        --max-a-chars ${params.threshold_mhc_max_allele_chars} \
        --out-predictions mutation_netmhcpan.tsv \
        --out-batches mutation_netmhcpan_hla_batches.tsv
    """
}

process MHC_PREDICT_SPLICING {
    tag "netmhcpan_splicing"
    label "prediction_medium"

    input:
    path peptides
    path hla

    output:
    path "splicing_netmhcpan.tsv", emit: predictions
    path "splicing_netmhcpan_hla_batches.tsv", emit: batches

    script:
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/predict_mhc_batches.py \
        --peptides ${peptides} \
        --hla ${hla} \
        --source-type splicing \
        --cmd ${params.tool_mhc_predictor_cmd} \
        --extra-args "${params.tool_mhc_predictor_extra_args}" \
        --max-a-chars ${params.threshold_mhc_max_allele_chars} \
        --out-predictions splicing_netmhcpan.tsv \
        --out-batches splicing_netmhcpan_hla_batches.tsv
    """
}

process MHC_MERGE_BATCHES {
    tag "hla_batches"
    label "utility_light"

    input:
    path mutation_batches
    path splicing_batches

    output:
    path "netmhcpan_hla_batches.tsv", emit: batches

    script:
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/merge_tsv.py --out netmhcpan_hla_batches.tsv ${mutation_batches} ${splicing_batches}
    """
}

process MHC_FILTER_CANDIDATES {
    tag "filter_pmhc"
    label "utility_light"

    input:
    path mutation_predictions
    path splicing_predictions

    output:
    path "filtered_candidate_pmHC.tsv", emit: pmhc

    script:
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/filter_pmhc_candidates.py \
        --mutation ${mutation_predictions} \
        --splicing ${splicing_predictions} \
        --affinity-max ${params.threshold_mhc_affinity_max} \
        --rank-el-max ${params.threshold_mhc_rank_el_max} \
        --ignore-bind-level \
        --out filtered_candidate_pmHC.tsv
    """
}
