process RESULT_INTEGRATE {
    tag "integrate"
    label "utility_light"

    input:
    path pmhc
    path iedb
    path pcidb
    path tabr
    path manifest_summary

    output:
    path "final_candidate_neoantigens.tsv", emit: final_tsv
    path "final_candidate_neoantigens_by_patient", emit: by_patient
    path "pipeline_summary.tsv", emit: summary

    script:
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/integrate_candidate_results.py \
        --pmhc ${pmhc} \
        --iedb ${iedb} \
        --pcidb ${pcidb} \
        --tabr ${tabr} \
        --manifest-summary ${manifest_summary} \
        --out final_candidate_neoantigens.tsv \
        --by-patient-dir final_candidate_neoantigens_by_patient \
        --summary pipeline_summary.tsv
    """
}

process CANDIDATE_INTEGRATE {
    tag "candidate_integrate"
    label "utility_light"

    input:
    path pmhc
    path mutation_maf
    path iedb
    path immunopeptidome
    path tabr
    path expression
    path manifest_summary

    output:
    path "integrated_candidate_neoantigens.tsv", emit: integrated
    path "integrated_candidate_neoantigens_by_patient", emit: by_patient
    path "pipeline_summary.tsv", emit: summary

    script:
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/integrate_candidate_results.py \
        --pmhc ${pmhc} \
        --mutation-maf ${mutation_maf} \
        --iedb ${iedb} \
        --immunopeptidome ${immunopeptidome} \
        --tabr ${tabr} \
        --expression ${expression} \
        --manifest-summary ${manifest_summary} \
        --out integrated_candidate_neoantigens.tsv \
        --by-patient-dir integrated_candidate_neoantigens_by_patient \
        --summary pipeline_summary.tsv
    """
}

process CANDIDATE_TOPSIS {
    tag "topsis"
    label "utility_light"

    input:
    path candidates

    output:
    path "final_candidate_neoantigens.tsv", emit: final_tsv

    script:
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/score_candidates_topsis.py \
        --candidates ${candidates} \
        --out final_candidate_neoantigens.tsv
    """
}
