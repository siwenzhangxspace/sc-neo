process PEPTIDE_GENERATE_MUTATION {
    tag "mutation_peptides"
    label "io_light"

    input:
    path maf

    output:
    path "mutation_peptides.tsv", emit: peptides

    script:
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/generate_mutation_peptides.py \
        --maf ${maf} \
        --protein-fasta ${params.ref_protein_fasta} \
        --lengths ${params.peptide_class_i_lengths} \
        --out mutation_peptides.tsv
    """
}

process PEPTIDE_GENERATE_SPLICING {
    tag "splicing_peptides"
    label "cpu_medium"

    input:
    path junctions
    path exon_table
    path ref_genome_fasta
    path genome_fasta_fai
    path reference_gtf
    path start_codon_table

    output:
    path "splicing_peptides.tsv", emit: peptides
    path "rejected_splicing_junctions.tsv", emit: rejected

    script:
    def strictPhase = params.enable_strict_splice_phase ? "--strict-phase" : ""
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/generate_splicing_peptides.py \
        --junctions ${junctions} \
        --exon-table ${exon_table} \
        --genome-fasta ${ref_genome_fasta} \
        --reference-gtf ${reference_gtf} \
        --start-codon-table ${start_codon_table} \
        ${strictPhase} \
        --lengths ${params.peptide_class_i_lengths} \
        --out splicing_peptides.tsv \
        --rejected-out rejected_splicing_junctions.tsv
    """
}

process PEPTIDE_COMPARE_IMMUNOPEPTIDOME {
    tag "immunopeptidome"
    label "cpu_medium"

    input:
    path mutation_peptides
    path splicing_peptides
    path healthy_reference
    path tumor_reference

    output:
    path "mutation_peptides.immunopeptidome_filtered.tsv", emit: mutation
    path "splicing_peptides.immunopeptidome_filtered.tsv", emit: splicing
    path "immunopeptidome_evidence.tsv", emit: evidence

    script:
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/compare_immunopeptidome.py \
        --mutation ${mutation_peptides} \
        --splicing ${splicing_peptides} \
        --healthy ${healthy_reference} \
        --tumor ${tumor_reference} \
        --makeblastdb-cmd ${params.tool_makeblastdb_cmd} \
        --blastp-cmd ${params.tool_blastp_cmd} \
        --extra-args "${params.tool_blastp_extra_args}" \
        --lengths ${params.peptide_class_i_lengths} \
        --out-mutation mutation_peptides.immunopeptidome_filtered.tsv \
        --out-splicing splicing_peptides.immunopeptidome_filtered.tsv \
        --out-evidence immunopeptidome_evidence.tsv
    """
}

process PEPTIDE_BYPASS_IMMUNOPEPTIDOME {
    tag "immunopeptidome_disabled"
    label "utility_light"

    input:
    path mutation_peptides
    path splicing_peptides

    output:
    path "mutation_peptides.immunopeptidome_filtered.tsv", emit: mutation
    path "splicing_peptides.immunopeptidome_filtered.tsv", emit: splicing
    path "immunopeptidome_evidence.tsv", emit: evidence

    script:
    """
    ${params.tool_python_cmd} ${workflow.projectDir}/bin/compare_immunopeptidome.py \
        --mutation ${mutation_peptides} \
        --splicing ${splicing_peptides} \
        --disabled \
        --lengths ${params.peptide_class_i_lengths} \
        --out-mutation mutation_peptides.immunopeptidome_filtered.tsv \
        --out-splicing splicing_peptides.immunopeptidome_filtered.tsv \
        --out-evidence immunopeptidome_evidence.tsv
    """
}
