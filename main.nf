nextflow.enable.dsl = 2

include { INPUT_VALIDATE_MANIFEST } from './modules/input'
include { READS_PREPARE; ALIGN_READS } from './modules/reads'
include { HLA_CALL; HLA_SUMMARIZE } from './modules/hla'
include { VARIANT_CALL; VARIANT_FILTER_QUALITY; VARIANT_SELECT_RECURRENT; VARIANT_TEST_ENRICHMENT } from './modules/variant_calling'
include { VARIANT_ANNOTATE; VARIANT_CONVERT_MAF; VARIANT_FILTER_SOMATIC } from './modules/variant_annotation'
include { SPLICE_FILTER_MANIFEST; SPLICE_BUILD_MATRIX; SPLICE_NORMALIZE_MATRIX; SPLICE_FILTER_BACKGROUND } from './modules/splicing'
include { PEPTIDE_GENERATE_MUTATION; PEPTIDE_GENERATE_SPLICING; PEPTIDE_COMPARE_IMMUNOPEPTIDOME; PEPTIDE_BYPASS_IMMUNOPEPTIDOME } from './modules/peptides'
include { MHC_PREDICT_MUTATION; MHC_PREDICT_SPLICING; MHC_MERGE_BATCHES; MHC_FILTER_CANDIDATES } from './modules/mhc'
include { TCR_ASSEMBLE; TCR_PARSE_CDR3 } from './modules/tcr'
include { EXPRESSION_QUANTIFY; EXPRESSION_NORMALIZE; CANDIDATE_EXPRESSION_SUPPORT } from './modules/expression'
include { PREDICT_IMMUNOGENICITY; PREDICT_TCR_RECOGNITION } from './modules/scoring'
include { CANDIDATE_INTEGRATE; CANDIDATE_TOPSIS } from './modules/results'


def requireParam(name) {
    if (!params[name]) {
        throw new IllegalArgumentException("Missing required parameter: --${name}")
    }
}


workflow {
    main:
    if (!params.input_manifest) {
        error "Missing required parameter: --input_manifest"
    }

    [
        "ref_star_index",
        "ref_genome_fasta",
        "ref_genome_gtf",
        "ref_protein_fasta",
        "ref_splice_gtex_h5ad",
        "ref_splice_exon_table",
        "ref_splice_transcript_table",
        "ref_splice_reference_gtf",
        "ref_splice_start_codon_table",
        "ref_vep_cache_dir",
        "ref_dbsnp_vcf",
        "ref_cosmic_vcf",
        "ref_trust4_sequence",
        "ref_trust4_imgt",
        "ref_immunogenicity_script",
    ].each { it -> requireParam(it) }

    if (params.enable_immunopeptidome_compare) {
        requireParam("ref_immunopeptidome_healthy_tsv")
        requireParam("ref_immunopeptidome_tumor_tsv")
    }

    manifest_input = channel.fromPath(params.input_manifest, checkIfExists: true)
    INPUT_VALIDATE_MANIFEST(manifest_input)
    manifest_valid = INPUT_VALIDATE_MANIFEST.out.manifest

    cells = manifest_valid
        .splitCsv(header: true, sep: "\t")
        .map { row ->
            def read_layout = (row.read_layout ?: (row.fastq_2 ? "paired" : "single")) as String
            def read_files = read_layout == "paired" ? [file(row.fastq_1), file(row.fastq_2)] : [file(row.fastq_1)]
            tuple(row.patient_id as String, row.cell_id as String, row.cell_type as String, read_layout, read_files)
        }

    READS_PREPARE(cells)
    reads = READS_PREPARE.out.reads

    ALIGN_READS(reads)
    aligned = ALIGN_READS.out.aligned
    alignment_bams = aligned.map { _patient_id, _cell_id, _cell_type, bam, bai, _sj, _log, _qc -> [bam, bai] }
    star_logs = aligned.map { _patient_id, _cell_id, _cell_type, _bam, _bai, _sj, log, _qc -> log }
    cell_qc = aligned.map { _patient_id, _cell_id, _cell_type, _bam, _bai, _sj, _log, qc -> qc }

    // Patient-level HLA typing: pool tumor/normal reads before one OptiType call per patient.
    hla_reads = reads.filter { _patient_id, _cell_id, cell_type, _read_layout, _read_files, _fastp_json ->
        cell_type == "tumor" || cell_type == "normal"
    }
    tumor_hla_input = hla_reads
        .map { patient_id, cell_id, _cell_type, read_layout, read_files, _fastp_json ->
            tuple(patient_id, cell_id, read_layout, read_files)
        }
        .groupTuple(by: 0)
        .map { patient_id, cell_ids, read_layouts, read_files ->
            tuple(patient_id, cell_ids, read_layouts, read_files.flatten())
        }
    HLA_CALL(tumor_hla_input)
    hla_files = HLA_CALL.out.hla.map { _patient_id, hla -> hla }
    HLA_SUMMARIZE(manifest_valid, hla_files.collect())

    // Recurrent tumor variants followed by VCF-based tumor/normal enrichment.
    variant_aligned = aligned.filter { _patient_id, _cell_id, cell_type, _bam, _bai, _sj, _log, _qc ->
        cell_type == "tumor" || cell_type == "normal"
    }
    VARIANT_CALL(variant_aligned)
    VARIANT_FILTER_QUALITY(VARIANT_CALL.out.vcf)
    scallele_vcfs = VARIANT_FILTER_QUALITY.out.vcf.map { _patient_id, _cell_id, _cell_type, vcf -> vcf }

    VARIANT_SELECT_RECURRENT(manifest_valid, scallele_vcfs.collect())
    VARIANT_TEST_ENRICHMENT(
        manifest_valid,
        VARIANT_SELECT_RECURRENT.out.tsv,
        scallele_vcfs.collect()
    )

    VARIANT_ANNOTATE(
        VARIANT_TEST_ENRICHMENT.out.vcf,
        file(params.ref_genome_fasta),
        file("${params.ref_genome_fasta}.fai"),
        file(params.ref_vep_cache_dir),
        file(params.ref_dbsnp_vcf),
        file(params.ref_dbsnp_vcf_index ?: "${params.ref_dbsnp_vcf}.tbi"),
        file(params.ref_cosmic_vcf),
        file(params.ref_cosmic_vcf_index ?: "${params.ref_cosmic_vcf}.tbi")
    )
    VARIANT_CONVERT_MAF(VARIANT_ANNOTATE.out.vcf)
    VARIANT_FILTER_SOMATIC(
        VARIANT_CONVERT_MAF.out.maf,
        VARIANT_TEST_ENRICHMENT.out.tsv
    )
    PEPTIDE_GENERATE_MUTATION(
        VARIANT_FILTER_SOMATIC.out.maf,
        file(params.ref_protein_fasta),
        file(params.ref_genome_fasta),
        file("${params.ref_genome_fasta}.fai"),
        file(params.ref_genome_gtf)
    )

    // Alternative-splicing candidates and junction-spanning peptides.
    SPLICE_FILTER_MANIFEST(manifest_valid)
    splice_manifest = SPLICE_FILTER_MANIFEST.out.manifest
    altanalyze_files = aligned
        .filter { _patient_id, _cell_id, cell_type, _bam, _bai, _sj, _log, _qc -> cell_type == "tumor" }
        .map { _patient_id, _cell_id, _cell_type, bam, bai, _sj, _log, _qc -> [bam, bai] }
        .flatten()
        .collect()
    SPLICE_BUILD_MATRIX(altanalyze_files)
    SPLICE_NORMALIZE_MATRIX(splice_manifest, SPLICE_BUILD_MATRIX.out.raw_matrix)
    SPLICE_FILTER_BACKGROUND(
        splice_manifest,
        SPLICE_NORMALIZE_MATRIX.out.matrix,
        file(params.ref_splice_gtex_h5ad),
        file(params.ref_splice_transcript_table),
        params.ref_splice_tcga_control_h5ad ? file(params.ref_splice_tcga_control_h5ad) : []
    )
    PEPTIDE_GENERATE_SPLICING(
        SPLICE_FILTER_BACKGROUND.out.junctions,
        file(params.ref_splice_exon_table),
        file(params.ref_genome_fasta),
        file("${params.ref_genome_fasta}.fai"),
        file(params.ref_splice_reference_gtf),
        file(params.ref_splice_start_codon_table)
    )

    // Optional immunopeptidome filtering/evidence. Same-length healthy BLAST hits with 100% identity are removed;
    // tumor-library matches are retained as evidence only.
    if (params.enable_immunopeptidome_compare) {
        PEPTIDE_COMPARE_IMMUNOPEPTIDOME(
            PEPTIDE_GENERATE_MUTATION.out.peptides,
            PEPTIDE_GENERATE_SPLICING.out.peptides,
            file(params.ref_immunopeptidome_healthy_tsv),
            file(params.ref_immunopeptidome_tumor_tsv)
        )
        mutation_peptides_for_mhc = PEPTIDE_COMPARE_IMMUNOPEPTIDOME.out.mutation
        splicing_peptides_for_mhc = PEPTIDE_COMPARE_IMMUNOPEPTIDOME.out.splicing
        immunopeptidome_evidence = PEPTIDE_COMPARE_IMMUNOPEPTIDOME.out.evidence
    } else {
        PEPTIDE_BYPASS_IMMUNOPEPTIDOME(
            PEPTIDE_GENERATE_MUTATION.out.peptides,
            PEPTIDE_GENERATE_SPLICING.out.peptides
        )
        mutation_peptides_for_mhc = PEPTIDE_BYPASS_IMMUNOPEPTIDOME.out.mutation
        splicing_peptides_for_mhc = PEPTIDE_BYPASS_IMMUNOPEPTIDOME.out.splicing
        immunopeptidome_evidence = PEPTIDE_BYPASS_IMMUNOPEPTIDOME.out.evidence
    }

    // Patient-specific netMHCpan; formal hard filter is rank_EL < 2 and Aff.nM < 150.
    MHC_PREDICT_MUTATION(mutation_peptides_for_mhc, HLA_SUMMARIZE.out.patient)
    MHC_PREDICT_SPLICING(splicing_peptides_for_mhc, HLA_SUMMARIZE.out.patient)
    MHC_MERGE_BATCHES(MHC_PREDICT_MUTATION.out.batches, MHC_PREDICT_SPLICING.out.batches)
    MHC_FILTER_CANDIDATES(
        MHC_PREDICT_MUTATION.out.predictions,
        MHC_PREDICT_SPLICING.out.predictions
    )

    // TCR beta reconstruction. parse_trust4.py keeps the highest-count valid TRB CDR3 per T cell.
    tcell_reads = reads.filter { _patient_id, _cell_id, cell_type, _read_layout, _read_files, _fastp_json ->
        cell_type == "tcell"
    }
    tcell_fastqs = tcell_reads
        .map { _patient_id, _cell_id, _cell_type, _read_layout, read_files, _fastp_json -> read_files }
        .collect()
    TCR_ASSEMBLE(
        tcell_fastqs,
        file(params.ref_trust4_sequence),
        file(params.ref_trust4_imgt)
    )
    trust_reports = TCR_ASSEMBLE.out.report
    trust4_raw = TCR_ASSEMBLE.out.raw.map { report, annot -> [report, annot] }
    TCR_PARSE_CDR3(manifest_valid, trust_reports.collect())

    // Tumor-cell gene expression from the coordinate-sorted BAMs.
    tumor_expression_input = aligned
        .filter { _patient_id, _cell_id, cell_type, _bam, _bai, _sj, _log, _qc -> cell_type == "tumor" }
        .map { patient_id, cell_id, _cell_type, bam, bai, _sj, _log, _qc ->
            tuple(patient_id, cell_id, bam, bai)
        }
    EXPRESSION_QUANTIFY(tumor_expression_input, file(params.ref_genome_gtf))
    EXPRESSION_NORMALIZE(EXPRESSION_QUANTIFY.out.counts.collect())
    CANDIDATE_EXPRESSION_SUPPORT(
        MHC_FILTER_CANDIDATES.out.pmhc,
        EXPRESSION_NORMALIZE.out.normalized
    )

    // Candidate-level scoring evidence.
    PREDICT_IMMUNOGENICITY(MHC_FILTER_CANDIDATES.out.pmhc, file(params.ref_immunogenicity_script))
    PREDICT_TCR_RECOGNITION(
        MHC_FILTER_CANDIDATES.out.pmhc,
        TCR_PARSE_CDR3.out.tcr
    )

    CANDIDATE_INTEGRATE(
        MHC_FILTER_CANDIDATES.out.pmhc,
        VARIANT_FILTER_SOMATIC.out.maf,
        PREDICT_IMMUNOGENICITY.out.netmhcpan,
        immunopeptidome_evidence,
        PREDICT_TCR_RECOGNITION.out.tabr,
        CANDIDATE_EXPRESSION_SUPPORT.out.support,
        INPUT_VALIDATE_MANIFEST.out.summary
    )
    CANDIDATE_TOPSIS(CANDIDATE_INTEGRATE.out.integrated)

    publish:
    manifest_validated = INPUT_VALIDATE_MANIFEST.out.manifest
    manifest_summary = INPUT_VALIDATE_MANIFEST.out.summary
    alignment_bams = alignment_bams
    star_logs = star_logs
    cell_qc = cell_qc

    hla_patient = HLA_SUMMARIZE.out.patient
    optitype_raw = HLA_CALL.out.raw.map { _patient_id, result, coverage_plot -> [result, coverage_plot] }

    scallele_vcfs = scallele_vcfs
    recurrent_tumor_variants_tsv = VARIANT_SELECT_RECURRENT.out.tsv
    recurrent_tumor_variants_vcf = VARIANT_SELECT_RECURRENT.out.vcf
    vep_annotated_vcf = VARIANT_ANNOTATE.out.vcf
    vcf2maf_raw_maf = VARIANT_CONVERT_MAF.out.maf
    somatic_mutation_candidates_maf = VARIANT_FILTER_SOMATIC.out.maf

    altanalyze_raw_junction_matrix = SPLICE_BUILD_MATRIX.out.raw_matrix
    altanalyze_groups = SPLICE_BUILD_MATRIX.out.groups
    altanalyze_comparisons = SPLICE_BUILD_MATRIX.out.comparisons
    altanalyze_junction_matrix = SPLICE_NORMALIZE_MATRIX.out.matrix
    splice_exon_table = file(params.ref_splice_exon_table)
    filtered_splice_junctions = SPLICE_FILTER_BACKGROUND.out.junctions

    mutation_peptides = PEPTIDE_GENERATE_MUTATION.out.peptides
    rejected_mutation_variants = PEPTIDE_GENERATE_MUTATION.out.rejected
    splicing_peptides = PEPTIDE_GENERATE_SPLICING.out.peptides
    rejected_splicing_junctions = PEPTIDE_GENERATE_SPLICING.out.rejected
    immunopeptidome_filtered_mutation_peptides = mutation_peptides_for_mhc
    immunopeptidome_filtered_splicing_peptides = splicing_peptides_for_mhc
    immunopeptidome_evidence_out = immunopeptidome_evidence

    mutation_netmhcpan = MHC_PREDICT_MUTATION.out.predictions
    mutation_netmhcpan_batches = MHC_PREDICT_MUTATION.out.batches
    splicing_netmhcpan = MHC_PREDICT_SPLICING.out.predictions
    splicing_netmhcpan_batches = MHC_PREDICT_SPLICING.out.batches
    netmhcpan_hla_batches = MHC_MERGE_BATCHES.out.batches
    filtered_candidate_pmhc = MHC_FILTER_CANDIDATES.out.pmhc
    netmhcpan_immunogenicity = PREDICT_IMMUNOGENICITY.out.netmhcpan

    trust4_raw = trust4_raw
    tcr_beta_cdr3 = TCR_PARSE_CDR3.out.tcr

    expression_featurecounts_raw = EXPRESSION_QUANTIFY.out.raw
    tumor_gene_counts = EXPRESSION_NORMALIZE.out.raw_matrix
    tumor_expression_normalized = EXPRESSION_NORMALIZE.out.normalized
    tumor_expression_matrix = EXPRESSION_NORMALIZE.out.normalized_matrix
    candidate_expression_support = CANDIDATE_EXPRESSION_SUPPORT.out.support

    tabr_bert_predictions = PREDICT_TCR_RECOGNITION.out.tabr
    integrated_candidate_neoantigens = CANDIDATE_INTEGRATE.out.integrated
    integrated_candidate_neoantigens_by_patient = CANDIDATE_INTEGRATE.out.by_patient
    pipeline_summary = CANDIDATE_INTEGRATE.out.summary
    final_candidate_neoantigens = CANDIDATE_TOPSIS.out.final_tsv
}


output {
    manifest_validated {
        path '01_qc/deliverables'
        mode params.output_publish_mode
        overwrite true
    }
    manifest_summary {
        path '01_qc/deliverables'
        mode params.output_publish_mode
        overwrite true
    }
    alignment_bams {
        path '01_qc/intermediate/alignments'
        mode params.output_publish_mode
        overwrite true
    }
    star_logs {
        path '01_qc/intermediate/star_logs'
        mode params.output_publish_mode
        overwrite true
    }
    cell_qc {
        path '01_qc/deliverables'
        mode params.output_publish_mode
        overwrite true
    }

    hla_patient {
        path '02_hla/deliverables'
        mode params.output_publish_mode
        overwrite true
    }
    optitype_raw {
        path '02_hla/intermediate/optitype'
        mode params.output_publish_mode
        overwrite true
    }

    scallele_vcfs {
        path '03_variants/intermediate/calls'
        mode params.output_publish_mode
        overwrite true
    }
    recurrent_tumor_variants_tsv {
        path '03_variants/intermediate/recurrent'
        mode params.output_publish_mode
        overwrite true
    }
    recurrent_tumor_variants_vcf {
        path '03_variants/intermediate/recurrent'
        mode params.output_publish_mode
        overwrite true
    }
    vcf2maf_raw_maf {
        path '03_variants/intermediate/maf'
        mode params.output_publish_mode
        overwrite true
    }
    vep_annotated_vcf {
        path '03_variants/intermediate/annotation'
        mode params.output_publish_mode
        overwrite true
    }
    somatic_mutation_candidates_maf {
        path '03_variants/deliverables'
        mode params.output_publish_mode
        overwrite true
    }

    altanalyze_raw_junction_matrix {
        path '04_splicing/intermediate/matrix'
        mode params.output_publish_mode
        overwrite true
    }
    altanalyze_groups {
        path '04_splicing/intermediate/matrix'
        mode params.output_publish_mode
        overwrite true
    }
    altanalyze_comparisons {
        path '04_splicing/intermediate/matrix'
        mode params.output_publish_mode
        overwrite true
    }
    altanalyze_junction_matrix {
        path '04_splicing/deliverables'
        mode params.output_publish_mode
        overwrite true
    }
    splice_exon_table {
        path '04_splicing/deliverables'
        mode params.output_publish_mode
        overwrite true
    }
    filtered_splice_junctions {
        path '04_splicing/deliverables'
        mode params.output_publish_mode
        overwrite true
    }

    mutation_peptides {
        path '05_peptides/intermediate'
        mode params.output_publish_mode
        overwrite true
    }
    rejected_mutation_variants {
        path '05_peptides/intermediate'
        mode params.output_publish_mode
        overwrite true
    }
    splicing_peptides {
        path '05_peptides/intermediate'
        mode params.output_publish_mode
        overwrite true
    }
    rejected_splicing_junctions {
        path '05_peptides/intermediate'
        mode params.output_publish_mode
        overwrite true
    }
    immunopeptidome_filtered_mutation_peptides {
        path '05_peptides/deliverables'
        mode params.output_publish_mode
        overwrite true
    }
    immunopeptidome_filtered_splicing_peptides {
        path '05_peptides/deliverables'
        mode params.output_publish_mode
        overwrite true
    }
    immunopeptidome_evidence_out {
        path '05_peptides/deliverables'
        mode params.output_publish_mode
        overwrite true
    }

    mutation_netmhcpan {
        path '06_mhc/deliverables'
        mode params.output_publish_mode
        overwrite true
    }
    mutation_netmhcpan_batches {
        path '06_mhc/intermediate'
        mode params.output_publish_mode
        overwrite true
    }
    splicing_netmhcpan {
        path '06_mhc/deliverables'
        mode params.output_publish_mode
        overwrite true
    }
    splicing_netmhcpan_batches {
        path '06_mhc/intermediate'
        mode params.output_publish_mode
        overwrite true
    }
    netmhcpan_hla_batches {
        path '06_mhc/intermediate'
        mode params.output_publish_mode
        overwrite true
    }
    filtered_candidate_pmhc {
        path '06_mhc/deliverables'
        mode params.output_publish_mode
        overwrite true
    }

    trust4_raw {
        path '07_tcr/intermediate/trust4'
        mode params.output_publish_mode
        overwrite true
    }
    tcr_beta_cdr3 {
        path '07_tcr/deliverables'
        mode params.output_publish_mode
        overwrite true
    }

    expression_featurecounts_raw {
        path '08_scoring/intermediate/expression'
        mode params.output_publish_mode
        overwrite true
    }
    tumor_gene_counts {
        path '08_scoring/deliverables/expression'
        mode params.output_publish_mode
        overwrite true
    }
    tumor_expression_normalized {
        path '08_scoring/deliverables/expression'
        mode params.output_publish_mode
        overwrite true
    }
    tumor_expression_matrix {
        path '08_scoring/deliverables/expression'
        mode params.output_publish_mode
        overwrite true
    }
    candidate_expression_support {
        path '08_scoring/deliverables'
        mode params.output_publish_mode
        overwrite true
    }
    netmhcpan_immunogenicity {
        path '08_scoring/deliverables'
        mode params.output_publish_mode
        overwrite true
    }
    tabr_bert_predictions {
        path '08_scoring/deliverables'
        mode params.output_publish_mode
        overwrite true
    }

    integrated_candidate_neoantigens {
        path '09_prioritization/intermediate'
        mode params.output_publish_mode
        overwrite true
    }
    integrated_candidate_neoantigens_by_patient {
        path '09_prioritization/intermediate'
        mode params.output_publish_mode
        overwrite true
    }
    pipeline_summary {
        path '09_prioritization/deliverables'
        mode params.output_publish_mode
        overwrite true
    }
    final_candidate_neoantigens {
        path '09_prioritization/deliverables'
        mode params.output_publish_mode
        overwrite true
    }
}
