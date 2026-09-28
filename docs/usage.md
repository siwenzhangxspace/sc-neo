# Workflow Usage Guide

This guide explains how to prepare, configure, and run `sc-neo`. The workflow is maintained for Nextflow `26.04.2` and supports the Apptainer, Singularity, Docker, and conda profiles.

## 1. Requirements

A complete run requires Java 11+, Nextflow, the appropriate containers or host tools, and all reference files. External tools include STAR, samtools, OptiType, scAllele, VEP for annotation, mafsmith for VCF-to-MAF conversion, netMHCpan, TRUST4, TABR-BERT, and BLAST+.

Create a local configuration from the template:

```bash
cp conf/reference_template.config conf/local.config
```

Set the local reference, tool, and container paths in `conf/local.config`. Do not commit this file, patient manifests, FASTQ/BAM files, `results/`, or `work/`.

## 2. Input manifest

The manifest is a tab-separated file with these five required columns:

```tsv
patient_id	cell_id	cell_type	fastq_1	fastq_2
P001	TumorCell_001	tumor	/path/to/TumorCell_001_R1.fastq.gz	/path/to/TumorCell_001_R2.fastq.gz
P001	NormalCell_001	normal	/path/to/NormalCell_001_R1.fastq.gz	/path/to/NormalCell_001_R2.fastq.gz
P001	TCell_001	tcell	/path/to/TCell_001_R1.fastq.gz	/path/to/TCell_001_R2.fastq.gz
```

- `cell_type` must be `tumor`, `normal`, or `tcell`.
- `cell_id` must be globally unique and must not contain `.`.
- `fastq_1` is required. By default, `allow_single_end_reads=true`, so `fastq_2` may be empty; otherwise, both FASTQ files are required.
- Each patient must have at least `threshold_min_cells_per_type` cells of each type. The default is `3`.
- The validator derives `read_layout` from whether `fastq_2` is empty. It ignores any `read_layout` column in the manifest.
- HLA calling requires tumor and normal reads from the same patient to use the same layout. TRUST4 batches also require all T-cell reads to use the same layout.

## 3. Required references

At startup, `main.nf` requires the following parameters to be non-empty:

| Parameter | Purpose |
| --- | --- |
| `ref_star_index` | STAR genome index. |
| `ref_genome_fasta` | Genome FASTA; the `.fai` file at the same path is also used. |
| `ref_genome_gtf` | GTF used by STAR/featureCounts. |
| `ref_protein_fasta` | Protein reference for mutation peptides; supports FASTA or a three-column TSV with gene, alias, and sequence. |
| `ref_splice_gtex_h5ad` | GTEx junction-control H5AD. |
| `ref_splice_exon_table` | AltAnalyze subexon coordinate table. |
| `ref_splice_transcript_table` | Transcript-exon table used to annotate documented junctions. |
| `ref_splice_reference_gtf` | GTF used to infer splicing translation phase. |
| `ref_splice_start_codon_table` | Start-codon table for splicing translation. |
| `ref_vep_cache_dir` | Offline VEP cache. |
| `ref_dbsnp_vcf`, `ref_cosmic_vcf` | VEP custom VCFs; each index defaults to the `.tbi` file next to its VCF. |
| `ref_trust4_sequence`, `ref_trust4_imgt` | TRUST4 references. |
| `ref_immunogenicity_script` | External immunogenicity prediction script. |

`ref_splice_tcga_control_h5ad` is an optional additional control. Immunopeptidome comparison is disabled by default (`enable_immunopeptidome_compare=false`). When enabled, it also requires `ref_immunopeptidome_healthy_tsv` and `ref_immunopeptidome_tumor_tsv`. Both TSVs need a peptide/sequence column; HLA and disease/tissue columns are optional metadata. All splicing references must use the same release as the AltAnalyze UID.

## 4. Run the workflow

Apptainer is the recommended profile:

```bash
nextflow run main.nf \
  -profile apptainer \
  -c conf/local.config \
  --input_manifest manifest.tsv \
  --output_dir results \
  -work-dir work \
  -resume
```

You can replace the profile with `singularity`, `docker`, or `conda`. To resume after an interruption, keep the same inputs and configuration and use `-resume` to reuse completed tasks.

## 5. Common parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `output_dir` | `results` | Root directory for published workflow outputs. |
| `output_publish_mode` | `copy` | Output publishing mode. |
| `enable_read_trimming` | `false` | Whether to run fastp. |
| `enable_immunopeptidome_compare` | `false` | Whether to compare against healthy and tumor immunopeptidome references. |
| `allow_single_end_reads` | `true` | Whether an empty `fastq_2` is allowed. |
| `threshold_min_cells_per_type` | `3` | Minimum number of cells of each type per patient. |
| `threshold_min_tumor_variant_cells` | `3` | Minimum number of tumor cells supporting a recurrent variant. |
| `threshold_variant_enrichment_fdr` | `0.1` | FDR threshold for variant enrichment. |
| `threshold_splice_gtex_n_max` | `3` | Upper bound for the control mean; the comparison is strictly less than this value. |
| `threshold_splice_gtex_t_min` | `20` | Lower bound for tumor count minus control mean; the comparison is strictly greater than this value. |
| `threshold_splice_min_reads` | `3` | Minimum number of cells supporting a junction. |
| `threshold_mhc_affinity_max` | `150` | `Aff.nM` upper bound; the comparison is strictly less than this value. |
| `threshold_mhc_rank_el_max` | `2` | `rank_el` upper bound; the comparison is strictly less than this value. |
| `enable_strict_splice_phase` | `false` | Whether to keep only splicing translation frames with phase evidence. |
| `enable_splice_exclude_documented_junctions` | `false` | Whether to remove junctions documented in the transcript table. |

Configure tool commands and extra arguments with `tool_*_cmd` and `tool_*_extra_args`. See `nextflow.config` for all defaults.

## 6. Container profiles

Apptainer and Singularity enable `autoMounts` by default. A non-empty `runtime_container` overrides the process containers globally; otherwise, each process uses the image selected by `process.withName` under `runtime_container_dir`. `runtime_container_dir` is unset by default, so set the image directory in the local configuration when using a container profile.

| Image | Processes |
| --- | --- |
| `sc-neo-core.sif` | Input validation, HLA summary, variant and splicing filters, peptide generation, MHC, expression, scoring, and integration helpers. |
| `sc-neo-align.sif` | `READS_PREPARE`, `ALIGN_READS`. |
| `sc-neo-hla.sif` | `HLA_CALL`. |
| `s5_latest.sif` | `VARIANT_CALL`. |
| `vep.sif` | `VARIANT_ANNOTATE`. |
| `sc-neo-mafsmith.sif` | `VARIANT_CONVERT_MAF`. |
| `sc-neo-altanalyze.sif` | `SPLICE_BUILD_MATRIX`. |
| `trust4-plus.sif` | `TCR_ASSEMBLE`. |
| `tabr-bert.sif` | `PREDICT_TCR_RECOGNITION`. |

For example, build the core and VEP images with:

```bash
DOCKER_BUILDKIT=1 docker build \
  -f containers/sc-neo-core.Dockerfile \
  -t sc-neo-core:local .
apptainer build sc-neo-core.sif docker-daemon://sc-neo-core:local
apptainer build vep.sif docker://ensemblorg/ensembl-vep:release_112.0
```

The VEP cache is not included in the image. It must match `ref_vep_species`, `ref_vep_assembly`, and `ref_vep_cache_version`.

## 7. Mount paths

Containers must be able to access the FASTQ files in the manifest, all `ref_*` paths, `runtime_container_dir`, `output_dir`, and the work directory. If data is stored outside automatically mounted locations, add site-specific bind mounts in the local configuration:

```groovy
apptainer {
    runOptions = '--bind /path/to/references,/path/to/results,/path/to/scratch'
}

singularity {
    runOptions = '--bind /path/to/references,/path/to/results,/path/to/scratch'
}
```

Replace the placeholder directories with paths available at your site. Prefer using the same absolute paths inside and outside the container in the manifest and configuration.

## 8. Resource configuration

| Label | Default resources | Typical workload |
| --- | --- | --- |
| `utility_light` | 1 CPU, 4 GB, `maxForks=64` | Validation, parsing, merging, and integration. |
| `io_light` | 2 CPU, 8 GB, `maxForks=48` | Read preparation, small conversions, and peptide generation. |
| `cpu_medium` | 8 CPU, 32 GB, `maxForks=32` | Medium CPU workloads. |
| `cpu_threaded` | 16 CPU, 48 GB, `maxForks=20` | Tools with explicit multithreading. |
| `mem_high_threaded` | 16 CPU, 96 GB, `maxForks=9` | STAR and AltAnalyze. |
| `mem_extreme_batch` | 96 CPU, 256 GB, `maxForks=1` | Aggregated TRUST4 batches. |
| `annotation_heavy` | 12 CPU, 96 GB, `maxForks=8` | VEP, models, and database predictions. |
| `prediction_medium` | 8 CPU, 32 GB, `maxForks=24` | MHC, IEDB, and similar prediction wrappers. |

For cluster deployments, override the executor, resources, or selectors in `conf/local.config`. The AltAnalyze wrapper creates an isolated runtime directory in the task work directory and limits multiprocessing to `task.cpus`; rebuild the image after modifying the wrapper.

## 9. Results

All outputs are published under `output_dir`. Each stage directory contains `deliverables/` for analysis and `intermediate/` for raw tool outputs and traceable intermediate files.

The final candidates are written to:

```text
09_prioritization/deliverables/final_candidate_neoantigens.tsv
```
