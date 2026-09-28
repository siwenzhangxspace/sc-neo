# sc-neo

`sc-neo` is a Nextflow DSL2 workflow for neoantigen candidate discovery and prioritization from Smart-seq2 single-cell data. It orchestrates FASTQ files from tumor, normal, and T cells, combining evidence from RNA variants, aberrant splicing, patient-specific HLA-I, MHC-I presentation, TCRs, expression, and immunogenicity to rank candidates within each patient.

## Quick start

```bash
cp conf/reference_template.config conf/local.config
# Set the site-specific reference, tool, and container paths in conf/local.config
nextflow run main.nf \
  -profile apptainer \
  -c conf/local.config \
  --input_manifest manifest.tsv \
  --output_dir results \
  -resume
```

Immunopeptidome comparison is disabled by default. To enable it, set `enable_immunopeptidome_compare = true` in `conf/local.config` and provide the corresponding reference TSVs. See the [workflow usage guide](docs/usage.md) for details.

## Repository structure

- `main.nf`: workflow connections and published outputs.
- `nextflow.config`: default parameters, profiles, resources, and process container selectors.
- `modules/`: DSL2 processes organized by business domain.
- `bin/`: Python helper CLIs called by processes.
- `conf/reference_template.config`: template for site-specific configuration.
- `containers/`: container build files and wrappers.

## Local checks

```bash
python3 -m py_compile bin/*.py
python3 -m unittest discover -s tests
nextflow run main.nf --help
```

A full run also requires Java 11+, Nextflow, site-specific reference resources, and external tools such as STAR, samtools, OptiType, scAllele, VEP for annotation, mafsmith for VCF-to-MAF conversion, netMHCpan, TRUST4, TABR-BERT, and BLAST+.
