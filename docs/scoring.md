# scNeo evidence scoring

The executable specification is `bin/score_candidates_topsis.py`. Candidate identity remains patient, gene, HLA, peptide length, peptide, event and protein/junction change (key7). Scoring does not change the upstream MHC hard filter (`rank_el < 2` and `Aff.nM < 150`).

## Evidence and transforms

| Module | Transform to [0,1] |
| --- | --- |
| Presentation P | Mean of 150/(150 + affinity_nM) and 0.5/(0.5 + rank_EL). |
| Tumor immunopeptidome C | 0.5 + 0.5 × support; support = identity/100 × query_coverage/(1 + E-value). Full-sequence exact tumor match sets support to 1. Completed search with no hit has support 0. |
| TCR T | Maximum valid TABR-BERT rank across patient-matched tumor CD8 TCRs, without cohort rescaling. |
| RNA R | min(1, log1p(mean_malignant_TPM) / log1p(198.738794605166)). |
| IEDB I | clip((official immunogenicity score + 1)/2, 0, 1). |

Healthy-library filtering removes only full-sequence exact matches. A 100% identity local HSP covering part of a peptide is insufficient. The complete tumor reference is searched with BLASTP-short without partitioning by peptide length. Coverage uses query coordinates, not alignment length including gaps. E-value, bit score, identity, alignment length and reference ID provide deterministic best-HSP ordering. Preserve the same reference release and records, including duplicates, when reproducing archived E-values.

IEDB uses the configured official script. Supported nine-residue peptide/HLA pairs use the script's allele-specific anchor mask; other lengths and unsupported alleles use its default mask. Integration joins on both peptide and normalized HLA, and records the mask rule.

## RNA input and missing values

From FASTQ, featureCounts gene counts and gene lengths yield per-cell TPM:
`TPM_g = (count_g / length_g) / sum(count / length) * 1e6`.
The RNA module averages source-gene TPM over **all sampled malignant cells of the patient**, including measured zeros; it does not restrict to variant/junction-supporting cells or take a maximum.

To reproduce an existing quantification, provide `--input_gene_tpm measured_tpm.tsv`:
`patient_id, cell_id, gene_name (or gene_id), tpm` as tab-separated columns. Include explicit zero rows for measured zeros and every sampled malignant cell for each scored source gene. Use unlogged TPM, not Seurat normalized data or log-TPM. Only malignant cells from the manifest belong in this file. If the source matrix is log-TPM, invert its documented log transform before import. The pipeline does not infer that transform.

The optional TPM input retains the original quantification and skips featureCounts. Raw-count and old normalized-expression companion matrices are blank in this mode; use the TPM matrix. A fresh featureCounts quantification can differ from TPM produced by another quantifier; identical formulas alone do not guarantee identical source TPM.

Missing gene/cell measurements, absent source-gene mapping and zero-total libraries are missing evidence, not zero expression. True measured zero is scored as RNA=0. The configured reference is shared across single-cell and external TPM inputs; do not recompute it independently within each cohort.

## Fixed weights and TOPSIS

| Evidence available | P | C | T | R | I |
| --- | --- | --- | --- | --- | --- |
| Five modules | 0.25 | 0.20 | 0.25 | 0.15 | 0.15 |
| TCR unavailable | 0.40 | 0.30 | 0 | 0.15 | 0.15 |
| PCI disabled, TCR available (optional mode) | 0.35 | 0 | 0.35 | 0.15 | 0.15 |

A missing TCR repertoire selects the predefined four-module rule; failed or incomplete predictions must be rerun. Disabling PCI is an explicit optional analysis, not a failed search interpreted as no hit. Missing both TCR and PCI cannot be scored. At least four modules and all active inputs are required; missing P, R or I leaves the row unscored.

With fixed negative and positive ideals 0 and 1,
`d0 = sqrt(sum((w*x)^2))`,
`d1 = sqrt(sum((w*(1-x))^2))`,
`TOPSIS_score = 100*d0/(d0+d1)`.
There is no cohort min-max or vector normalization, no completeness multiplier, and no final min-max scaling. Weights are applied before Euclidean distances. The diagnostic `evidence_completeness` is not a multiplier; the deprecated `topsis_completeness_adjusted_score` is left blank.

## Selection and interpretation

Rank all evaluable key7 candidates in the submitted cohort by descending score, resolving ties by the lexical joined key7 string. Select exactly ceil(0.10 × number evaluable). This is a cohort-wide allocation, not 10% within each patient. Outputs include `cohort_rank`, `top_decile`, `score_eligible`, `weight_policy` and transformed module values. `priority_class=High` now means the selected top decile; other evaluable rows are `Other`, and ineligible rows are `Not scored`. The historical ≥70 rule is no longer the selection rule.

Scores are deterministic for fixed evidence; cohort composition affects the selection budget and boundary, not an individual's score. For the external benchmark, the evaluation unit is patient–peptide: average scores over that unit's HLA records before ranking and computing metrics. Do not equate top-decile pHLA rows with top-decile patient–peptide units.

These fixed weights are the current analysis configuration. This code does not establish that the configuration was prespecified, or turn development data into independent validation.

## Reproducibility

Run `python bin/test_scoring.py`. After upgrading an existing workflow, use new work/output directories for the first run. Old intermediate tables lack the TPM scale, search-completion/coverage fields or HLA-aware IEDB scores and must be regenerated. Legacy missing fields fail explicitly where interpreting them would change the model. Keep the same reference files, official IEDB script, netMHCpan/TABR models and source TPM when comparing scores across versions.
