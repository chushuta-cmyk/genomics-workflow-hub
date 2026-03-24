# Pipeline Overview

This repository contains reusable scripts for genomic analysis workflows.

## Modules

### core_gwas

Quality control, PCA, covariate preparation, and association testing.
Number of scripts: 109

### fine_mapping

Fine-mapping: regional plotting, LD heatmaps, annotation, and haplotype effect estimation.
Number of scripts: 2

### post_gwas

Post-GWAS analysis: significant SNP extraction, locus merging, cross-trait analysis, LD analysis, haplotype, coloc, and comparison.
Number of scripts: 20

### unknown

Miscellaneous scripts.
Number of scripts: 9

## Execution Order

1. **Data Preparation**: Use scripts in `core_gwas/qc/` for quality control and filtering.
2. **Covariate Calculation**: Use `core_gwas/covariates/` for PCA and kinship matrix generation.
3. **Association Testing**: Use `core_gwas/association/` for GWAS using GEMMA or other tools.
4. **Post-GWAS**: Extract significant SNPs, merge loci, perform cross-trait analysis, LD analysis, haplotype block detection, and coloc.
5. **Fine-mapping**: Regional association plots, LD heatmaps, candidate gene annotation, haplotype effect estimation.
6. **Workflow Orchestration**: Use `workflows/` for pipeline automation.

## Example Use Case: Soybean GWAS

These scripts were generalized from soybean GWAS projects. They can be adapted to other plant or animal genomes by modifying path placeholders and configuration parameters.
