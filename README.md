# Genomics Workflow Hub

A curated collection of reusable GWAS, post-GWAS, fine-mapping, and workflow scripts generalized from soybean genomics projects.

## Purpose

This repository provides portable, well-documented scripts for genomic analysis workflows, with hard-coded paths replaced by generic placeholders to facilitate adaptation to new datasets.

## Curation Summary

- **Total scripts curated**: 144 (Python: 101, Shell: 22, YAML: 4, Markdown: 17)
- **Data files excluded**: 78
- **Duplicate/debug scripts archived**: 11
- **Documentation files created**: 4

All scripts have been generalized by replacing project-specific absolute paths with portable placeholders. No analysis code was executed during curation.

## Folder Structure

- `core_gwas/`: Quality control, PCA, covariate preparation, and association testing scripts
- `post_gwas/`: Significant SNP extraction, locus merging, cross-trait analysis, LD analysis, haplotype, coloc, and comparison scripts
- `fine_mapping/`: Regional plotting, LD heatmaps, annotation, and haplotype effect scripts
- `workflows/`: Workflow runners and orchestration scripts
- `agent_workflows/`: Agent prompts, skills, and report generation templates
- `configs/`: Configuration file templates
- `docs/`: Documentation and pipeline overview
- `archive/`: Project-specific examples and duplicate scripts

## Curation Notes

- **Included**: Reusable Python, R, shell, YAML, and workflow documentation files
- **Excluded**: Raw data, intermediate results, tables, figures, and output files
- **Generalized**: All absolute paths replaced with generic placeholders (e.g., `data/input/`, `data/reference/`, `results/`)

## How to Adapt Scripts

1. Update path placeholders in scripts to match your local directory structure
2. Review configuration templates in `configs/`
3. Consult `docs/path_generalization_rules.md` for mapping details
4. Use soybean project examples in `archive/` as reference implementations

## Warning

No data or analysis results are bundled with this repository. Scripts are provided as-is and may require modification for your specific analysis pipeline.# genomics-workflow-hub
