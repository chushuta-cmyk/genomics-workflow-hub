# Path Generalization Rules

During curation, hard-coded absolute paths in scripts were replaced with generic placeholders to improve portability.

## Mapping Table

| Original Path Pattern | Generic Replacement | Notes |
|-----------------------|---------------------|-------|
| `/data03/karama/soybean_gwas_filtered/` | `data/input/` | Input data directory |
| `/data03/karama/soybean_gwas/` | `data/input/` | Input data directory |
| `/data03/karama/projects/soybean_analysis/` | `results/` | Analysis results directory |
| `/data03/karama/workflows/` | `workflows/` | Workflow scripts directory |
| `/data03/karama/documents/` | `docs/` | Documentation directory |
| `/data03/karama/` | `data/` | Generic data root |

## How to Adapt

1. Replace generic placeholders in scripts with your actual directory paths.
2. Update configuration files in `configs/` to match your dataset.
3. Ensure required input files are placed in the appropriate directories.

## Example

Original line in script:
```
/data03/karama/soybean_gwas_filtered/genotypes.vcf
```
Generalized to:
```
data/input/genotypes.vcf
```
