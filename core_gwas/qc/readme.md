# GWAS Post-Analysis Workflows

This directory contains reusable workflows and scripts for GWAS post-analysis, including significant SNP extraction, cross-trait effect comparison, locus merging, and Bayesian colocalization.

## Contents

### 1. Workflow Documentation
- `gwas_post_analysis_workflow.md`: Complete workflow documentation with 7-step pipeline
- `config_analysis_template.yaml`: Configuration template for custom analyses

### 2. Analysis Scripts (in `scripts/` directory)
- `extract_gwas_tables.py`: Extract essential columns from GWAS results
- `extract_significant_snps.py`: Filter genome-wide significant SNPs (p < 5e-8)
- `cross_trait_effects.py`: Extract SNP effects across multiple traits
- `merge_loci.py`: Merge nearby SNPs into loci (±250kb window)
- `run_coloc.R`: Bayesian colocalization analysis (coloc and eCAVIAR)
- `generate_summary.py`: Generate comprehensive analysis report

### 3. Pipeline Orchestration
- `run_pipeline.py`: Main pipeline runner with configuration support

## Quick Start

### Option 1: Quick Execution
```bash
python run_pipeline.py \
    --gwas-files size.assoc.txt protein.assoc.txt oil.assoc.txt \
    --trait-names size protein oil \
    --output-dir ./results
```

### Option 2: Configuration-Based Execution
1. Copy and edit configuration template:
   ```bash
   cp config_analysis_template.yaml config_analysis.yaml
   # Edit config_analysis.yaml with your file paths and parameters
   ```

2. Run pipeline:
   ```bash
   python run_pipeline.py --config config_analysis.yaml
   ```

## Pipeline Steps

The workflow consists of 7 standardized steps:

1. **GWAS Result Indexing**: Create index of GWAS files
2. **GWAS Table Extraction**: Extract essential columns for efficient processing
3. **Significant SNP Extraction**: Filter SNPs with p_wald < 5e-8
4. **Cross-Trait Effect Extraction**: Retrieve SNP effects across all traits
5. **Locus Merging**: Group SNPs within ±250kb into loci
6. **Colocalization Analysis**: Test for shared causal variants between trait pairs
7. **Report Generation**: Create comprehensive summary report

## Output Structure

```
results/
├── analysis_index/
│   └── gwas_files.txt          # Step 1: GWAS file index
├── SNP_analyze/
│   ├── gwas_*.tsv              # Step 2: Extracted GWAS tables
│   ├── sig_*.tsv               # Step 3: Significant SNPs
│   ├── cross_trait_effect.tsv  # Step 4: Cross-trait effects
│   ├── locus_summary.tsv       # Step 5: Merged loci
│   ├── coloc_results.tsv       # Step 6: Colocalization results
│   └── analysis_summary.md     # Step 7: Summary report
├── logs/                       # Execution logs
└── scripts/                    # Copy of analysis scripts
```

## Customization

### Adding New Traits
1. Add trait to configuration file:
   ```yaml
   traits:
     names: ["size", "protein", "oil", "new_trait"]
   dataset:
     gwas_files:
       new_trait: "/path/to/new_trait.assoc.txt"
   ```

2. Rerun pipeline

### Adjusting Parameters
- **Significance threshold**: Modify `p_threshold` in configuration
- **Locus window**: Adjust `locus.window` (default: 250000 bp)
- **Colocalization threshold**: Change `coloc.pp4_threshold` (default: 0.8)

## Requirements

### Python Dependencies
- pandas >= 1.3.0
- numpy >= 1.21.0
- pyyaml >= 6.0

### R Dependencies
- R >= 4.0.0
- coloc package
- dplyr package
- data.table package

## Validation

Test the pipeline with a small subset:
```bash
# Create test files
head -n 10000 size.assoc.txt > test_size.assoc.txt
head -n 10000 protein.assoc.txt > test_protein.assoc.txt
head -n 10000 oil.assoc.txt > test_oil.assoc.txt

# Run quick test
python run_pipeline.py \
    --gwas-files test_size.assoc.txt test_protein.assoc.txt test_oil.assoc.txt \
    --trait-names size protein oil \
    --output-dir ./test_results \
    --verbose
```

## Troubleshooting

### Common Issues

1. **Memory errors with large GWAS files**
   - Use `--chunksize` parameter in `extract_gwas_tables.py`
   - Process traits sequentially

2. **Missing R packages**
   ```r
   install.packages(c("coloc", "dplyr", "data.table"))
   ```

3. **No significant SNPs found**
   - Check GWAS file format and column names
   - Adjust p-value threshold if appropriate

### Log Files
Check `logs/` directory for detailed error messages from each step.

## Applications

This workflow is designed for:

1. **Wild soybean analysis**: Adapt for specific population characteristics
2. **Cultivated soybean analysis**: Adjust for different genetic backgrounds
3. **Combined analyses**: Pool data with appropriate covariates
4. **Other crop species**: Modify configuration for different genomes

## Citation

If using this workflow in research, please cite:

```
GWAS Post-Analysis Pipeline v1.0
Bioinformatics Team, Soybean Genomics Research Center
```

## Support

For questions or issues:
1. Check the workflow documentation
2. Review example configurations
3. Examine log files for error messages
4. Contact: bioinformatics-pipeline@example.org

---

*Last updated: March 6, 2025*
*Pipeline version: 1.0*