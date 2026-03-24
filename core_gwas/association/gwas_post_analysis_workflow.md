# GWAS Post-Analysis Workflow: Multi-Trait SNP Integration and Colocalization

## Overview
This workflow provides a standardized pipeline for post-GWAS analysis, including significant SNP extraction, cross-trait effect comparison, locus merging, and Bayesian colocalization. The workflow is designed to be reusable for wild, cultivated, or combined soybean datasets.

**Version:** 1.0  
**Last Updated:** March 6, 2025  
**Author:** Bioinformatics Pipeline Team

---

## 1. Workflow Structure

```
├── analysis_index/                 # Step 1: GWAS file indexing
│   └── gwas_files.txt
├── SNP_analyze/                    # Steps 2-7: Main analysis outputs
│   ├── gwas_size.tsv              # Step 2: Extracted GWAS tables
│   ├── gwas_protein.tsv
│   ├── gwas_oil.tsv
│   ├── sig_size.tsv               # Step 3: Significant SNPs
│   ├── sig_protein.tsv
│   ├── sig_oil.tsv
│   ├── sig_all_traits.tsv
│   ├── cross_trait_effect.tsv     # Step 4: Cross-trait effects
│   ├── locus_summary.tsv          # Step 5: Merged loci
│   ├── coloc_results.tsv          # Step 6: Colocalization results
│   └── analysis_summary.md        # Step 7: Final report
└── scripts/                       # Reusable analysis scripts
    ├── extract_gwas_tables.py
    ├── extract_significant_snps.py
    ├── cross_trait_effects.py
    ├── merge_loci.py
    ├── run_coloc.R
    └── generate_summary.py
```

---

## 2. Step-by-Step Workflow

### Step 1: GWAS Result Index
**Purpose:** Locate GWAS result files and create a structured index.

**Input Files:**
- User-provided paths to GWAS result files (`.assoc.txt` format from GEMMA)
- Example: `/path/to/size_analysis.assoc.txt`, `/path/to/protein_analysis.assoc.txt`, `/path/to/oil_analysis.assoc.txt`

**Commands:**
```bash
# Create analysis directory structure
mkdir -p analysis_index SNP_analyze logs

# Index GWAS files (manual or automated)
echo "size:/path/to/size_analysis.assoc.txt" > analysis_index/gwas_files.txt
echo "protein:/path/to/protein_analysis.assoc.txt" >> analysis_index/gwas_files.txt
echo "oil:/path/to/oil_analysis.assoc.txt" >> analysis_index/gwas_files.txt
```

**Output Files:**
- `analysis_index/gwas_files.txt` - Index of GWAS files with trait labels

---

### Step 2: Extract GWAS Tables
**Purpose:** Extract essential columns from GWAS result files for efficient processing.

**Input Files:**
- GWAS result files indexed in `analysis_index/gwas_files.txt`

**Commands:**
```bash
# Using the Python script
python scripts/extract_gwas_tables.py \
    --index analysis_index/gwas_files.txt \
    --output-dir SNP_analyze \
    --columns chr rs ps allele1 allele0 af beta se p_wald
```

**Script:** `extract_gwas_tables.py`  
**Parameters:**
- `--index`: Path to GWAS file index
- `--output-dir`: Directory for output files
- `--columns`: Columns to extract (default: essential GWAS columns)

**Output Files:**
- `SNP_analyze/gwas_size.tsv`
- `SNP_analyze/gwas_protein.tsv`
- `SNP_analyze/gwas_oil.tsv`

**File Format:**
```
chr     rs      ps      allele1 allele0 af      beta    se      p_wald
1       1:114:A:C       114     C       A       0.500   0.9471  0.9626  0.3258
...
```

---

### Step 3: Significant SNP Extraction
**Purpose:** Filter SNPs with genome-wide significance threshold.

**Input Files:**
- Extracted GWAS tables from Step 2

**Threshold:** `p_wald < 5e-8`

**Commands:**
```bash
# Using the Python script
python scripts/extract_significant_snps.py \
    --gwas-tables SNP_analyze/gwas_size.tsv SNP_analyze/gwas_protein.tsv SNP_analyze/gwas_oil.tsv \
    --trait-names size protein oil \
    --p-threshold 5e-8 \
    --output-dir SNP_analyze
```

**Script:** `extract_significant_snps.py`  
**Parameters:**
- `--gwas-tables`: List of GWAS table files
- `--trait-names`: Corresponding trait names
- `--p-threshold`: Significance threshold (default: 5e-8)
- `--output-dir`: Output directory

**Output Files:**
- `SNP_analyze/sig_size.tsv` - Significant SNPs for size trait
- `SNP_analyze/sig_protein.tsv` - Significant SNPs for protein trait
- `SNP_analyze/sig_oil.tsv` - Significant SNPs for oil trait
- `SNP_analyze/sig_all_traits.tsv` - All significant SNPs with trait annotation

**File Format (sig_all_traits.tsv):**
```
rs      chr     pos     trait   beta    p_wald  allele1 allele0 af
1:46003733:C:T  1       46003733        size    -13.99011       2.44532e-08    T       C       0.499
5:5567603:G:A   5       5567603 protein 3.190022       4.189179e-08    A       G       0.493
...
```

---

### Step 4: Cross-Trait Effect Extraction
**Purpose:** For each significant SNP, retrieve effect estimates from all GWAS results.

**Input Files:**
- All GWAS tables (Step 2)
- All significant SNPs (Step 3)

**Commands:**
```bash
# Using the Python script
python scripts/cross_trait_effects.py \
    --gwas-tables SNP_analyze/gwas_size.tsv SNP_analyze/gwas_protein.tsv SNP_analyze/gwas_oil.tsv \
    --trait-names size protein oil \
    --significant-snps SNP_analyze/sig_all_traits.tsv \
    --output SNP_analyze/cross_trait_effect.tsv
```

**Script:** `cross_trait_effects.py`  
**Parameters:**
- `--gwas-tables`: GWAS table files for all traits
- `--trait-names`: Corresponding trait names
- `--significant-snps`: File with all significant SNPs
- `--output`: Output file path

**Output Files:**
- `SNP_analyze/cross_trait_effect.tsv`

**File Format:**
```
rs      chr     pos     trait_significant  beta_size  p_wald_size  beta_protein  p_wald_protein  beta_oil  p_wald_oil
1:46003733:C:T  1       46003733        size    -13.99011   2.44532e-08    1.769797     0.1222849    3.916591   0.1459349
5:5567603:G:A   5       5567603        protein -3.216118   0.01158952     3.190022     4.189179e-08  0.332493   0.8308358
...
```

---

### Step 5: Locus Merge
**Purpose:** Merge nearby SNPs into loci using a genomic window.

**Window:** ±250kb (500kb total)

**Input Files:**
- Cross-trait effect table from Step 4

**Commands:**
```bash
# Using the Python script
python scripts/merge_loci.py \
    --snp-table SNP_analyze/cross_trait_effect.tsv \
    --window 250000 \
    --output SNP_analyze/locus_summary.tsv
```

**Script:** `merge_loci.py`  
**Parameters:**
- `--snp-table`: SNP effect table from Step 4
- `--window`: Window size in base pairs (default: 250000)
- `--output`: Output file path

**Output Files:**
- `SNP_analyze/locus_summary.tsv`

**File Format:**
```
locus_id        chr     start   end     num_snps        lead_snp        lead_p  traits_in_locus
locus_1         1       45950000        46200000        3       1:46003733:C:T  2.44532e-08    size
locus_2         5       5500000 6000000 2       5:5567603:G:A   4.189179e-08    protein
...
```

**Algorithm:**
1. Sort SNPs by chromosome and position
2. Merge SNPs within ±250kb into same locus
3. Identify lead SNP (lowest p-value) for each locus
4. Summarize trait associations in each locus

---

### Step 6: Colocalization Analysis
**Purpose:** Run Bayesian colocalization analysis for trait pairs.

**Trait Pairs:**
1. Oil vs Protein
2. Oil vs Size  
3. Protein vs Size

**Methods:**
- **coloc**: Bayesian colocalization (Giambartolomei et al. 2014)
- **eCAVIAR**: Colocalization across multiple traits (Hormozdiari et al. 2016)

**Input Files:**
- GWAS tables for all traits (Step 2)
- Locus summary (Step 5)

**Commands:**
```bash
# Using the R script
Rscript scripts/run_coloc.R \
    --gwas-files SNP_analyze/gwas_size.tsv SNP_analyze/gwas_protein.tsv SNP_analyze/gwas_oil.tsv \
    --trait-names size protein oil \
    --loci SNP_analyze/locus_summary.tsv \
    --output SNP_analyze/coloc_results.tsv
```

**Script:** `run_coloc.R`  
**Parameters:**
- `--gwas-files`: GWAS table files
- `--trait-names`: Corresponding trait names
- `--loci`: Locus summary file
- `--output`: Output file path

**Output Files:**
- `SNP_analyze/coloc_results.tsv`

**File Format:**
```
locus_id        trait_pair       coloc_pp4     coloc_pp3     ecaviar_clpp    ecaviar_set
locus_1         oil_vs_protein   0.95          0.03          0.92            {1:46003733:C:T}
locus_1         oil_vs_size      0.87          0.10          0.85            {1:46003733:C:T,1:46012345:A:G}
...
```

**Interpretation:**
- `coloc_pp4`: Posterior probability of colocalization (PP4 > 0.8 suggests colocalization)
- `coloc_pp3`: Posterior probability of distinct causal variants (PP3 > 0.8 suggests distinct variants)
- `ecaviar_clpp`: Combined likelihood of colocalization from eCAVIAR
- `ecaviar_set`: Set of candidate causal SNPs

---

### Step 7: Report Export
**Purpose:** Generate final analysis summary report.

**Input Files:**
- All intermediate results (Steps 1-6)

**Commands:**
```bash
# Using the Python script
python scripts/generate_summary.py \
    --sig-snps SNP_analyze/sig_all_traits.tsv \
    --cross-trait SNP_analyze/cross_trait_effect.tsv \
    --loci SNP_analyze/locus_summary.tsv \
    --coloc SNP_analyze/coloc_results.tsv \
    --output SNP_analyze/analysis_summary.md
```

**Script:** `generate_summary.py`  
**Parameters:**
- `--sig-snps`: Significant SNPs file
- `--cross-trait`: Cross-trait effects file
- `--loci`: Locus summary file
- `--coloc`: Colocalization results file
- `--output`: Output markdown file path

**Output Files:**
- `SNP_analyze/analysis_summary.md`

**Report Contents:**
1. Executive summary
2. Significant SNP counts by trait
3. Cross-trait effect patterns
4. Locus summary table
5. Colocalization results
6. Key findings and biological interpretation

---

## 3. Configuration and Customization

### 3.1 Dataset Configuration
The workflow supports three dataset types:
- **Wild soybean**: Use wild-specific allele frequencies and population structure
- **Cultivated soybean**: Use cultivated-specific parameters
- **Combined**: Pooled analysis with population as covariate

**Configuration file template (`config_analysis.yaml`):**
```yaml
dataset:
  type: "wild"  # wild, cultivated, or combined
  name: "soybean_gwas"
  
traits:
  names: ["size", "protein", "oil"]
  files:
    size: "/path/to/size_analysis.assoc.txt"
    protein: "/path/to/protein_analysis.assoc.txt"
    oil: "/path/to/oil_analysis.assoc.txt"
  
thresholds:
  p_significant: 5e-8
  coloc_pp4: 0.8
  locus_window: 250000
  
output:
  directory: "./results"
  formats: ["tsv", "md", "log"]
```

### 3.2 Parameter Adjustments

| Parameter | Default | Description | Adjustment Guidelines |
|-----------|---------|-------------|----------------------|
| `p_threshold` | 5e-8 | Genome-wide significance | Increase for less stringent (1e-6), decrease for more stringent (1e-9) |
| `locus_window` | 250000 | Locus merging window | Decrease for fine-mapping (100000), increase for broad regions (500000) |
| `coloc_pp4_threshold` | 0.8 | Colocalization confidence | Decrease for sensitive detection (0.7), increase for stringent (0.9) |

### 3.3 Adding New Traits
To add new traits to the analysis:

1. **Update configuration:**
```yaml
traits:
  names: ["size", "protein", "oil", "new_trait"]
  files:
    new_trait: "/path/to/new_trait.assoc.txt"
```

2. **Rerun workflow:**
```bash
# Update index
echo "new_trait:/path/to/new_trait.assoc.txt" >> analysis_index/gwas_files.txt

# Rerun from Step 2
python scripts/extract_gwas_tables.py --index analysis_index/gwas_files.txt --output-dir SNP_analyze
```

---

## 4. Quality Control and Validation

### 4.1 Input Validation
Each script includes input validation:
- File existence and format checks
- Column name verification
- Data type validation
- Missing value handling

### 4.2 Output Verification
- All output files include header rows
- Consistent column naming across files
- Log files record processing steps and warnings
- Summary statistics reported in final output

### 4.3 Error Handling
- Graceful failure with informative error messages
- Checkpoint files allow restart from failed steps
- Log files capture all warnings and errors

---

## 5. Performance Considerations

### 5.1 Memory Management
- **Large GWAS files**: Process in chunks using pandas `chunksize` parameter
- **Cross-trait lookup**: Use dictionary indexing (O(1) lookup)
- **Locus merging**: Sort-then-merge algorithm (O(n log n))

### 5.2 Parallel Processing
**For large datasets (>1M SNPs):**
```python
# Parallel trait processing example
from concurrent.futures import ProcessPoolExecutor

with ProcessPoolExecutor(max_workers=4) as executor:
    futures = []
    for trait, filepath in gwas_files.items():
        future = executor.submit(process_trait, trait, filepath)
        futures.append(future)
```

### 5.3 Disk Usage
- **Intermediate files**: Can be compressed (`.tsv.gz`) to save space
- **Cleanup option**: Optional removal of intermediate files
- **Checkpointing**: Save progress to allow incremental processing

---

## 6. Example Execution

### 6.1 Complete Pipeline Run
```bash
#!/bin/bash
# Complete workflow execution script

set -e  # Exit on error

# Step 1: Index GWAS files
echo "=== Step 1: Indexing GWAS files ==="
mkdir -p analysis_index SNP_analyze logs
echo "size:/data/gwas/size_analysis.assoc.txt" > analysis_index/gwas_files.txt
echo "protein:/data/gwas/protein_analysis.assoc.txt" >> analysis_index/gwas_files.txt
echo "oil:/data/gwas/oil_analysis.assoc.txt" >> analysis_index/gwas_files.txt

# Step 2: Extract GWAS tables
echo "=== Step 2: Extracting GWAS tables ==="
python scripts/extract_gwas_tables.py \
    --index analysis_index/gwas_files.txt \
    --output-dir SNP_analyze \
    > logs/step2_extract.log 2>&1

# Step 3: Extract significant SNPs
echo "=== Step 3: Extracting significant SNPs ==="
python scripts/extract_significant_snps.py \
    --gwas-tables SNP_analyze/gwas_size.tsv SNP_analyze/gwas_protein.tsv SNP_analyze/gwas_oil.tsv \
    --trait-names size protein oil \
    --p-threshold 5e-8 \
    --output-dir SNP_analyze \
    > logs/step3_significant.log 2>&1

# Step 4: Cross-trait effects
echo "=== Step 4: Calculating cross-trait effects ==="
python scripts/cross_trait_effects.py \
    --gwas-tables SNP_analyze/gwas_size.tsv SNP_analyze/gwas_protein.tsv SNP_analyze/gwas_oil.tsv \
    --trait-names size protein oil \
    --significant-snps SNP_analyze/sig_all_traits.tsv \
    --output SNP_analyze/cross_trait_effect.tsv \
    > logs/step4_cross_trait.log 2>&1

# Step 5: Merge loci
echo "=== Step 5: Merging loci ==="
python scripts/merge_loci.py \
    --snp-table SNP_analyze/cross_trait_effect.tsv \
    --window 250000 \
    --output SNP_analyze/locus_summary.tsv \
    > logs/step5_loci.log 2>&1

# Step 6: Colocalization analysis
echo "=== Step 6: Running colocalization ==="
Rscript scripts/run_coloc.R \
    --gwas-files SNP_analyze/gwas_size.tsv SNP_analyze/gwas_protein.tsv SNP_analyze/gwas_oil.tsv \
    --trait-names size protein oil \
    --loci SNP_analyze/locus_summary.tsv \
    --output SNP_analyze/coloc_results.tsv \
    > logs/step6_coloc.log 2>&1

# Step 7: Generate summary
echo "=== Step 7: Generating summary report ==="
python scripts/generate_summary.py \
    --sig-snps SNP_analyze/sig_all_traits.tsv \
    --cross-trait SNP_analyze/cross_trait_effect.tsv \
    --loci SNP_analyze/locus_summary.tsv \
    --coloc SNP_analyze/coloc_results.tsv \
    --output SNP_analyze/analysis_summary.md \
    > logs/step7_summary.log 2>&1

echo "=== Pipeline completed successfully ==="
echo "Results saved in: SNP_analyze/"
echo "Logs saved in: logs/"
```

### 6.2 Quick Test Run
For testing with small datasets:
```bash
# Test with subset of data
head -n 10000 /data/gwas/size_analysis.assoc.txt > test_size.assoc.txt
# Update index and run minimal workflow
```

---

## 7. Troubleshooting Guide

### 7.1 Common Issues

| Issue | Solution |
|-------|----------|
| Missing column error | Check GWAS file format, ensure required columns exist |
| Memory error | Use `--chunksize` parameter, process in batches |
| R package installation error | Install required R packages: `install.packages(c('coloc', 'dplyr'))` |
| No significant SNPs | Adjust p-threshold or check GWAS quality |
| Locus merging too aggressive | Decrease `--window` parameter |

### 7.2 Log Files
- Check `logs/` directory for detailed error messages
- Each step generates timestamped log file
- Logs include parameter settings and processing statistics

### 7.3 Debug Mode
Enable verbose output:
```bash
python scripts/extract_gwas_tables.py --verbose --debug ...
```

---

## 8. References and Dependencies

### 8.1 Software Dependencies
- **Python 3.8+** with packages: pandas, numpy, argparse, pyyaml
- **R 4.0+** with packages: coloc, dplyr, data.table, argparse
- **Bash** for workflow orchestration

### 8.2 Methodological References
1. **GWAS significance threshold**: p < 5e-8 (genome-wide significance)
2. **Locus definition**: ±250kb window (Wellcome Trust Case Control Consortium)
3. **Colocalization**: Giambartolomei et al. (2014) PLoS Genet
4. **eCAVIAR**: Hormozdiari et al. (2016) Am J Hum Genet

### 8.3 Related Workflows
- **GWAS QC pipeline**: `gwas_skills_flow.md`
- **MR analysis pipeline**: `gwas_mr_report.md`
- **Phenotype transformation**: `pheno_transform_template.md`

---

## 9. Maintenance and Updates

### 9.1 Version History
- **v1.0 (2025-03-06)**: Initial release with 7-step workflow

### 9.2 Contributing
To extend or modify this workflow:
1. Fork the workflow repository
2. Create feature branch
3. Update scripts and documentation
4. Test with sample data
5. Submit pull request

### 9.3 Support
For questions or issues:
- Check troubleshooting guide (Section 7)
- Review example execution (Section 6)
- Contact: bioinformatics-pipeline@example.org

---

*Workflow designed for soybean GWAS post-analysis but adaptable to other crop species.*  
*Maintained by: Soybean Genomics Bioinformatics Team*  
*Last validated: March 6, 2025*