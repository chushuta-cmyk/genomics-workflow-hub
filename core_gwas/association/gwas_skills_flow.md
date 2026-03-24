# GWAS and MR Analysis Skills Flow (Parameterized Template)

## Overview
This document provides abstract, parameterized templates for all skills in the GWAS and Mendelian Randomization pipeline. All paths, trait names, and specific values are replaced with placeholders for maximum reusability.

---

## Skill_01_PhenoTransform: Phenotype Transformation and Ratio Calculation

### Purpose
Transform raw phenotype data by calculating ratios and applying normalization transformations.

### Abstract Function Signature
```python
def transform_phenotypes(
    input_df: pd.DataFrame,
    numerator_col: str,
    denominator_col: str,
    ratio_col_prefix: str = "ratio",
    log_transform: bool = True,
    normalize: bool = False
) -> pd.DataFrame:
    """
    Calculate ratio traits and apply transformations.
    
    Parameters:
    -----------
    input_df : DataFrame
        Input phenotype DataFrame
    numerator_col : str
        Column name for numerator trait
    denominator_col : str
        Column name for denominator trait  
    ratio_col_prefix : str
        Prefix for ratio column names
    log_transform : bool
        Whether to apply log transformation
    normalize : bool
        Whether to z-score normalize traits
    
    Returns:
    --------
    DataFrame with added ratio and transformed columns
    """
```

### Parameterized Workflow
```python
# 1. Calculate raw ratio
${RATIO_COL}_raw = ${NUMERATOR_COL} / ${DENOMINATOR_COL}

# 2. Apply log transformation (optional)
if ${LOG_TRANSFORM}:
    ${RATIO_COL}_log = np.log(${RATIO_COL}_raw)
    
# 3. Normalize (optional)  
if ${NORMALIZE}:
    for col in ${TRAIT_COLUMNS}:
        df[col] = (df[col] - df[col].mean()) / df[col].std()
```

### Usage Example
```bash
# Calculate Oil:Protein ratio with log transformation
python phenotype_transform.py \
    --input ${PHENO_INPUT_TSV} \
    --numerator Oil \
    --denominator Protein \
    --output ${PHENO_OUTPUT_TSV} \
    --log-transform \
    --ratio-prefix Oil_Prot
```

---

## Skill_02_StrictQC: Genotype Quality Control

### Purpose
Apply standardized quality control filters to genotype data using PLINK2.

### Abstract Command Template
```bash
plink2 \
    --pfile ${INPUT_PFILE_PREFIX} \
    --allow-extra-chr \
    --maf ${MAF_THRESHOLD} \
    --geno ${GENO_THRESHOLD} \
    --max-alleles ${MAX_ALLELES} \
    --snps-only \
    ${ADDITIONAL_FILTERS} \
    --make-pgen \
    --out ${OUTPUT_PREFIX}
```

### Parameter Definitions
| Parameter | Default Value | Description |
|-----------|---------------|-------------|
| `${INPUT_PFILE_PREFIX}` | Required | Input PLINK2 pfile prefix |
| `${MAF_THRESHOLD}` | 0.05 | Minor allele frequency threshold |
| `${GENO_THRESHOLD}` | 0.1 | Maximum missing genotype rate |
| `${MAX_ALLELES}` | 2 | Maximum number of alleles per variant |
| `${ADDITIONAL_FILTERS}` | (Optional) | Additional PLINK2 filters |
| `${OUTPUT_PREFIX}` | Required | Output file prefix |

### Standard QC Pipeline
```bash
# Step 1: Basic QC filtering
plink2 --pfile ${GENO_PREFIX} \
       --maf ${MAF_THRESHOLD} \
       --geno ${GENO_THRESHOLD} \
       --max-alleles ${MAX_ALLELES} \
       --make-pgen --out ${QC_PREFIX}

# Step 2: Sample QC (kinship filtering)
plink2 --pfile ${QC_PREFIX} \
       --king-cutoff ${KING_CUTOFF} \
       --out ${SAMPLE_QC_PREFIX}

# Step 3: LD pruning
plink2 --pfile ${QC_PREFIX} \
       --indep-pairwise ${WINDOW_SIZE} ${STEP_SIZE} ${R2_THRESHOLD} \
       --out ${LD_PREFIX}
```

---

## Skill_03_GWAS_LMM: Multi-trait GWAS with Linear Mixed Models

### Purpose
Perform genome-wide association analysis for multiple traits using GEMMA with linear mixed models.

### Abstract Workflow Template

#### 1. PCA Generation
```bash
# Convert to PLINK1 format for GEMMA
plink2 --pfile ${QC_PREFIX} --make-bed --out ${BED_PREFIX}

# Run PCA
plink --bfile ${BED_PREFIX} \
      --allow-extra-chr \
      --pca ${NUM_PCS} \
      --out ${PCA_PREFIX}
```

#### 2. Prepare GEMMA Input Files
```python
def prepare_gemma_input(
    pheno_df: pd.DataFrame,
    eigenvec_file: str,
    traits: List[str],
    id_col: str = "IID",
    num_pcs: int = 3
) -> Tuple[str, str]:
    """
    Prepare phenotype and covariate files for GEMMA.
    
    Returns:
    --------
    (pheno_file_path, cov_file_path)
    """
    # Read PCA results
    pca = pd.read_csv(eigenvec_file, sep=r"\s+", header=None)
    pca.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, num_pcs+1)]
    
    # Merge with phenotypes
    merged = pca.merge(pheno_df, left_on="IID", right_on=id_col, how="inner")
    
    # Save phenotype file (NA -> -9 for GEMMA)
    pheno_data = merged[traits].fillna(-9)
    pheno_data.to_csv(${PHENO_OUTPUT}, sep="\t", index=False, header=False)
    
    # Save covariate file
    merged["Intercept"] = 1
    cov_data = merged[["Intercept"] + [f"PC{i}" for i in range(1, num_pcs+1)]]
    cov_data.to_csv(${COV_OUTPUT}, sep="\t", index=False, header=False)
```

#### 3. Kinship Matrix Calculation
```bash
gemma -bfile ${BED_PREFIX} \
      -p ${PHENO_FILE} \
      -gk 1 \
      -outdir ${OUTPUT_DIR} \
      -o ${KINSHIP_PREFIX}
```

#### 4. Multi-trait GWAS Loop
```bash
# Loop through traits
for TRAIT_INDEX in {1..${NUM_TRAITS}}; do
    TRAIT_NAME=$(get_trait_name ${TRAIT_INDEX})  # Abstract trait name retrieval
    
    echo "Running GWAS for trait ${TRAIT_INDEX}: ${TRAIT_NAME}"
    
    gemma -bfile ${BED_PREFIX} \
          -k ${KINSHIP_FILE} \
          -c ${COV_FILE} \
          -p ${PHENO_FILE} \
          -n ${TRAIT_INDEX} \
          -lmm 4 \
          -o ${OUTPUT_PREFIX}_${TRAIT_NAME}
done
```

#### 5. Results Extraction Template
```python
def extract_gwas_results(
    assoc_files: Dict[str, str],
    p_threshold: float = 5e-8,
    output_file: str = "gwas_summary.tsv"
) -> pd.DataFrame:
    """
    Extract and summarize GWAS results from multiple traits.
    """
    summaries = []
    for trait, filepath in assoc_files.items():
        df = pd.read_csv(filepath, sep="\t")
        df["trait"] = trait
        df["-log10P"] = -np.log10(df["p_wald"])
        summaries.append(df)
    
    return pd.concat(summaries)
```

---

## Skill_04_Asymmetric_MR: Bidirectional Mendelian Randomization

### Purpose
Perform bidirectional Mendelian Randomization to test causal relationships between traits, with asymmetry detection.

### Abstract Algorithm

#### 1. Instrument Variable Selection
```r
# R function template for IV selection
select_ivs <- function(gwas_results, 
                       p_threshold = 5e-8,
                       clump_r2 = 0.001,
                       clump_kb = 10000) {
  
  # Filter significant SNPs
  ivs <- gwas_results %>%
    filter(pval < p_threshold) %>%
    
    # LD clumping
    ieugwasr::ld_clump(
      clump_r2 = clump_r2,
      clump_kb = clump_kb,
      pop = "EUR"  # ${POPULATION} placeholder
    )
  
  return(ivs)
}
```

#### 2. Bidirectional MR Analysis
```r
# Template for bidirectional MR
run_bidirectional_mr <- function(exposure_gwas, 
                                 outcome_gwas,
                                 exposure_name,
                                 outcome_name) {
  
  # Harmonize effect alleles
  dat <- harmonise_data(
    exposure_dat = exposure_gwas,
    outcome_dat = outcome_gwas
  )
  
  # Run MR methods
  mr_results <- mr(dat, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger"))
  
  # Calculate asymmetry
  forward_result <- mr_results %>% filter(method == "mr_ivw")
  # ... reverse analysis would be run separately
  
  return(list(
    mr_results = mr_results,
    harmonised_data = dat,
    asymmetry_ratio = NA  # Calculated when both directions available
  ))
}
```

#### 3. Asymmetry Detection Logic
```python
def detect_asymmetry(forward_beta: float, 
                     forward_se: float,
                     reverse_beta: float, 
                     reverse_se: float,
                     threshold: float = 1.5) -> Dict:
    """
    Detect asymmetric causal effects.
    
    Parameters:
    -----------
    threshold : float
        Ratio threshold for declaring asymmetry (e.g., 1.5 = 50% difference)
    
    Returns:
    --------
    Dictionary with asymmetry metrics
    """
    ratio = abs(forward_beta) / abs(reverse_beta)
    
    return {
        "forward_beta": forward_beta,
        "reverse_beta": reverse_beta,
        "asymmetry_ratio": ratio,
        "is_asymmetric": ratio > threshold or ratio < 1/threshold,
        "stronger_direction": "forward" if abs(forward_beta) > abs(reverse_beta) else "reverse"
    }
```

#### 4. Sensitivity Analysis Template
```r
# MR-Egger regression for pleiotropy assessment
mr_egger <- mr_egger_regression(b_exp = dat$beta.exposure,
                                b_out = dat$beta.outcome,
                                se_exp = dat$se.exposure,
                                se_out = dat$se.outcome)

# Weighted median for robust estimation  
mr_weighted_median <- mr_weighted_median(b_exp = dat$beta.exposure,
                                         b_out = dat$beta.outcome,
                                         se_exp = dat$se.exposure,
                                         se_out = dat$se.outcome)

# Leave-one-out analysis
loo_results <- mr_leaveoneout(dat)
```

---

## Skill_05_Genetic_Coupling: Genetic Architecture Analysis

### Purpose
Distinguish between genetic coupling (shared variants) vs. causal chains using conditional genetic correlation analysis.

### Mathematical Framework

#### 1. Marginal vs. Conditional Genetic Correlation
```
Let:
  R²_marginal(X,Y) = genetic correlation between traits X and Y
  R²_conditional(X|Y) = genetic correlation of X with Y after conditioning on Y
  R²_conditional(Y|X) = genetic correlation of Y with X after conditioning on X

Decision rules:
  1. If R²_marginal is HIGH and both R²_conditional are ≈0:
     → Traits are GENETICALLY COUPLED (same variants)
     
  2. If R²_marginal is HIGH and R²_conditional asymmetric:
     → Potential CAUSAL CHAIN (one trait drives the other)
     
  3. If R²_marginal is LOW:
     → Traits are GENETICALLY INDEPENDENT
```

#### 2. Implementation Template
```python
def analyze_genetic_coupling(
    trait1_gwas: pd.DataFrame,
    trait2_gwas: pd.DataFrame,
    genotype_matrix: np.ndarray,
    pve_threshold: float = 0.01
) -> Dict:
    """
    Analyze genetic coupling between two traits.
    
    Steps:
    1. Calculate marginal genetic correlation
    2. Calculate conditional genetic correlations
    3. Apply decision rules
    """
    
    # 1. Marginal genetic correlation
    r2_marginal = calculate_genetic_correlation(trait1_gwas, trait2_gwas)
    
    # 2. Conditional correlations
    # Regress trait2 SNPs out of trait1
    trait1_residual = regress_out(trait1_gwas, trait2_gwas)
    r2_cond_1 = calculate_genetic_correlation(trait1_residual, trait2_gwas)
    
    # Regress trait1 SNPs out of trait2
    trait2_residual = regress_out(trait2_gwas, trait1_gwas)
    r2_cond_2 = calculate_genetic_correlation(trait2_residual, trait1_gwas)
    
    # 3. Decision rules
    if r2_marginal > 0.1 and r2_cond_1 < pve_threshold and r2_cond_2 < pve_threshold:
        coupling_type = "GENETIC_COUPLING"
    elif r2_marginal > 0.1 and abs(r2_cond_1 - r2_cond_2) > 0.05:
        coupling_type = "CAUSAL_CHAIN"
    else:
        coupling_type = "INDEPENDENT"
    
    return {
        "r2_marginal": r2_marginal,
        "r2_conditional_1|2": r2_cond_1,
        "r2_conditional_2|1": r2_cond_2,
        "coupling_type": coupling_type,
        "asymmetry": abs(r2_cond_1 - r2_cond_2)
    }
```

#### 3. Visualization Template
```r
# R template for genetic coupling visualization
plot_genetic_coupling <- function(r2_matrix, trait_names) {
  
  # Create correlation heatmap
  heatmap <- ggplot(r2_matrix, aes(x = Trait1, y = Trait2, fill = R2)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0.5, limits = c(0, 1)) +
    theme_minimal()
  
  # Create conditional vs marginal plot
  scatter <- ggplot(coupling_results, 
                    aes(x = R2_marginal, y = R2_conditional, color = Coupling_Type)) +
    geom_point(size = 3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(x = "Marginal R²", y = "Conditional R²")
  
  return(list(heatmap = heatmap, scatter = scatter))
}
```

### Interpretation Guidelines

| Pattern | R²_marginal | R²_conditional | Interpretation |
|---------|-------------|----------------|----------------|
| Coupling | High (>0.2) | Both near-zero (<0.01) | Shared genetic variants |
| Causal chain | High (>0.2) | Asymmetric (one >0.1, one <0.01) | Directional causality |
| Independent | Low (<0.1) | Both low | Genetically independent |
| Partial overlap | Moderate (0.1-0.3) | Moderate (0.05-0.1) | Some shared, some unique variants |

---

## Skill_06_Visualization_Export: Abstract Plot Generation

### Purpose
Generate standardized visualizations for GWAS and MR results with parameterized templates.

### Plot Generation Template
```python
def generate_visualizations(
    gwas_results: Dict[str, pd.DataFrame],
    mr_results: Dict[str, Dict],
    coupling_results: Dict,
    output_dir: str,
    trait_names: List[str]
) -> List[str]:
    """
    Generate all visualization files.
    
    Returns:
    --------
    List of generated file paths
    """
    generated_files = []
    
    # 1. Manhattan plots (one per trait)
    for trait in trait_names:
        if trait in gwas_results:
            filepath = f"{output_dir}/gwas_manhattan_{trait}.png"
            plot_manhattan(gwas_results[trait], filepath)
            generated_files.append(filepath)
    
    # 2. MR plots
    mr_files = [
        f"{output_dir}/mr_forest_bidirectional.png",
        f"{output_dir}/mr_scatter_{mr_results['forward']['exposure']}_to_{mr_results['forward']['outcome']}.png",
        f"{output_dir}/mr_scatter_{mr_results['reverse']['exposure']}_to_{mr_results['reverse']['outcome']}.png"
    ]
    
    for filepath in mr_files:
        # Generate plot (abstracted)
        generated_files.append(filepath)
    
    # 3. Genetic coupling plots
    coupling_files = [
        f"{output_dir}/genetic_correlation_matrix.png",
        f"{output_dir}/conditional_r2_plot.png"
    ]
    
    return generated_files
```

### File Naming Convention
```
# GWAS plots
gwas_manhattan_${TRAIT_NAME}.png
gwas_qq_${TRAIT_NAME}.png

# MR plots  
mr_forest_${DIRECTION}.png  # e.g., bidirectional, forward, reverse
mr_scatter_${EXPOSURE}_to_${OUTCOME}.png
mr_funnel_${EXPOSURE}_to_${OUTCOME}.png

# Genetic architecture plots
genetic_correlation_matrix.png
conditional_r2_plot.png
coupling_diagram.png

# Summary plots
summary_effect_sizes.png
asymmetry_comparison.png
```

---

## Skill_07_Pipeline_Orchestration: Abstract Workflow

### Purpose
Orchestrate the complete analysis pipeline with parameterized configuration.

### Configuration Template (YAML)
```yaml
# config_pipeline.yaml
pipeline:
  name: "gwas_mr_analysis"
  version: "2.0"
  
input:
  genotypes:
    format: "pfile"
    prefix: "${GENO_PREFIX}"
  phenotypes:
    file: "${PHENO_FILE}"
    id_column: "IID"
    traits: ["${TRAIT_1}", "${TRAIT_2}", "${TRAIT_3}"]
  
qc:
  maf_threshold: 0.05
  geno_threshold: 0.1
  max_alleles: 2
  
analysis:
  gwas:
    method: "GEMMA_LMM"
    covariates: ["Intercept", "PC1", "PC2", "PC3"]
    kinship: true
  mr:
    p_threshold: 5e-8
    methods: ["IVW", "Weighted_median", "MR_Egger"]
    bidirectional: true
  
output:
  directory: "${OUTPUT_DIR}"
  formats: ["tsv", "png", "html"]
  reports: true
```

### Pipeline Execution Template
```python
def execute_pipeline(config: Dict) -> Dict:
    """
    Execute complete GWAS-MR pipeline.
    """
    results = {}
    
    # 1. QC
    results['qc'] = execute_skill_02_strictqc(config)
    
    # 2. GWAS
    results['gwas'] = execute_skill_03_gwas_lmm(config, results['qc'])
    
    # 3. MR
    results['mr'] = execute_skill_04_asymmetric_mr(config, results['gwas'])
    
    # 4. Genetic coupling
    results['coupling'] = execute_skill_05_genetic_coupling(config, results['gwas'])
    
    # 5. Visualization
    results['visualization'] = execute_skill_06_visualization_export(
        config, results
    )
    
    return results
```

---

## Implementation Notes

### Placeholder Convention
- `${UPPER_CASE_WITH_UNDERSCORES}`: Required parameters
- `[OPTIONAL_PARAM]`: Optional parameters with default values
- `#{COMMENT}`: Implementation notes

### Extension Points
1. **New trait transformations**: Add to `Skill_01_PhenoTransform`
2. **Additional QC filters**: Extend `Skill_02_StrictQC`
3. **Alternative association methods**: Modify `Skill_03_GWAS_LMM`
4. **New MR methods**: Add to `Skill_04_Asymmetric_MR`
5. **Additional genetic architecture tests**: Extend `Skill_05_Genetic_Coupling`

### Validation Rules
- All placeholders must be replaced before execution
- Parameter ranges should be validated (e.g., 0 < MAF < 0.5)
- File existence should be checked for all inputs
- Results should be saved with version metadata

---

## Example Instantiation

```python
# Concrete instantiation for soybean analysis
config = {
    "GENO_PREFIX": "/data/soybean/genotypes",
    "PHENO_FILE": "/data/soybean/phenotypes.tsv",
    "TRAITS": ["100SW", "Protein", "Oil", "Oil_Prot_ratio"],
    "OUTPUT_DIR": "/results/soybean_gwas_mr",
    "MAF_THRESHOLD": 0.05,
    "NUM_PCS": 3
}

# Execute pipeline
results = execute_pipeline(config)
```

---

*Document Version: 2.0 - Abstract Parameterized Template*
*Last Updated: March 15, 2025*
*Designed for: Reusable GWAS-MR analysis pipeline*