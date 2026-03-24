# Next Steps for Wild Soybean GWAS Analysis

## Immediate Tasks (Post‑GWAS Finishing)
1. **LD decay summary**  
   - Parse `wild_ld_decay_full.ld` to compute average \(r^2\) vs. distance bins.  
   - Generate LD decay curve plot (distance vs. \(r^2\)).  
   - Save summary table (mean \(r^2\) per 10‑kb bin).  

2. **Haplotype block detection**  
   - Re‑run PLINK with `--blocks no-pheno-req` to obtain haplotype blocks.  
   - Summarize block sizes, number of blocks per chromosome.  
   - Overlap GWAS loci with haplotype blocks.  

3. **GWAS loci overlap with LD blocks**  
   - Intersect locus intervals (`locus_summary_*.tsv`) with haplotype blocks.  
   - Identify loci that fall within the same LD block (potential co‑inherited signals).  

## Medium‑Term Integrative Analyses
4. **Bayesian colocalization summary**  
   - Obtain eQTL summary statistics from soybean seed tissues (public databases).  
   - Run coloc for each GWAS locus overlapping an eQTL region.  
   - Generate a summary table of posterior probabilities for colocalization.  

5. **Selection scan integration**  
   - Obtain selection scan results (XP‑EHH, F<sub>ST</sub>) between wild and cultivated soybean.  
   - Test for enrichment of GWAS signals in regions under selection.  
   - Prioritize loci that are both associated with a trait and show signatures of selection.  

6. **Gene annotation & functional enrichment**  
   - Annotate each locus with genes within 50 kb (using soybean reference annotation).  
   - Perform Gene Ontology (GO) and KEGG pathway enrichment for each trait.  
   - Use tools such as `clusterProfiler` or `g:Profiler`.  

7. **Pathway analysis (MAGMA)**  
   - Convert GWAS summary statistics to MAGMA format.  
   - Run gene‑based and pathway‑based tests.  
   - Compare results across traits.  

8. **Comparison with cultivated soybean GWAS**  
   - Collect published GWAS results for the same traits in cultivated soybean.  
   - Compare lead SNP positions, effect directions, and allele frequencies.  
   - Identify conserved vs. lineage‑specific associations.  

## Long‑Term / Advanced
9. **Multi‑trait meta‑analysis**  
   - Perform multivariate GWAS (e.g., MTAG) to increase power and detect shared genetic components.  

10. **Transcriptome‑wide association study (TWAS)**  
    - Integrate with soybean expression reference panels (if available).  
    - Impute gene expression in wild soybean and test for association with traits.  

11. **Network‑based analysis**  
    - Build gene‑interaction networks for candidate genes.  
    - Identify key hub genes and subnetworks enriched for trait associations.  

## Data & Code Management
12. **Documentation**  
    - Update README with final output descriptions.  
    - Archive raw GWAS summary statistics (gwas_tables/) and significant SNP lists.  

13. **Reproducibility**  
    - Create a Snakemake/Nextflow pipeline that chains all steps from GWAS to downstream reporting.  
    - Version‑control scripts and configuration files.  

---
**Last Updated:** 2026‑03‑09  
**Priority:** 1–3 are high priority for immediate follow‑up.  
**Estimated Effort:** 2–3 weeks for tasks 1–5.