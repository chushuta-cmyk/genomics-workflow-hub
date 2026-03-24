#!/bin/bash

echo "=== Checking PLINK .fam file ==="
echo "First 5 lines of .fam file:"
head -5 01_plink/05_final/soybean_final_plink1.fam
echo ""
echo "Total individuals in .fam:"
wc -l 01_plink/05_final/soybean_final_plink1.fam
echo ""

echo "=== Checking phenotype file ==="
echo "First 5 lines of phenotype file:"
head -5 03_gwas_input/gwas_ready.pheno
echo ""
echo "Total lines in phenotype file:"
wc -l 03_gwas_input/gwas_ready.pheno
echo ""

echo "=== Checking for ID matches ==="
# Extract IDs from .fam (FID and IID, columns 1 and 2)
awk '{print $1"\t"$2}' 01_plink/05_final/soybean_final_plink1.fam | sort > /tmp/fam_ids.txt
echo "Sample .fam IDs:"
head -3 /tmp/fam_ids.txt
echo ""

# Extract IDs from phenotype file (assuming first two columns are FID and IID)
awk '{print $1"\t"$2}' 03_gwas_input/gwas_ready.pheno | sort > /tmp/pheno_ids.txt
echo "Sample phenotype IDs:"
head -3 /tmp/pheno_ids.txt
echo ""

# Find common IDs
echo "Number of matching IDs:"
comm -12 /tmp/fam_ids.txt /tmp/pheno_ids.txt | wc -l
echo ""

echo "=== Checking phenotype file structure ==="
# Check number of columns
echo "Number of columns in phenotype file:"
head -1 03_gwas_input/gwas_ready.pheno | awk '{print NF}'
echo ""

# Check for missing values (typically -9 or NA)
echo "Count of -9 values in phenotype column(s):"
awk '{for(i=3;i<=NF;i++) if($i==-9) count++} END{print count}' 03_gwas_input/gwas_ready.pheno
echo ""
echo "Count of NA values:"
grep -o "NA" 03_gwas_input/gwas_ready.pheno | wc -l
