
awk '{$6=1; print $0}' 01_plink/01_convert/soybean_final.fam > temp.fam
mv temp.fam 01_plink/01_convert/soybean_final.fam
gemma -bfile 01_plink/01_convert/soybean_final -gk 1 -o kinship