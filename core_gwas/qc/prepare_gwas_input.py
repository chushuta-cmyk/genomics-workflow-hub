import pandas as pd

def prepare_data():
    pheno_file = 'phenotype_long.tsv'
    pca_file = '01_plink/04_pca/soybean_pca.eigenvec'
    
    # 1. 读取表型：ID 在第二列
    df = pd.read_csv(pheno_file, sep='\t')
    df['ID'] = df['ID'].astype(str).str.strip()
    
    # 计算两年均值
    traits = ['Protein', 'Oil', '100SW']
    pheno_avg = df.groupby('ID')[traits].mean().reset_index()

    # 2. 读取 PCA：第一列是 ID (此时 header 是 #IID)
    # 使用 sep='\s+' 自动处理空格或 Tab
    pca = pd.read_csv(pca_file, sep='\s+')
    
    # 根据你提供的 head，ID 列名现在是 '#IID'
    # 提取 ID 和前三个主成分
    covariates = pca[['#IID', 'PC1', 'PC2', 'PC3']].copy()
    covariates.columns = ['ID', 'PC1', 'PC2', 'PC3']
    covariates['ID'] = covariates['ID'].astype(str).str.strip()

    # 3. 合并：以 PCA 的顺序为基准 (非常重要，必须对齐基因型顺序)
    combined = pd.merge(covariates, pheno_avg, on='ID', how='left')

    # 4. 生成 GEMMA 格式文件 (无表头)
    # 表型文件
    combined[traits].to_csv('gwas_ready.pheno', sep='\t', index=False, header=False, na_rep='NA')
    
    # 协变量文件 (第一列加 Intercept)
    combined['Intercept'] = 1
    combined[['Intercept', 'PC1', 'PC2', 'PC3']].to_csv('gwas_ready.cov', sep='\t', index=False, header=False)
    
    print(f"成功！匹配样本总数: {len(combined)}")
    print("前 3 行表型预览:\n", combined[traits].head(3))
    print("前 3 行协变量预览:\n", combined[['Intercept', 'PC1', 'PC2', 'PC3']].head(3))

if __name__ == "__main__":
    prepare_data()