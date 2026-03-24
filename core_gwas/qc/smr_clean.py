import pandas as pd

def clean_for_smr(file_path, output_path):
    # 1. 加载数据
    df = pd.read_csv(file_path, sep='\t')
    
    # 2. 过滤主染色体 (1-20)
    df['chr'] = df['chr'].astype(str)
    main_chrs = [str(i) for i in range(1, 21)]
    df = df[df['chr'].isin(main_chrs)]
    
    # 3. 处理“离群值”：
    # SMR 主要是看 Top SNP，但如果 Beta 值大得离谱（如之前看到的 -18），
    # 可能是单样本偏差。通常保留 p_wald < 1e-5 的点作为候选。
    # 同时剔除等位基因频率 (af) 极低的点 (MAF < 0.05)
    df = df[df['af'] > 0.05]
    df = df[df['af'] < 0.95]
    
    # 4. 转换格式以符合 SMR/GCTA 要求 (假设你需要: SNP, CHR, BP, A1, A2, freq, beta, se, p)
    # 注意：A1 是效应等位基因
    df_smr = df[['rs', 'chr', 'ps', 'allele1', 'allele0', 'af', 'beta', 'se', 'p_wald']]
    df_smr.columns = ['SNP', 'CHR', 'BP', 'A1', 'A2', 'freq', 'b', 'se', 'p']
    
    df_smr.to_csv(output_path, sep='\t', index=False)
    print(f"清洗完成: {output_path}")

# 执行清洗
for trait in ['oil', 'protein', 'size']:
    clean_for_smr(f'gwas_{trait}_results.assoc.txt', f'{trait}_smr_ready.txt')