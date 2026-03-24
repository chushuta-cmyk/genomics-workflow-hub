import pandas as pd

traits = ['size', 'protein', 'oil']

for t in traits:
    file_path = f'05_gsmr_correlation/{t}.ma'
    # 使用 sep='\s+' 自动处理各种空格/制表符混乱
    df = pd.read_csv(file_path, sep='\s+')
    
    # 强制统一列名，确保没有多余空格
    df.columns = ['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']
    
    # 强制将 SNP ID 转为字符串并去除两端空格
    df['SNP'] = df['SNP'].astype(str).str.strip()
    
    # 保存为标准的 Tab 分隔文件，quoting=3 确保不带引号
    output_path = f'05_gsmr_correlation/{t}_pure.ma'
    df.to_csv(output_path, sep='\t', index=False, quoting=3)
    print(f"✅ {t}_pure.ma 已生成")