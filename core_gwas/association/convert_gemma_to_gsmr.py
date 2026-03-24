
import pandas as pd
import sys
import os

def convert_gemma_to_gsmr(input_file, output_file, sample_size):
    print(f"Reading: {input_file}")
    df = pd.read_csv(input_file, sep='\t')
    
    # 映射 GEMMA 列名到 GSMR 标准格式
    gsmr_df = pd.DataFrame({
        'SNP': df['rs'].astype(str),
        'A1': df['allele1'].astype(str),
        'A2': df['allele0'].astype(str),
        'freq': df['af'].astype(float),
        'beta': df['beta'].astype(float),
        'se': df['se'].astype(float),
        'p': df['p_wald'].astype(float),
        'n': int(sample_size)
    })
    
    # 极简过滤：移除空值和异常 SE/P
    gsmr_df = gsmr_df.dropna()
    gsmr_df = gsmr_df[(gsmr_df['se'] > 0) & (gsmr_df['p'] >= 0) & (gsmr_df['p'] <= 1)]
    
    gsmr_df.to_csv(output_file, sep='\t', index=False)
    print(f"✓ Saved to: {output_file} (Rows: {len(gsmr_df)})")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 convert_gemma_to_gsmr.py <input> <output> <sample_size>")
        sys.exit(1)
    
    # 关键修正：添加索引 [1], [2], [3]
    input_f, output_f, n = sys.argv[1], sys.argv[2], sys.argv[3]
    convert_gemma_to_gsmr(input_f, output_f, n)