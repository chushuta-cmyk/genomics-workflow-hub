import pandas as pd
import os

files = ['gwas_protein_results.assoc.txt', 'gwas_oil_results.assoc.txt', 'gwas_size_results.assoc.txt']
output_suffix = '.clean.txt'

# 定义大豆标准的 20 条染色体
main_chrs = [str(i) for i in range(1, 21)]

for f in files:
    if not os.path.exists(f):
        continue
    print("正在清理文件: {} ...".format(f))
    
    # 使用 chunksize 处理大文件，节省内存
    reader = pd.read_csv(f, sep='\t', chunksize=100000)
    first_chunk = True
    
    out_file = f.replace('.assoc.txt', output_suffix)
    
    for chunk in reader:
        # 强制转换第一列为字符串并过滤
        chunk['chr'] = chunk['chr'].astype(str)
        clean_chunk = chunk[chunk['chr'].isin(main_chrs)]
        
        # 写入新文件
        mode = 'w' if first_chunk else 'a'
        header = True if first_chunk else False
        clean_chunk.to_csv(out_file, sep='\t', index=False, mode=mode, header=header)
        first_chunk = False
        
    print("清理完成，已保存至: {}".format(out_file))

print("\n现在你可以使用这些 .clean.txt 文件重新跑曼哈顿图，给董老师看一个清爽的结果。")