import pandas as pd
import matplotlib
matplotlib.use('Agg') # 关键：强制不使用显示器窗口，只写文件
import matplotlib.pyplot as plt

# 确保读取当前目录下的文件
pca_data = pd.read_csv('soybean_pca.eigenvec', sep='\s+')

plt.figure(figsize=(10, 7))
plt.scatter(pca_data.iloc[:, 1], pca_data.iloc[:, 2], s=15, alpha=0.7, c='navy')
plt.title('PCA Plot - Soybean Population Structure')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.grid(True, linestyle='--', alpha=0.6)

plt.savefig('soybean_pca_result.png', dpi=300)
print("成功！请查看当前目录下的 soybean_pca_result.png")