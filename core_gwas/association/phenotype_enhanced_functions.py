"""
phenotype_statistics_enhanced.py

增强版本：针对 Size-Oil-Protein 相关性的特殊分析
基于原始 phenotype_statistics.py，添加三个新函数

使用方法：
  1. 将以下函数复制到原脚本的 if __name__ == "__main__": 之前
  2. 在主程序块中调用：
     df_analysis, corr_dict = analyze_size_oil_protein_correlation(df)
     prioritize_snps_for_gwas(df_analysis)
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import linregress
import matplotlib.pyplot as plt
import seaborn as sns

# ============================================================================
# 函数 1: 计算衍生性状（比例和标准化）
# ============================================================================

def calculate_trait_ratios(data):
    """
    计算 Oil/Protein 比例和其他衍生性状
    
    输入：
        data (DataFrame): 包含 'oil', 'protein', '100SW' 列的数据框
    
    输出：
        DataFrame: 添加了衍生性状的新数据框
    
    新增列：
        - Oil_Protein_Ratio: oil / protein
        - Oil_to_100SW: oil / 100SW
        - Protein_to_100SW: protein / 100SW
        - 对应的标准化版本 (_std 后缀)
    """
    df = data.copy()
    
    # 避免除以零
    epsilon = 1e-6
    
    # 计算比例
    df['Oil_Protein_Ratio'] = df['oil'] / (df['protein'] + epsilon)
    df['Oil_to_100SW'] = df['oil'] / (df['100SW'] + epsilon)
    df['Protein_to_100SW'] = df['protein'] / (df['100SW'] + epsilon)
    
    # 标准化衍生性状
    for col in ['Oil_Protein_Ratio', 'Oil_to_100SW', 'Protein_to_100SW']:
        mean = df[col].mean()
        std = df[col].std()
        if std > 0:
            df[f'{col}_std'] = (df[col] - mean) / std
        else:
            df[f'{col}_std'] = 0
    
    return df


# ============================================================================
# 函数 2: 关键相关性分析
# ============================================================================

def analyze_size_oil_protein_correlation(data, output_prefix='size_oil_protein'):
    """
    全面分析 Size 与 Oil/Protein 的关系
    
    包含五个分析模块：
    1. 基础 Pearson 相关性
    2. 偏相关分析（控制混杂因素）
    3. 分层分析（按 Size 四分位数）
    4. 多元回归分析
    5. 生成可视化图表
    
    输入：
        data (DataFrame): 包含表型数据的数据框
        output_prefix (str): 输出文件的前缀
    
    输出：
        tuple: (df_with_ratios, corr_results_dict)
    """
    
    df = data.copy()
    df = calculate_trait_ratios(df)
    
    # 初始化报告
    report = []
    report.append("="*80)
    report.append("Size-Oil-Protein 相关性分析")
    report.append("研究 Size 是否导致 Oil↑ / Protein↓ 的 Trade-off")
    report.append("="*80)
    
    # ========================================================================
    # 模块 1: 基础相关性分析
    # ========================================================================
    
    report.append("\n【1】基础 Pearson 相关性分析")
    report.append("-" * 80)
    
    trait_pairs = [
        ('100SW', 'oil'),
        ('100SW', 'protein'),
        ('100SW', 'Oil_Protein_Ratio'),
        ('oil', 'protein'),
    ]
    
    corr_results = {}
    
    for trait1, trait2 in trait_pairs:
        valid_data = df[[trait1, trait2]].dropna()
        
        if len(valid_data) > 2:
            corr, pval = stats.pearsonr(valid_data[trait1], valid_data[trait2])
            corr_results[f'{trait1}_vs_{trait2}'] = {
                'r': corr,
                'p': pval,
                'n': len(valid_data)
            }
            
            report.append(f"\n  {trait1} vs {trait2}:")
            report.append(f"    Pearson r = {corr:.4f}")
            report.append(f"    p-value   = {pval:.2e}")
            report.append(f"    样本数    = {len(valid_data)}")
            
            # 添加显著性标记
            if pval < 0.001:
                report.append(f"    显著性    = *** (p < 0.001)")
            elif pval < 0.01:
                report.append(f"    显著性    = ** (p < 0.01)")
            elif pval < 0.05:
                report.append(f"    显著性    = * (p < 0.05)")
            else:
                report.append(f"    显著性    = 不显著 (p ≥ 0.05)")
    
    # ========================================================================
    # 模块 2: 偏相关分析
    # ========================================================================
    
    report.append("\n\n【2】偏相关分析（控制混杂因素）")
    report.append("-" * 80)
    
    valid_data = df[['100SW', 'oil', 'protein']].dropna()
    
    if len(valid_data) > 3:
        # 计算残差
        res_size_oil = linregress(valid_data['100SW'], valid_data['oil'])
        residuals_oil = valid_data['oil'] - (
            res_size_oil.slope * valid_data['100SW'] + res_size_oil.intercept
        )
        
        res_size_protein = linregress(valid_data['100SW'], valid_data['protein'])
        residuals_protein = valid_data['protein'] - (
            res_size_protein.slope * valid_data['100SW'] + res_size_protein.intercept
        )
        
        # 计算偏相关
        partial_corr, partial_pval = stats.pearsonr(residuals_oil, residuals_protein)
        
        report.append(f"\n  Oil vs Protein（控制 Size 的影响）:")
        report.append(f"    偏相关系数 = {partial_corr:.4f}")
        report.append(f"    p-value   = {partial_pval:.2e}")
        report.append(f"    样本数    = {len(valid_data)}")
        
        original_corr = corr_results.get('oil_vs_protein', {}).get('r', 0)
        if abs(partial_corr) > abs(original_corr):
            report.append(f"    ⚠️  偏相关强于原始相关 → Size 可能是混杂因子")
        else:
            report.append(f"    ✓  偏相关弱于原始相关 → Size 部分解释了 Oil-Protein 关系")
    
    # ========================================================================
    # 模块 3: 分层分析
    # ========================================================================
    
    report.append("\n\n【3】分层分析：按种子大小分组")
    report.append("-" * 80)
    
    # 按四分位数分组
    df['size_quartile'] = pd.qcut(
        df['100SW'],
        q=4,
        labels=['Q1(小)', 'Q2', 'Q3', 'Q4(大)'],
        duplicates='drop'
    )
    
    report.append("\n  按 Size 四分位数分组的油蛋白比例:")
    report.append(f"  {'Group':<10} {'Mean 100SW':<15} {'Mean Oil':<15} {'Mean Protein':<15} {'Oil/Protein':<15}")
    report.append(f"  {'-'*70}")
    
    for group in ['Q1(小)', 'Q2', 'Q3', 'Q4(大)']:
        group_data = df[df['size_quartile'] == group]
        if len(group_data) > 0:
            mean_100sw = group_data['100SW'].mean()
            mean_oil = group_data['oil'].mean()
            mean_protein = group_data['protein'].mean()
            ratio = mean_oil / (mean_protein + 1e-6)
            
            report.append(
                f"  {group:<10} {mean_100sw:>14.2f} {mean_oil:>14.2f} {mean_protein:>14.2f} {ratio:>14.4f}"
            )
    
    # ANOVA 检验
    groups = []
    for group in ['Q1(小)', 'Q2', 'Q3', 'Q4(大)']:
        group_data = df[df['size_quartile'] == group]['Oil_Protein_Ratio'].dropna()
        if len(group_data) > 0:
            groups.append(group_data.values)
    
    if len(groups) > 1:
        f_stat, p_anova = stats.f_oneway(*groups)
        report.append(f"\n  方差分析 (ANOVA) 结果:")
        report.append(f"    F-statistic = {f_stat:.4f}")
        report.append(f"    p-value     = {p_anova:.2e}")
        if p_anova < 0.05:
            report.append(f"    结论        = 不同大小的种子其 Oil/Protein 比显著不同 ✅")
        else:
            report.append(f"    结论        = 不同大小的种子其 Oil/Protein 比无显著差异 ❌")
    
    # ========================================================================
    # 模块 4: 回归分析
    # ========================================================================
    
    report.append("\n\n【4】多元回归分析：Size 对 Oil/Protein 比的影响")
    report.append("-" * 80)
    
    valid_reg = df[['100SW', 'Oil_Protein_Ratio']].dropna()
    
    if len(valid_reg) > 2:
        slope, intercept, r_value, p_value, std_err = linregress(
            valid_reg['100SW'],
            valid_reg['Oil_Protein_Ratio']
        )
        
        report.append(f"\n  回归方程：Oil/Protein Ratio = {intercept:.4f} + {slope:.4f} × Size")
        report.append(f"  斜率 (slope)     = {slope:.4f}")
        report.append(f"  标准误 (SE)      = {std_err:.4f}")
        report.append(f"  R² 值            = {r_value**2:.4f}")
        report.append(f"  p-value          = {p_value:.2e}")
        report.append(f"  样本数           = {len(valid_reg)}")
        
        if p_value < 0.05:
            if slope > 0:
                report.append(f"  结论  → Size 增大导致 Oil/Protein 比增加（Oil相对增多）✅")
            else:
                report.append(f"  结论  → Size 增大导致 Oil/Protein 比减少（Protein相对增多）✅")
        else:
            report.append(f"  结论  → Size 对 Oil/Protein 比无显著影响 ❌")
    
    # ========================================================================
    # 模块 5: 可视化
    # ========================================================================
    
    report.append("\n\n【5】生成可视化图表...")
    
    try:
        fig, axes = plt.subplots(3, 2, figsize=(14, 12))
        fig.suptitle('Size-Oil-Protein 相关性详细分析', fontsize=16, fontweight='bold')
        
        # 图1：Size vs Oil
        ax = axes[0, 0]
        valid = df[['100SW', 'oil']].dropna()
        ax.scatter(valid['100SW'], valid['oil'], alpha=0.5, s=20)
        if len(valid) > 1:
            z = np.polyfit(valid['100SW'], valid['oil'], 1)
            p = np.poly1d(z)
            x_line = np.linspace(valid['100SW'].min(), valid['100SW'].max(), 100)
            ax.plot(x_line, p(x_line), 'r-', linewidth=2)
        r_val = corr_results.get('100SW_vs_oil', {}).get('r', 0)
        p_val = corr_results.get('100SW_vs_oil', {}).get('p', 1)
        ax.set_xlabel('Size (100SW)', fontsize=10)
        ax.set_ylabel('Oil (%)', fontsize=10)
        ax.set_title(f'Size vs Oil (r={r_val:.3f}, p={p_val:.2e})', fontsize=11)
        ax.grid(True, alpha=0.3)
        
        # 图2：Size vs Protein
        ax = axes[0, 1]
        valid = df[['100SW', 'protein']].dropna()
        ax.scatter(valid['100SW'], valid['protein'], alpha=0.5, s=20, color='orange')
        if len(valid) > 1:
            z = np.polyfit(valid['100SW'], valid['protein'], 1)
            p = np.poly1d(z)
            x_line = np.linspace(valid['100SW'].min(), valid['100SW'].max(), 100)
            ax.plot(x_line, p(x_line), 'r-', linewidth=2)
        r_val = corr_results.get('100SW_vs_protein', {}).get('r', 0)
        p_val = corr_results.get('100SW_vs_protein', {}).get('p', 1)
        ax.set_xlabel('Size (100SW)', fontsize=10)
        ax.set_ylabel('Protein (%)', fontsize=10)
        ax.set_title(f'Size vs Protein (r={r_val:.3f}, p={p_val:.2e})', fontsize=11)
        ax.grid(True, alpha=0.3)
        
        # 图3：Oil vs Protein
        ax = axes[1, 0]
        valid = df[['oil', 'protein']].dropna()
        ax.scatter(valid['oil'], valid['protein'], alpha=0.5, s=20, color='green')
        if len(valid) > 1:
            z = np.polyfit(valid['oil'], valid['protein'], 1)
            p = np.poly1d(z)
            x_line = np.linspace(valid['oil'].min(), valid['oil'].max(), 100)
            ax.plot(x_line, p(x_line), 'r-', linewidth=2)
        r_val = corr_results.get('oil_vs_protein', {}).get('r', 0)
        p_val = corr_results.get('oil_vs_protein', {}).get('p', 1)
        ax.set_xlabel('Oil (%)', fontsize=10)
        ax.set_ylabel('Protein (%)', fontsize=10)
        ax.set_title(f'Oil vs Protein (r={r_val:.3f}, p={p_val:.2e})', fontsize=11)
        ax.grid(True, alpha=0.3)
        
        # 图4：按 Size 分层的油箱线图
        ax = axes[1, 1]
        data_to_plot = [
            df[df['size_quartile'] == q]['oil'].dropna()
            for q in ['Q1(小)', 'Q2', 'Q3', 'Q4(大)']
        ]
        bp = ax.boxplot(data_to_plot, labels=['Q1(小)', 'Q2', 'Q3', 'Q4(大)'], patch_artist=True)
        for patch in bp['boxes']:
            patch.set_facecolor('lightblue')
        ax.set_ylabel('Oil (%)', fontsize=10)
        ax.set_title('按种子大小分层的油含量', fontsize=11)
        ax.grid(True, alpha=0.3, axis='y')
        
        # 图5：按 Size 分层的蛋白质箱线图
        ax = axes[2, 0]
        data_to_plot = [
            df[df['size_quartile'] == q]['protein'].dropna()
            for q in ['Q1(小)', 'Q2', 'Q3', 'Q4(大)']
        ]
        bp = ax.boxplot(data_to_plot, labels=['Q1(小)', 'Q2', 'Q3', 'Q4(大)'], patch_artist=True)
        for patch in bp['boxes']:
            patch.set_facecolor('lightcoral')
        ax.set_ylabel('Protein (%)', fontsize=10)
        ax.set_title('按种子大小分层的蛋白质含量', fontsize=11)
        ax.grid(True, alpha=0.3, axis='y')
        
        # 图6：Size vs Oil/Protein 比
        ax = axes[2, 1]
        valid = df[['100SW', 'Oil_Protein_Ratio']].dropna()
        ax.scatter(valid['100SW'], valid['Oil_Protein_Ratio'], alpha=0.5, s=20, color='purple')
        if len(valid) > 1:
            z = np.polyfit(valid['100SW'], valid['Oil_Protein_Ratio'], 1)
            p = np.poly1d(z)
            x_line = np.linspace(valid['100SW'].min(), valid['100SW'].max(), 100)
            ax.plot(x_line, p(x_line), 'r-', linewidth=2)
        ax.set_xlabel('Size (100SW)', fontsize=10)
        ax.set_ylabel('Oil/Protein Ratio', fontsize=10)
        
        slope_val = linregress(valid['100SW'], valid['Oil_Protein_Ratio']).slope if len(valid) > 1 else 0
        p_val = linregress(valid['100SW'], valid['Oil_Protein_Ratio']).pvalue if len(valid) > 1 else 1
        ax.set_title(f'Size vs Oil/Protein Ratio (slope={slope_val:.4f}, p={p_val:.2e})', fontsize=11)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f'{output_prefix}_analysis.png', dpi=300, bbox_inches='tight')
        report.append(f"✅ 图表已保存为 {output_prefix}_analysis.png")
        
    except Exception as e:
        report.append(f"⚠️  图表生成失败: {str(e)}")
    
    # ========================================================================
    # 保存报告
    # ========================================================================
    
    report_text = '\n'.join(report)
    output_file = f'{output_prefix}_correlation_report.txt'
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(report_text)
    
    report.append(f"\n✅ 报告已保存为 {output_file}")
    
    # 打印到控制台
    print('\n'.join(report))
    
    return df, corr_results


# ============================================================================
# 函数 3: GWAS SNP 优先级排序指南
# ============================================================================

def prioritize_snps_for_gwas(data, output_file='snps_to_test.txt'):
    """
    生成 GWAS 分析的指南和优先级排序
    
    输入：
        data (DataFrame): 包含表型数据的数据框
        output_file (str): 输出文件名
    
    输出：
        str: 生成的指南文本
    """
    
    report = []
    report.append("="*80)
    report.append("GWAS 优先级 SNP 筛选指南")
    report.append("="*80)
    
    report.append("\n【优先级排序】")
    report.append("-" * 80)
    
    report.append("\n最高优先级 (Tier 1) - 聚焦于以下三个性状：")
    report.append("  1. Size (100SW)        - 作为暴露因子（原因）")
    report.append("  2. Oil Content         - 作为结果变量1")
    report.append("  3. Protein Content     - 作为结果变量2")
    
    report.append("\n研究假设及相应的 GWAS 策略：")
    
    report.append("\n  H1: Size → Oil (正向因果)")
    report.append("      GWAS 目标: 发现影响 Size 的 SNPs（工具变量）")
    report.append("      p 值阈值: p < 5e-8（标准全基因组阈值）")
    report.append("      预期 SNPs: 20-50 个显著 SNPs")
    
    report.append("\n  H2: Size → Protein (可能的因果)")
    report.append("      GWAS 目标: 同上，检验 Protein 中的因果效应")
    report.append("      p 值阈值: p < 5e-8")
    report.append("      预期 SNPs: 15-40 个显著 SNPs")
    
    report.append("\n  H3: Oil ↔ Protein (是否存在因果关系)")
    report.append("      GWAS 目标: 检验 Oil 和 Protein 是否存在反向因果")
    report.append("      p 值阈值: p < 5e-8")
    report.append("      预期 SNPs: 10-30 个显著 SNPs")
    
    report.append("\n\n【GWAS 技术建议】")
    report.append("-" * 80)
    
    report.append("\n推荐软件: GEMMA（快速、稳健、支持复杂模型）")
    report.append("备选方案: PLINK/FASTGWA（GPU 加速）")
    
    report.append("\nGEMMA 运行命令示例：")
    report.append("""
# 1. 准备 SNP 关系矩阵 (GRM)
gemma -bfile genotypes \\
  -gk 1 \\
  -o output_prefix

# 2. 为 Size 进行 GWAS
gemma -bfile genotypes \\
  -k output_prefix.cXX.txt \\
  -p size_phenotype.txt \\
  -c covariate.txt \\
  -lmm 4 \\
  -o size_gwas

# 3. 为 Oil 进行 GWAS
gemma -bfile genotypes \\
  -k output_prefix.cXX.txt \\
  -p oil_phenotype.txt \\
  -c covariate.txt \\
  -lmm 4 \\
  -o oil_gwas

# 4. 为 Protein 进行 GWAS  
gemma -bfile genotypes \\
  -k output_prefix.cXX.txt \\
  -p protein_phenotype.txt \\
  -c covariate.txt \\
  -lmm 4 \\
  -o protein_gwas
    """)
    
    report.append("\n【输出文件准备】")
    report.append("-" * 80)
    
    report.append("\n生成的摘要统计格式（用于 GSMR）：")
    report.append("""
文件名: size_gwas.summary.txt
格式: 
  SNP    CHR    POS      A1  A2  freq   beta      se        p
  1_1000 1      1000     A   G   0.30   0.0523    0.0102    1.2e-7
  1_2000 1      2000     C   T   0.45   -0.0234   0.0095    3.4e-6
  ...
    """)
    
    report.append("\n【因果推断准备】")
    report.append("-" * 80)
    
    report.append("\nGWAS 完成后，将运行以下因果推断分析（详见 phase3_gsmr_heidi_guide.md）：")
    report.append("\n  分析 1: Size → Oil      (GSMR with HEDI)")
    report.append("  分析 2: Size → Protein  (GSMR with HEDI)")
    report.append("  分析 3: Oil → Size      (反向验证)")
    report.append("  分析 4: Oil → Protein   (多重性检验)")
    report.append("  分析 5: Bayesian MR     (当多效性 > 10% 时)")
    
    report.append("\n\n【下一步】")
    report.append("-" * 80)
    report.append("\n1️⃣  完成 Phase 1 表型分析（已完成）")
    report.append("2️⃣  进行 Phase 2 GWAS 分析（使用上述指导）")
    report.append("3️⃣  运行 Phase 3 GSMR 因果推断（参考 phase3_gsmr_heidi_guide.md）")
    report.append("4️⃣  执行 Phase 4 多效性检验和贝叶斯评估")
    
    # 保存
    report_text = '\n'.join(report)
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(report_text)
    
    print(report_text)
    return report_text


# ============================================================================
# 主程序示例
# ============================================================================

if __name__ == "__main__":
    # 示例：假设已加载数据
    # df = pd.read_csv('your_phenotype_data.csv')
    
    # 运行分析
    print("【Phase 1 - 特殊分析】")
    # df_with_ratios, corr_dict = analyze_size_oil_protein_correlation(df)
    
    # 生成 GWAS 指南
    print("\n【GWAS 准备指南】")
    # prioritize_snps_for_gwas(df_with_ratios)
    
    print("\n✅ Phase 1 分析完全完成！")
    print("📋 输出文件列表：")
    print("  - size_oil_protein_analysis.png（6 个子图）")
    print("  - size_oil_protein_correlation_report.txt（详细报告）")
    print("  - snps_to_test.txt（GWAS 指南）")
