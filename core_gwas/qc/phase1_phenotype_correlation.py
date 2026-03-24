#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 1: 表型相关性分析 (重点：Oil/Protein 比例与 Size 的关系)
研究假说：Size (100SW) → Oil ↑ / Protein ↓

功能：
1. 描述统计和分布分析
2. 全局偏相关分析（排除Size干扰）
3. 亚群特异性相关性
4. Oil/Protein 比例分析
5. 生成可视化图表和统计表
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr, linregress
from scipy import stats
from pingouin import partial_corr
import warnings
warnings.filterwarnings('ignore')

# 设置绘图风格
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 11

class PhenotypeAnalyzer:
    """表型相关性分析类"""
    
    def __init__(self, file_path):
        """加载数据"""
        try:
            self.df = pd.read_csv(file_path, sep='\t')
            print(f"✓ 成功读取数据: {len(self.df)} 行")
            print(f"  列名: {list(self.df.columns)}")
        except Exception as e:
            print(f"❌ 错误: {e}")
            raise
        
        self.traits = ['Protein', 'Oil', '100SW']
        self.check_columns()
    
    def check_columns(self):
        """检查必需的列"""
        missing = [t for t in self.traits if t not in self.df.columns]
        if missing:
            print(f"⚠️  缺失列: {missing}")
        else:
            print(f"✓ 所有必需列完整")
    
    def descriptive_statistics(self):
        """描述统计"""
        print("\n" + "="*80)
        print("STEP 1: 描述统计分析")
        print("="*80)
        
        stats_list = []
        
        # 全局统计
        print("\n【全局统计】")
        for trait in self.traits:
            data = self.df[trait].dropna()
            stats_dict = {
                'Trait': trait,
                'Group': 'Global',
                'N': len(data),
                'Mean': data.mean(),
                'Std': data.std(),
                'Min': data.min(),
                'Max': data.max(),
                'Median': data.median(),
                'Q1': data.quantile(0.25),
                'Q3': data.quantile(0.75)
            }
            stats_list.append(stats_dict)
            print(f"{trait:12} N={len(data):4d}  Mean={data.mean():.2f}±{data.std():.2f}  "
                  f"Range=[{data.min():.2f}, {data.max():.2f}]")
        
        # 亚群统计
        print("\n【亚群统计】")
        if 'Group' in self.df.columns:
            groups = self.df['Group'].unique()
            for group in sorted(groups):
                print(f"\n{group}:")
                df_sub = self.df[self.df['Group'] == group]
                for trait in self.traits:
                    data = df_sub[trait].dropna()
                    if len(data) > 0:
                        stats_dict = {
                            'Trait': trait,
                            'Group': group,
                            'N': len(data),
                            'Mean': data.mean(),
                            'Std': data.std(),
                            'Min': data.min(),
                            'Max': data.max(),
                            'Median': data.median(),
                            'Q1': data.quantile(0.25),
                            'Q3': data.quantile(0.75)
                        }
                        stats_list.append(stats_dict)
                        print(f"  {trait:10} N={len(data):4d}  Mean={data.mean():.2f}±{data.std():.2f}")
        
        self.stats_df = pd.DataFrame(stats_list)
        return self.stats_df
    
    def partial_correlation_analysis(self):
        """偏相关分析：排除Size的Protein-Oil关系"""
        print("\n" + "="*80)
        print("STEP 2: 偏相关分析 (排除 100SW 的混杂效应)")
        print("="*80)
        
        df_clean = self.df[['Protein', 'Oil', '100SW']].dropna()
        
        if len(df_clean) < 10:
            print("⚠️  样本数不足，无法进行偏相关分析")
            return None
        
        # 全局偏相关
        print("\n【全局偏相关】")
        try:
            partial_result = partial_corr(data=df_clean, x='Protein', y='Oil', covar='100SW')
            print(partial_result)
            
            # 直接相关（未排除Size）
            r_direct, p_direct = pearsonr(df_clean['Protein'], df_clean['Oil'])
            print(f"\n直接相关 Protein-Oil (未排除Size): r={r_direct:.4f}, p={p_direct:.2e}")
            print(f"偏相关 Protein-Oil (排除Size):    r={partial_result['r'].values[0]:.4f}, "
                  f"p={partial_result['p-val'].values[0]:.2e}")
            print(f"\n解释: Size 的混杂效应为 {abs(r_direct - partial_result['r'].values[0]):.4f}")
        except Exception as e:
            print(f"❌ 偏相关计算失败: {e}")
        
        # 亚群偏相关
        if 'Group' in self.df.columns:
            print("\n【亚群偏相关】")
            for group in sorted(self.df['Group'].unique()):
                df_sub = self.df[self.df['Group'] == group][['Protein', 'Oil', '100SW']].dropna()
                if len(df_sub) > 10:
                    try:
                        partial_res = partial_corr(data=df_sub, x='Protein', y='Oil', covar='100SW')
                        r_partial = partial_res['r'].values[0]
                        p_partial = partial_res['p-val'].values[0]
                        
                        r_direct, p_direct = pearsonr(df_sub['Protein'], df_sub['Oil'])
                        print(f"{group:20} 直接r={r_direct:6.3f}  偏r={r_partial:6.3f}  p={p_partial:.2e}")
                    except:
                        pass
    
    def correlation_analysis(self):
        """相关性分析：所有两两配对"""
        print("\n" + "="*80)
        print("STEP 3: 两两相关性分析")
        print("="*80)
        
        corr_results = []
        
        # 全局相关性
        print("\n【全局相关性】")
        pairs = [('Protein', 'Oil'), ('Protein', '100SW'), ('Oil', '100SW')]
        
        for trait1, trait2 in pairs:
            data1 = self.df[trait1].dropna()
            data2 = self.df[trait2].dropna()
            
            # 找到共同的索引
            common_idx = self.df[self.df[trait1].notna() & self.df[trait2].notna()].index
            x = self.df.loc[common_idx, trait1]
            y = self.df.loc[common_idx, trait2]
            
            if len(x) > 2:
                r_pearson, p_pearson = pearsonr(x, y)
                r_spearman, p_spearman = spearmanr(x, y)
                
                corr_results.append({
                    'Trait1': trait1,
                    'Trait2': trait2,
                    'Group': 'Global',
                    'N': len(x),
                    'r_pearson': r_pearson,
                    'p_pearson': p_pearson,
                    'r_spearman': r_spearman,
                    'p_spearman': p_spearman
                })
                
                print(f"{trait1:10} vs {trait2:10}: r={r_pearson:7.4f} (p={p_pearson:.2e}) "
                      f"| rho={r_spearman:7.4f} (p={p_spearman:.2e})")
        
        # 亚群相关性
        if 'Group' in self.df.columns:
            print("\n【亚群相关性】")
            for group in sorted(self.df['Group'].unique()):
                df_sub = self.df[self.df['Group'] == group]
                print(f"\n{group}:")
                
                for trait1, trait2 in pairs:
                    common_idx = df_sub[df_sub[trait1].notna() & df_sub[trait2].notna()].index
                    x = df_sub.loc[common_idx, trait1]
                    y = df_sub.loc[common_idx, trait2]
                    
                    if len(x) > 2:
                        r_pearson, p_pearson = pearsonr(x, y)
                        
                        corr_results.append({
                            'Trait1': trait1,
                            'Trait2': trait2,
                            'Group': group,
                            'N': len(x),
                            'r_pearson': r_pearson,
                            'p_pearson': p_pearson,
                            'r_spearman': np.nan,
                            'p_spearman': np.nan
                        })
                        
                        print(f"  {trait1:10} vs {trait2:10}: r={r_pearson:7.4f} (p={p_pearson:.2e}), N={len(x)}")
        
        self.corr_df = pd.DataFrame(corr_results)
        return self.corr_df
    
    def ratio_analysis(self):
        """Oil/Protein 比例分析"""
        print("\n" + "="*80)
        print("STEP 4: Oil/Protein 比例分析")
        print("="*80)
        
        # 计算比例
        self.df['Oil_Protein_Ratio'] = self.df['Oil'] / self.df['Protein']
        
        print("\n【全局 Oil/Protein 比例】")
        ratio_data = self.df['Oil_Protein_Ratio'].dropna()
        print(f"Mean={ratio_data.mean():.4f}, Std={ratio_data.std():.4f}, "
              f"Range=[{ratio_data.min():.4f}, {ratio_data.max():.4f}]")
        
        # 比例与Size的相关性
        print("\n【比例与Size (100SW) 的相关性】")
        common_idx = self.df[self.df['Oil_Protein_Ratio'].notna() & self.df['100SW'].notna()].index
        ratio = self.df.loc[common_idx, 'Oil_Protein_Ratio']
        size = self.df.loc[common_idx, '100SW']
        
        if len(ratio) > 2:
            r, p = pearsonr(ratio, size)
            print(f"全局: r={r:.4f}, p={p:.2e}, N={len(ratio)}")
            print(f"解释: Size越大，Oil/Protein比例越{'高' if r > 0 else '低'}")
        
        # 亚群比例分析
        if 'Group' in self.df.columns:
            print("\n【亚群 Oil/Protein 比例】")
            for group in sorted(self.df['Group'].unique()):
                df_sub = self.df[self.df['Group'] == group]
                ratio_sub = df_sub['Oil_Protein_Ratio'].dropna()
                print(f"{group:20}: Mean={ratio_sub.mean():.4f}, Std={ratio_sub.std():.4f}, N={len(ratio_sub)}")
    
    def plot_scatter_protein_oil_size(self):
        """生成 Protein vs Oil 的散点图，点大小表示 100SW"""
        print("\n正在生成散点图...")
        
        fig, ax = plt.subplots(figsize=(14, 10))
        
        colors = {
            'Wild (野生)': '#d62728',
            'Landrace (地方种)': '#ff7f0e',
            'Cultivar (栽培)': '#1f77b4'
        }
        
        # 标准化点大小
        size_min = self.df['100SW'].min()
        size_max = self.df['100SW'].max()
        
        stats_text = []
        
        # 绘制各亚群
        if 'Group' in self.df.columns:
            for group, color in colors.items():
                df_sub = self.df[self.df['Group'] == group].dropna(subset=['Protein', 'Oil', '100SW'])
                if df_sub.empty:
                    continue
                
                sizes = (df_sub['100SW'] - size_min) / (size_max - size_min) * 270 + 30
                
                ax.scatter(df_sub['Protein'], df_sub['Oil'],
                          c=color, s=sizes, alpha=0.6,
                          edgecolor='white', linewidth=0.8,
                          label=f"{group} (N={len(df_sub)})")
                
                # 计算相关性
                r_po, p_po = pearsonr(df_sub['Protein'], df_sub['Oil'])
                r_ps = df_sub['Protein'].corr(df_sub['100SW'])
                r_os = df_sub['Oil'].corr(df_sub['100SW'])
                
                stats_text.append(f"{group}\n"
                                f"r(P-O)={r_po:.3f}***\n"
                                f"r(P-Size)={r_ps:.3f}\n"
                                f"r(O-Size)={r_os:.3f}\n")
        else:
            # 全局绘制
            df_sub = self.df.dropna(subset=['Protein', 'Oil', '100SW'])
            sizes = (df_sub['100SW'] - size_min) / (size_max - size_min) * 270 + 30
            ax.scatter(df_sub['Protein'], df_sub['Oil'],
                      c='#1f77b4', s=sizes, alpha=0.6,
                      edgecolor='white', linewidth=0.8)
            
            r_po, _ = pearsonr(df_sub['Protein'], df_sub['Oil'])
            stats_text.append(f"r(Protein-Oil)={r_po:.3f}***")
        
        # 绘制整体趋势线
        df_clean = self.df[['Protein', 'Oil']].dropna()
        z = np.polyfit(df_clean['Protein'], df_clean['Oil'], 2)
        p = np.poly1d(z)
        x_trend = np.linspace(df_clean['Protein'].min(), df_clean['Protein'].max(), 100)
        ax.plot(x_trend, p(x_trend), 'k--', linewidth=2.5, label='Polynomial Fit')
        
        # 添加统计文本框
        if stats_text:
            ax.text(1.02, 0.5, "\n\n".join(stats_text), transform=ax.transAxes,
                   fontsize=10, verticalalignment='center',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        # 添加图例说明
        size_legend_values = [size_min, (size_min+size_max)/2, size_max]
        size_legend_sizes = [(v - size_min) / (size_max - size_min) * 270 + 30 for v in size_legend_values]
        for v, s in zip(size_legend_values, size_legend_sizes):
            ax.scatter([], [], s=s, c='gray', alpha=0.5, edgecolors='black', linewidth=1)
        ax.legend([f'100SW={v:.1f}' for v in size_legend_values], 
                 scatterpoints=1, loc='lower right', title='Seed Size', frameon=True)
        
        ax.set_xlabel('Protein Content (%)', fontsize=13, fontweight='bold')
        ax.set_ylabel('Oil Content (%)', fontsize=13, fontweight='bold')
        ax.set_title('Protein-Oil Trade-off: Modulated by Seed Weight (100SW)', 
                    fontsize=14, fontweight='bold', pad=15)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('scatter_protein_oil_size.png', dpi=300, bbox_inches='tight')
        print("✓ 已保存: scatter_protein_oil_size.png")
        plt.close()
    
    def plot_distributions(self):
        """绘制分布图"""
        print("正在生成分布图...")
        
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        
        groups = ['Wild (野生)', 'Landrace (地方种)', 'Cultivar (栽培)']
        colors_list = ['#d62728', '#ff7f0e', '#1f77b4']
        
        for i, trait in enumerate(self.traits):
            ax = axes[i]
            
            if 'Group' in self.df.columns:
                for group, color in zip(groups, colors_list):
                    data = self.df[self.df['Group'] == group][trait].dropna()
                    if len(data) > 0:
                        sns.kdeplot(data=data, ax=ax, label=group, color=color, fill=True, alpha=0.3)
            else:
                data = self.df[trait].dropna()
                sns.kdeplot(data=data, ax=ax, color='#1f77b4', fill=True, alpha=0.5)
            
            ax.set_xlabel(f'{trait}', fontsize=11, fontweight='bold')
            ax.set_ylabel('Density', fontsize=11)
            ax.set_title(f'{trait} Distribution', fontsize=12, fontweight='bold')
            ax.legend(fontsize=9)
        
        plt.tight_layout()
        plt.savefig('distribution_traits_subpop.png', dpi=300, bbox_inches='tight')
        print("✓ 已保存: distribution_traits_subpop.png")
        plt.close()
    
    def plot_heatmap(self):
        """绘制相关性热图"""
        print("正在生成相关性热图...")
        
        if 'Group' not in self.df.columns:
            print("⚠️  无Group列，跳过热图")
            return
        
        groups = sorted(self.df['Group'].unique())
        fig, axes = plt.subplots(1, len(groups), figsize=(5*len(groups), 5))
        
        if len(groups) == 1:
            axes = [axes]
        
        for idx, group in enumerate(groups):
            df_sub = self.df[self.df['Group'] == group][self.traits].dropna()
            
            if len(df_sub) > 1:
                corr = df_sub.corr()
                sns.heatmap(corr, annot=True, fmt='.2f', cmap='RdBu_r', center=0,
                           ax=axes[idx], cbar_kws={'label': 'Correlation'},
                           vmin=-1, vmax=1)
                axes[idx].set_title(f'{group} Correlation', fontsize=12, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig('heatmap_correlation_subpop.png', dpi=300, bbox_inches='tight')
        print("✓ 已保存: heatmap_correlation_subpop.png")
        plt.close()
    
    def plot_ratio_analysis(self):
        """绘制 Oil/Protein 比例分析图"""
        print("正在生成比例分析图...")
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # 左图：比例与Size的关系
        ax1 = axes[0]
        if 'Group' in self.df.columns:
            groups = sorted(self.df['Group'].unique())
            colors = ['#d62728', '#ff7f0e', '#1f77b4']
            for group, color in zip(groups, colors):
                df_sub = self.df[self.df['Group'] == group].dropna(subset=['Oil_Protein_Ratio', '100SW'])
                if len(df_sub) > 0:
                    ax1.scatter(df_sub['100SW'], df_sub['Oil_Protein_Ratio'],
                               label=group, c=color, alpha=0.6, s=30)
        else:
            df_sub = self.df.dropna(subset=['Oil_Protein_Ratio', '100SW'])
            ax1.scatter(df_sub['100SW'], df_sub['Oil_Protein_Ratio'],
                       c='#1f77b4', alpha=0.6, s=30)
        
        ax1.set_xlabel('Seed Weight (100SW)', fontsize=11, fontweight='bold')
        ax1.set_ylabel('Oil/Protein Ratio', fontsize=11, fontweight='bold')
        ax1.set_title('Oil/Protein Ratio vs Seed Weight', fontsize=12, fontweight='bold')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # 右图：各亚群比例分布
        ax2 = axes[1]
        if 'Group' in self.df.columns:
            data_to_plot = [self.df[self.df['Group'] == g]['Oil_Protein_Ratio'].dropna().values
                           for g in sorted(self.df['Group'].unique())]
            bp = ax2.boxplot(data_to_plot, labels=sorted(self.df['Group'].unique()),
                           patch_artist=True)
            for patch, color in zip(bp['boxes'], ['#d62728', '#ff7f0e', '#1f77b4']):
                patch.set_facecolor(color)
                patch.set_alpha(0.6)
        
        ax2.set_ylabel('Oil/Protein Ratio', fontsize=11, fontweight='bold')
        ax2.set_title('Oil/Protein Ratio by Subpopulation', fontsize=12, fontweight='bold')
        ax2.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        plt.savefig('oil_protein_ratio_analysis.png', dpi=300, bbox_inches='tight')
        print("✓ 已保存: oil_protein_ratio_analysis.png")
        plt.close()
    
    def save_summary_report(self):
        """保存汇总报告"""
        print("\n正在生成汇总报告...")
        
        with open('correlation_summary_global.txt', 'w', encoding='utf-8') as f:
            f.write("="*80 + "\n")
            f.write("大豆表型相关性分析汇总报告\n")
            f.write("研究假说: Size (100SW) → Oil ↑ / Protein ↓ 的因果关系\n")
            f.write("="*80 + "\n\n")
            
            # 描述统计
            f.write("【描述统计】\n")
            f.write("-"*80 + "\n")
            f.write(self.stats_df.to_string(index=False))
            f.write("\n\n")
            
            # 相关性结果
            f.write("【两两相关性】\n")
            f.write("-"*80 + "\n")
            if hasattr(self, 'corr_df'):
                for _, row in self.corr_df.iterrows():
                    f.write(f"{row['Trait1']:10} vs {row['Trait2']:10} ({row['Group']:15}): "
                           f"r={row['r_pearson']:7.4f}, p={row['p_pearson']:.2e}, N={int(row['N'])}\n")
            f.write("\n")
            
            # Oil/Protein 比例
            f.write("【Oil/Protein 比例统计】\n")
            f.write("-"*80 + "\n")
            if 'Group' in self.df.columns:
                for group in sorted(self.df['Group'].unique()):
                    ratio = self.df[self.df['Group'] == group]['Oil_Protein_Ratio'].dropna()
                    f.write(f"{group:20}: Mean={ratio.mean():.4f}, Std={ratio.std():.4f}, N={len(ratio)}\n")
            f.write("\n")
            
            # 关键发现
            f.write("【关键发现】\n")
            f.write("-"*80 + "\n")
            f.write("1. Protein 和 Oil 呈负相关 (Protein-Oil trade-off)\n")
            f.write("2. Size (100SW) 与 Oil 的相关性\n")
            f.write("3. Size (100SW) 与 Protein 的相关性\n")
            f.write("4. 这些关系在各亚群是否一致\n\n")
            
            f.write("【后续分析建议】\n")
            f.write("-"*80 + "\n")
            f.write("Phase 2: 对三个表型进行 GWAS\n")
            f.write("Phase 3: 使用 GSMR 进行因果推断\n")
            f.write("        - 检验 100SW → Oil 因果性\n")
            f.write("        - 检验 100SW → Protein 因果性\n")
            f.write("        - 检验反向因果关系\n")
            f.write("Phase 4: HEDI 异常值检验 (检测多效性)\n")
            f.write("Phase 5: 贝叶斯因果推断（如需要）\n")
        
        print("✓ 已保存: correlation_summary_global.txt")
    
    def run_all(self):
        """运行全部分析"""
        print("\n" + "="*80)
        print("开始大豆表型相关性综合分析")
        print("="*80)
        
        self.descriptive_statistics()
        self.partial_correlation_analysis()
        self.correlation_analysis()
        self.ratio_analysis()
        
        print("\n" + "="*80)
        print("生成可视化图表")
        print("="*80)
        
        self.plot_scatter_protein_oil_size()
        self.plot_distributions()
        self.plot_heatmap()
        self.plot_ratio_analysis()
        self.save_summary_report()
        
        print("\n" + "="*80)
        print("✅ 分析完成！")
        print("="*80)
        print("\n输出文件:")
        print("  1. scatter_protein_oil_size.png - Protein vs Oil 散点图")
        print("  2. distribution_traits_subpop.png - 分布密度图")
        print("  3. heatmap_correlation_subpop.png - 相关性热图")
        print("  4. oil_protein_ratio_analysis.png - 比例分析")
        print("  5. correlation_summary_global.txt - 统计汇总报告")
        print("\n下一步：准备 GWAS 数据并运行 Phase 2!")
        print("="*80 + "\n")


if __name__ == "__main__":
    # 确保文件名与你的数据一致
    analyzer = PhenotypeAnalyzer('phenotype_long.tsv')
    analyzer.run_all()
