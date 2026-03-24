#!/usr/bin/env python3
"""
GSMR结果可视化脚本
功能: 生成因果效应图和综合报告
作者: GWAS Pipeline
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from pathlib import Path

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 300

class GSMRVisualizer:
    def __init__(self, result_dir='05_gsmr_correlation/gsmr_results'):
        self.result_dir = Path(result_dir)
        self.output_dir = self.result_dir / 'plots'
        self.output_dir.mkdir(exist_ok=True, parents=True)
        
    def parse_gsmr_results(self):
        """解析GSMR结果文件"""
        results = []
        
        # 分析配置
        analyses = [
            ('size_to_protein', 'Size → Protein'),
            ('size_to_oil', 'Size → Oil')
        ]
        
        for file_prefix, label in analyses:
            gsmr_file = self.result_dir / f'{file_prefix}.gsmr'
            
            if not gsmr_file.exists():
                print(f"⚠️  未找到: {gsmr_file}")
                continue
                
            try:
                df = pd.read_csv(gsmr_file, sep=r'\s+')
                if not df.empty:
                    row = df.iloc[0]
                    results.append({
                        'Analysis': label,
                        'Exposure': 'Size',
                        'Outcome': label.split(' → ')[1],
                        'N_SNPs': row.get('nsnp', row.get('n_SNPs', 0)),
                        'Beta': float(row.get('bxy', 0)),
                        'SE': float(row.get('se', 0)),
                        'P_value': float(row.get('p', 1))
                    })
            except Exception as e:
                print(f"❌ 解析{file_prefix}失败: {e}")
                
        return pd.DataFrame(results)
    
    def plot_causal_effects(self, df):
        """绘制因果效应森林图"""
        if df.empty:
            print("无有效数据，跳过可视化")
            return
            
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # 计算置信区间
        df['CI_lower'] = df['Beta'] - 1.96 * df['SE']
        df['CI_upper'] = df['Beta'] + 1.96 * df['SE']
        
        # 绘制误差条
        y_pos = np.arange(len(df))
        colors = ['#2E86AB' if p < 0.05 else '#A23B72' for p in df['P_value']]
        
        ax.errorbar(df['Beta'], y_pos, 
                   xerr=1.96 * df['SE'],
                   fmt='o', 
                   markersize=10,
                   capsize=8,
                   capthick=2,
                   ecolor=colors,
                   markerfacecolor=colors,
                   markeredgecolor='black',
                   markeredgewidth=1.5,
                   linewidth=2,
                   alpha=0.85)
        
        # 添加零线
        ax.axvline(x=0, color='#F18F01', linestyle='--', linewidth=2, alpha=0.7)
        
        # 设置坐标轴
        ax.set_yticks(y_pos)
        ax.set_yticklabels(df['Analysis'], fontsize=12, fontweight='bold')
        ax.set_xlabel('Causal Effect (β ± 95% CI)', fontsize=13, fontweight='bold')
        ax.set_title('GSMR Causal Inference Results\nSize → Protein/Oil', 
                    fontsize=15, fontweight='bold', pad=20)
        
        # 添加P值标注
        for i, row in df.iterrows():
            p_str = f"p={row['P_value']:.2e}" if row['P_value'] < 0.001 else f"p={row['P_value']:.3f}"
            sig_marker = '***' if row['P_value'] < 0.001 else ('**' if row['P_value'] < 0.01 else ('*' if row['P_value'] < 0.05 else 'NS'))
            
            ax.text(row['Beta'] + row['SE'] * 2.5, i, 
                   f"{p_str}\n{sig_marker}", 
                   va='center', 
                   fontsize=10,
                   bbox=dict(boxstyle='round,pad=0.5', 
                            facecolor='white', 
                            edgecolor='gray', 
                            alpha=0.8))
        
        # 网格优化
        ax.grid(True, axis='x', linestyle=':', alpha=0.5, linewidth=1.2)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        
        # 保存
        output_file = self.output_dir / 'GSMR_Causal_Effects.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ 森林图已保存: {output_file}")
        plt.close()
    
    def plot_snp_usage(self, df):
        """绘制工具变量SNP数量"""
        if df.empty:
            return
            
        fig, ax = plt.subplots(figsize=(8, 5))
        
        colors = sns.color_palette("husl", len(df))
        bars = ax.bar(df['Outcome'], df['N_SNPs'], 
                     color=colors, 
                     edgecolor='black', 
                     linewidth=1.5,
                     alpha=0.8)
        
        # 添加数值标签
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{int(height)}',
                   ha='center', va='bottom',
                   fontsize=12, fontweight='bold')
        
        ax.set_ylabel('Number of Instrumental SNPs', fontsize=12, fontweight='bold')
        ax.set_xlabel('Outcome Trait', fontsize=12, fontweight='bold')
        ax.set_title('GSMR Instrumental SNP Usage', fontsize=14, fontweight='bold')
        ax.grid(axis='y', linestyle='--', alpha=0.4)
        
        plt.tight_layout()
        output_file = self.output_dir / 'GSMR_SNP_Usage.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ SNP使用量图已保存: {output_file}")
        plt.close()
    
    def generate_report(self, df):
        """生成Markdown格式报告"""
        report = []
        report.append("# GSMR 因果推断分析报告\n")
        report.append("## 分析概述\n")
        report.append("**暴露性状**: Size (百粒重)")
        report.append("**结局性状**: Protein (蛋白质含量), Oil (脂肪含量)\n")
        report.append("---\n")
        
        report.append("## 主要结果\n")
        
        if not df.empty:
            report.append("| 分析 | 工具SNP数 | 因果效应 (β) | 标准误 (SE) | P值 | 显著性 |")
            report.append("|------|-----------|--------------|-------------|-----|--------|")
            
            for _, row in df.iterrows():
                sig = '***' if row['P_value'] < 0.001 else ('**' if row['P_value'] < 0.01 else ('*' if row['P_value'] < 0.05 else 'NS'))
                report.append(
                    f"| {row['Analysis']} | {int(row['N_SNPs'])} | "
                    f"{row['Beta']:.4f} | {row['SE']:.4f} | "
                    f"{row['P_value']:.2e} | {sig} |"
                )
        else:
            report.append("⚠️ 无有效结果数据\n")
        
        report.append("\n---\n")
        report.append("## 结果解读\n")
        
        for _, row in df.iterrows():
            report.append(f"### {row['Analysis']}")
            
            if row['P_value'] < 0.05:
                direction = "正向" if row['Beta'] > 0 else "负向"
                report.append(f"- ✓ **显著因果关系** (p={row['P_value']:.2e})")
                report.append(f"- 方向: {direction} ({row['Beta']:.4f} ± {row['SE']:.4f})")
                report.append(f"- 工具变量: {int(row['N_SNPs'])} 个有效SNP\n")
            else:
                report.append(f"- ⚠️ **无显著因果关系** (p={row['P_value']:.3f})")
                report.append(f"- 效应估计: {row['Beta']:.4f} ± {row['SE']:.4f}")
                report.append(f"- 工具变量: {int(row['N_SNPs'])} 个有效SNP\n")
        
        report.append("---\n")
        report.append("## 生物学意义\n")
        report.append("- **P < 0.05**: 暴露性状对结局性状存在因果效应")
        report.append("- **β > 0**: 暴露性状增加导致结局性状增加")
        report.append("- **β < 0**: 暴露性状增加导致结局性状减少\n")
        
        report_text = "\n".join(report)
        
        # 保存Markdown报告
        report_file = self.result_dir / 'GSMR_Analysis_Report.md'
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write(report_text)
        
        print(f"✓ Markdown报告已保存: {report_file}")
        
        # 同时打印到终端
        print("\n" + "="*70)
        print(report_text)
        print("="*70)
    
    def run(self):
        """执行完整可视化流程"""
        print("\n" + "="*70)
        print("GSMR结果可视化工具")
        print("="*70 + "\n")
        
        # 解析结果
        print("[1/4] 解析GSMR结果...")
        df = self.parse_gsmr_results()
        
        if df.empty:
            print("❌ 未找到有效的GSMR结果，请先运行 run_gsmr_analysis.sh")
            return
        
        print(f"✓ 成功解析 {len(df)} 个分析结果\n")
        
        # 绘制图表
        print("[2/4] 生成因果效应森林图...")
        self.plot_causal_effects(df)
        
        print("[3/4] 生成SNP使用量图...")
        self.plot_snp_usage(df)
        
        # 生成报告
        print("[4/4] 生成综合报告...")
        self.generate_report(df)
        
        print("\n" + "="*70)
        print("✓ 可视化完成!")
        print(f"📁 输出目录: {self.output_dir}")
        print("="*70 + "\n")

if __name__ == '__main__':
    visualizer = GSMRVisualizer()
    visualizer.run()
