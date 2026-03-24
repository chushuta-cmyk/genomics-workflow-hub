# Generalized from soybean project-specific path layout.
#!/usr/bin/env python3
"""
Generate comparison summary plots.
"""

import numpy as np
import matplotlib.pyplot as plt
import os

def create_comparison_plot1(output_file):
    """Create multi-panel summary plot (Format 1)."""
    
    # Create figure with 2x2 grid
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('GWAS and MR Analysis Summary - Format 1', fontsize=16, y=0.98)
    
    # Panel 1: Manhattan plot highlights (simplified)
    ax1 = axes[0, 0]
    traits = ['100SW', 'Protein', 'Oil', 'Ratio']
    max_logp = [12.4, 8.7, 9.2, 11.8]  # From trace log
    
    bars = ax1.bar(traits, max_logp, color=['#4C72B0', '#55A868', '#C44E52', '#8172B3'])
    ax1.axhline(y=-np.log10(5e-8), color='red', linestyle='--', label='Genome-wide sig.')
    ax1.set_ylabel('Max -log10(P)', fontsize=11)
    ax1.set_title('GWAS Peak Associations', fontsize=12)
    ax1.legend(fontsize=9)
    
    # Add value labels
    for bar, val in zip(bars, max_logp):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                f'{val:.1f}', ha='center', va='bottom', fontsize=9)
    
    # Panel 2: MR effect sizes
    ax2 = axes[0, 1]
    comparisons = ['Ratio→Size', 'Size→Ratio', 'Size→Protein', 'Size→Oil']
    beta_values = [2.0956, 0.3788, -0.124, 0.087]  # From trace log
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4']
    
    bars2 = ax2.bar(comparisons, beta_values, color=colors)
    ax2.axhline(y=0, color='black', linewidth=0.8)
    ax2.set_ylabel('MR Effect Size (β)', fontsize=11)
    ax2.set_title('Mendelian Randomization Results', fontsize=12)
    ax2.tick_params(axis='x', rotation=45)
    
    # Add value labels
    for bar, val in zip(bars2, beta_values):
        height = bar.get_height()
        va = 'bottom' if height >= 0 else 'top'
        y_offset = 0.02 if height >= 0 else -0.02
        ax2.text(bar.get_x() + bar.get_width()/2., height + y_offset,
                f'{val:.3f}', ha='center', va=va, fontsize=9)
    
    # Panel 3: Genetic correlation heatmap (simplified)
    ax3 = axes[1, 0]
    traits_small = ['Size', 'Ratio']
    corr_matrix = np.array([
        [1.0, 0.577],  # sqrt(0.334) = 0.577
        [0.577, 1.0]
    ])
    
    im = ax3.imshow(corr_matrix, cmap='RdBu_r', vmin=0, vmax=1)
    ax3.set_xticks([0, 1])
    ax3.set_yticks([0, 1])
    ax3.set_xticklabels(traits_small)
    ax3.set_yticklabels(traits_small)
    ax3.set_title('Genetic Correlation (Size vs Ratio)', fontsize=12)
    
    # Add correlation values
    for i in range(2):
        for j in range(2):
            ax3.text(j, i, f'{corr_matrix[i, j]:.3f}',
                    ha='center', va='center', color='white', fontsize=10)
    
    # Panel 4: Asymmetry visualization
    ax4 = axes[1, 1]
    directions = ['Forward\n(Ratio→Size)', 'Reverse\n(Size→Ratio)']
    effect_sizes = [2.0956, 0.3788]
    asymmetry_ratio = effect_sizes[0] / effect_sizes[1]
    
    bars4 = ax4.bar(directions, effect_sizes, color=['#FF6B6B', '#4ECDC4'])
    ax4.set_ylabel('Effect Size (β)', fontsize=11)
    ax4.set_title(f'Causal Asymmetry\nRatio = {asymmetry_ratio:.1f}x', fontsize=12)
    
    # Add value labels and asymmetry indicator
    for bar, val in zip(bars4, effect_sizes):
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., height + 0.05,
                f'{val:.3f}', ha='center', va='bottom', fontsize=9)
    
    # Add asymmetry arrow
    ax4.annotate('', xy=(1, effect_sizes[1]), xytext=(0, effect_sizes[0]),
                arrowprops=dict(arrowstyle='<->', color='black', lw=1.5))
    ax4.text(0.5, (effect_sizes[0] + effect_sizes[1])/2, 
            f'{asymmetry_ratio:.1f}x', ha='center', va='center',
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved comparison plot 1: {output_file}")

def create_comparison_plot2(output_file):
    """Create MR results comparison plot (Format 2)."""
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('MR Method Consistency Comparison - Format 2', fontsize=16, y=0.98)
    
    # Left panel: Ratio → Size comparison
    methods = ['IVW', 'Weighted\nMedian', 'MR-Egger']
    beta_ratio_to_size = [2.0956, 2.1023, 2.0887]  # From trace log
    se_ratio_to_size = [0.0571, 0.0623, 0.0854]
    
    x = np.arange(len(methods))
    ax1.bar(x, beta_ratio_to_size, yerr=se_ratio_to_size, 
            capsize=5, color='#4C72B0', alpha=0.7)
    ax1.set_xticks(x)
    ax1.set_xticklabels(methods)
    ax1.set_ylabel('Effect Size (β)', fontsize=11)
    ax1.set_title('Ratio → 100SW\n(Bidirectional MR - Forward)', fontsize=12)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Add value labels
    for i, (beta, se) in enumerate(zip(beta_ratio_to_size, se_ratio_to_size)):
        ax1.text(i, beta + se + 0.05, f'{beta:.3f}±{se:.3f}', 
                ha='center', va='bottom', fontsize=9)
    
    # Right panel: Size → Ratio comparison
    beta_size_to_ratio = [0.3788, 0.3821, 0.3754]
    se_size_to_ratio = [0.0098, 0.0112, 0.0156]
    
    ax2.bar(x, beta_size_to_ratio, yerr=se_size_to_ratio,
            capsize=5, color='#55A868', alpha=0.7)
    ax2.set_xticks(x)
    ax2.set_xticklabels(methods)
    ax2.set_ylabel('Effect Size (β)', fontsize=11)
    ax2.set_title('100SW → Ratio\n(Bidirectional MR - Reverse)', fontsize=12)
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Add value labels
    for i, (beta, se) in enumerate(zip(beta_size_to_ratio, se_size_to_ratio)):
        ax2.text(i, beta + se + 0.005, f'{beta:.3f}±{se:.3f}', 
                ha='center', va='bottom', fontsize=9)
    
    # Add overall comparison text
    asymmetry = beta_ratio_to_size[0] / beta_size_to_ratio[0]
    fig.text(0.5, 0.02, 
             f'Asymmetry Ratio (Forward/Reverse): {asymmetry:.1f}x\n' +
             'Interpretation: Ratio has stronger causal effect on Size than vice versa',
             ha='center', fontsize=10, style='italic',
             bbox=dict(boxstyle="round,pad=0.5", facecolor="lightyellow", alpha=0.7))
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved comparison plot 2: {output_file}")

def create_comparison_plot3(output_file):
    """Create genetic architecture visualization (Format 3)."""
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Create network visualization
    traits = ['100SW', 'Protein', 'Oil', 'Oil:Protein\nRatio']
    positions = {
        '100SW': (0, 1),
        'Protein': (-1, 0),
        'Oil': (1, 0),
        'Oil:Protein\nRatio': (0, -1)
    }
    
    # Draw nodes
    node_colors = {'100SW': '#4C72B0', 'Protein': '#55A868', 
                   'Oil': '#C44E52', 'Oil:Protein\nRatio': '#8172B3'}
    
    for trait, (x, y) in positions.items():
        ax.scatter(x, y, s=1200, color=node_colors[trait], alpha=0.7, edgecolors='black', linewidth=2)
        ax.text(x, y, trait, ha='center', va='center', fontsize=11, fontweight='bold')
    
    # Draw edges with weights based on effect sizes
    edges = [
        ('Oil:Protein\nRatio', '100SW', 2.0956, '#FF6B6B', 'bold'),  # Strong effect
        ('100SW', 'Oil:Protein\nRatio', 0.3788, '#4ECDC4', 'normal'),  # Weaker effect
        ('100SW', 'Protein', -0.124, '#45B7D1', 'normal'),  # Negative effect
        ('100SW', 'Oil', 0.087, '#96CEB4', 'normal'),  # Positive effect
    ]
    
    for start, end, weight, color, style in edges:
        x1, y1 = positions[start]
        x2, y2 = positions[end]
        
        # Draw arrow (direction matters)
        dx, dy = x2 - x1, y2 - y1
        arrow_length = np.sqrt(dx**2 + dy**2)
        
        # Normalize for arrowhead
        if arrow_length > 0:
            dx_norm, dy_norm = dx/arrow_length, dy/arrow_length
            
            # Adjust start and end points to avoid overlapping nodes
            start_adj = (x1 + dx_norm * 0.15, y1 + dy_norm * 0.15)
            end_adj = (x2 - dx_norm * 0.15, y2 - dy_norm * 0.15)
            
            linewidth = 2.5 if style == 'bold' else 1.5
            alpha = 0.8 if style == 'bold' else 0.5
            
            ax.annotate('', xy=end_adj, xytext=start_adj,
                       arrowprops=dict(arrowstyle='->', color=color, 
                                      lw=linewidth, alpha=alpha, shrinkA=5, shrinkB=5))
        
        # Add weight label at midpoint
        mid_x, mid_y = (x1 + x2)/2, (y1 + y2)/2
        offset_x, offset_y = -dy*0.1, dx*0.1  # Perpendicular offset
        
        ax.text(mid_x + offset_x, mid_y + offset_y, f'β={weight:.3f}',
                ha='center', va='center', fontsize=9,
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8))
    
    # Add title and interpretation
    ax.set_title('Genetic Architecture Network\n(Edge weights = MR β values)', fontsize=14, pad=20)
    
    # Add legend for edge strengths
    ax.text(-1.5, -1.5, 'Strong effect: β=2.096 (Ratio→Size)\n' +
                       'Weak effect: β=0.379 (Size→Ratio)\n' +
                       'Negative: β=-0.124 (Size→Protein)\n' +
                       'Positive: β=0.087 (Size→Oil)',
            bbox=dict(boxstyle="round,pad=0.5", facecolor="lightyellow", alpha=0.7),
            fontsize=9, ha='left')
    
    # Set axis limits and remove ticks
    ax.set_xlim(-1.8, 1.8)
    ax.set_ylim(-1.8, 1.8)
    ax.set_aspect('equal')
    ax.axis('off')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved comparison plot 3: {output_file}")

def main():
    """Main function."""
    output_dir = "docs/pics"
    os.makedirs(output_dir, exist_ok=True)
    
    # Create comparison plots
    create_comparison_plot1(os.path.join(output_dir, "compare_plot1.png"))
    create_comparison_plot2(os.path.join(output_dir, "compare_plot2.png"))
    create_comparison_plot3(os.path.join(output_dir, "compare_plot3.png"))

if __name__ == "__main__":
    main()