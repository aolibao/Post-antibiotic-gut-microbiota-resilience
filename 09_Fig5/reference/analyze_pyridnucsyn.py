
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import matplotlib.ticker as ticker

# Set font type for editable text in PDF
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

# Configuration
meta_path = "species/kraken.P1.P4/endpoint_sample_metadata.csv"
ko_path = "genefamily/ko_cpm_unstratified.csv"
output_plot = "result/Fig5/PYRIDNUCSYN_genes_subplots.pdf"
output_stats = "result/Fig5/PYRIDNUCSYN_genes_stats.csv"

# Gene Metadata for PYRIDNUCSYN-PWY (NAD biosynthesis I)
gene_info = {
    'K00278': {'gene': 'nadB', 'ec': 'EC 1.4.3.16', 'desc': 'L-aspartate oxidase'},
    'K03704': {'gene': 'nadA', 'ec': 'EC 2.5.1.72', 'desc': 'Quinolinate synthase'},
    'K00763': {'gene': 'nadC', 'ec': 'EC 2.4.2.19', 'desc': 'Quinolinate phosphoribosyltransferase'},
    'K00963': {'gene': 'nadD', 'ec': 'EC 2.7.7.18', 'desc': 'NaMN adenylyltransferase'},
    'K01951': {'gene': 'nadE (Gln)', 'ec': 'EC 6.3.5.1', 'desc': 'NAD+ synthase (Gln)'},
    'K01916': {'gene': 'nadE (NH3)', 'ec': 'EC 6.3.1.5', 'desc': 'NAD+ synthase (NH3)'}
}

# Load data
print("Loading metadata...")
meta = pd.read_csv(meta_path)
meta = meta[meta['Group'].isin(['P1', 'P4'])]
print(f"Loaded {len(meta)} samples in P1/P4")

print("Loading KO table...")
ko_df = pd.read_csv(ko_path, index_col=0)

# Filter for target KOs
available_kos = [ko for ko in gene_info.keys() if ko in ko_df.columns]
print(f"Found KOs: {available_kos}")

# Prepare plot data
data = ko_df[available_kos].merge(meta[['sample_id', 'Group']], left_index=True, right_on='sample_id')

# Statistical Analysis
stats_results = []
for ko in available_kos:
    p1_vals = data[data['Group'] == 'P1'][ko]
    p4_vals = data[data['Group'] == 'P4'][ko]
    
    stat, pval = mannwhitneyu(p1_vals, p4_vals)
    info = gene_info[ko]
    stats_results.append({
        'KO': ko,
        'Gene': info['gene'],
        'P1_mean': p1_vals.mean(),
        'P4_mean': p4_vals.mean(),
        'P_value': pval
    })

stats_df = pd.DataFrame(stats_results)
stats_df.to_csv(output_stats, index=False)
print("\nStatistical Results:")
print(stats_df)

# Plotting
n_genes = len(available_kos)
n_cols = 3
n_rows = (n_genes + n_cols - 1) // n_cols

# Set figsize to make subplots square (e.g., 4x4 inches per subplot)
fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows))
axes = axes.flatten()

palette = {'P1': '#E64B35', 'P4': '#4DBBD5'}

for i, ko in enumerate(available_kos):
    ax = axes[i]
    info = gene_info[ko]
    
    # Data for this gene
    plot_df = data[['Group', ko]].dropna()
    p1 = plot_df[plot_df['Group'] == 'P1'][ko]
    p4 = plot_df[plot_df['Group'] == 'P4'][ko]
    
    # Stats (re-calculate or use stored)
    stat, pval = mannwhitneyu(p1, p4)
    
    # Boxplot
    sns.boxplot(data=plot_df, x='Group', y=ko, hue='Group', palette=palette, ax=ax, dodge=False)
    sns.stripplot(data=plot_df, x='Group', y=ko, color='black', alpha=0.5, size=3, ax=ax)
    
    # Title and Labels
    title = f"{info['gene']}\n({info['ec']})"
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('Abundance (CPM)')
    
    # Add P-value
    y_max = plot_df[ko].max()
    y_min = plot_df[ko].min()
    
    # Scientific notation for p-value
    if pval < 0.001:
        p_text = f"p = {pval:.2e}"
    else:
        p_text = f"p = {pval:.4f}"
        
    ax.text(0.5, 0.9, p_text, transform=ax.transAxes, ha='center', va='center', 
            fontsize=12, bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

    # Remove legend inside subplots
    if ax.get_legend():
        ax.get_legend().remove()

# Hide empty subplots
for j in range(i + 1, len(axes)):
    axes[j].axis('off')

plt.tight_layout()
plt.savefig(output_plot)
print(f"Saved plot to {output_plot}")
