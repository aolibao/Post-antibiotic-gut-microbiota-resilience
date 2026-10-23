
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
output_plot = "result/Fig5/PWY-7761_genes_subplots.pdf"

# Gene Metadata
gene_info = {
    'K01951': {'gene': 'nadE', 'ec': 'EC 6.3.5.1', 'desc': 'NAD+ synthase (Gln)'},
    'K01916': {'gene': 'nadE', 'ec': 'EC 6.3.1.5', 'desc': 'NAD+ synthase (NH3)'},
    'K00963': {'gene': 'nadD', 'ec': 'EC 2.7.7.18', 'desc': 'NaMN adenylyltransferase'},
    'K03426': {'gene': 'pncC', 'ec': 'EC 3.5.1.42', 'desc': 'Nicotinamidase'},
    'K00615': {'gene': 'nudC', 'ec': 'EC 3.6.1.22', 'desc': 'NAD+ pyrophosphatase'},
    # K01501 was excluded as it had near zero abundance in previous check
}

# Load data
print("Loading metadata...")
meta = pd.read_csv(meta_path)
meta = meta[meta['Group'].isin(['P1', 'P4'])]

print("Loading KO table...")
ko_df = pd.read_csv(ko_path, index_col=0)

# Filter for target KOs
available_kos = [ko for ko in gene_info.keys() if ko in ko_df.columns]
print(f"Found KOs: {available_kos}")

# Prepare plot data
data = ko_df[available_kos].merge(meta[['sample_id', 'Group']], left_index=True, right_on='sample_id')

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
    
    # Stats
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
    # Get y-axis range to position text
    y_max = plot_df[ko].max()
    y_min = plot_df[ko].min()
    y_range = y_max - y_min
    
    # Scientific notation for p-value
    if pval < 0.001:
        p_text = f"p = {pval:.2e}"
    else:
        p_text = f"p = {pval:.4f}"
        
    ax.text(0.5, 0.9, p_text, transform=ax.transAxes, ha='center', va='center', 
            fontsize=12, bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

    # Remove legend inside subplots (can add a global one if needed, but groups are clear)
    if ax.get_legend():
        ax.get_legend().remove()

# Hide empty subplots
for j in range(i + 1, len(axes)):
    axes[j].axis('off')

plt.tight_layout()
plt.savefig(output_plot)
print(f"Saved plot to {output_plot}")
