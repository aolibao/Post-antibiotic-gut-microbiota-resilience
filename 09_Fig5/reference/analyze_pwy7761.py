
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

# Configuration
meta_path = "species/kraken.P1.P4/endpoint_sample_metadata.csv"
ko_path = "genefamily/ko_cpm_unstratified.csv"
output_plot = "result/Fig5/PWY-7761_genes_boxplot.pdf"
output_stats = "result/Fig5/PWY-7761_genes_stats.csv"

# Target KOs for PWY-7761 (NAD salvage pathway II)
# Based on MetaCyc and KEGG mapping
target_kos = {
    'K01951': 'nadE (Gln)',
    'K01916': 'nadE (NH3)',
    'K00963': 'nadD',
    'K03426': 'pncC',
    'K01501': 'pnuE/nudC',
    'K00615': 'nudC'
}

# Load data
print("Loading metadata...")
meta = pd.read_csv(meta_path)
meta = meta[meta['Group'].isin(['P1', 'P4'])]
print(f"Loaded {len(meta)} samples in P1/P4")

print("Loading KO table...")
ko_df = pd.read_csv(ko_path, index_col=0)
print(f"Loaded KO table with {ko_df.shape[1]} KOs and {ko_df.shape[0]} samples")

# Filter for target KOs
available_kos = [ko for ko in target_kos.keys() if ko in ko_df.columns]
print(f"Found {len(available_kos)} target KOs: {available_kos}")

if not available_kos:
    print("No target KOs found in the dataset.")
    exit()

# Merge data
data = ko_df[available_kos].merge(meta[['sample_id', 'Group']], left_index=True, right_on='sample_id')

# Statistical Analysis
stats_results = []
for ko in available_kos:
    p1_vals = data[data['Group'] == 'P1'][ko]
    p4_vals = data[data['Group'] == 'P4'][ko]
    
    stat, pval = mannwhitneyu(p1_vals, p4_vals)
    gene_name = target_kos[ko]
    stats_results.append({
        'KO': ko,
        'Gene': gene_name,
        'P1_mean': p1_vals.mean(),
        'P4_mean': p4_vals.mean(),
        'P_value': pval
    })

stats_df = pd.DataFrame(stats_results)
stats_df.to_csv(output_stats, index=False)
print("\nStatistical Results:")
print(stats_df)

# Plotting
plt.figure(figsize=(12, 6))
# Melt for seaborn
plot_data = data.melt(id_vars=['sample_id', 'Group'], value_vars=available_kos, 
                      var_name='KO', value_name='CPM')
# Add gene names
plot_data['Gene'] = plot_data['KO'].map(target_kos)
plot_data['Label'] = plot_data['KO'] + '\n' + plot_data['Gene']

sns.boxplot(data=plot_data, x='Label', y='CPM', hue='Group', palette={'P1': '#E64B35', 'P4': '#4DBBD5'})
sns.stripplot(data=plot_data, x='Label', y='CPM', hue='Group', dodge=True, color='black', alpha=0.5, size=3)

plt.title('Abundance of PWY-7761 Functional Genes (P1 vs P4)')
plt.xticks(rotation=45)
plt.tight_layout()

plt.savefig(output_plot)
print(f"\nPlot saved to {output_plot}")
