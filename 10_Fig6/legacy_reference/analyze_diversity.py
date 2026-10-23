import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import zipfile
import os
import io

# Set style
sns.set(style="whitegrid")

# Load metadata
metadata = pd.read_csv('sample-metadata.tsv', sep='\t', skiprows=[1])
metadata = metadata.set_index('sample-id')

# Create treatment group column (assuming first letter of group is the treatment)
# Check if 'group' column exists
if 'group' in metadata.columns:
    metadata['treatment_group'] = metadata['group'].apply(lambda x: x.split('-')[0] if '-' in str(x) else str(x)[0])
else:
    print("Error: 'group' column not found in metadata.")
    exit(1)

# Function to extract data from qza
def extract_qza_data(qza_path, data_filename):
    with zipfile.ZipFile(qza_path, 'r') as z:
        for filename in z.namelist():
            if filename.endswith(f'data/{data_filename}'):
                with z.open(filename) as f:
                    return pd.read_csv(f, sep='\t', index_col=0)
    return None

# Alpha diversity metrics
alpha_metrics = {
    'Shannon': 'shannon_vector.qza',
    'Observed Features': 'observed_features_vector.qza',
    'Faith PD': 'faith_pd_vector.qza',
    'Evenness': 'evenness_vector.qza'
}

alpha_results = pd.DataFrame()

# Extract and merge alpha diversity data
for metric, filename in alpha_metrics.items():
    path = os.path.join('core-metrics-results', filename)
    if os.path.exists(path):
        df = extract_qza_data(path, 'alpha-diversity.tsv')
        if df is not None:
            df.columns = [metric]
            metadata = metadata.join(df)
            alpha_results[metric] = metadata[metric]

# Statistical Analysis for Alpha Diversity
print("Alpha Diversity Statistical Analysis")
print("====================================")

for metric in alpha_results.columns:
    print(f"\nMetric: {metric}")
    
    # 1. Compare time points within each group
    print("\n  Comparison of time points within each group (Kruskal-Wallis):")
    for group in sorted(metadata['treatment_group'].unique()):
        subset = metadata[metadata['treatment_group'] == group]
        groups_data = [subset[subset['days'] == day][metric].dropna() for day in sorted(subset['days'].unique())]
        groups_data = [g for g in groups_data if len(g) > 1] # Only include if > 1 sample
        
        if len(groups_data) > 1:
            stat, p = stats.kruskal(*groups_data)
            print(f"    Group {group}: statistic={stat:.3f}, p-value={p:.3e}")
        else:
            print(f"    Group {group}: Not enough data for comparison")

    # 2. Compare groups at each time point
    print("\n  Comparison of groups at each time point (Kruskal-Wallis):")
    for day in sorted(metadata['days'].unique()):
        subset = metadata[metadata['days'] == day]
        groups_data = [subset[subset['treatment_group'] == group][metric].dropna() for group in sorted(subset['treatment_group'].unique())]
        groups_data = [g for g in groups_data if len(g) > 1]
        
        if len(groups_data) > 1:
            stat, p = stats.kruskal(*groups_data)
            print(f"    Day {day}: statistic={stat:.3f}, p-value={p:.3e}")
        else:
            print(f"    Day {day}: Not enough data for comparison")

# Plotting Alpha Diversity
print("\nGenerating Alpha Diversity Plots...")
for metric in alpha_results.columns:
    plt.figure(figsize=(12, 6))
    sns.boxplot(x='days', y=metric, hue='treatment_group', data=metadata)
    plt.title(f'{metric} by Group and Time')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(f'{metric.replace(" ", "_")}_boxplot.png')
    plt.close()

# Beta Diversity Analysis (PCoA)
print("\nBeta Diversity Analysis")
print("=======================")

beta_metrics = {
    'Bray-Curtis': 'bray_curtis_pcoa_results.qza',
    'Weighted UniFrac': 'weighted_unifrac_pcoa_results.qza',
    'Unweighted UniFrac': 'unweighted_unifrac_pcoa_results.qza',
    'Jaccard': 'jaccard_pcoa_results.qza'
}

def extract_pcoa_data(qza_path):
    with zipfile.ZipFile(qza_path, 'r') as z:
        for filename in z.namelist():
            if filename.endswith('data/ordination.txt'):
                with z.open(filename) as f:
                    # Read the ordination file manually
                    lines = f.read().decode('utf-8').split('\n')
                    
                    # Find the start of site data
                    site_start = 0
                    for i, line in enumerate(lines):
                        if line.startswith('Site'):
                            site_start = i
                            break
                    
                    # Read site data (coordinates)
                    # The number of samples is in a line before 'Site' usually, or we can just read until empty line
                    # But pandas read_csv might be tricky with variable header length
                    # Let's try to parse the 'Site' block
                    
                    data_lines = []
                    for i in range(site_start + 1, len(lines)):
                        line = lines[i].strip()
                        if not line or line.startswith('Biplot') or line.startswith('Proportion explained'):
                            break
                        data_lines.append(line.split('\t'))
                    
                    if not data_lines:
                        return None
                        
                    df = pd.DataFrame(data_lines)
                    df = df.set_index(0) # Sample ID
                    df.index.name = 'sample-id'
                    # Take first 3 PC columns
                    df = df.iloc[:, 0:3]
                    df.columns = ['PC1', 'PC2', 'PC3']
                    df = df.astype(float)
                    
                    # Extract proportion explained
                    prop_explained = []
                    for i, line in enumerate(lines):
                        if line.startswith('Proportion explained'):
                            parts = lines[i+1].strip().split('\t')
                            prop_explained = [float(x) for x in parts[:3]]
                            break
                            
                    return df, prop_explained
    return None, None

for metric, filename in beta_metrics.items():
    path = os.path.join('core-metrics-results', filename)
    if os.path.exists(path):
        pcoa_df, props = extract_pcoa_data(path)
        if pcoa_df is not None:
            # Merge with metadata
            plot_data = metadata.join(pcoa_df, how='inner')
            
            # Plot PC1 vs PC2
            plt.figure(figsize=(10, 8))
            sns.scatterplot(x='PC1', y='PC2', hue='treatment_group', style='days', data=plot_data, s=100)
            
            pc1_label = f"PC1 ({props[0]*100:.1f}%)" if props else "PC1"
            pc2_label = f"PC2 ({props[1]*100:.1f}%)" if props else "PC2"
            
            plt.xlabel(pc1_label)
            plt.ylabel(pc2_label)
            plt.title(f'{metric} PCoA')
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            plt.tight_layout()
            plt.savefig(f'{metric.replace(" ", "_")}_PCoA.png')
            plt.close()
            print(f"Generated PCoA plot for {metric}")

print("\nAnalysis Complete.")
