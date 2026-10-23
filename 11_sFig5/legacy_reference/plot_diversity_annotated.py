import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import zipfile
import os
import io
import skbio
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms

# Set style
sns.set(style="whitegrid")

# Define Group Mapping and Colors (User provided)
# A: ABX-Trp; B: ABX-Trp-During; C: ABX-NMN-During; D: ABX-NMN; E: ABX-Con
group_mapping = {
    'A': 'ABX-Trp',
    'B': 'ABX-Trp-During',
    'C': 'ABX-NMN-During',
    'D': 'ABX-NMN',
    'E': 'ABX-Con'
}

# User image legend colors
color_mapping = {
    'ABX-Con': '#1f77b4',       # Blue
    'ABX-Trp': '#d62728',       # Red
    'ABX-NMN': '#ff7f0e',       # Orange
    'ABX-Trp-During': '#9467bd',# Purple
    'ABX-NMN-During': '#2ca02c' # Green
}

def draw_confidence_ellipse(x, y, ax, confidence=0.3, **kwargs):
    """
    Draw a confidence ellipse of *confidence* level.
    """
    if x.size < 3:
        return

    cov = np.cov(x, y)
    
    # Check if cov is singular
    if np.linalg.det(cov) == 0:
        return
        
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    
    # Eigenvalues of covariance matrix
    try:
        lambda_, v = np.linalg.eig(cov)
    except:
        return
        
    lambda_ = np.sqrt(lambda_)
    
    # Calculate scale factor for confidence level
    # Chi-square with 2 degrees of freedom
    scale = np.sqrt(stats.chi2.ppf(confidence, 2))
    
    ell_radius_x = lambda_[0] * scale
    ell_radius_y = lambda_[1] * scale
    
    # Angle of the ellipse
    angle = np.degrees(np.arctan2(v[1, 0], v[0, 0]))
    
    ell = Ellipse(xy=(np.mean(x), np.mean(y)),
                  width=ell_radius_x * 2,
                  height=ell_radius_y * 2,
                  angle=angle,
                  **kwargs)
    
    ax.add_patch(ell)

# Load Metadata
metadata_file = 'metas2.csv'
print(f"Loading metadata from {metadata_file}...")
metadata = pd.read_csv(metadata_file, index_col=0)

# Create 'treatment_group' column based on 'group' (assuming first letter)
# But wait, 'metas2.csv' has 'group' column with A-1, A-2, etc.
# Let's check 'group' column again
if 'group' in metadata.columns:
    # Extract A, B, C, D, E from group ID
    metadata['group_code'] = metadata['group'].apply(lambda x: str(x).split('-')[0] if '-' in str(x) else str(x)[0])
    # Map to full names
    metadata['treatment_group'] = metadata['group_code'].map(group_mapping)
else:
    print("Error: 'group' column not found.")
    exit(1)

# Function to extract data from qza
def extract_qza_data(qza_path, data_filename):
    with zipfile.ZipFile(qza_path, 'r') as z:
        for filename in z.namelist():
            if filename.endswith(f'data/{data_filename}'):
                with z.open(filename) as f:
                    return pd.read_csv(f, sep='\t', index_col=0)
    return None

# ==========================================
# 1. Shannon Boxplot with Significance
# ==========================================
print("\nGenerating Shannon Boxplot...")
shannon_path = 'core-metrics-results/shannon_vector.qza'
if os.path.exists(shannon_path):
    df_shannon = extract_qza_data(shannon_path, 'alpha-diversity.tsv')
    if df_shannon is not None:
        df_shannon.columns = ['Shannon']
        merged_shannon = metadata.join(df_shannon, how='inner')
        
        # Calculate Stats (Mann-Whitney U test vs ABX-Con) per Timepoint
        stats_data = []
        for day in sorted(merged_shannon['days'].unique()):
            day_data = merged_shannon[merged_shannon['days'] == day]
            
            # Get control data
            con_data = day_data[day_data['treatment_group'] == 'ABX-Con']['Shannon']
            
            if len(con_data) > 0:
                for group in sorted(day_data['treatment_group'].unique()):
                    if group == 'ABX-Con': continue
                    
                    treat_data = day_data[day_data['treatment_group'] == group]['Shannon']
                    
                    if len(treat_data) > 0:
                        try:
                            u_stat, p_val = stats.mannwhitneyu(con_data, treat_data, alternative='two-sided')
                            stats_data.append({
                                'days': day,
                                'group': group,
                                'p_value': p_val
                            })
                        except:
                            pass

        # Plot
        plt.figure(figsize=(14, 8))
        
        # Boxplot with custom palette
        # Sort groups for consistent order: ABX-Con first? Or alphabetical?
        # Let's put ABX-Con first if possible, or just standard sort
        order = ['ABX-Con', 'ABX-Trp', 'ABX-NMN', 'ABX-Trp-During', 'ABX-NMN-During']
        # Filter order to only include present groups
        present_groups = [g for g in order if g in merged_shannon['treatment_group'].unique()]
        
        ax = sns.boxplot(x='days', y='Shannon', hue='treatment_group', data=merged_shannon, 
                         palette=color_mapping, hue_order=present_groups, showfliers=False)
        
        # Add swarmplot for points? Maybe too crowded. Let's skip swarmplot to keep it clean like user provided.
        # But adding significance stars is key.
        
        # Significance Annotation Logic
        # We need to find the position of the boxes to place stars.
        # Box positions:
        # Each x-tick (day) has clusters of boxes.
        # If there are N groups, the center is at x (0, 1, 2...).
        # The width of cluster is usually 0.8.
        # Each box width is 0.8 / N.
        # The offset for group i (0..N-1) is (i - (N-1)/2) * width.
        
        n_groups = len(present_groups)
        width = 0.8 / n_groups
        
        y_max = merged_shannon['Shannon'].max()
        y_range = y_max - merged_shannon['Shannon'].min()
        y_offset = y_range * 0.05
        
        # Iterate stats and add stars
        for stat in stats_data:
            day = stat['days']
            group = stat['group']
            p = stat['p_value']
            
            if p < 0.05:
                sig = '*'
                if p < 0.01: sig = '**'
                if p < 0.001: sig = '***'
                
                # Find x position
                day_idx = sorted(merged_shannon['days'].unique()).index(day)
                group_idx = present_groups.index(group)
                
                x_pos = day_idx + (group_idx - (n_groups - 1) / 2) * width
                
                # Find y position (max of that box + offset)
                box_data = merged_shannon[(merged_shannon['days'] == day) & (merged_shannon['treatment_group'] == group)]['Shannon']
                y_pos = box_data.max() + y_offset
                
                plt.text(x_pos, y_pos, sig, ha='center', va='bottom', fontsize=12, color='black')

        plt.title('Shannon Diversity by Group and Time')
        plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0., title='Group')
        plt.tight_layout()
        plt.savefig('Shannon_boxplot_meta2_annotated.pdf', format='pdf')
        plt.savefig('Shannon_boxplot_meta2_annotated.png', format='png')
        plt.close()

# ==========================================
# 2. Weighted UniFrac PCoA with Significance
# ==========================================
print("\nGenerating Weighted UniFrac PCoA...")
pcoa_path = 'core-metrics-results/weighted_unifrac_pcoa_results.qza'
dm_path = 'core-metrics-results/weighted_unifrac_distance_matrix.qza'

if os.path.exists(pcoa_path) and os.path.exists(dm_path):
    # Extract Distance Matrix
    print("  Extracting Distance Matrix...")
    dm = extract_qza_data(dm_path, 'distance-matrix.tsv')
    
    # Extract PCoA Coords
    print("  Extracting PCoA Coords...")
    # (Using manual parsing logic from previous script for robustness)
    with zipfile.ZipFile(pcoa_path, 'r') as z:
        for filename in z.namelist():
            if filename.endswith('data/ordination.txt'):
                with z.open(filename) as f:
                    lines = f.read().decode('utf-8').split('\n')
                    site_start = 0
                    for i, line in enumerate(lines):
                        if line.startswith('Site'):
                            site_start = i
                            break
                    data_lines = []
                    for i in range(site_start + 1, len(lines)):
                        line = lines[i].strip()
                        if not line or line.startswith('Biplot') or line.startswith('Proportion explained'):
                            break
                        data_lines.append(line.split('\t'))
                    
                    df_pcoa = pd.DataFrame(data_lines)
                    df_pcoa = df_pcoa.set_index(0)
                    df_pcoa.index.name = 'sample-id'
                    df_pcoa = df_pcoa.iloc[:, 0:2] # PC1, PC2
                    df_pcoa.columns = ['PC1', 'PC2']
                    df_pcoa = df_pcoa.astype(float)
                    
                    # Proportion explained
                    prop_explained = []
                    for i, line in enumerate(lines):
                        if line.startswith('Proportion explained'):
                            parts = lines[i+1].strip().split('\t')
                            prop_explained = [float(x) for x in parts[:2]]
                            break
    
    # Merge
    merged_pcoa = metadata.join(df_pcoa, how='inner')
    
    # Filter DM to match samples in metadata
    # Ensure indices match
    # Note: DM index might have duplicates if metadata has duplicates? 
    # Or DM might be smaller.
    
    # Check for duplicates in index
    if merged_pcoa.index.duplicated().any():
        print("Warning: Duplicate sample IDs in merged metadata. Dropping duplicates.")
        merged_pcoa = merged_pcoa[~merged_pcoa.index.duplicated()]
        
    common_samples = merged_pcoa.index.intersection(dm.index)
    
    # Ensure DM has unique index too
    if dm.index.duplicated().any():
        dm = dm[~dm.index.duplicated()]
        dm = dm.loc[:, ~dm.columns.duplicated()]

    # Reorder DM to match common_samples
    dm_filtered = dm.loc[common_samples, common_samples]
    metadata_filtered = merged_pcoa.loc[common_samples]
    
    # PERMANOVA (Global Group Effect)
    print("  Calculating PERMANOVA...")
    try:
        # Create DistanceMatrix object (Ensure C-contiguous)
        dist_values = np.ascontiguousarray(dm_filtered.values)
        dist_matrix = skbio.DistanceMatrix(dist_values, ids=dm_filtered.index)
        
        # Run PERMANOVA
        permanova_results = skbio.stats.distance.permanova(dist_matrix, metadata_filtered, column='treatment_group', permutations=999)
        
        p_val_global = permanova_results['p-value']
        f_stat = permanova_results['test statistic']
        stats_text = f"PERMANOVA (Group): F={f_stat:.2f}, p={p_val_global:.3f}"
        print(f"  {stats_text}")
        
    except Exception as e:
        print(f"PERMANOVA Error: {e}")
        stats_text = ""

    # Plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Scatterplot
    sns.scatterplot(x='PC1', y='PC2', hue='treatment_group', style='days', data=merged_pcoa, 
                    palette=color_mapping, s=100, hue_order=present_groups, ax=ax)
    
    # Add Confidence Ellipses (0.3)
    # Loop through each group and day combination
    for group in present_groups:
        for day in sorted(merged_pcoa['days'].unique()):
            subset = merged_pcoa[(merged_pcoa['treatment_group'] == group) & (merged_pcoa['days'] == day)]
            
            # Need at least 3 points to draw a meaningful ellipse
            if len(subset) >= 3:
                # Use same color as group, maybe with some transparency
                color = color_mapping[group]
                
                # Calculate and draw ellipse
                # Use a helper function or simpler matplotlib ellipse
                # Here we use the custom function defined above
                try:
                    # Confidence 0.3 means 30% confidence interval
                    draw_confidence_ellipse(subset['PC1'], subset['PC2'], ax, confidence=0.3, 
                                            edgecolor=color, facecolor=color, alpha=0.2, linestyle='--')
                    # Draw outline stronger?
                    draw_confidence_ellipse(subset['PC1'], subset['PC2'], ax, confidence=0.3, 
                                            edgecolor=color, facecolor='none', alpha=0.8, linestyle='-', linewidth=1.5)
                except Exception as e:
                    print(f"Error drawing ellipse for {group} Day {day}: {e}")

    pc1_label = f"PC1 ({prop_explained[0]*100:.1f}%)" if prop_explained else "PC1"
    pc2_label = f"PC2 ({prop_explained[1]*100:.1f}%)" if prop_explained else "PC2"
    
    plt.xlabel(pc1_label)
    plt.ylabel(pc2_label)
    
    title = 'Weighted UniFrac PCoA'
    if stats_text:
        title += f'\n{stats_text}'
        
    plt.title(title)
    plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0., title='Group')
    plt.tight_layout()
    plt.savefig('Weighted_UniFrac_PCoA_meta2_annotated.pdf', format='pdf')
    plt.savefig('Weighted_UniFrac_PCoA_meta2_annotated.png', format='png')
    plt.close()

print("\nDone.")
