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

# Load Metadata
metadata_file = 'metas2.csv'
print(f"Loading metadata from {metadata_file}...")
metadata = pd.read_csv(metadata_file, index_col=0)

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
# Evenness Boxplot with Significance
# ==========================================
print("\nGenerating Evenness Boxplot...")
evenness_path = 'core-metrics-results/evenness_vector.qza'

if os.path.exists(evenness_path):
    df_evenness = extract_qza_data(evenness_path, 'alpha-diversity.tsv')
    if df_evenness is not None:
        # Check column name, usually 'pielou_evenness' or similar
        # But extract_qza_data reads it. Let's print columns if needed.
        # Usually it is the only column.
        # Let's rename it to 'Evenness'
        df_evenness.columns = ['Evenness']
        
        merged_evenness = metadata.join(df_evenness, how='inner')
        
        # Calculate Stats (Mann-Whitney U test vs ABX-Con) per Timepoint
        stats_data = []
        for day in sorted(merged_evenness['days'].unique()):
            day_data = merged_evenness[merged_evenness['days'] == day]
            
            # Get control data
            con_data = day_data[day_data['treatment_group'] == 'ABX-Con']['Evenness']
            
            if len(con_data) > 0:
                for group in sorted(day_data['treatment_group'].unique()):
                    if group == 'ABX-Con': continue
                    
                    treat_data = day_data[day_data['treatment_group'] == group]['Evenness']
                    
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
        order = ['ABX-Con', 'ABX-Trp', 'ABX-NMN', 'ABX-Trp-During', 'ABX-NMN-During']
        # Filter order to only include present groups
        present_groups = [g for g in order if g in merged_evenness['treatment_group'].unique()]
        
        ax = sns.boxplot(x='days', y='Evenness', hue='treatment_group', data=merged_evenness, 
                         palette=color_mapping, hue_order=present_groups, showfliers=False)
        
        # Significance Annotation Logic
        n_groups = len(present_groups)
        width = 0.8 / n_groups
        
        y_max = merged_evenness['Evenness'].max()
        y_range = y_max - merged_evenness['Evenness'].min()
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
                day_idx = sorted(merged_evenness['days'].unique()).index(day)
                group_idx = present_groups.index(group)
                
                x_pos = day_idx + (group_idx - (n_groups - 1) / 2) * width
                
                # Find y position (max of that box + offset)
                box_data = merged_evenness[(merged_evenness['days'] == day) & (merged_evenness['treatment_group'] == group)]['Evenness']
                y_pos = box_data.max() + y_offset
                
                plt.text(x_pos, y_pos, sig, ha='center', va='bottom', fontsize=12, color='black')

        plt.title('Evenness by Group and Time')
        plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0., title='Group')
        plt.tight_layout()
        plt.savefig('Evenness_boxplot_meta2_annotated.pdf', format='pdf')
        plt.savefig('Evenness_boxplot_meta2_annotated.png', format='png')
        plt.close()
        print("Evenness plot saved.")

print("\nDone.")
