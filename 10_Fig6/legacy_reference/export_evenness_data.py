import pandas as pd
import zipfile
import os

# Define Group Mapping
# A: ABX-Trp; B: ABX-Trp-During; C: ABX-NMN-During; D: ABX-NMN; E: ABX-Con
group_mapping = {
    'A': 'ABX-Trp',
    'B': 'ABX-Trp-During',
    'C': 'ABX-NMN-During',
    'D': 'ABX-NMN',
    'E': 'ABX-Con'
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

# Extract Evenness Data
evenness_path = 'core-metrics-results/evenness_vector.qza'

if os.path.exists(evenness_path):
    print("Extracting Evenness data...")
    df_evenness = extract_qza_data(evenness_path, 'alpha-diversity.tsv')
    if df_evenness is not None:
        df_evenness.columns = ['Evenness']
        
        # Merge with metadata
        merged_evenness = metadata.join(df_evenness, how='inner')
        
        # Select and rename columns for clarity
        # We want: Sample ID (index), days, treatment_group, Evenness
        output_df = merged_evenness[['days', 'treatment_group', 'Evenness']].copy()
        output_df.index.name = 'SampleID'
        
        # Sort by day and group
        output_df = output_df.sort_values(by=['days', 'treatment_group'])
        
        # Save to CSV
        output_filename = 'Evenness_data_annotated.csv'
        output_df.to_csv(output_filename)
        print(f"Successfully exported data to {output_filename}")
        print(output_df.head())
    else:
        print("Failed to extract data from QZA.")
else:
    print(f"File not found: {evenness_path}")
