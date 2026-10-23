import pandas as pd
import zipfile
import os

# Metadata file path
metadata_file = 'sample-metadata.tsv'

# Alpha diversity files
alpha_files = {
    'Shannon': 'core-metrics-results/shannon_vector.qza',
    'Observed_Features': 'core-metrics-results/observed_features_vector.qza',
    'Faith_PD': 'core-metrics-results/faith_pd_vector.qza',
    'Evenness': 'core-metrics-results/evenness_vector.qza'
}

# Target time points
target_days = [-10, 7, 14, 28]

def load_qza_data(qza_path, data_filename='alpha-diversity.tsv'):
    """Extracts data from a QIIME 2 artifact (.qza)"""
    if not os.path.exists(qza_path):
        print(f"Warning: File not found: {qza_path}")
        return None
        
    with zipfile.ZipFile(qza_path, 'r') as z:
        for filename in z.namelist():
            if filename.endswith(f'data/{data_filename}'):
                with z.open(filename) as f:
                    return pd.read_csv(f, sep='\t', index_col=0)
    return None

def main():
    # 1. Load Metadata
    print("Loading metadata...")
    try:
        metadata = pd.read_csv(metadata_file, sep='\t', skiprows=[1], index_col=0)
    except Exception as e:
        print(f"Error loading metadata: {e}")
        return

    # Filter metadata for target days
    filtered_metadata = metadata[metadata['days'].isin(target_days)].copy()
    print(f"Filtered samples (Days {target_days}): {len(filtered_metadata)}")

    # 2. Load and Merge Alpha Diversity Data
    merged_data = filtered_metadata[['days', 'group']].copy()

    for metric_name, qza_path in alpha_files.items():
        print(f"Loading {metric_name} from {qza_path}...")
        df = load_qza_data(qza_path)
        
        if df is not None:
            # Rename column to metric name
            df.columns = [metric_name]
            # Merge with metadata (left join to keep only filtered samples that have diversity data)
            merged_data = merged_data.join(df, how='inner')
        else:
            print(f"Failed to load {metric_name}")

    print(f"Final dataset shape: {merged_data.shape}")

    # 3. Export to CSV
    output_filename = 'alpha_diversity_selected_days.csv'
    merged_data.to_csv(output_filename)
    print(f"Successfully exported to {output_filename}")
    
    # Preview
    print("\nFirst 5 rows:")
    print(merged_data.head())

if __name__ == "__main__":
    main()
