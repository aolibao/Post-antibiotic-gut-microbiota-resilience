
import pandas as pd
from scipy.stats import mannwhitneyu

# Load KO file
ko_path = "genefamily/ko_cpm_unstratified.csv"
meta_path = "species/kraken.P1.P4/endpoint_sample_metadata.csv"

# Load metadata
meta = pd.read_csv(meta_path)
meta = meta[meta['Group'].isin(['P1', 'P4'])]
p1_samples = meta[meta['Group'] == 'P1']['sample_id'].tolist()
p4_samples = meta[meta['Group'] == 'P4']['sample_id'].tolist()

# Load KO data
ko_df = pd.read_csv(ko_path, index_col=0)

# Target KOs for PYRIDNUCSYN-PWY
target_kos = {
    'K00278': 'nadB (Aspartate oxidase)',
    'K03704': 'nadA (Quinolinate synthase)',
    'K00763': 'nadC (Quinolinate PRT)',
    'K00963': 'nadD (NaMN adenylyltransferase)',
    'K01951': 'nadE (NAD synthase Gln)',
    'K01916': 'nadE (NAD synthase NH3)'
}

results = []

for ko, gene in target_kos.items():
    if ko in ko_df.columns:
        p1_vals = ko_df.loc[p1_samples, ko]
        p4_vals = ko_df.loc[p4_samples, ko]
        
        stat, pval = mannwhitneyu(p1_vals, p4_vals)
        mean_p1 = p1_vals.mean()
        mean_p4 = p4_vals.mean()
        
        direction = "NS"
        if pval < 0.05:
            direction = "P4" if mean_p4 > mean_p1 else "P1"
        
        results.append({
            'KO': ko,
            'Gene': gene,
            'P1_mean': round(mean_p1, 2),
            'P4_mean': round(mean_p4, 2),
            'P_value': pval,
            'Enriched': direction
        })

res_df = pd.DataFrame(results)
print(res_df.to_string())
