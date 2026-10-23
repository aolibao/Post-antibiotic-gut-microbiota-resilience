
import pandas as pd
import plotly.graph_objects as go
import re

# Configuration
meta_path = "species/kraken.P1.P4/endpoint_sample_metadata.csv"
stratified_path = "genefamily/ko_cpm_stratified.tsv"
output_html = "result/Fig5/sankey_NAD_pathways.html"
output_pdf = "result/Fig5/sankey_NAD_pathways.pdf"

# Pathway Definitions
# Genes for each pathway
pathways = {
    'PWY-7761': {
        'name': 'NAD salvage II',
        'kos': {
            'K01951': 'nadE',
            'K00963': 'nadD',
            'K03426': 'pncC',
            'K00615': 'nudC',
            'K01501': 'pnuE'
        }
    },
    'PYRIDNUCSYN-PWY': {
        'name': 'NAD biosynthesis I',
        'kos': {
            'K00278': 'nadB',
            'K03704': 'nadA',
            'K00763': 'nadC',
            'K00963': 'nadD', # Shared
            'K01951': 'nadE'  # Shared
        }
    }
}

# Collect all relevant KOs
all_target_kos = set()
for p in pathways.values():
    all_target_kos.update(p['kos'].keys())

print(f"Target KOs: {all_target_kos}")

# 1. Load Metadata and Identify P4 Samples
print("Loading metadata...")
meta = pd.read_csv(meta_path)
p4_samples = meta[meta['Group'] == 'P4']['sample_id'].tolist()
print(f"Found {len(p4_samples)} P4 samples.")

# 2. Process Stratified File
print("Processing stratified file...")
# We will read line by line or chunks because the file might be large? 
# Or just read it if it's not huge. The file size is unknown but usually these are big.
# However, we only need rows matching our KOs.
# Let's try reading with pandas and filtering immediately. 
# The column names in the file likely have suffixes like "-CPM".
# We need to map sample IDs to column names.

# Read header first
with open(stratified_path, 'r') as f:
    header = f.readline().strip().split('\t')

# Map P4 sample IDs to columns
# Sample ID "SRR..." usually matches column "SRR...-CPM" or just "SRR..."
# Based on previous `head`, columns are like "SRR24704728-CPM".
p4_cols = []
for s in p4_samples:
    # Try exact match or with suffix
    if s in header:
        p4_cols.append(s)
    elif f"{s}-CPM" in header:
        p4_cols.append(f"{s}-CPM")
    # else: print(f"Warning: Sample {s} not found in header")

print(f"Matched {len(p4_cols)} columns for P4 samples.")

# Load only relevant rows
# We can't easily filter rows with read_csv without reading all.
# So we will read the file and filter.
data_records = []

# Read file line by line to save memory
with open(stratified_path, 'r') as f:
    header = f.readline().strip().split('\t')
    col_indices = [header.index(c) for c in p4_cols]
    
    for line in f:
        parts = line.strip().split('\t')
        gene_family = parts[0]
        
        # Check if row matches any target KO
        # Format: K01951|g__...
        ko_match = None
        for ko in all_target_kos:
            if gene_family.startswith(ko + "|"):
                ko_match = ko
                break
        
        if ko_match:
            # Extract Species
            # gene_family looks like "K01951|g__Absiella.s__Absiella_dolichum"
            try:
                species_part = gene_family.split('|')[1]
                # Clean species name: remove "g__Genus.s__" prefix if standard
                # Or keep it readable. "s__Absiella_dolichum" -> "Absiella dolichum"
                if "s__" in species_part:
                    species_name = species_part.split("s__")[1].replace("_", " ")
                else:
                    species_name = species_part # Keep as is (e.g. unclassified)
                
                # Calculate mean abundance for P4
                values = [float(parts[i]) for i in col_indices]
                mean_abundance = sum(values) / len(values)
                
                if mean_abundance > 0:
                    data_records.append({
                        'KO': ko_match,
                        'Species': species_name,
                        'Abundance': mean_abundance
                    })
            except IndexError:
                continue

df = pd.DataFrame(data_records)
print(f"Loaded {len(df)} records.")
print(df.head())

# 3. Aggregate Data for Sankey
# Filter Top Species
# Calculate total abundance per species across all relevant KOs
species_totals = df.groupby('Species')['Abundance'].sum().sort_values(ascending=False)
top_n = 15 # Top 15 species
top_species = species_totals.head(top_n).index.tolist()

print(f"Top species: {top_species}")

# Replace other species with "Other"
df['Species_Grouped'] = df['Species'].apply(lambda x: x if x in top_species else 'Other')

# Re-aggregate with grouped species
df_grouped = df.groupby(['KO', 'Species_Grouped'])['Abundance'].sum().reset_index()

# Construct Nodes and Links
# Nodes: Pathways, Genes (per pathway), Species
nodes = []
node_map = {} # Name -> Index

def get_node_index(name):
    if name not in node_map:
        node_map[name] = len(nodes)
        nodes.append(name)
    return node_map[name]

links = {'source': [], 'target': [], 'value': [], 'color': []}

# Colors
# Pathways
pwy_colors = {'NAD salvage II': '#E64B35', 'NAD biosynthesis I': '#4DBBD5'}
# Genes
gene_colors = '#00A087'
# Species
species_colors = '#3C5488'
link_color_alpha = 'rgba(200, 200, 200, 0.5)'

# Build flow: Pathway -> Gene -> Species
for pwy_code, pwy_data in pathways.items():
    pwy_name = pwy_data['name']
    pwy_node = get_node_index(pwy_name)
    
    for ko, gene_name in pwy_data['kos'].items():
        # Gene Node Name (unique per pathway to split flow)
        # e.g. "nadD (Salvage)"
        gene_node_label = f"{gene_name}" 
        # Actually, for Sankey, if we want to show shared genes distinctly, we need unique IDs but maybe same label?
        # Plotly allows duplicate labels if we manage indices.
        # But let's use unique labels for clarity in this specific chart structure.
        # Or better: "nadD" node.
        # If we use a single "nadD" node, flows from both pathways merge into it.
        # Then flow from "nadD" goes to species.
        # This is actually a very good representation:
        # PWY1 -> nadD
        # PWY2 -> nadD
        # nadD -> Species
        # But wait, the flow magnitude from PWY1 to nadD is the Total Abundance of nadD?
        # No.
        # The logic is: Species contribute to Genes. Genes constitute Pathways.
        # So: Species -> Gene -> Pathway.
        # Let's do Species -> Gene -> Pathway.
        
        # Calculate total abundance of this gene (from grouped species)
        gene_records = df_grouped[df_grouped['KO'] == ko]
        
        if gene_records.empty:
            continue
            
        # Gene Node
        # We use a unique node for the gene to allow merging from species
        # But if the gene belongs to multiple pathways, we need to handle the output flow.
        # If "nadD" goes to PWY1 and PWY2, we have to split the abundance.
        # Since we can't split biologically (it's the same pool), we will DUPLICATE the gene node for each pathway.
        # Species -> nadD (Salvage) -> Salvage Pathway
        # Species -> nadD (De novo) -> De novo Pathway
        # This shows the contribution of species to that specific functional context.
        
        gene_display_name = f"{gene_name} " # Add space to differentiate if needed or use hidden ID
        # Let's use "Gene (Pathway)" format for the node ID, but display just "Gene"
        gene_node_id = f"{gene_name} ({pwy_name})"
        gene_node_idx = get_node_index(gene_node_id)
        
        # Link: Gene -> Pathway
        total_abund = gene_records['Abundance'].sum()
        links['source'].append(gene_node_idx)
        links['target'].append(pwy_node)
        links['value'].append(total_abund)
        links['color'].append(link_color_alpha)
        
        # Link: Species -> Gene
        for _, row in gene_records.iterrows():
            sp_name = row['Species_Grouped']
            sp_node = get_node_index(sp_name)
            
            links['source'].append(sp_node)
            links['target'].append(gene_node_idx)
            links['value'].append(row['Abundance'])
            links['color'].append(link_color_alpha)

# Create Plot
# Clean labels (remove pathway suffix from genes for display)
node_labels = [n.split(' (')[0] if '(' in n and 'NAD' not in n else n for n in nodes]

# Save data for user
df_grouped.to_csv("result/Fig5/sankey_NAD_data.csv", index=False)
print("Data saved to result/Fig5/sankey_NAD_data.csv")

fig = go.Figure(data=[go.Sankey(
    node=dict(
        pad=15,
        thickness=20,
        line=dict(color="black", width=0.5),
        label=node_labels,
        color="blue" # Default, can customize
    ),
    link=dict(
        source=links['source'],
        target=links['target'],
        value=links['value'],
        color=links['color']
    )
)])

fig.update_layout(title_text="Functional Contribution to NAD Pathways (P4 Group)", font_size=12, width=1000, height=800)
fig.write_html(output_html)
# fig.write_image(output_pdf) # Requires kaleido

print(f"Sankey saved to {output_html}")

# For PDF, we need kaleido. If not available, we rely on user opening HTML.
# I'll check if I can save as PDF.
try:
    fig.write_image(output_pdf)
    print(f"PDF saved to {output_pdf}")
except Exception as e:
    print(f"Could not save PDF: {e}")

