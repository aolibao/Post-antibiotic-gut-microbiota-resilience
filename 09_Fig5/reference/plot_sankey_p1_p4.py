
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio

# Configuration
meta_path = "species/kraken.P1.P4/endpoint_sample_metadata.csv"
stratified_path = "genefamily/ko_cpm_stratified.tsv"
output_html = "result/Fig5/sankey_NAD_P1_P4.html"
output_pdf = "result/Fig5/sankey_NAD_P1_P4.pdf"

# 1. Pathway Definitions
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
all_target_kos = set()
for p in pathways.values():
    all_target_kos.update(p['kos'].keys())

# 2. Taxonomy Mapping (Genus -> Phylum)
# Simplified mapping based on common gut genera
genus_phylum_map = {
    # Bacteroidetes
    'Bacteroides': 'Bacteroidetes', 'Prevotella': 'Bacteroidetes', 'Parabacteroides': 'Bacteroidetes', 
    'Alistipes': 'Bacteroidetes', 'Phocaeicola': 'Bacteroidetes', 'Muribaculum': 'Bacteroidetes',
    # Firmicutes
    'Faecalibacterium': 'Firmicutes', 'Ruminococcus': 'Firmicutes', 'Eubacterium': 'Firmicutes',
    'Roseburia': 'Firmicutes', 'Blautia': 'Firmicutes', 'Coprococcus': 'Firmicutes',
    'Lachnospira': 'Firmicutes', 'Clostridium': 'Firmicutes', 'Anaerostipes': 'Firmicutes',
    'Lactobacillus': 'Firmicutes', 'Streptococcus': 'Firmicutes', 'Enterococcus': 'Firmicutes',
    'Veillonella': 'Firmicutes', 'Megasphaera': 'Firmicutes', 'Dorea': 'Firmicutes',
    'Subdoligranulum': 'Firmicutes', 'Oscillibacter': 'Firmicutes', 'Butyrivibrio': 'Firmicutes',
    # Proteobacteria
    'Escherichia': 'Proteobacteria', 'Klebsiella': 'Proteobacteria', 'Sutterella': 'Proteobacteria',
    'Parasutterella': 'Proteobacteria', 'Desulfovibrio': 'Proteobacteria', 'Bilophila': 'Proteobacteria',
    'Enterobacter': 'Proteobacteria', 'Citrobacter': 'Proteobacteria', 'Shigella': 'Proteobacteria',
    # Actinobacteria
    'Bifidobacterium': 'Actinobacteria', 'Collinsella': 'Actinobacteria', 'Adlercreutzia': 'Actinobacteria',
    # Verrucomicrobia
    'Akkermansia': 'Verrucomicrobia',
    # Other
    'Other': 'Other'
}

# Colors
phylum_colors = {
    'Bacteroidetes': '#4DBBD5', # Blue
    'Firmicutes': '#E64B35',    # Red
    'Proteobacteria': '#00A087',# Green
    'Actinobacteria': '#3C5488',# Dark Blue
    'Verrucomicrobia': '#F39B7F',# Orange
    'Other': '#84919E'          # Grey
}

pathway_colors = {
    'NAD salvage II': '#E64B35',      # Reddish
    'NAD biosynthesis I': '#4DBBD5'   # Blueish
}

# Gene colors (shared scheme)
# We will use a neutral or mixed color, or color by pathway if possible.
# But since genes are shared, we can use a distinct color.
gene_node_color = '#7E6148' # Brownish

# 3. Load Metadata
print("Loading metadata...")
meta = pd.read_csv(meta_path)
p1_samples = meta[meta['Group'] == 'P1']['sample_id'].tolist()
p4_samples = meta[meta['Group'] == 'P4']['sample_id'].tolist()

# 4. Load Stratified Data
print("Processing stratified file...")
data_records = []
with open(stratified_path, 'r') as f:
    header = f.readline().strip().split('\t')
    
    # Map columns
    p1_cols = [i for i, h in enumerate(header) if any(s in h for s in p1_samples)]
    p4_cols = [i for i, h in enumerate(header) if any(s in h for s in p4_samples)]
    
    for line in f:
        parts = line.strip().split('\t')
        gene_family = parts[0]
        
        ko_match = None
        for ko in all_target_kos:
            if gene_family.startswith(ko + "|"):
                ko_match = ko
                break
        
        if ko_match:
            try:
                # Extract Species and Genus
                # Format: K01951|g__Bacteroides.s__Bacteroides_vulgatus
                taxa_part = gene_family.split('|')[1]
                
                # Extract Genus for coloring
                if "g__" in taxa_part:
                    genus = taxa_part.split("g__")[1].split(".")[0]
                else:
                    genus = "Other"
                
                # Extract Species Name
                if "s__" in taxa_part:
                    species_name = taxa_part.split("s__")[1].replace("_", " ")
                else:
                    species_name = taxa_part
                
                # Abundance
                vals_p1 = [float(parts[i]) for i in p1_cols]
                vals_p4 = [float(parts[i]) for i in p4_cols]
                
                mean_p1 = sum(vals_p1) / len(vals_p1) if vals_p1 else 0
                mean_p4 = sum(vals_p4) / len(vals_p4) if vals_p4 else 0
                
                if mean_p1 > 0 or mean_p4 > 0:
                    data_records.append({
                        'KO': ko_match,
                        'Species': species_name,
                        'Genus': genus,
                        'Abundance_P1': mean_p1,
                        'Abundance_P4': mean_p4
                    })
            except:
                continue

df = pd.DataFrame(data_records)

# 5. Build Sankey Function
def build_sankey_data(group_col, title):
    # Filter Top Species for this group
    total_abund = df.groupby('Species')[group_col].sum().sort_values(ascending=False)
    top_species = total_abund.head(15).index.tolist()
    
    df_g = df.copy()
    df_g['Species_Grouped'] = df_g['Species'].apply(lambda x: x if x in top_species else 'Other')
    df_g['Genus_Grouped'] = df_g.apply(lambda row: row['Genus'] if row['Species'] in top_species else 'Other', axis=1)
    
    # Aggregate
    agg = df_g.groupby(['KO', 'Species_Grouped', 'Genus_Grouped'])[group_col].sum().reset_index()
    agg = agg[agg[group_col] > 0]
    
    nodes = []
    node_map = {}
    node_colors = []
    
    def get_node(name, color):
        if name not in node_map:
            node_map[name] = len(nodes)
            nodes.append(name)
            node_colors.append(color)
        return node_map[name]
    
    links = {'source': [], 'target': [], 'value': [], 'color': []}
    
    # Logic: Species -> Gene -> Pathway
    # To merge genes, we use ONE node for each gene (e.g. "nadD")
    # But we need to split flow from Gene to multiple Pathways if applicable.
    
    # 1. Species -> Gene links
    # We iterate over aggregated data
    # But wait, we need to know which Gene connects to which Pathway FIRST to define Gene nodes?
    # No, Gene nodes are defined by KO.
    
    # Pre-calculate gene-to-pathway mapping
    # Some KOs are in multiple pathways.
    # We will split the gene abundance equally to pathways for visualization?
    # Or keep gene node as a sink for species, and source for pathways.
    
    ko_to_pathways = {}
    for p_code, p_data in pathways.items():
        for ko in p_data['kos']:
            if ko not in ko_to_pathways:
                ko_to_pathways[ko] = []
            ko_to_pathways[ko].append(p_data['name'])
            
    # Iterate KOs present in data
    present_kos = agg['KO'].unique()
    
    for ko in present_kos:
        # Get gene name
        # Find first matching name from definitions
        gene_name = None
        for p in pathways.values():
            if ko in p['kos']:
                gene_name = p['kos'][ko]
                break
        
        # Create Gene Node
        gene_node_idx = get_node(gene_name, gene_node_color)
        
        # Species -> Gene Links
        ko_data = agg[agg['KO'] == ko]
        for _, row in ko_data.iterrows():
            sp = row['Species_Grouped']
            genus = row['Genus_Grouped']
            
            # Get Phylum Color
            phylum = genus_phylum_map.get(genus, 'Other')
            sp_color = phylum_colors.get(phylum, phylum_colors['Other'])
            
            sp_node_idx = get_node(sp, sp_color)
            
            links['source'].append(sp_node_idx)
            links['target'].append(gene_node_idx)
            links['value'].append(row[group_col])
            # Link color same as source (species) but transparent
            import matplotlib.colors as mcolors
            c_rgb = mcolors.to_rgb(sp_color)
            c_rgba = f"rgba({c_rgb[0]*255},{c_rgb[1]*255},{c_rgb[2]*255},0.4)"
            links['color'].append(c_rgba)
            
        # Gene -> Pathway Links
        # Split total abundance of this gene to its pathways
        total_gene_abund = ko_data[group_col].sum()
        target_pwys = ko_to_pathways.get(ko, [])
        
        if target_pwys:
            # Split equally
            split_val = total_gene_abund / len(target_pwys)
            for pwy in target_pwys:
                pwy_color = pathway_colors.get(pwy, '#888888')
                pwy_node_idx = get_node(pwy, pwy_color)
                
                links['source'].append(gene_node_idx)
                links['target'].append(pwy_node_idx)
                links['value'].append(split_val)
                # Link color gradient? Or just gene color
                links['color'].append("rgba(150,150,150,0.3)")

    return dict(
        type='sankey',
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=nodes,
            color=node_colors
        ),
        link=links,
        domain=dict(x=[0, 1], y=[0, 1]) # Will be updated for subplots
    )

# Create Subplots
trace_p1 = build_sankey_data('Abundance_P1', 'P1 Group')
trace_p4 = build_sankey_data('Abundance_P4', 'P4 Group')

# Update domains for side-by-side
trace_p1['domain'] = dict(x=[0, 0.46], y=[0, 1])
trace_p4['domain'] = dict(x=[0.54, 1], y=[0, 1])

fig = go.Figure(data=[trace_p1, trace_p4])

fig.update_layout(
    title_text="Functional Contribution to NAD Pathways (P1 vs P4)",
    font_size=12,
    width=1200,
    height=800,
    annotations=[
        dict(x=0.23, y=1.05, showarrow=False, text="P1 Group (Pathological)", font=dict(size=16, color="black")),
        dict(x=0.77, y=1.05, showarrow=False, text="P4 Group (Recovery)", font=dict(size=16, color="black"))
    ]
)

# Save
fig.write_html(output_html)
print(f"HTML saved to {output_html}")

try:
    fig.write_image(output_pdf)
    print(f"PDF saved to {output_pdf}")
except Exception as e:
    print(f"PDF Error: {e}")

