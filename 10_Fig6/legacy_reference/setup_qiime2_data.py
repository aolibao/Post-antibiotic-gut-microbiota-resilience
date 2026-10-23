import os
import csv
import glob

# Paths
base_dir = "/home/lihui/qiaohan/16s测序--paa处理"
raw_data_dir = os.path.join(base_dir, "G2260120810/rawdata")
metadata_file = os.path.join(base_dir, "metadata.csv")
output_manifest = os.path.join(base_dir, "manifest.tsv")
output_metadata = os.path.join(base_dir, "sample-metadata.tsv")

def clean_id(s):
    return s.strip().replace(" ", "_").replace("+", "Pos").replace("-", "Neg").replace("(", "").replace(")", "")

def main():
    # Read metadata
    samples = []
    try:
        # Try reading with GBK encoding which is common for Chinese Windows files
        with open(metadata_file, 'r', encoding='gbk') as f:
            reader = csv.reader(f)
            header = next(reader)
            # Identify columns
            # Column 0: Directory Name / Run ID (e.g., G22601316254-1)
            # Column 1: Sample File Name (e.g., samlpe1)
            # Column 2: days
            # Column 3: group
            # Column 4: DNA (Unique ID?)
            
            print(f"Header: {header}")
            
            for row in reader:
                if not row: continue
                run_id = row[0].strip()
                sample_file_name = row[1].strip().replace("samlpe", "sample") # Fix typo if exists in filename mapping
                days = row[2].strip()
                group = row[3].strip()
                dna_id = row[4].strip()
                
                # Construct Sample ID
                # Use Group_Day as SampleID (e.g., A-1_Day-10)
                # Ensure it is unique and QIIME2 compliant
                # sample_id = f"{group}_Day{days}".replace("-", "_").replace("+", "plus")
                # Actually, group contains hyphens (A-1). 
                # Let's clean it up.
                clean_group = clean_id(group)
                clean_days = days.replace("-", "Neg").replace("+", "Pos")
                if clean_days.startswith("Neg"):
                     # Handle negative numbers nicely if needed, e.g. Neg10
                     pass
                
                sample_id = f"{clean_group}_Day_{clean_days}"
                
                # Verify file existence
                # Expected path: raw_data_dir / run_id / sample_file_name + ".fastq"
                # Check if .fastq or .fq
                file_path = os.path.join(raw_data_dir, run_id, sample_file_name + ".fastq")
                
                if not os.path.exists(file_path):
                    # Try checking directory listing to match
                    # Maybe it is .fq?
                    if os.path.exists(file_path.replace(".fastq", ".fq")):
                        file_path = file_path.replace(".fastq", ".fq")
                    else:
                        print(f"Warning: File not found for {sample_id}: {file_path}")
                        # Try to find case-insensitive or partial match?
                        # For now, skip or mark as missing
                        continue
                
                samples.append({
                    'sample-id': sample_id,
                    'absolute-filepath': file_path,
                    'days': days,
                    'group': group,
                    'dna_id': dna_id
                })

    except Exception as e:
        print(f"Error reading metadata: {e}")
        return

    # Write Manifest
    # sample-id, absolute-filepath
    with open(output_manifest, 'w') as f:
        f.write("sample-id\tabsolute-filepath\n")
        for s in samples:
            f.write(f"{s['sample-id']}\t{s['absolute-filepath']}\n")
    
    print(f"Manifest written to {output_manifest} with {len(samples)} samples.")

    # Write Sample Metadata
    # sample-id, days, group, dna_id
    with open(output_metadata, 'w') as f:
        f.write("sample-id\tdays\tgroup\tdna_id\n")
        # Define types for QIIME2 metadata (optional but good)
        f.write("#q2:types\tnumeric\tcategorical\tcategorical\n")
        for s in samples:
            f.write(f"{s['sample-id']}\t{s['days']}\t{s['group']}\t{s['dna_id']}\n")
            
    print(f"Metadata written to {output_metadata}.")

if __name__ == "__main__":
    main()
