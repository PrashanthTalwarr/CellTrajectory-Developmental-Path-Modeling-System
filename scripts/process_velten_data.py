# scripts/process_velten_compressed.py
import os
import glob
import gzip
import pandas as pd
import numpy as np
import scanpy as sc
from anndata import AnnData
from tqdm import tqdm

# Set up paths
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(base_dir, "Data compressed")
processed_dir = os.path.join(base_dir, "data", "processed")
os.makedirs(processed_dir, exist_ok=True)

print(f"Looking for compressed CSV files in: {data_dir}")

# Get all compressed CSV files
csv_gz_files = sorted(glob.glob(os.path.join(data_dir, "*.csv.gz")))
print(f"Found {len(csv_gz_files)} compressed CSV files")

# Extract metadata from filename
def extract_metadata(filename):
    """Extract metadata from filename"""
    base_name = os.path.basename(filename)
    parts = base_name.split('_')
    
    metadata = {
        'sample_id': parts[0],  # GSM ID
    }
    
    # Extract plate and well information
    for i, part in enumerate(parts):
        if 'plate' in part:
            metadata['plate'] = part
            if i+1 < len(parts):
                metadata['well_letter'] = parts[i+1]
            if i+2 < len(parts):
                metadata['well_number'] = parts[i+2].split('.')[0]  # Remove file extension
    
    return metadata

# Sample a few files to determine the structure
print("Examining sample files to determine structure...")
sample_size = 3
sample_dfs = []

for i, gz_file in enumerate(csv_gz_files[:sample_size]):
    file_name = os.path.basename(gz_file)
    print(f"\nSample {i+1}: {file_name}")
    
    # Read the compressed file
    try:
        with gzip.open(gz_file, 'rt') as f:
            df = pd.read_csv(f)
            sample_dfs.append(df)
            print(f"Shape: {df.shape}")
            print(f"Columns: {df.columns.tolist()}")
            print(f"First few rows:")
            print(df.head(2))
    except Exception as e:
        print(f"Error reading file: {e}")
        try:
            # Try as tab-delimited
            with gzip.open(gz_file, 'rt') as f:
                df = pd.read_csv(f, sep='\t')
                sample_dfs.append(df)
                print(f"Shape (tab-delimited): {df.shape}")
                print(f"Columns: {df.columns.tolist()}")
                print(f"First few rows:")
                print(df.head(2))
        except Exception as e:
            print(f"Error reading as tab-delimited: {e}")

if not sample_dfs:
    print("Failed to read any sample files. Exiting.")
    exit(1)

# Determine data structure based on samples
first_df = sample_dfs[0]

# Determine if genes are in rows or columns
genes_in_rows = True
if len(first_df.columns) > 2:  # More than two columns suggests genes are in columns
    genes_in_rows = False
    
print(f"\nBased on sample analysis, genes appear to be in {'rows' if genes_in_rows else 'columns'}")

# Process files in batches to conserve memory
def process_files(files, batch_size=50, genes_in_rows=True):
    """Process compressed CSV files in batches"""
    all_data = []
    all_metadata = []
    unique_genes = set()
    
    total_batches = (len(files) + batch_size - 1) // batch_size
    for batch_idx in range(total_batches):
        start_idx = batch_idx * batch_size
        end_idx = min(start_idx + batch_size, len(files))
        batch_files = files[start_idx:end_idx]
        
        print(f"Processing batch {batch_idx + 1}/{total_batches} ({len(batch_files)} files)")
        
        batch_data = []
        batch_metadata = []
        
        for file_path in tqdm(batch_files):
            try:
                # Extract metadata
                meta = extract_metadata(file_path)
                
                # Read the compressed file
                with gzip.open(file_path, 'rt') as f:
                    try:
                        df = pd.read_csv(f)
                    except:
                        try:
                            # Try tab-delimited
                            f.seek(0)  # Go back to start of file
                            df = pd.read_csv(f, sep='\t')
                        except:
                            print(f"Error: Could not read {os.path.basename(file_path)}")
                            continue
                
                # Extract expression data
                if genes_in_rows:
                    # Assuming gene IDs in first column, expression in second
                    gene_col = df.columns[0]
                    expr_col = df.columns[1]
                    expr = df.set_index(gene_col)[expr_col]
                else:
                    # Assuming genes are columns, one row of expression
                    expr = pd.Series(df.iloc[0].values, index=df.columns)
                
                # Update unique genes
                unique_genes.update(expr.index)
                
                # Store data
                batch_data.append(expr)
                batch_metadata.append(meta)
                
            except Exception as e:
                print(f"Error processing {file_path}: {e}")
        
        # Add batch results to overall results
        all_data.extend(batch_data)
        all_metadata.extend(batch_metadata)
        
        print(f"Batch {batch_idx + 1} complete: {len(batch_data)} files processed successfully")
    
    return all_data, all_metadata, unique_genes

# Process all files
print("\nProcessing all files...")
all_series, all_metadata, all_genes = process_files(csv_gz_files, genes_in_rows=genes_in_rows)

print(f"Successfully processed {len(all_series)} files")
print(f"Found {len(all_genes)} unique genes")

if not all_series:
    print("No data was processed successfully. Exiting.")
    exit(1)

# Create combined expression matrix
print("Creating combined expression matrix...")
all_genes_list = sorted(list(all_genes))

# Create a matrix with all genes for all samples
expr_matrix = pd.DataFrame(0, index=range(len(all_series)), columns=all_genes_list)

# Fill in the expression values
print("Populating expression matrix...")
for i, expr in enumerate(tqdm(all_series)):
    expr_matrix.loc[i, expr.index] = expr.values

# Create metadata DataFrame
metadata_df = pd.DataFrame(all_metadata)

# Create cell IDs combining plate and well info
cell_ids = []
for i, meta in enumerate(all_metadata):
    plate = meta.get('plate', 'unknown')
    well_letter = meta.get('well_letter', '')
    well_number = meta.get('well_number', '')
    cell_id = f"{meta['sample_id']}_{plate}_{well_letter}_{well_number}"
    cell_ids.append(cell_id)

# Set cell IDs as indices
expr_matrix.index = metadata_df.index = cell_ids

# Save the expression matrix and metadata
expr_file = os.path.join(processed_dir, "velten_expression_matrix.csv")
metadata_file = os.path.join(processed_dir, "velten_metadata.csv")

print(f"Saving expression matrix to: {expr_file}")
expr_matrix.to_csv(expr_file)

print(f"Saving metadata to: {metadata_file}")
metadata_df.to_csv(metadata_file)

# Create AnnData object
print("Creating AnnData object...")
adata = AnnData(X=expr_matrix.values, 
                obs=metadata_df, 
                var=pd.DataFrame(index=all_genes_list))

# Save raw AnnData object
adata_file = os.path.join(processed_dir, "velten_hematopoiesis_raw.h5ad")
print(f"Saving AnnData object to: {adata_file}")
adata.write(adata_file)

print("Processing complete!")