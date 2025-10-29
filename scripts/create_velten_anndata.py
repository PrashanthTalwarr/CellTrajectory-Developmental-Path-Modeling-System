import os
import glob
import pandas as pd
import numpy as np
import scanpy as sc
from anndata import AnnData
from tqdm import tqdm

# Set up paths
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(base_dir, "data")
raw_dir = os.path.join(data_dir, "raw", "GSE75478_extracted")
processed_dir = os.path.join(data_dir, "processed")
os.makedirs(processed_dir, exist_ok=True)

print(f"Looking for CSV files in: {raw_dir}")

# Get all CSV files
csv_files = sorted(glob.glob(os.path.join(raw_dir, "*.csv")))
if not csv_files:
    print(f"No CSV files found in {raw_dir}")
    exit(1)

print(f"Found {len(csv_files)} CSV files")

# Function to read a CSV file with flexible delimiter detection
def read_file(file_path):
    """Read a CSV file with automatic delimiter detection"""
    try:
        return pd.read_csv(file_path)
    except:
        try:
            return pd.read_csv(file_path, sep="\t")
        except:
            try:
                return pd.read_csv(file_path, sep=None, engine="python")
            except Exception as e:
                print(f"Error reading {file_path}: {e}")
                return None

# Sample a few files to determine the structure
print("Sampling files to determine structure...")
sample_files = csv_files[:min(5, len(csv_files))]
sample_dfs = []

for file in sample_files:
    df = read_file(file)
    if df is not None:
        sample_dfs.append(df)
        print(f"File: {os.path.basename(file)}, Shape: {df.shape}")
        print(f"Columns: {df.columns.tolist()}")
        print(f"First few rows:\n{df.head(2)}\n")

if not sample_dfs:
    print("Failed to read any sample files.")
    exit(1)

# Extract metadata from filename
def extract_metadata(filename):
    """Extract metadata from filename"""
    base_name = os.path.basename(filename)
    parts = base_name.split("_")
    
    metadata = {
        "sample_id": parts[0],  # e.g., GSM1955701
    }
    
    # Extract plate and well information if available
    if len(parts) > 2 and parts[2].startswith("plate"):
        metadata["plate"] = parts[2]
        
    if len(parts) > 3:
        metadata["well"] = "_".join(parts[3:]).split(".")[0]  # Remove file extension
        
    return metadata

# Determine data structure based on samples
first_df = sample_dfs[0]
col_names = first_df.columns.tolist()

# Check if genes are in rows or columns
if len(first_df.columns) < len(first_df):
    # More rows than columns - genes likely in rows
    print("Detected data format: Genes in rows, samples in columns")
    genes_in_rows = True
else:
    # More columns than rows - genes likely in columns
    print("Detected data format: Genes in columns, samples in rows")
    genes_in_rows = False

# Process all files
print("Processing all files...")
all_data = []
all_metadata = []

for file_path in tqdm(csv_files):
    try:
        # Read the file
        df = read_file(file_path)
        if df is None:
            continue
            
        # Extract metadata
        meta = extract_metadata(file_path)
        
        # Extract expression data
        if genes_in_rows:
            # Assuming first column contains gene identifiers
            if len(df.columns) >= 2:  # Make sure there's data beyond the gene column
                # Determine which column is gene ID and which is expression
                gene_col = df.columns[0]
                expr_col = df.columns[1]
                
                # Set gene IDs as index and get expression values
                expr = df.set_index(gene_col)[expr_col]
                
                # Store as Series
                all_data.append(expr)
                all_metadata.append(meta)
        else:
            # Genes are columns
            # Assuming one row per sample
            expr = df.iloc[0] if df.shape[0] > 0 else None
            
            if expr is not None:
                all_data.append(expr)
                all_metadata.append(meta)
    except Exception as e:
        print(f"Error processing {file_path}: {e}")

print(f"Successfully processed {len(all_data)} files")

if not all_data:
    print("No data was processed successfully.")
    exit(1)

# Create a combined expression matrix
print("Creating combined expression matrix...")

# Determine all unique genes
all_genes = set()
for expr in all_data:
    all_genes.update(expr.index)
    
all_genes = sorted(list(all_genes))
print(f"Found {len(all_genes)} unique genes")

# Create a matrix with all genes for all samples
expr_matrix = pd.DataFrame(0, index=range(len(all_data)), columns=all_genes)

# Fill in the expression values
for i, expr in enumerate(tqdm(all_data)):
    expr_matrix.loc[i, expr.index] = expr.values
    
# Create metadata DataFrame
metadata_df = pd.DataFrame(all_metadata)

# Ensure matching indices
expr_matrix.index = metadata_df.index = [f"cell_{i}" for i in range(len(all_metadata))]

# Save the expression matrix and metadata
expr_file = os.path.join(processed_dir, "velten_expression_matrix.csv")
metadata_file = os.path.join(processed_dir, "velten_metadata.csv")

print(f"Saving expression matrix to: {expr_file}")
expr_matrix.to_csv(expr_file)

print(f"Saving metadata to: {metadata_file}")
metadata_df.to_csv(metadata_file)

# Create AnnData object
print("Creating AnnData object...")
adata = AnnData(X=expr_matrix.values, obs=metadata_df, var=pd.DataFrame(index=all_genes))

# Save raw AnnData object
adata_file = os.path.join(processed_dir, "velten_hematopoiesis_raw.h5ad")
print(f"Saving AnnData object to: {adata_file}")
adata.write(adata_file)

print("Processing complete!")
