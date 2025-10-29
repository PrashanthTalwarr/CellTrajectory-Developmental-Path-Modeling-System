# scripts/scan_velten_data.py
import os
import glob
import pandas as pd
import numpy as np

# Set up paths
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(base_dir, "Data compressed")

print(f"Scanning directory: {data_dir}")

# Get all files
all_files = sorted(glob.glob(os.path.join(data_dir, "*")))

# Group files by size
size_groups = {}
for file_path in all_files:
    file_size = os.path.getsize(file_path)
    if file_size not in size_groups:
        size_groups[file_size] = []
    size_groups[file_size].append(file_path)

print("\nFile size distribution:")
for size, files in sorted(size_groups.items()):
    print(f"{size} bytes: {len(files)} files")
    # Print examples of files in each size group
    if len(files) > 0:
        print(f"  Examples: {os.path.basename(files[0])}")
        if len(files) > 1:
            print(f"           {os.path.basename(files[1])}")
        print()

# Look for common file extensions
extensions = {}
for file_path in all_files:
    ext = os.path.splitext(file_path)[1].lower()
    if ext not in extensions:
        extensions[ext] = []
    extensions[ext].append(file_path)

print("\nFile extensions found:")
for ext, files in sorted(extensions.items()):
    print(f"{ext}: {len(files)} files")
    
# Look for specific file patterns
patterns = {
    "count_matrix": ["count", "matrix", "expression"],
    "metadata": ["meta", "annotation", "phenotype"],
    "barcodes": ["barcode", "cell"],
    "features": ["feature", "gene"],
}

pattern_matches = {}
for pattern_name, keywords in patterns.items():
    pattern_matches