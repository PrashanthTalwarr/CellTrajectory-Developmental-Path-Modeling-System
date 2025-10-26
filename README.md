# CellTrajectory-Developmental-Path-Modeling-System


[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

A computational framework for modeling cellular differentiation trajectories in single-cell RNA-seq data, with emphasis on rare population identification and developmental path prediction.

## Overview

CellTrajectory is a Python-based analytical framework designed to address key challenges in understanding cellular differentiation:

- **Tracking Development**: Reconstructs continuous developmental trajectories from stem cells to mature lineages
- **Rare Population Analysis**: Identifies and characterizes rare cell populations (<1% of cells) 
- **Path Prediction**: Maps developmental origins of specialized cell types back to progenitor populations
- **Trajectory Visualization**: Provides publication-quality visualizations of differentiation paths

This framework integrates graph-based trajectory inference with isolation-based rare population detection to provide a comprehensive analysis of cellular development.

![CellTrajectory Overview](https://raw.githubusercontent.com/yourusername/CellTrajectory/main/data/processed/results/figure8_combined_visualization.png)

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/CellTrajectory.git
cd CellTrajectory

# Create and activate a conda environment
conda create -n celltrajectory python=3.9
conda activate celltrajectory

# Install dependencies
pip install -r requirements.txt
```

## Dependencies

- Python 3.9+
- scanpy
- anndata
- numpy
- pandas
- matplotlib
- seaborn
- leidenalg

## Usage

CellTrajectory is structured as a modular analysis pipeline:

### 1. Data Preprocessing

```python
import scanpy as sc
from celltrajectory.preprocessing import preprocess_data

# Load data
adata = sc.read("path/to/your/data.h5ad")

# Preprocess data
adata = preprocess_data(adata, min_genes=200, min_cells=3)
```

### 2. Trajectory Analysis

```python
from celltrajectory.trajectory import compute_pseudotime, compute_trajectory_graph

# Calculate pseudotime
adata = compute_pseudotime(adata, root_cluster='HSC')

# Create trajectory graph
adata = compute_trajectory_graph(adata)
```

### 3. Rare Population Detection

```python
from celltrajectory.trajectory import identify_rare_populations, predict_developmental_paths

# Find rare populations
rare_pops, metrics = identify_rare_populations(
    adata, 
    min_cells=10, 
    isolation_score_threshold=0.8
)

# Predict developmental path for a rare population
path_indices = predict_developmental_paths(
    adata, 
    rare_cluster=rare_pops[0], 
    root_cluster='HSC'
)
```

### 4. Visualization

```python
import matplotlib.pyplot as plt
import scanpy as sc

# Visualize cell types
sc.pl.umap(adata, color='predicted_cell_type')

# Visualize developmental pseudotime
sc.pl.umap(adata, color='dpt_pseudotime', cmap='viridis')

# Visualize trajectory
sc.pl.paga(adata, color='leiden', threshold=0.03)

# Visualize rare populations
sc.pl.umap(adata, color='is_rare', palette=['lightgray', 'red'])
```

## Example Notebooks

The repository includes example Jupyter notebooks demonstrating the complete workflow:

1. [`01_preprocessing.ipynb`](notebooks/01_preprocessing.ipynb) - Data loading, QC, and preprocessing
2. [`02_trajectory_inference.ipynb`](notebooks/02_trajectory_inference.ipynb) - Trajectory analysis and rare population detection
3. [`03_visualization.ipynb`](notebooks/03_visualization.ipynb) - Results visualization and report generation

## Results

Applied to the Velten et al. (2017) human hematopoiesis dataset, CellTrajectory:

- Identified 8 distinct cell clusters in the hematopoietic hierarchy
- Successfully reconstructed developmental trajectories through pseudotime analysis
- Identified rare cell populations (cluster 7) comprising 1% of total cells
- Traced developmental origins of specialized cell types back to progenitor cells

| Cluster | Count | Percentage |
|---------|-------|------------|
| 0 | 440 | 19.32% |
| 1 | 437 | 19.19% |
| 2 | 347 | 15.24% |
| 3 | 296 | 13.00% |
| 4 | 258 | 11.33% |
| 5 | 250 | 10.98% |
| 6 | 226 | 9.93% |
| 7 | 23 | 1.01% |

## Project Structure

```
CellTrajectory/
├── data/                       # Data directory
│   ├── raw/                    # Raw data files
│   └── processed/              # Processed data and results
├── notebooks/                  # Example Jupyter notebooks
├── celltrajectory/             # Main package
│   ├── preprocessing/          # Data preprocessing modules
│   ├── trajectory/             # Trajectory inference modules
│   └── visualization/          # Visualization utilities
├── scripts/                    # Utility scripts
├── tests/                      # Test suite
├── LICENSE                     # License file
├── README.md                   # This file
└── requirements.txt            # Package dependencies
```

## Methods

### Rare Population Detection

CellTrajectory identifies rare populations through a two-step process:
1. **Size-based filtering**: Identifies clusters with fewer than a threshold number of cells
2. **Isolation scoring**: Calculates how isolated each small cluster is from other clusters in the neighborhood graph

### Developmental Path Prediction

To predict developmental paths:
1. **Pseudotime calculation**: Orders cells along a continuous progression
2. **Root and target identification**: Selects source (stem cell) and target (differentiated/rare) populations
3. **Waypoint selection**: Identifies representative cells along the developmental continuum
4. **Path construction**: Creates a connected path through the high-dimensional space

## Citation

If you use CellTrajectory in your research, please cite:



## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Human hematopoiesis dataset from Velten et al., *Nature Cell Biology* (2017)
- Based on methodologies from the field of computational single-cell genomics
- Inspired by the need to better understand rare cell populations in development
