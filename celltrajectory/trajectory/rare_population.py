# celltrajectory/trajectory/rare_population.py

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from typing import Optional, List, Tuple, Dict

def identify_rare_populations(
    adata: AnnData,
    min_cells: int = 10,
    isolation_score_threshold: float = 0.8,
    cluster_key: str = 'leiden'
) -> Tuple[List[str], pd.DataFrame]:
    """
    Identify rare cell populations based on size and isolation metrics.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    min_cells : int, default=10
        Minimum number of cells to be considered as a population.
    isolation_score_threshold : float, default=0.8
        Threshold for isolation score to identify rare populations.
    cluster_key : str, default='leiden'
        Key in adata.obs for cluster labels.
        
    Returns
    -------
    rare_populations : List[str]
        List of cluster IDs identified as rare populations.
    rare_pop_metrics : pd.DataFrame
        DataFrame with metrics for all small clusters.
    """
    # Verify cluster key exists
    if cluster_key not in adata.obs:
        raise ValueError(f"Cluster key '{cluster_key}' not found in adata.obs")
    
    # Calculate cluster sizes
    cluster_sizes = adata.obs[cluster_key].value_counts()
    small_clusters = cluster_sizes[cluster_sizes < min_cells].index.tolist()
    
    if not small_clusters:
        print("No small clusters found")
        return [], pd.DataFrame()
    
    # Make sure we have neighborhood graph
    if 'neighbors' not in adata.uns:
        sc.pp.neighbors(adata)
    
    # Get connectivity matrix
    connectivities = adata.obsp['connectivities']
    
    # Calculate isolation metrics for each small cluster
    isolation_scores = {}
    
    for cluster in small_clusters:
        # Get cells in this cluster
        cluster_mask = adata.obs[cluster_key] == cluster
        
        if sum(cluster_mask) == 0:
            continue
            
        # Get indices of cells in this cluster
        cluster_indices = np.where(cluster_mask)[0]
        
        # Calculate connections within cluster
        within_connections = connectivities[np.ix_(cluster_indices, cluster_indices)].sum()
        
        # Calculate total connections from this cluster
        total_connections = connectivities[cluster_indices].sum()
        
        # Calculate isolation score
        if total_connections > 0:
            # Lower score means more isolated
            isolation_scores[cluster] = 1 - (within_connections / total_connections)
        else:
            isolation_scores[cluster] = 1.0
    
    # Identify rare and isolated populations
    rare_populations = [
        cluster for cluster in small_clusters
        if isolation_scores.get(cluster, 0) > isolation_score_threshold
    ]
    
    # Create results DataFrame
    results = pd.DataFrame({
        'cluster': small_clusters,
        'size': [cluster_sizes[c] for c in small_clusters],
        'isolation_score': [isolation_scores.get(c, 0) for c in small_clusters],
        'is_rare': [c in rare_populations for c in small_clusters]
    })
    
    return rare_populations, results


def predict_developmental_paths(
    adata: AnnData,
    rare_cluster: str,
    root_cluster: str = 'HSC',
    n_waypoints: int = 10,
    pseudotime_key: str = 'dpt_pseudotime',
    cluster_key: str = 'leiden'
) -> List[int]:
    """
    Predict developmental path from root cluster to rare population.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    rare_cluster : str
        Cluster ID of the rare population.
    root_cluster : str, default='HSC'
        Cluster ID of the root population.
    n_waypoints : int, default=10
        Number of cells to include in the path.
    pseudotime_key : str, default='dpt_pseudotime'
        Key for pseudotime values in adata.obs.
    cluster_key : str, default='leiden'
        Key for cluster labels in adata.obs.
        
    Returns
    -------
    path : List[int]
        Indices of cells forming the developmental path.
    """
    # Verify keys exist
    if cluster_key not in adata.obs:
        raise ValueError(f"Cluster key '{cluster_key}' not found in adata.obs")
    
    # Get cells in root and target clusters
    root_indices = np.where(adata.obs[cluster_key] == root_cluster)[0]
    target_indices = np.where(adata.obs[cluster_key] == rare_cluster)[0]
    
    if len(root_indices) == 0:
        raise ValueError(f"Root cluster '{root_cluster}' has no cells")
    if len(target_indices) == 0:
        raise ValueError(f"Target cluster '{rare_cluster}' has no cells")
    
    # Set root for pseudotime calculation
    adata.uns['iroot'] = root_indices[0]
    
    # Calculate pseudotime if not already present
    if pseudotime_key not in adata.obs or np.isnan(adata.obs[pseudotime_key]).all():
        sc.tl.dpt(adata)
    
    # Get pseudotime values
    pseudotime = adata.obs[pseudotime_key].values
    
    # Find a representative path between root and target
    # Use the cell with median pseudotime in each cluster
    root_pt = np.median(pseudotime[root_indices])
    target_pt = np.median(pseudotime[target_indices])
    
    # Create evenly spaced pseudotime values between root and target
    path_pts = np.linspace(root_pt, target_pt, n_waypoints)
    
    # Find cells closest to these pseudotime values
    path_indices = []
    for pt in path_pts:
        # Find the cell closest to this pseudotime
        closest_cell = np.argmin(np.abs(pseudotime - pt))
        path_indices.append(closest_cell)
    
    return path_indices
