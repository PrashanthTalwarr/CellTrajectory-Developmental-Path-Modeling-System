# celltrajectory/trajectory/pseudotime.py

import numpy as np
import scanpy as sc
from anndata import AnnData
from typing import Optional, List, Tuple, Dict

def compute_pseudotime(
    adata: AnnData,
    root_cells: Optional[np.ndarray] = None,
    root_cluster: Optional[str] = None,
    n_dcs: int = 15,
    use_leiden: bool = True
) -> AnnData:
    """
    Compute pseudotime ordering of cells.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    root_cells : np.ndarray, optional
        Boolean mask or indices of root cells.
    root_cluster : str, optional
        Cluster label to use as root cells.
    n_dcs : int, default=15
        Number of diffusion components to use.
    use_leiden : bool, default=True
        Whether to use leiden clusters if root_cluster is specified.
        
    Returns
    -------
    adata : AnnData
        Annotated data matrix with pseudotime information added.
    """
    # Make a copy to avoid modifying original
    adata = adata.copy()
    
    # Calculate diffusion map
    sc.tl.diffmap(adata, n_comps=n_dcs)
    
    # Determine root cells
    if root_cells is None and root_cluster is not None:
        if use_leiden:
            cluster_key = 'leiden'
        else:
            cluster_key = 'louvain'
            
        if cluster_key not in adata.obs:
            raise ValueError(f"Cluster key '{cluster_key}' not found in adata.obs")
            
        root_cells = adata.obs[cluster_key] == root_cluster
        
    if root_cells is None:
        # Use cell with lowest diffusion component 1 as root
        root_idx = np.argmin(adata.obsm['X_diffmap'][:, 1])
        adata.uns['iroot'] = root_idx
    elif isinstance(root_cells, np.ndarray) and root_cells.dtype == bool:
        # Convert boolean mask to index
        root_idx = np.where(root_cells)[0][0]
        adata.uns['iroot'] = root_idx
    else:
        # Assume root_cells is an index
        adata.uns['iroot'] = root_cells
    
    # Compute DPT pseudotime
    sc.tl.dpt(adata)
    
    return adata


def identify_branch_points(
    adata: AnnData,
    n_branches: int = 2,
    pseudotime_key: str = 'dpt_pseudotime'
) -> Tuple[Dict[str, np.ndarray], AnnData]:
    """
    Identify branch points in cellular trajectories.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with pseudotime information.
    n_branches : int, default=2
        Number of branches to identify.
    pseudotime_key : str, default='dpt_pseudotime'
        Key for pseudotime values in adata.obs.
        
    Returns
    -------
    branches : Dict[str, np.ndarray]
        Dictionary of branch indices.
    adata : AnnData
        Annotated data matrix with branch information added.
    """
    # Make sure we have pseudotime
    if pseudotime_key not in adata.obs:
        raise ValueError(f"Pseudotime key '{pseudotime_key}' not found in adata.obs")
    
    # Make sure we have diffusion map
    if 'X_diffmap' not in adata.obsm:
        sc.tl.diffmap(adata)
    
    # Compute DPT with branching
    sc.tl.dpt(adata, n_branchings=n_branches)
    
    # Extract branch information
    branches = {}
    
    # Get DPT branch indices
    if 'dpt_groups' in adata.obs:
        for branch_id in adata.obs['dpt_groups'].unique():
            branches[f'branch_{branch_id}'] = np.where(adata.obs['dpt_groups'] == branch_id)[0]
    
    return branches, adata


def compute_trajectory_graph(
    adata: AnnData,
    use_paga: bool = True,
    use_umap: bool = True,
    root_cluster: Optional[str] = None
) -> AnnData:
    """
    Compute a graph representation of the trajectory.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    use_paga : bool, default=True
        Whether to use PAGA for trajectory graph.
    use_umap : bool, default=True
        Whether to use UMAP for visualization.
    root_cluster : str, optional
        Cluster to use as root.
        
    Returns
    -------
    adata : AnnData
        Annotated data matrix with trajectory graph information.
    """
    # Make sure we have clusters
    if 'leiden' not in adata.obs:
        sc.tl.leiden(adata)
    
    # Compute PAGA graph
    if use_paga:
        sc.tl.paga(adata, groups='leiden')
        
        # Root the graph if specified
        if root_cluster is not None:
            root_idx = np.where(adata.obs['leiden'] == root_cluster)[0][0]
            adata.uns['iroot'] = root_idx
            
            # Compute pseudotime
            sc.tl.dpt(adata)
    
    # Compute force-directed layout for visualization
    if use_umap:
        if 'X_umap' not in adata.obsm:
            sc.tl.umap(adata)
    
    return adata