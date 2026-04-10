from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def compute_interaction_matrices(adata, cell_types = None, graph_key='radius_32'):
    assert adata.n_obs > 15, "Anndata object has too few observations to compute interactions."
    if cell_types is None:
        cell_types = adata.obs['label'].unique().tolist()
        cell_types.sort()
    interaction_matrix = pd.DataFrame(0, index=cell_types, columns=cell_types)
    # compute interaction matrix for this niche
    adj = adata.obsp[graph_key].toarray()
    np.fill_diagonal(adj, 0)  # remove self-loops
    from scipy.sparse import csr_matrix
    sparse_adj = csr_matrix(adj)
    adata.obsp[f"{graph_key}_no_self_loops"] = sparse_adj
    from athena.utils.general import get_nx_graph_from_anndata
    g = get_nx_graph_from_anndata(adata, key=f"{graph_key}_no_self_loops")
    edges = list(g.edges()) 
    sample_interaction = interaction_matrix.copy()
    for edge in edges:
        (c1, c2) = edge

        (type1, type2) = (adata.obs.loc[c1]['label'], adata.obs.loc[c2]['label'])
        #print(f"Edge between cell {c1} (type: {type1}) and cell {c2} (type: {type2})")
        ## sort type 1 and type 2 alphabetically to make the matrix symmetric
        if type1 > type2:
            (type1, type2) = (type2, type1)
        sample_interaction.loc[type1, type2] += 1
    num_edges = len(edges)
    print(f"Total number of edges in graph: {num_edges}")
    sample_interaction_freqs = sample_interaction / num_edges
    # set lower triangle to nan for better visualization
    # assuming sample_interaction_freqs is a pandas DataFrame
    sample_interaction_freqs = sample_interaction_freqs.T
    mask = np.triu(np.ones_like(sample_interaction_freqs, dtype=bool), k=1)
    sample_interaction_freqs = sample_interaction_freqs.mask(mask)

    return sample_interaction_freqs


def plot_interaction_matrix(interaction_freqs, title=None):
    plt.figure(figsize=(16, 12))
    sns.heatmap(interaction_freqs, annot=True, fmt=".2f", cmap="Reds", cbar_kws={'label': 'Interaction Frequency'})
    plt.title(title)
    plt.xlabel('Cell Type 2')
    plt.ylabel('Cell Type 1')
    plt.tight_layout()
    return plt.gcf()  # return the current figure



