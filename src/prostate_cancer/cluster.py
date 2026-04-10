# %%
import igraph as ig
import leidenalg as la
import numpy as np
import pandas as pd
import scipy.sparse
import umap
from ai4bmr_core.utils.logging import get_logger
from numba import cuda

logger = get_logger("02.0_clustering", verbose=1)

from anndata import AnnData

if cuda.is_available():
    logger.info("GPU available")
    import rapids_singlecell as sc
else:
    import scanpy as sc


def cluster(
    data: pd.DataFrame,
    resolution: float = 1,
    graph_method: str = "umap",
    n_neighbors: int = 15,
    cluster_method: str = "leiden",
):

    logger.info(f"cluster data with {len(data)} observations")

    # %% create AnnData
    obs = data.index.to_frame()
    tmp_ = ["-".join(i) for i in zip(obs["sample_name"], obs["object_id"].astype(str))]
    obs = obs.assign(id=tmp_).set_index("id")
    x = data.reset_index(drop=True).assign(id=tmp_).set_index("id")
    ad = AnnData(X=x, obs=obs)

    if cuda.is_available():
        sc.get.anndata_to_GPU(ad)

    # %%
    logger.info(f"compute graph with method `{graph_method}`")

    if graph_method == "umap":
        min_dist = 0.1
        metric = "euclidean"
        n_components = 2

        mapper = umap.UMAP(
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            metric=metric,
            n_components=n_components,
            n_jobs=-1,
        )
        mapper.fit(data)
        csr = mapper.graph_
    elif graph_method == "knn":
        raise NotImplementedError()
        # csr = kneighbors_graph(data, n_neighbors=n_neighbors, mode='distance', metric=metric, include_self=False)
        # assert False, 'we need to convert distance to weights'
    elif graph_method == "scanpy":
        sc.pp.neighbors(ad, n_neighbors=n_neighbors, use_rep="X")
        csr = ad.obsp["connectivities"]

    def csr_to_ig(csr: scipy.sparse.csr_matrix, directed=True, weighted=True):
        rows, cols = csr.nonzero()

        g = ig.Graph(directed=directed)
        g.add_vertices(csr.shape[0])
        g.add_edges(list(zip(rows, cols)))

        if weighted:
            weights = csr.data
        else:
            weights = np.ones_like(csr.data)

        g.es["weight"] = weights
        return g

    G = csr_to_ig(csr)

    # %% 02.0_clustering

    logger.info(f"compute 02.0_clustering with method `{cluster_method}`")

    def get_leiden_membership(G, resolution: float = 1):
        partition = la.find_partition(
            G,
            la.RBConfigurationVertexPartition,
            weights=G.es["weight"],
            resolution_parameter=resolution,
        )
        return partition.membership

    def get_igraph_membership(G, resolution: float = 1):
        return G.community_leiden(
            objective_function="modularity", resolution=resolution, weights="weight"
        ).membership

    def get_kmeans_membership(data, n_clusters: int = 20):
        from sklearn.cluster import KMeans

        kmeans = KMeans(n_clusters=n_clusters)
        kmeans.fit(data)
        return kmeans.labels_

    if cluster_method == "leiden":
        membership = get_leiden_membership(G, resolution=resolution)
    elif cluster_method == "igraph":
        membership = get_igraph_membership(G, resolution=resolution)
    elif cluster_method == "kmeans":
        n_clusters = 20
        membership = get_kmeans_membership(data, n_clusters=n_clusters)
    elif cluster_method == "scanpy":
        assert graph_method == "scanpy"
        sc.tl.leiden(ad, resolution=resolution)
        # sc.tl.leiden(ad, resolution=resolution, flavor="igraph", n_iterations=2)
        membership = ad.obs["leiden"].astype(str).values

    membership = pd.Categorical([str(i) for i in membership])

    num_cluster = np.unique(membership).size
    data = data.assign(membership=membership).set_index("membership", append=True)

    # %%
    # from sklearn.metrics import silhouette_score
    # logger.info('compute silhouette score')
    #
    # score = silhouette_score(data, data.index.get_level_values('membership'))
    #
    # result = dict(score=score, n_cluster=n_cluster)
    result = dict(num_obs=len(data), num_cluster=num_cluster)

    # %% embedding
    logger.info("compute embedding")

    # n_neighbors = 15
    # min_dist = 0.1
    # metric = 'euclidean'
    # n_components = 2

    if "neighbors" not in ad.uns.keys():
        sc.pp.neighbors(ad, use_rep="X")
    sc.tl.umap(ad)

    if cuda.is_available():
        sc.get.anndata_to_CPU(ad)

    embedding = ad.obsm["X_umap"]

    # mapper = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, metric=metric, n_components=n_components)
    # mapper.fit(data)

    # %%
    logger.info("completed")
    return data, embedding, result
