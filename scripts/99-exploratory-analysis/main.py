# %%
from ai4bmr_datasets import PCa
from pathlib import Path
from dotenv import load_dotenv
import os
assert load_dotenv(override=True)

# %%
# check documentation: https://github.com/AI4SCR/ai4bmr-datasets
ds = PCa(base_dir=Path(os.environ['BASE_DIR']), image_version='filtered', mask_version='annotated',
         load_metadata=True, load_intensity=True, align=True)
ds.setup(engine='pyarrow')

sample_id = ds.sample_ids[0]
len(ds.sample_ids)

ds.images[sample_id]  # lazy image obj
ds.images[sample_id].data  # actual image
ds.panel  # panel
ds.clinical # clinical metadata

# %% SAMPLE-LEVEL AnnDatas
import athena as ath
from anndata import AnnData

X = ds.intensity.loc[sample_id].values

obs = ds.metadata.loc[sample_id]
obs = obs.astype('category')
obs.index = obs.index.astype(str)

ad = AnnData(X=X, obs=obs, var=ds.panel)
ad.uns['mask'] = ds.masks[sample_id].data

ath.pp.compute_centroids(ad)
ath.graph.build_graph(ad, topology='radius', radius=32, include_self=False, graph_key='radius32')
# ROI-level score
ath.metrics.shannon(ad, attr='label', local=False)
ad.uns['shannon_label']  # sample-level scores are stored in uns

# Cell-level score
ath.metrics.shannon(ad, attr='label', graph_key='radius32', local=True)
ad.obs.shannon_label_radius32  # cell-level scores are stored in obs
ax = ath.pl.spatial(ad, attr='shannon_label_radius32')