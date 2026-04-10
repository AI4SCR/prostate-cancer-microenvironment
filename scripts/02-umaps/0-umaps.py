from ai4bmr_datasets.datasets import PCa
from ai4bmr_learn.plotting.umap import run_umap
from ai4bmr_core.utils.sampling import sample_min_per_group_then_uniform
import pandas as pd
from prostate_cancer.utils import normalize
from pathlib import Path
import pickle
from loguru import logger
from prostate_cancer.utils import get_colormap_dict
from itertools import product
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.cm import ScalarMappable
import matplotlib
from matplotlib.colors import Normalize

save_dir = Path('/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/0-paper/2-umaps/0-all-cells')
# params = list(product(['cuml', 'umap-learn'], [15, 30, 50], [0.1, 0.3, 0.5]))
# params = list(product(
#     ['umap-learn'],
#     [15, 30, 50],
#     [0.1, 0.3, 0.5],
#     [None, ['dna1', 'dna2', 'icsk1', 'icsk2', 'icsk3', 'fap']]
# ))
params = list(product(['umap-learn'], [50], [0.1], [['dna1', 'dna2', 'icsk1', 'icsk2', 'icsk3', 'fap']]))

ds = PCa(image_version='filtered', mask_version='annotated',
         load_metadata=True, load_intensity=True)
ds.setup()

assert set(ds.intensity.index.get_level_values('sample_id')) - set(ds.clinical.index) == set()
assert set(ds.metadata.index.get_level_values('sample_id')) - set(ds.clinical.index) == set()
assert len(ds.intensity) == len(ds.metadata)

data = ds.intensity.copy()
# data = data.sample(10_000)

metadata, data = ds.metadata.align(data, axis=0, join='inner')
data = normalize(data, exclude_zeros=True)

sid_to_pid = ds.clinical['pat_id'].to_dict()


# metadata['pat_id'] = metadata.index.get_level_values('sample_id').map(sid_to_pid)

def plot_points(*, reducer, subset_points=None,
                values=None, cmap=None, cbar: bool = False,
                labels=None, color_key=None,
                all_points: bool = False,
                shuffle: bool = True, order=None):
    fig, ax = plt.subplots()
    if all_points:
        x = reducer.embedding_[:, 0]
        y = reducer.embedding_[:, 1]
        c = to_rgba('#D3D3D3', alpha=.1)
        ax.scatter(x=x, y=y, s=0.1, color=c)

    if order is not None:
        x = reducer.embedding_[order, 0]
        y = reducer.embedding_[order, 1]
    else:
        x = reducer.embedding_[:, 0]
        y = reducer.embedding_[:, 1]

    if subset_points is not None:
        x = x[subset_points]
        y = y[subset_points]

    if labels is not None:
        l = labels[subset_points]
        c = np.array([to_rgba(color_key[i], alpha=1) for i in l])
    elif values is not None:
        v = values[subset_points]
        c = cmap(v)
    else:
        raise ValueError('Must provide labels or values')

    if shuffle:
        order = np.arange(len(x))
        np.random.shuffle(order)

        x, y = x[order], y[order]
        c = c[order]

    ax.scatter(x=x, y=y, s=1, c=c)
    ax.set_axis_off()

    if values is not None and cmap is not None and cbar:
        sm = ScalarMappable(cmap=cmap)
        sm.set_array(values[subset_points])  # ensures proper range
        fig.colorbar(sm, ax=ax, orientation='vertical', fraction=0.046, pad=0.04)

    return ax


def compute_umap(data: pd.DataFrame, save_path: Path, n_neighbors: int, min_dist: float, engine: str):
    reducer = run_umap(data=data, n_neighbors=n_neighbors, min_dist=min_dist, engine=engine)

    container = {'index': data.index,
                 'reducer': reducer,
                 'engine': engine,
                 'min_dist': min_dist,
                 'n_neighbors': n_neighbors}

    save_path.parent.mkdir(parents=True, exist_ok=True)
    with open(save_path, 'wb') as f:
        pickle.dump(container, f)

    return container


# %% COMPUTE UMAPs
def get_reducer_path(save_dir, n_neighbors: int, min_dist: float, engine: str, exclude_markers):
    excl = '-excl_markers=' + '_'.join(exclude_markers) if exclude_markers is not None else ''
    name = f'n_neighbors={n_neighbors}-min_dist={min_dist}-engine={engine}' + excl
    return save_dir / name / f'reducer.pkl'


for (engine, n_neighbors, min_dist, exclude_markers) in params:
    # break
    logger.info(
        f'Computing UMAP for {len(data)} data points with n_neighbors {n_neighbors}, min_dist {min_dist} and engine {engine}. Excluding: {exclude_markers}')

    reducer_path = get_reducer_path(save_dir=save_dir,
                                    n_neighbors=n_neighbors, min_dist=min_dist, engine=engine,
                                    exclude_markers=exclude_markers)
    if reducer_path.exists():
        continue
    else:
        if exclude_markers is not None:
            tmp = data.loc[:, ~data.columns.isin(exclude_markers)].copy()
        else:
            tmp = data.copy()

        container = compute_umap(data=tmp, save_path=reducer_path, n_neighbors=n_neighbors, min_dist=min_dist,
                                 engine=engine)

# %%
for label in ['main_group', 'label', 'pat_id']:
    for (engine, n_neighbors, min_dist, exclude_markers) in params:
        # break
        logger.info(
            f'UMAP for {label} with {len(data)} data points with n_neighbors {n_neighbors}, min_dist {min_dist} and engine {engine}. Excluding: {exclude_markers}')

        reducer_path = get_reducer_path(save_dir=save_dir,
                                        n_neighbors=n_neighbors, min_dist=min_dist, engine=engine,
                                        exclude_markers=exclude_markers)
        with open(reducer_path, 'rb') as f:
            container = pickle.load(f)

        index = container['index']
        reducer = container['reducer']

        colormap_dict = get_colormap_dict(name=label)

        metadata = ds.metadata.loc[index].copy()
        metadata['pat_id'] = metadata.index.get_level_values('sample_id').map(sid_to_pid)

        labels = metadata[label].values
        values = None

        # subset data
        grouped = metadata.groupby(label)
        num_samples = 100_000
        min_per_group = int(num_samples / grouped.ngroups)
        min_per_group = min(min_per_group, grouped.size().min())
        metadata = sample_min_per_group_then_uniform(grouped=grouped, n=num_samples, min_per_group=min_per_group,
                                                     random_state=0)

        subset_points = np.zeros(len(index), dtype=bool)
        keep_idc = index.isin(metadata.index)
        subset_points[keep_idc] = True

        save_path = reducer_path.parent / f'label={label}.pdf'
        ax = plot_points(reducer=reducer, subset_points=subset_points, labels=labels, color_key=colormap_dict)
        ax.figure.tight_layout()
        ax.figure.savefig(save_path, transparent=True)

# %% SUBSET
for main_group in ['immune', 'epithelial', 'endothelial', 'stromal']:
    for (engine, n_neighbors, min_dist, exclude_markers) in params:
        logger.info(
            f'UMAP for {main_group} with {len(data)} data points with n_neighbors {n_neighbors}, min_dist {min_dist} and engine {engine}. Excluding: {exclude_markers}')

        reducer_path = get_reducer_path(save_dir=save_dir,
                                        n_neighbors=n_neighbors, min_dist=min_dist, engine=engine,
                                        exclude_markers=exclude_markers)

        with open(reducer_path, 'rb') as f:
            container = pickle.load(f)

        index = container['index']
        reducer = container['reducer']

        colormap_dict = get_colormap_dict(name='label')

        metadata = ds.metadata.loc[index].copy()
        labels = metadata['label'].values

        metadata = metadata[metadata['main_group'] == main_group]
        uniq_labels = metadata.label.unique()

        # subset data
        grouped = metadata.groupby('label')
        num_samples = 100_000
        min_per_group = int(num_samples / grouped.ngroups)
        min_per_group = min(min_per_group, grouped.size().min())
        metadata = sample_min_per_group_then_uniform(grouped=grouped, n=num_samples, min_per_group=min_per_group,
                                                     random_state=0)

        subset_points = np.zeros(len(index), dtype=bool)
        keep_idc = index.isin(metadata.index)
        subset_points[keep_idc] = True

        save_path = reducer_path.parent / f'label={main_group}.pdf'
        ax = plot_points(reducer=reducer, subset_points=subset_points, labels=labels, color_key=colormap_dict,
                         all_points=True)
        ax.figure.tight_layout()
        ax.figure.savefig(save_path, transparent=True)
        # ax.figure.show()

# %% INTENSITIES
for value in data.columns:
    for (engine, n_neighbors, min_dist, exclude_markers) in params:
        logger.info(
            f'UMAP for {value} with {len(data)} data points with n_neighbors {n_neighbors}, min_dist {min_dist} and engine {engine}. Excluding: {exclude_markers}')

        reducer_path = get_reducer_path(save_dir=save_dir,
                                        n_neighbors=n_neighbors, min_dist=min_dist, engine=engine,
                                        exclude_markers=exclude_markers)

        with open(reducer_path, 'rb') as f:
            container = pickle.load(f)

        index = container['index']
        reducer = container['reducer']

        colormap_dict = get_colormap_dict(name=label)

        metadata = ds.metadata.loc[index].copy()
        metadata['pat_id'] = metadata.index.get_level_values('sample_id').map(sid_to_pid)

        labels = None
        values = data[value].values
        order = np.argsort(values)
        order = None
        cmap = matplotlib.colormaps['Reds']
        norm = Normalize(vmin=values.min(), vmax=values.max())

        # subset data
        grouped = metadata.groupby('label')
        num_samples = 100_000
        min_per_group = int(num_samples / grouped.ngroups)
        min_per_group = min(min_per_group, grouped.size().min())
        metadata = sample_min_per_group_then_uniform(grouped=grouped, n=num_samples, min_per_group=min_per_group,
                                                     random_state=0)

        subset_points = np.zeros(len(index), dtype=bool)
        keep_idc = index.isin(metadata.index)
        subset_points[keep_idc] = True

        save_path = reducer_path.parent / f'value={value}.pdf'
        ax = plot_points(reducer=reducer, subset_points=subset_points,
                         values=norm(values), cmap=cmap, cbar=True,
                         shuffle=False, order=order)
        ax.figure.tight_layout()
        ax.figure.savefig(save_path, transparent=True)
        # ax.figure.show()
