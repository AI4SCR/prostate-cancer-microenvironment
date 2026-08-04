# %%
"""Port a legacy `reducer.pkl` (umap-learn `UMAP` object) to a portable parquet.

The original pipeline's `reducer.pkl` files were pickled by a umap-learn/
numba/pynndescent version stack from ~September 2025. Unpickling them with
this repo's own environment (or any modern numba/umap-learn) fails with
`ModuleNotFoundError: No module named 'numba.core.types.old_scalars'` --
that module was removed upstream in later numba releases. This script must
be run with `pixi run` inside `scripts/port/pixi.toml`'s pinned environment
(numba==0.61, numpy==2.0, pandas>=3, umap-learn>=0.5.5, pyarrow), NOT this
repo's default `.venv`.

Even in that older-but-still-incompatible env, unpickling needs one more
fix: a `pynndescent.NNDescent.__setstate__` compatibility patch, since the
pickled `NNDescent` state predates the `_min_distance`/`quantization`
attributes newer pynndescent expects. Patched here before loading.

`UMAP.fit()` is never seeded anywhere in the legacy pipeline, so these
embeddings are NOT reproducible by re-running the fit -- each fit gives a
geometrically different (if topologically similar) embedding. The ported
parquet files ARE the ground truth for the published figure; this script
refuses to overwrite an existing output file unless `--overwrite` is passed
explicitly.

Usage (from this env's own pixi project, not the repo's .venv):
    pixi run python scripts/port/port_umap_reducer.py <pkl_path> <out_path>
"""
import pickle
from pathlib import Path

import pandas as pd
from jsonargparse import CLI


def _patch_pynndescent_setstate() -> None:
    from pynndescent import pynndescent_

    original_setstate = pynndescent_.NNDescent.__setstate__

    def patched_setstate(self, state):
        if "_min_distance" not in state:
            state["_min_distance"] = 0.0
        if "quantization" not in state:
            state["quantization"] = None
        original_setstate(self, state)

    pynndescent_.NNDescent.__setstate__ = patched_setstate


def main(pkl_path: Path, out_path: Path, overwrite: bool = False):
    assert pkl_path.exists(), f"{pkl_path} does not exist"
    assert overwrite or not out_path.exists(), (
        f"{out_path} already exists -- this may be the ground-truth ported embedding "
        "for a figure with no seeded UMAP fit to regenerate it from. Pass --overwrite "
        "explicitly if you are certain you want to replace it."
    )

    _patch_pynndescent_setstate()

    with open(pkl_path, "rb") as f:
        container = pickle.load(f)
    reducer = container["reducer"] if isinstance(container, dict) else container
    embedding = reducer.embedding_

    index = container.get("index") if isinstance(container, dict) else None
    df = pd.DataFrame(embedding, columns=["umap_1", "umap_2"])
    if isinstance(index, pd.MultiIndex):
        for i, name in enumerate(index.names):
            df.insert(i, name, index.get_level_values(i).to_numpy())
    elif index is not None:
        df.insert(0, index.name or "index", pd.Index(index).to_numpy())

    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_path, index=False)
    print(f"saved {len(df)} rows to {out_path}")


if __name__ == "__main__":
    CLI(main)
