# %%
from pathlib import Path

import pandas as pd
from jsonargparse import CLI


# %%


def main(data_dir: Path):
    # data_dir = Path(
    #     '/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/datasets/PCa/')
    data_dir = data_dir.expanduser()
    path_annotations = Path(data_dir / "clustering/annotations.parquet")

    # %%
    annotations = pd.read_parquet(path_annotations, engine="pyarrow")
    levels_to_drop = set(annotations.index.names) - {"sample_name", "object_id"}
    index = annotations.index.droplevel(level=list(levels_to_drop))
    annotations.index = index

    # %%
    # roi_240223_012 = annotations.loc[('240223_012',), :]
    # roi_240223_012.group_name.value_counts()
    # annotations = annotations.drop(labels='240223_012', level=0)
    assert "240223_012" not in set(annotations.index.get_level_values("sample_name"))

    # filter out specific `group_name` based on assessment of Francesco
    filter_ = annotations.group_name == "epithelial-Synaptophysin_pos-p53_pos"
    assert filter_.sum() == 56
    annotations = annotations[~filter_]

    # %%
    path_labels = Path(data_dir / "metadata/01_raw/label-names.xlsx")
    labels = pd.read_excel(path_labels, sheet_name="label-names", skiprows=0)

    cols = ["name_in_annotations", "main_group", "label"]
    filter_ = labels.columns.str.contains("meta_group")
    meta_group_cols = labels.columns[filter_]
    cols = cols + meta_group_cols.tolist()
    labels = labels[cols]

    assert labels.isna().any().any() == False
    labels = labels.set_index("name_in_annotations")

    assert set(labels.index) - set(annotations.group_name) == {
        "epithelial-Synaptophysin_pos-p53_pos"
    }

    # %% reassign parent group names according to Francescos labels
    main_group = labels["main_group"].to_dict()
    annotations = annotations.assign(main_group=annotations.group_name.map(main_group))

    for meta_group_col in meta_group_cols:
        meta_group = labels[meta_group_col].to_dict()
        annotations[meta_group_col] = annotations.group_name.map(meta_group)

    label = labels["label"].to_dict()
    annotations = annotations.assign(label=annotations.group_name.map(label))

    assert annotations.isna().any().any() == False

    # %% create ids
    main_group_id_dict = sorted(annotations.main_group.unique())
    main_group_id_dict = {k: v for v, k in enumerate(main_group_id_dict)}
    annotations["main_group_id"] = annotations.main_group.map(main_group_id_dict)

    for meta_group_col in meta_group_cols:
        meta_group_id_dict = sorted(annotations[meta_group_col].unique())
        meta_group_id_dict = {k: v for v, k in enumerate(meta_group_id_dict)}
        annotations[f"{meta_group_col}_id"] = annotations[meta_group_col].map(
            meta_group_id_dict
        )

    label_id_dict = sorted(annotations.label.unique())
    label_id_dict = {k: v for v, k in enumerate(label_id_dict)}
    annotations["label_id"] = annotations.label.map(label_id_dict)

    assert annotations.isna().any().any() == False

    # %%
    id_cols = [col for col in annotations.columns if col.endswith("_id")]
    for id_col in id_cols:
        annotations.loc[:, id_col] = annotations[id_col].astype(int)

    cols = ["main_group", "label"]
    cell_counts = annotations[cols].groupby(cols).value_counts()
    cell_counts.to_csv(data_dir / "metadata" / "02_processed" / "labels-counts.csv")

    # %%
    annotations.to_parquet(data_dir / "metadata" / "02_processed" / "labels.parquet")


if __name__ == "__main__":
    CLI(main, as_positional=False)
