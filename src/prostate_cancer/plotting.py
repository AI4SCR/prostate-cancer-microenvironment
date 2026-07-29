def legend_from_dict(label_to_color: dict):
    """Create legend elements from a dictionary mapping labels to colors."""
    from matplotlib.patches import Patch

    return [
        Patch(facecolor=color, label=label)
        for label, color in label_to_color.items()
    ]


def plot_embedding_categorical(cells, col: str, save_path, title: str, x="umap_1", y="umap_2", point_size: float = 2.0):
    """Scatter an embedding (e.g. UMAP) colored by a categorical column, with a legend when tractable."""
    import matplotlib.pyplot as plt

    from prostate_cancer.utils import create_color_maps

    color_map = create_color_maps(cells[[col]].astype({col: "category"}))[col]

    fig, ax = plt.subplots(figsize=(9, 8))
    ax.scatter(cells[x], cells[y], c=[color_map[v] for v in cells[col]], s=point_size, alpha=0.5, linewidths=0)
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])

    # too many categories (e.g. 195 patients) to legend legibly -- skip it
    if len(color_map) <= 40:
        handles = legend_from_dict(color_map)
        ax.legend(handles=handles, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=6, ncols=2 if len(handles) > 15 else 1)

    fig.tight_layout()
    fig.savefig(save_path, dpi=200)
    plt.close(fig)


def plot_embedding_marker(cells, marker: str, save_path, title: str, x="umap_1", y="umap_2", point_size: float = 2.0):
    """Scatter an embedding (e.g. UMAP) colored by a continuous marker intensity, with a colorbar."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 7))
    sca = ax.scatter(cells[x], cells[y], c=cells[marker], s=point_size, alpha=0.6, cmap="viridis", linewidths=0)
    fig.colorbar(sca, ax=ax, label=marker)
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(save_path, dpi=200)
    plt.close(fig)

