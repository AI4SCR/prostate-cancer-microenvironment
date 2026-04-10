def legend_from_dict(label_to_color: dict):
    """Create legend elements from a dictionary mapping labels to colors."""
    from matplotlib.patches import Patch

    return [
        Patch(facecolor=color, label=label)
        for label, color in label_to_color.items()
    ]

