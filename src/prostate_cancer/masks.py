"""Compare instance-segmentation mask TIFFs across pipeline stages.

Object count = number of unique nonzero labels in a mask TIFF (0 = background,
each other integer = one segmented object).
"""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import tifffile
from skimage.color import label2rgb


def n_objects(path: Path) -> int:
    img = tifffile.imread(path)
    return int(np.unique(img[img != 0]).size)


def compare_stages(stage_dirs: dict[str, Path], max_workers: int = 8) -> list[dict]:
    """Count objects per stem across an ordered set of stage directories.

    `stage_dirs` maps stage name -> directory of `<stem>.tiff` files, in pipeline
    order (e.g. `{"deepcell": ..., "filtered": ..., "annotated": ...}`). Returns one
    row per stem seen in any stage, with `n_<stage>` (None if absent) and
    `missing_<stage_a>_to_<stage_b>` for each consecutive stage pair.
    """
    stems_per_stage = {stage: {p.stem for p in d.glob("*.tiff")} for stage, d in stage_dirs.items()}
    all_stems = set.union(*stems_per_stage.values())

    paths = {
        (stage, stem): stage_dirs[stage] / f"{stem}.tiff"
        for stage, stems in stems_per_stage.items()
        for stem in stems
    }
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        counts = dict(zip(paths, pool.map(n_objects, paths.values())))

    stage_names = list(stage_dirs)
    rows = []
    for stem in sorted(all_stems):
        row = {"stem": stem}
        for stage in stage_names:
            row[f"n_{stage}"] = counts.get((stage, stem))
        for a, b in zip(stage_names, stage_names[1:]):
            row[f"missing_{a}_to_{b}"] = (
                row[f"n_{a}"] - row[f"n_{b}"] if row[f"n_{a}"] is not None and row[f"n_{b}"] is not None else None
            )
        rows.append(row)
    return rows


def write_comparison_report(rows: list[dict], stage_names: list[str], title: str, out_path: Path) -> None:
    transitions = list(zip(stage_names, stage_names[1:]))
    common = {(a, b): [r for r in rows if r[f"n_{a}"] is not None and r[f"n_{b}"] is not None] for a, b in transitions}

    lines = [f"# {title}\n", f"Total samples seen (union of all stages): {len(rows)}"]
    all_present = sum(1 for r in rows if all(r[f"n_{s}"] is not None for s in stage_names))
    lines.append(f"Samples present in ALL {len(stage_names)} versions: {all_present}\n")

    lines.append("## Sample-level (whole file) discrepancies\n")
    for a, b in transitions:
        missing = [r for r in rows if r[f"n_{a}"] is not None and r[f"n_{b}"] is None]
        lines.append(f"- Samples in {a} but missing from {b}: {len(missing)}")
        for r in missing:
            lines.append(f"  - `{r['stem']}`")
    lines.append("")

    lines.append("## Object (mask) counts\n")
    for stage in stage_names:
        total = sum(r[f"n_{stage}"] for r in rows if r[f"n_{stage}"] is not None)
        lines.append(f"- Total objects, {stage}: {total}")
    for a, b in transitions:
        removed = sum(r[f"missing_{a}_to_{b}"] for r in common[(a, b)])
        lines.append(f"- Objects removed, {a} -> {b} (common samples only, n={len(common[(a, b)])}): {removed}")
    lines.append("")

    lines.append("## Per-sample table\n")
    header = ["stem"] + [f"n_{s}" for s in stage_names] + [f"missing ({a}->{b})" for a, b in transitions]
    lines.append("| " + " | ".join(header) + " |")
    lines.append("|" + "---|" * len(header))
    for r in rows:
        cells = [r["stem"]] + [r[f"n_{s}"] for s in stage_names] + [r[f"missing_{a}_to_{b}"] for a, b in transitions]
        lines.append("| " + " | ".join("—" if c is None else str(c) for c in cells) + " |")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n")


def plot_removed_objects(before_path: Path, after_path: Path, out_path: Path, before_label: str, after_label: str) -> None:
    """Two-panel figure: full `before` mask (left) vs. only the objects absent from `after` (right)."""
    before = tifffile.imread(before_path)
    after = tifffile.imread(after_path)

    before_labels = set(np.unique(before)) - {0}
    after_labels = set(np.unique(after)) - {0}
    removed_labels = before_labels - after_labels
    assert removed_labels, f"no removed objects for {before_path.stem}"

    removed_mask = np.where(np.isin(before, list(removed_labels)), before, 0)

    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(14, 7))
    ax_left.imshow(label2rgb(before, bg_label=0))
    ax_left.set_title(f"{before_label} (all {len(before_labels)} objects)")
    ax_left.axis("off")

    ax_right.imshow(label2rgb(removed_mask, bg_label=0))
    ax_right.set_title(f"removed ({len(removed_labels)} objects, not in {after_label})")
    ax_right.axis("off")

    fig.suptitle(before_path.stem)
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
