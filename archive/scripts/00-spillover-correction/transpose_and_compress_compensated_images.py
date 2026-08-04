from pathlib import Path

import numpy as np
from loguru import logger as base_logger
from jsonargparse import CLI
from skimage.io import imread, imsave


def main(data_dir: Path, verbose: int = 0):
    logger = base_logger.bind(task="transpose-compress", verbose=verbose)
    for img_path in Path(data_dir).rglob("*.tif"):
        logger.info(f"processing {img_path}")
        i = imread(img_path)
        imsave(
            img_path.with_suffix(".tiff"),
            np.transpose(i),
            plugin="tifffile",
            compression="deflate",
        )
        img_path.unlink()


if __name__ == "__main__":
    CLI(main, as_positional=False)
