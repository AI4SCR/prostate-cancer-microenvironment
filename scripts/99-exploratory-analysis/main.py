# %%
import os

from ai4bmr_datasets import PCa
from pathlib import Path
from dotenv import load_dotenv
import os
assert load_dotenv(override=True)

# %%
# check documentation: https://github.com/AI4SCR/ai4bmr-datasets
ds = PCa(base_dir=Path(os.environ['BASE_DIR']), image_version='filtered', mask_version='annotated',
         load_metadata=True, align=True)
ds.setup(engine='pyarrow')

sample_id = ds.sample_ids[0]
len(ds.sample_ids)

ds.images[sample_id]  # lazy image obj
ds.images[sample_id].data  # actual image
ds.panel  # panel

# %% Normalize an image

