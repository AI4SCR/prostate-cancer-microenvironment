import numpy as np
from ai4bmr_datasets import PCa
import numpy as np
import pickle
import pandas as pd
from tqdm import tqdm
from pathlib import Path

ds = PCa(
    base_dir=Path('/work/FAC/FBM/DBC/mrapsoma/prometex/data/datasets/PCa'),
    image_version='filtered', mask_version='annotated',
    load_metadata=True, load_intensity=False, load_spatial=False)
ds.setup()

clinical = ds.clinical.copy()
vars = [
    'os_status',
    'gleason_grp',
    'gs_grp',
    'last_fu',
    'psa_progr', 'psa_progr_time',
    'clinical_progr', 'clinical_progr_time',
    'disease_progr', 'disease_progr_time',
]

clinical = clinical.groupby('pat_id')[vars].max()
clinical['os_event'] = clinical.os_status.map({'alive': 0, 'dead': 1})
clinical = clinical.drop(columns='os_status')
clinical['gs_grp'] = clinical.gs_grp.map(lambda x: {'nan': np.nan}.get(x, x))
clinical = clinical.astype(float).convert_dtypes()

save_path = Path('/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/0-paper/0-export/survival-gleason.parquet')
clinical.to_parquet(save_path, engine='fastparquet')

