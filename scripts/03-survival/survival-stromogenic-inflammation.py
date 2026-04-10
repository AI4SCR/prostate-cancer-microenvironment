from ai4bmr_datasets import PCa
import pickle
import pandas as pd
from tqdm import tqdm
from pathlib import Path

ds = PCa(
    base_dir=Path('/work/FAC/FBM/DBC/mrapsoma/prometex/data/datasets/PCa'),
    image_version='filtered', mask_version='annotated',
    load_metadata=True, load_intensity=False, load_spatial=False)
ds.setup()

metadata = ds.metadata.copy()
clinical = ds.clinical.copy()
clinical.groupby('pat_id')
clinical['has_inflamed'] = clinical['inflammation'].map({'yes': 1, 'no': 0})
clinical['has_stromogenic'] = clinical['stromogenic_smc_loss_reactive_stroma_present'].map({'yes': 1, 'no': 0})
vars = [
    'os_status',
    'cause_of_death',
    'last_fu',
    'psa_progr', 'psa_progr_time',
    'clinical_progr', 'clinical_progr_time',
    'disease_progr', 'disease_progr_time',
    'has_inflamed', 'has_stromogenic'
]
clinical = clinical.groupby('pat_id')[vars].max()

save_path = Path('/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/0-paper/0-export/survival-stromogenic-inflammation.parquet')
clinical.to_parquet(save_path, engine='fastparquet')

