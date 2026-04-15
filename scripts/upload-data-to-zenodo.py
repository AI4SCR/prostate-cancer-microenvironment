import os

from dotenv import load_dotenv
assert load_dotenv()
import requests
from pathlib import Path

# test token
ACCESS_TOKEN = os.environ['ZENODO_API_KEY']
headers = {'Authorization': f'Bearer {ACCESS_TOKEN}'}
r = requests.get('https://zenodo.org/api/deposit/depositions',
                  headers=headers)
r.raise_for_status()
r.json()

# get bucket url
DEPOSITION_ID = 19552665
r = requests.get(
    f"https://zenodo.org/api/deposit/depositions/{DEPOSITION_ID}",
    headers=headers,
)
r.raise_for_status()

data = r.json()
bucket_url = data["links"]["bucket"]

print("Bucket URL:", bucket_url)

# upload
BASE_DIR = os.environ['BASE_DIR']

fname = 'PCa.tar.gz'
path = Path(BASE_DIR).parent / fname
assert path.exists()
assert path.is_file()

# path = '/users/amarti51/projects/prostate-cancer-microenvironment/test.txt'
# fname = 'test.txt'

headers = {'Authorization': f'Bearer {ACCESS_TOKEN}'}
url = "%s/%s" % (bucket_url, fname)
print("Uploading to", url)

with open(path, "rb") as fp:
    r = requests.put(
        url,
        data=fp,
        headers=headers,
    )

r.raise_for_status()
r.json()

