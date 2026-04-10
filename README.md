# PCa

## Set-up

### Python
Install the regular Python dependencies from [`pyproject.toml`](/Users/adrianomartinelli/projects/prostate-cancer-microenvironment/pyproject.toml).

Some project dependencies are installed manually and are not managed through `pyproject.toml`:

```bash
uv pip install -e ~/projects/ai4bmr-datasets
uv pip install -e ~/projects/ai4bmr-learn
```

Install ATHENA manually from its repository as well:

- [ATHENA](https://github.com/AI4SCR/ATHENA)

The editable installs are required because this codebase imports `ai4bmr_datasets` and `ai4bmr_learn`, and some scripts also depend on `athena`.

### Git Configuration
Add task specific ignores in the corresponding `.gitignore` of the relevant sub-folders.
I recommend to configure your global ignores with this [gist](https://gist.github.com/adrianomartinelli/7471ce6be2b6dd6a93ab9838fa201b3e).

### Branch Management
Create branches for the specific task prefixed with a personal identifier (`<IDENTIFIER>/<FEATURE>`) like `art/feature-1`. 
`main` branch is protected, you need to open a PR to merge with `main`.

### Environment
Your scripts should access the data in the machine specific folder. Use `.env` in this folder with the following content:
To make the environment variable available to your scripts, you need to load the `.env` file.

```R
# with R
library(dotenv)
PATH_BOX = Sys.getenv("PATH_BOX")
```

```python
# with python
from dotenv import load_dotenv
load_dotenv()

import os
PATH_BOX = os.getenv('PATH_BOX')
```

## Steinbock

### Setup
Install [docker dektop](https://docs.docker.com/desktop/install/mac-install/) and make sure that the daemon is running.
Add the following function to your `.zshrc` or `.bashrc` to make the steinbock docker image available as a command.

```bash
# in your terminal
vim .zshrc
```

```
steinbock0.16.1(){
docker run \
-v "$(pwd)":/data \
ghcr.io/bodenmillergroup/steinbock:0.16.1 \
"$@"
}
```

## Set-up
### R
Do not install R with brew on Mac!

Add CRAN Repository and Key
https://cran.r-project.org/bin/linux/ubuntu/fullREADME.html

```R
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")

Available CPU Processors
echo "CPU threads: $(grep -c processor /proc/cpuinfo)"
update.packages(lib.loc="/usr/local/lib/R/site-library", ask = FALSE,
  checkBuilt = TRUE, Ncpus = 28)

install.packages(
'rlang', 'ggrastr', 'ragg')
)

Warning messages:
1: In install.packages(...) :
  installation of package ‘textshaping’ had non-zero exit status
2: In install.packages(...) :
  installation of package ‘ragg’ had non-zero exit status
3: In install.packages(...) :
  installation of package ‘ggrastr’ had non-zero exit status
4: In install.packages(...) :
  installation of package ‘scater’ had non-zero exit status
> client_loop: send disconnect: Broken pipe
```

# R studio server
enbale incoming conncetions
```bash
sudo ufw allow 8787
```

### Commands

```bash
mkdir -p steinbock/raw && cd steinbock
find <PATH_TO_RAW_DATA> -name "*.mcd" -exec cp raw/{} \; # copy .mcd files
ln -s <PATH_TO_PROJECT>/03_spatial/03_spatial/process-steinbock-panel.py .
steinbock0.16.1 preprocess imc panel
cut -d "," -f 2 panel.csv  # show available markers

conda activate PCa
python process-steinbock-panel.py panel.csv

steinbock0.16.1 preprocess imc images --hpf 50 --imgout img

mkdir img_store && mv img/* img_store
find img_store -regex ".*220811_TMA 17_15_III_A1_Big_slide 2_00[1-9].tiff" -exec mv {} img/ \;
find img_store -regex ".*PCa\ TMA\ 15_02_IVB\ ROI\ 1\ to\ 8_00[1-9].tiff" -exec mv {} img/ \;

steinbock0.16.1 segment deepcell --minmax --type nuclear && \
steinbock0.16.1 measure intensities && \
steinbock0.16.1 measure regionprops && \
steinbock0.16.1 measure neighbors --type centroids --kmax 5
```

# Troubleshooting

- When using symlinks make sure that the location to which the symlink points is available in the docker image, i.e. mounted
