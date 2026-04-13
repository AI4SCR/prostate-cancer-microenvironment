#!/usr/bin/env bash

rsync -avhP \
  --exclude="*.zip" \
  --exclude="*.mcd" \
  --exclude="*.tar" \
  --exclude="*.tar.gz" \
  --exclude="images/raw/" \
  unil:/work/FAC/FBM/DBC/mrapsoma/prometex/data/datasets/PCa/02_processed \
  $HOME/projects/prostate-cancer-microenvironment/data/PCa