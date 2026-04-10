
library(arrow)
library(tidyverse)
library(survival)
library(ggsurvfit)
library(gtsummary)
library(compositions)
library(coxme)
library(survival)
library(survminer)
library(dplyr)

save_dir = file.path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/publications/PCa/outputs/7-survival')
clinical.path = file.path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/publications/PCa/0-export/survival-gleason.parquet')
clinical = read_parquet(clinical.path)
num.patients = clinical$pat_id |> n_distinct()

gs_grp = clinical[!is.na(clinical$gs_grp),]
gleason_grp = clinical[!is.na(clinical$gleason_grp),]

fit.gs_grp <- survfit(Surv(last_fu, os_event) ~ gs_grp, data = gs_grp)

g = ggsurvplot(
  fit.gs_grp,
  data = gs_grp,
  conf.int = TRUE,
  risk.table = TRUE,
  pval = TRUE,
  pval.method = TRUE,
  xlab = "Time",
  ylab = "Overall survival probability",
  legend.title = "Patient-level Gleason score"
)
g
ggsave(filename = file.path(save_dir, 'survival-gleason-patient.pdf'), 
       plot = g$plot, device = 'pdf', 
       width = 10, height = 4)


fit.gleason_grp <- survfit(Surv(last_fu, os_event) ~ gleason_grp, data = gleason_grp)

g = ggsurvplot(
  fit.gleason_grp,
  data = gleason_grp,
  conf.int = TRUE,
  risk.table = TRUE,
  pval = TRUE,
  pval.method = TRUE,
  xlab = "Time",
  ylab = "Overall survival probability",
  legend.title = "Max-pooled TMA-level Gleason score"
)
g
ggsave(filename = file.path(save_dir, 'survival-gleason-tma.pdf'), plot = g$plot, device = 'pdf', 
       width = 10, height = 4)

