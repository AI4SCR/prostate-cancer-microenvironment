
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
clinical.path = file.path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/publications/PCa/0-export/survival-cell-freq-groups.parquet')
clinical = read_parquet(clinical.path)
num.patients = clinical$pat_id |> n_distinct()

fit.os <- survfit(Surv(last_fu, os_event) ~ cluster_id, data = clinical)

g = ggsurvplot(
  fit.os,
  data = clinical,
  conf.int = TRUE,
  risk.table = TRUE,
  pval = TRUE,
  pval.method = TRUE,
  xlab = "Time",
  ylab = "Overall survival probability",
  legend.title = "Cell Frequency Groups"
)
g
ggsave(filename = file.path(save_dir, 'survival-cell-freq-groups-os.pdf'), 
       plot = g$plot, device = 'pdf', 
       width = 10, height = 4)


fit.progr <- survfit(Surv(disease_progr_time, disease_progr) ~ cluster_id, data = clinical)

g = ggsurvplot(
  fit.progr,
  data = clinical,
  conf.int = TRUE,
  risk.table = TRUE,
  pval = TRUE,
  pval.method = TRUE,
  xlab = "Time",
  ylab = "Disease progression probability",
  legend.title = "Cell Frequency Groups"
)
g
ggsave(filename = file.path(save_dir, 'survival-cell-freq-groups-progr.pdf'), 
       plot = g$plot, device = 'pdf', 
       width = 10, height = 4)

