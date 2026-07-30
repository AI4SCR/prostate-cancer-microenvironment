
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

clinical.path = file.path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/publications/PCa/0-export/survival-stromogenic-inflammation.parquet')
clinical = read_parquet(clinical.path)
num.patients = clinical$pat_id |> n_distinct()

clinical <- clinical |>
  mutate(event = os_status == "dead")

# clinical <- clinical |>
#   mutate(event = cause_of_death == "dead")

inflammed = clinical[!is.na(clinical$has_inflamed),]
stromogenic = clinical[!is.na(clinical$has_stromogenic),]


time_col = 'last_fu'
status_col = 'event'
group_col = 'has_inflamed'

fit <- survfit(Surv(last_fu, event) ~ has_inflamed, data = inflammed)

g = ggsurvplot(
  fit,
  data = inflammed,
  conf.int = TRUE,
  risk.table = TRUE,
  pval = TRUE,
  pval.method = TRUE,
  xlab = "Time",
  ylab = "Overall survival probability",
  legend.title = "Inflamed"
)
g
ggsave(filename = file.path(save_dir, 'survival-inflammation.pdf'), plot = g$plot, device = 'pdf', 
       width = 10, height = 4)


fit <- survfit(Surv(last_fu, event) ~ has_stromogenic, data = stromogenic)
g = ggsurvplot(
  fit,
  data = stromogenic,
  conf.int = TRUE,
  risk.table = TRUE,
  pval = TRUE,
  pval.method = TRUE,
  xlab = "Time",
  ylab = "Overall survival probability",
  legend.title = "Stromogenic"
)
g
ggsave(filename = file.path(save_dir, 'survival-stromogenic.pdf'), plot = g$plot, device = 'pdf', 
       width = 10, height = 4)
