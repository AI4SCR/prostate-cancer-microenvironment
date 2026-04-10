
library(arrow)
library(tidyverse)
library(survival)
library(broom)
library(dplyr)
library(ggplot2)
library(ggsurvfit)
library(gtsummary)
library(compositions)

scores.path = file.path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/publications/PCa/0-export/scores-v2.parquet')
clinical.path = file.path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/publications/PCa/0-export/clinical.parquet')

scores = read_parquet(scores.path)
clinical = read_parquet(clinical.path)
num.patients = clinical$pat_id |> n_distinct()

event.name = 'os_status'
event.value = 'dead'
time.name = 'last_fu'
tumors_only = TRUE
aggregation = 'max'
graph_key = 'radius_32'

score_type.parent = 'interaction_label_classic_diff'
clr_transform = FALSE

# score_type.parent = 'interaction_meta_label_classic_diff'
# clr_transform = FALSE
# 
# score_type.parent = 'interaction_label_proportion_diff'
# clr_transform = TRUE

# score_type.parent = 'interaction_meta_label_proportion_diff'
# clr_transform = TRUE

score_type = paste0(score_type.parent, '_', graph_key)

save.dir = file.path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/publications/PCa/outputs/7-survival/', score_type.parent, graph_key, event.name)
dir.create(save.dir, recursive = TRUE, showWarnings = FALSE)

data = scores[scores$score_type == score_type, ]
score_names = unique(data$score_name)

pat_ids.valid = intersect(scores$pat_id, clinical$pat_id)

clinical.names = c('sample_id', 'pat_id', 'tma_id',
                   'last_fu',
                   'os_status',  # "dead"  "alive"
                   'disease_progr',  # 1 0 NA
                   'disease_progr_time'
)
clinical = clinical |> select(all_of(clinical.names))

cols.surv = c('pat_id', event.name, time.name)
clinical.surv = clinical |> select(all_of(cols.surv)) |> distinct()
stopifnot(nrow(clinical.surv) == num.patients)

clinical.surv['event'] = clinical.surv[[event.name]] == event.value
clinical.surv['time'] = clinical.surv[[time.name]]

if (tumors_only){
  data = data[!is.na(data$is_tumor) & data$is_tumor == 'yes',]
}
stopifnot(all(!is.na(data$score)))

# pooling
if(aggregation == 'max'){
  data = data |> 
    group_by(pat_id, score_name) |> 
    summarise(score = max(score)) |> 
    ungroup()
}else if(aggregation == 'mean'){
  data = data |> 
    group_by(pat_id, score_name) |> 
    summarise(score = mean(score)) |> 
    ungroup()
}

data = data |> inner_join(clinical.surv, by = 'pat_id')
data = data |> pivot_wider(names_from = score_name, values_from = score)

# CLR transform of proportions
if(clr_transform){
  data[, 6:ncol(data)] = clr(data[, 6:ncol(data)])
}

data$pat_id |> is.na() %>% sum()
data$pat_id |> n_distinct()
data$event |> sum()

score_name = "endothelial-blood-vessel(ERG+)->epithelial-basal"

cox.fit = function(data, score_name){
  cols <- c("time","event", score_name)
  
  df <- data |> select(all_of(cols))
  df['present'] = !is.na(df[[score_name]])
  
  fml <- as.formula(paste0("Surv(time, event) ~ `", score_name, "` + present"))
  fit = coxph(fml, data = df)
  s = summary(fit)
  coefs <- s$coefficients
  
  hr <- unname(s$coefficients[1, "exp(coef)"])
  hr.lower = s$conf.int[1, "lower .95"]
  hr.upper = s$conf.int[1, "upper .95"]
  p_value <- unname(s$coefficients[1,"Pr(>|z|)"])
  
  out <- tibble(
    score_name = score_name,
    hr = hr,
    hr_lower = hr.lower,
    hr_upper = hr.upper,
    p_value = p_value,
  )
}

results = map_dfr(score_names, ~cox.fit(data=data, score_name=.x)) |>
  arrange(p_value)
results$score_type = score_type

results = results |>
  mutate(p_adj = p.adjust(p_value, method = "fdr"))

save.path = file.path(save.dir, paste0(score_type,'-aggregation=', aggregation, '-tumors_only=', tumors_only, '.csv'))
write_csv(results, save.path)