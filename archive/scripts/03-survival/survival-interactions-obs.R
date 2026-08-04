
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
scores$score_type |> unique()

# event.name = 'disease_progr'
event.name = 'os_status'

if(event.name == 'os_status'){
  event.value = 'dead'
  time.name = 'last_fu'
} else if(event.name == 'disease_progr'){
  event.value = 1
  time.name = 'disease_progr_time'
}

tumors_only = TRUE
# tumors_only = FALSE

aggregation = 'max'
# aggregation = 'mean'

# graph_key = 'radius_12'
# graph_key = 'radius_16'
# graph_key = 'radius_32'
# graph_key = 'radius_48'
graph_key = 'radius_64'

# score_type.parent = 'interaction_meta_label_proportion_observation'
score_type.parent = 'interaction_label_proportion_observation'
clr_transform = TRUE
log_transform = FALSE

# score_type.parent = 'interaction_meta_label_classic_observation'
# score_type.parent = 'interaction_label_classic_observation'
# log_transform = TRUE
# clr_transform = FALSE

stopifnot(clr_transform + log_transform <= 1)

score_type = paste0(score_type.parent, '_', graph_key)

data = scores[scores$score_type == score_type, ]
score_names = unique(data$score_name)
stopifnot(nrow(data) > 0)

save.dir = file.path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/publications/PCa/outputs/7-survival/', score_type.parent, graph_key, event.name)
dir.create(save.dir, recursive = TRUE, showWarnings = FALSE)

pat_ids.valid = intersect(scores$pat_id, clinical$pat_id)
stopifnot(length(pat_ids.valid) == 195)

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
if(tumors_only){
  stopifnot(data$pat_id |> unique() |> length() == 190)
}else{
  stopifnot(data$pat_id |> unique() |> length() == 195)
}

data = data |> 
  pivot_wider(names_from = score_name, values_from = score)

data[is.na(data)] = 0
# rowSums(data[, 6:ncol(data)])

# CLR transform of proportions
if(clr_transform){
  data[, 6:ncol(data)] = clr(data[, 6:ncol(data)])
}

if(log_transform){
  mask = data[, 6:ncol(data)] == 0
  thres = min(data[, 6:ncol(data)][!mask]) / 2
  data[, 6:ncol(data)][mask] = thres
  data[, 6:ncol(data)] = log1p(data[, 6:ncol(data)])
  # hist(unlist(data[, 6:ncol(data)]), breaks=100)
}

stopifnot(data$pat_id |> is.na() %>% sum() == 0)
data$pat_id |> n_distinct()
data$event |> sum()

score_name = "endothelial-blood-vessel(ERG+)->epithelial-basal"

cox.fit = function(data, score_name){
  cols <- c("time","event", score_name)
  
  df <- data |> select(all_of(cols))
  
  fml <- as.formula(paste0("Surv(time, event) ~ `", score_name, "`"))
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

save.path = file.path(save.dir, paste0(score_type,'-aggregation=', aggregation, '-tumors_only=', tumors_only, '-event=', event.name,'.csv'))
write_csv(results, save.path)