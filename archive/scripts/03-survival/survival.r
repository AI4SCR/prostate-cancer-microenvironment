
library(arrow)
library(tidyverse)
library(survival)
library(ggsurvfit)
library(gtsummary)
library(compositions)
library(coxme)

scores.path = file.path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/publications/PCa/0-export/scores.parquet')
metadata.path = file.path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/publications/PCa/0-export/metadata.parquet')
clinical.path = file.path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/publications/PCa/0-export/clinical.parquet')

scores = read_parquet(scores.path) |> select(-all_of('__index_level_0__'))
metadata = read_parquet(metadata.path)
clinical = read_parquet(clinical.path)
num.patients = clinical$pat_id |> n_distinct()

pat_ids.valid = intersect(scores$pat_id, clinical$pat_id)

agg = 'max'

clinical.names = c('pat_id', 'tma_id', 'age_at_surgery',
                   "gs_grp",  # 4  1  3  2  5 NA
                   'last_fu',
                   'os_status',  # "dead"  "alive"
                   'cause_of_death',  # "PCa_death", "non-PCa_death", "alive" 
                   'psa_progr',  # 1 0 NA
                   'clinical_progr',  # 1 0 NA
                   'disease_progr',  # 1 0 NA
                   'gleason_grp',  # 4  1  3  2  5 NA
                   'inflammation',  # "no"  "yes" NA
                   'cribriform'  # "yes" "no"  NA
)
clinical = clinical |> select(all_of(clinical.names)) |> distinct()

scores = scores[!is.na(scores$is_tumor) & scores$is_tumor == 'yes',]

# NOTE: we need to filter since for the interactions we have the same score_name,
#   but different score_type (e.g. proportion, observed)
score_types = unique(scores$score_type)
score_types = c('proportion')

# score_type = 'proportion'
# score_name = "epithelial-luminal"
# 
# s1 <- survfit(Surv(time=time, event=event) ~ 1, data=data)
# summary(s1)
#
# s1 |>
#   ggsurvfit() +
#   labs(
#     x = "Months",
#     y = "Overall survival probability"
#   ) + 
#   add_confidence_interval() +
#   add_risktable()


cox.fit = function(table, score_name){
  # print(score_name)
  
  cols <- c("time","event", score_name)
  df <- table |> select(all_of(cols))
  
  # cols.mixed <- c("time","event", "pat_id", score_name)
  # df.mixed = table.mixed |> select(all_of(cols.mixed))
  
  # NOTE: we introduce an indicator for NaN values
  #   and set the score to 0 for NaN values
  #   reason: for interactions, NaN is not equivalent to 0
  present = !is.na(df[[score_name]])
  if(all(present)){
    # fml <- as.formula(paste0("Surv(time, event) ~ `", score_name, "` + age_at_surgery"))
    fml <- as.formula(paste0("Surv(time, event) ~ `", score_name, "`"))
    
    # mixed-effect model
    # fml.mixed <- as.formula(paste0("Surv(time, event) ~ (1 | pat_id) + `", score_name, "`"))
  }else{
    df$present = 1
    df[!present, score_name] = 0
    
    fml <- as.formula(paste0("Surv(time, event) ~ `", score_name, "` + present"))
    
    # mixed-effect model, compare performance to baseline model
    # fml <- as.formula(paste0("Surv(time, event) ~ `", score_name, "` + present + (1 | pat_id)"))
  }
  
  fit = coxph(fml, data = df)
  # fit.mixed = coxme(fml.mixed, data = df.mixed)
  
  s = summary(fit)
  # tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  
  hr <- unname(s$coefficients[1, "exp(coef)"])
  hr.lower = s$conf.int[1, "lower .95"]
  hr.upper = s$conf.int[1, "upper .95"]
  p_value <- unname(s$coefficients[1,"Pr(>|z|)"])
  
  # Global tests
  p_lrt   <- unname(s$logtest["pvalue"])
  p_wald  <- unname(s$waldtest["pvalue"])
  p_score <- unname(s$sctest["pvalue"])
  
  out <- tibble(
    score_name = score_name,
    hr = hr,
    hr_lower = hr.lower,
    hr_upper = hr.upper,
    p_value = p_value,
    p_lrt = p_lrt,
    p_wald = p_wald,
    p_score = p_score
    )
  
  if(!all(present)) {
    out$hr_present <- unname(exp(s$coefficients["present", "coef"]))
    out$p_value_present <- unname(s$coefficients["present", "Pr(>|z|)"])
    out$hr_lower_present <- unname(s$conf.int["present", "lower .95"])
    out$hr_upper_present <- unname(s$conf.int["present", "upper .95"])
  }
  
  out
}

event_names = c('os_status', 'cause_of_death', 'clinical_progr', 'disease_progr')
event_dict <- c(
  os_status = "dead",
  cause_of_death = "PCa_death",
  clinical_progr = 1,
  disease_progr = 1
)

get_table_for_mixed = function(data, clinical){
  data = data |> inner_join(clinical, by = 'pat_id')
  data = data |> pivot_wider(names_from = score_name, values_from = score)
  data$pat_id |> is.na() %>% sum()
  data$pat_id |> n_distinct()
  data
}

get_table_for_baseline = function(data, clinical, agg='mean'){
  if(agg == 'mean'){
    data = data |> 
      group_by(pat_id, score_name) |> 
      summarise(score = mean(score)) |> 
      ungroup()
  }else if(agg == 'max'){
    data = data |> 
      group_by(pat_id, score_name) |> 
      summarise(score = max(score)) |> 
      ungroup()
  }
  
  data = data |> inner_join(clinical, by = 'pat_id')
  data = data |> pivot_wider(names_from = score_name, values_from = score)
  data$pat_id |> is.na() %>% sum()
  data$pat_id |> n_distinct()
  
  data
}

set_events = function(data, event_name){
  data$time = data$last_fu
  event = event_dict[[event_name]]
  data$event = data[[event_name]] == event
  data
}

for(event_name in event_names){
  
  # event_name = event_names[1]
  save.dir = file.path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/publications/PCa/outputs/7-survival/', 'without-proportions', event_name)
  dir.create(save.dir, recursive = TRUE, showWarnings = FALSE)
  
  cols.surv = c('pat_id', 'last_fu', event_name) # 'age_at_surgery'
  clinical.surv = clinical |> select(all_of(cols.surv)) |> distinct()
  stopifnot(nrow(clinical.surv) == num.patients)
  
  # score_type = 'interaction_label_classic_observation_radius_32'
  # score_type = "proportion"
  for(score_type in score_types){
    
    data = scores[scores$score_type == score_type, ]
    stopifnot(0 == is.na(data$score) |> sum())
    
    score_names = unique(data$score_name)
    
    data = data |> select(pat_id, tma_id, score_name, score)
    table = get_table_for_baseline(data=data, clinical=clinical.surv, agg=agg)
    # table.mixed = get_table_for_mixed(data=data, clinical=clinical.surv)
    
    # CLR transform of proportions
    if('proportion' %in% score_type){
      table[, 4:ncol(table)] = clr(table[, 4:ncol(table)])
      # table.mixed[, 5:ncol(table.mixed)] = clr(table.mixed[, 5:ncol(table.mixed)])
    }
    
    table = set_events(data=table, event_name=event_name)
    # table.mixed = set_events(data=table.mixed, event_name=event_name)
    
    # score_name = score_names[1]
    # score_name = "endothelial-blood-vessel(ERG+)->endothelial-blood-vessel(ERG-)"
    results = map_dfr(score_names, ~cox.fit(table=table, score_name=.x)) |>
      arrange(p_value)
    results$score_type = score_type
    
    results = results |>
      mutate(p_adj = p.adjust(p_value, method = "fdr"))
    
    if("p_value_present" %in% colnames(results)){
      results = results |>
        mutate(p_adj_present = p.adjust(p_value_present, method = "fdr"))
    }
    
    save.path = file.path(save.dir, paste0(score_type,'-', agg, '.csv'))
    write_csv(results, save.path)
    
  }
}



