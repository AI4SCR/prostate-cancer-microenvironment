
library(arrow)
library(tidyverse)
library(survival)
library(broom)
library(dplyr)
library(ggplot2)
library(ggsurvfit)
library(survminer)
library(gtsummary)
library(compositions)

scores.path = file.path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/publications/PCa/0-export/scores-v2.parquet')
clinical.path = file.path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/publications/PCa/0-export/clinical.parquet')

scores = read_parquet(scores.path)
clinical = read_parquet(clinical.path)
num.patients = clinical$pat_id |> n_distinct()

# checks
table(clinical$os_status, clinical$disease_progr)
has_inflammation = !is.na(clinical$inflammation)
clinical[has_inflammation,]$pat_id |> n_distinct() # 190

event.name = 'disease_progr'
event.name = 'os_status'

if(event.name == 'os_status'){
  event.value = 'dead'
  time.name = 'last_fu'
} else if(event.name == 'disease_progr'){
  event.value = 1
  time.name = 'disease_progr_time'
}

# score_type = 'proportion_tma'
# aggregation = 'max'
# tumors_only = TRUE

# score_type = 'proportion_pat'
# aggregation = ''
# tumors_only = FALSE

# score_type = 'proportion_pat_tumor'
# aggregation = ''
# tumors_only = FALSE

save.dir = file.path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/publications/PCa/outputs/7-survival/', score_type, event.name)
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

# 515 TMA_IDs before filtering for tumor cores only
if (tumors_only){
  data = data[!is.na(data$is_tumor) & data$is_tumor == 'yes',]
}
stopifnot(all(!is.na(data$score)))
data.tma = data |> 
  select(tma_id, score, score_name) |>
  pivot_wider(names_from = score_name, values_from = score)
data.tma = data.tma[order(data.tma$tma_id), order(colnames(data.tma))]
head(data.tma)

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
}else{
  print('no aggregation')
}

data = data |> 
  select(pat_id, score, score_name) |>
  pivot_wider(names_from = score_name, values_from = score)

# CLR transform of proportions
data[, 2:ncol(data)] = clr(data[, 2:ncol(data)])

data = data |> inner_join(clinical.surv, by = 'pat_id')

data$pat_id |> is.na() %>% sum()
data$pat_id |> n_distinct()
data$event |> sum()

stopifnot(!(is.na(data$event) %>% any()))
stopifnot(!(is.na(data$time) %>% any()))

# score_name = "epithelial-luminal(ERG+p53+)"

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

save.path = file.path(save.dir, paste0('aggregation=', aggregation, '-tumors_only=', tumors_only, '.csv'))
write_csv(results, save.path)

pdat = results |> select(score_name, score_type, hr, hr_lower, hr_upper, p_value, p_adj)
# order = order(pdat$p_value)
# pdat = pdat[order, ]

pdat <- pdat %>%
  arrange(p_value) %>%
  mutate(
    score_name = factor(score_name, levels = score_name),
    signf = p_adj < 0.05
  )

p <- ggplot(pdat, aes(y = score_name)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  geom_errorbarh(
    aes(xmin = hr_lower, xmax = hr_upper),
    height = 0.2, linewidth = 0.6, color = "black", alpha = 0.8
  ) +
  geom_point(
    aes(x = hr, fill = signf),
    size = 2.5, shape = 21,
  ) +
  labs(
    x = "Hazard ratio",
    y = NULL,
    title = paste("Cox -", event.name, "~", score_type),
    subtitle = "Points = HR; bars = 95% CI"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(size = 13, face = "bold")
  )
p
ggsave(
  filename = sub('.csv', '.pdf', save.path),
  plot = p,
  width = 8,
  height = 6
  # dpi = 300
)


# Kaplan-Meier Curves
# proportions_tma
# for os_status q = 0.90 | p-value = 0.01085
# disease_progr nothing is significant

# proportions_pat
# for os_status q = 0.85 | p-value = 0.007453
# disease_progr q = 0.20 | p-value = 0.03538

# proportion_pat_tumor
# for os_status q = 0.85 | p-value = 0.01027
# disease_progr q = 0.20 | p-value = 0.03538

qs = seq(.1, 0.95, length.out = 18)
for(q in qs){
  
  thres = quantile(data$`stromal-CAF1(CD105+)`, q)
  data['has_risk'] = data$`stromal-CAF1(CD105+)` > thres
  fit <- survfit(Surv(time, event) ~ has_risk, data = data)
  
  p.val = surv_pvalue(fit, data=data)
  cat(sprintf("q = %.2f | p-value = %.4g\n", q, p.val$pval))
  
  g = ggsurvplot(
    fit,
    data = data,
    conf.int = TRUE,
    risk.table = TRUE,
    pval = TRUE,
    pval.method = TRUE,
    xlab = "Time",
    ylab = "Overall survival probability",
    legend.title = "has_risk"
  )
  g
}
