library(arrow)
library(tidyverse)
library(survival)
library(broom)
library(dplyr)
library(ggplot2)
library(ggsurvfit)
library(gtsummary)

surv.dir = file.path('/Users/adrianomartinelli/Library/CloudStorage/OneDrive-ETHZurich/oneDrive-documents/data/publications/PCa/outputs/7-survival')
model_name = 'without-proportions'
# event_name = 'os_status'
# event_name = 'disease_progr'
event_name = 'clinical_progr'
score_type = 'interaction_label_classic_diff_radius_32'
score_type = 'proportion'
agg = 'mean'

event_names = c('os_status', 'cause_of_death', 'clinical_progr', 'disease_progr')
event_name = event_names[1]
for(event_name in event_names){
  fname = paste0(score_type, '-', agg, '.csv')
  scores.path = file.path(surv.dir, model_name, event_name, fname)
  
  dat = read_csv(scores.path)
  
  if(model_name == 'with-proportions'){
    
  select_ = dat$p_value_target > 0.05 | is.na(dat$p_value_target)
  pdat = dat[select_, ]
  
  select_ = pdat$p_value_source > 0.05 | is.na(pdat$p_value_source)
  pdat = pdat[select_, ]
  
  pdat = pdat |> select(score_name, score_type, hr, p_value, p_adj)
  } else {
    select_ = dat$p_value > 0.05 | is.na(dat$p_value)
    pdat = dat[select_, ]
    
    pdat = pdat |> select(score_name, score_type, hr, p_value, p_adj)
  }
  
  # Prepare and order your data
  order = order(dat$p_value)
  dat = dat[order, ]
  
  df <- dat %>%
    arrange(p_value) %>%
    mutate(score_name = factor(score_name, levels = dat$score_name))
  
  # Compute reasonable x-axis limits (log scale)
  lims <- range(c(df$hr_lower, df$hr_upper), na.rm = TRUE)
  lims <- c(max(lims[1], 1e-3), min(lims[2], 1e3))  # clamp extremes if needed
  
  df$hr = pmin(pmax(df$hr, lims[1]), lims[2])
  df$hr_upper = pmin(df$hr_upper, lims[2])
  df$hr_lower = pmax(df$hr_lower, lims[1])
  df$signf = df$p_adj < 0.05
  
  # Forest plot
  p <- ggplot(df, aes(y = score_name)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    geom_errorbarh(
      aes(xmin = hr_lower, xmax = hr_upper),
      height = 0.2, linewidth = 0.6, color = "black", alpha = 0.8
    ) +
    geom_point(
      aes(x = hr, fill = signf),
      size = 2.5, shape = 21,
    ) +
    scale_x_log10(limits = lims) +
    labs(
      x = "Hazard ratio (log scale)",
      y = NULL,
      title = paste("Cox —", event_name, "~", score_type),
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
    filename = sub('.csv', '.pdf', scores.path),
    plot = p,
    width = 6,
    height = 6
    # dpi = 300
  )
}




