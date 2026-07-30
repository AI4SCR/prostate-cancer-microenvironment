
clinical.path = file.path('/Users/me3312/Documents/Paper_PCa/0-paper/0-export/clinical.parquet')
old_clinical = read_parquet(clinical.path)
new_clinical_path <- "/Users/me3312/Desktop/check_clinical/0-export/clinical.parquet"
new_clinical <- read_parquet(new_clinical_path)

## check if equal
all.equal(old_clinical, new_clinical)
table(new_clinical$disease_progr, new_clinical$clinical_progr)

clinical = read_parquet(clinical.path)

result_dir = "/Users/me3312/Documents/Paper_PCa/5-niches/barplot_data/"
path_patient <- "/Users/me3312/Documents/Paper_PCa/5-niches/barplot_data/metadata_clustered_pat_id_label.csv"
path_tma <- "/Users/me3312/Documents/Paper_PCa/5-niches/barplot_data/metadata_clustered_tma_id_label.csv"

df_patient <- read.csv(path_patient)
df_tma <- read.csv(path_tma)

path_groups <- "/Users/me3312/Desktop/check_clinical/0-export/survival-cell-freq-groups.parquet"
path_gleason <- "/Users/me3312/Desktop/check_clinical/0-export/survival-gleason.parquet"
path_histo <- "/Users/me3312/Desktop/check_clinical/0-export/survival-stromogenic-inflammation.parquet"

groups <- read_parquet(path_groups)
gleason <- read_parquet(path_gleason)
histo <- read_parquet(path_histo)

df_check <- df_patient %>%
  select(pat_id, cluster_group) %>%
  filter(cluster_group == "C4") 
df_check_2 <- groups %>%
  select(pat_id, cluster_id) %>%
  filter(cluster_id == "5")


final_path <- "/Users/me3312/Documents/Paper_PCa/5-niches/barplot_data/metadata_with_dendrogram_colors_label_pat_id.parquet"
df_patient <- read_parquet(final_path)

library(dplyr)

df <- df_patient %>%
  filter(leaf_color_group != "black")# %>%
  # group_by(leaf_color_group) %>%
  # filter(n_distinct(pat_id) >= 5) %>%
  # ungroup()
  

df_colors = df %>%
  select(leaf_color_group, leaf_color) %>%
  distinct() %>%
  arrange(leaf_color_group)

custom_palette <- df_colors$leaf_color
names(custom_palette) <- df_colors$leaf_color_group

## merge with disease_progr, last_fu
clinical_time <- clinical %>%
  select(pat_id, disease_progr_time, last_fu) %>%
  distinct()

df <- df %>%
  left_join(clinical_time, by = "pat_id")

library(survival)
library(survminer)
df$cluster_group <- factor(df$leaf_color_group, levels = names(custom_palette))
fit <- survfit(Surv(disease_progr_time, disease_progr) ~ cluster_group, data = df)

# Extract order used internally by survfit
strata_order <- names(fit$strata)
strata_order <- gsub("cluster_group=", "", strata_order)
names(custom_palette) <- names(fit$strata)
p_prog <- ggsurvplot(
  fit,
  data = df,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = custom_palette,
  xlab = "Time",
  ylab = "Progression-free survival probability",
  legend.title = "Group",
  risk.table.height = 0.25
)
plot_name <- "progr_patient_cluster_group_final_all.pdf"
plot_path <- file.path(result_dir, plot_name)
#ggsave(plot_path, p_prog$plot, width = 10, height = 6)
p_prog$plot

df$overall_survival <- ifelse(df$os_status == 'alive', 0, 1)
fit <- survfit(Surv(last_fu, overall_survival) ~ cluster_group, data = df)

# Extract order used internally by survfit
strata_order <- names(fit$strata)
strata_order <- gsub("cluster_group=", "", strata_order)
names(custom_palette) <- names(fit$strata)
p_survival <- ggsurvplot(
  fit,
  data = df,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = custom_palette,
  xlab = "Time",
  ylab = "Survival probability",
  legend.title = "Group",
  risk.table.height = 0.25
)
plot_name <- "survival_patient_cluster_group_final_all.pdf"
plot_path <- file.path(result_dir, plot_name)
#ggsave(plot_path, p_survival$plot, width = 10, height = 6)
p_survival$plot


########################## INFLAMMATION #############
library(dplyr)
dir_inflam <- "/Users/me3312/Documents/Paper_PCa/5-niches/kaplan_meier/inflammation"
df_inflam <- clinical %>%
  select(pat_id, inflammation) %>%
  filter(!is.na(inflammation)) %>%
  mutate(
    inflammation_bin = case_when(
      inflammation == "yes" ~ 1,
      inflammation == "no"  ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  group_by(pat_id) %>%
  summarise(
    inflammation = max(inflammation_bin, na.rm = TRUE)
  ) %>%
  ungroup()

df_outcome <- clinical %>%
  select(pat_id, disease_progr, disease_progr_time, last_fu, os_status) %>%
  distinct()

df_final <- df_inflam %>%
  left_join(df_outcome, by = "pat_id")
df_final$overall_survival <- ifelse(df_final$os_status == 'alive', 0, 1)

fit <- survfit(Surv(disease_progr_time, disease_progr) ~ inflammation, data = df_final)
p_prog_inflam <- ggsurvplot(
  fit,
  data = df_final,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = c("blue", "red"),
  xlab = "Time",
  ylab = "Progression-free survival probability",
  legend.title = "Inflammation",
  risk.table.height = 0.25
)
plot_name <- "progr_inflammation.pdf"
plot_path <- file.path(dir_inflam, plot_name)
#ggsave(plot_path, p_prog_inflam$plot, width = 6, height = 4)
p_prog_inflam$plot
fit <- survfit(Surv(last_fu, overall_survival) ~ inflammation, data = df_final)
p_survival_inflam <- ggsurvplot(
  fit,
  data = df_final,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = c("blue", "red"),
  xlab = "Time",
  ylab = "Survival probability",
  legend.title = "Inflammation",
  risk.table.height = 0.25
)
plot_name <- "survival_inflammation.pdf"
plot_path <- file.path(dir_inflam, plot_name)
#ggsave(plot_path, p_survival_inflam$plot, width = 6, height = 4)
p_survival_inflam$plot


################## STROMOGENIC #############

dir_stromo <- "/Users/me3312/Documents/Paper_PCa/5-niches/kaplan_meier/stromogenic"
df_stromo <- clinical %>%
  select(pat_id, stromogenic_smc_loss_reactive_stroma_present) %>%
  filter(!is.na(stromogenic_smc_loss_reactive_stroma_present)) %>%
  mutate(
    stromogenic_bin = case_when(
      stromogenic_smc_loss_reactive_stroma_present == "yes" ~ 1,
      stromogenic_smc_loss_reactive_stroma_present == "no"  ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  group_by(pat_id) %>%
  summarise(
    stromogenic = max(stromogenic_bin, na.rm = TRUE)
  ) %>%
  ungroup()

df_outcome <- clinical %>%
  select(pat_id, disease_progr, disease_progr_time, last_fu, os_status) %>%
  distinct()

df_final <- df_stromo %>%
  left_join(df_outcome, by = "pat_id")
df_final$overall_survival <- ifelse(df_final$os_status == 'alive', 0, 1)

fit <- survfit(Surv(disease_progr_time, disease_progr) ~ stromogenic, data = df_final)
p_prog_stromo <- ggsurvplot(
  fit,
  data = df_final,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = c("blue", "red"),
  xlab = "Time",
  ylab = "Progression-free survival probability",
  legend.title = "stromogenic",
  risk.table.height = 0.25
)
plot_name <- "progr_stromogenic.pdf"
plot_path <- file.path(dir_stromo, plot_name)
#ggsave(plot_path, p_prog_stromo$plot, width = 6, height = 4)
p_prog_stromo$plot
fit <- survfit(Surv(last_fu, overall_survival) ~ stromogenic, data = df_final)
p_survival_stromo <- ggsurvplot(
  fit,
  data = df_final,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = c("blue", "red"),
  xlab = "Time",
  ylab = "Survival probability",
  legend.title = "stromogenic",
  risk.table.height = 0.25
)
plot_name <- "survival_stromogenic.pdf"
plot_path <- file.path(dir_stromo, plot_name)
#ggsave(plot_path, p_survival_stromo$plot, width = 6, height = 4)
p_survival_stromo$plot


################### GLEASON #############
dir_gleason <- "/Users/me3312/Documents/Paper_PCa/5-niches/kaplan_meier/gleason"
dir.create(dir_gleason, showWarnings = FALSE)


df_final <- clinical %>%
  select(pat_id, disease_progr, disease_progr_time, last_fu, os_status, gs_grp) %>%
  filter(!is.na(gs_grp)) %>%
  filter(gs_grp != "nan") %>%
  distinct()


df_final$overall_survival <- ifelse(df_final$os_status == 'alive', 0, 1)
df_final$gs_grp <- factor(df_final$gs_grp)

fit <- survfit(Surv(disease_progr_time, disease_progr) ~ gs_grp, data = df_final)
p_prog_gleason <- ggsurvplot(
  fit,
  data = df_final,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  #palette = c("blue", "red"),
  xlab = "Time",
  ylab = "Progression-free survival probability",
  legend.title = "stromogenic",
  risk.table.height = 0.25
)
plot_name <- "progr_gleason.pdf"
plot_path <- file.path(dir_gleason, plot_name)
#ggsave(plot_path, p_prog_gleason$plot, width = 6, height = 4)
p_prog_gleason$plot
fit <- survfit(Surv(last_fu, overall_survival) ~ gs_grp, data = df_final)
p_survival_gleason <- ggsurvplot(
  fit,
  data = df_final,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  #palette = c("blue", "red"),
  xlab = "Time",
  ylab = "Survival probability",
  legend.title = "stromogenic",
  risk.table.height = 0.25
)
plot_name <- "survival_gleason.pdf"
plot_path <- file.path(dir_gleason, plot_name)
#ggsave(plot_path, p_survival_gleason$plot, width = 6, height = 4)
p_survival_gleason$plot



dir_gleason <- "/Users/me3312/Documents/Paper_PCa/5-niches/kaplan_meier/gleason"
dir.create(dir_gleason, showWarnings = FALSE)


df_hetero <- clinical %>%
  select(pat_id, tma_id, gs_grp, gleason_grp) %>%
  filter(!is.na(gs_grp)) %>%
  filter(gs_grp != "nan") %>%
  filter(!is.na(gleason_grp)) %>%
  filter(gleason_grp != "nan") %>%
  distinct()


clinical[order(clinical$sample_id), ]

filter_ = duplicated(clinical$tma_id)
tmas = clinical[!filter_, ]

filter_ = duplicated(clinical$pat_id)
pats = data = clinical[!filter_, ]

# %%
tma.cat_names = c(
  "gs_pat_1",
  "gs_pat_2",
  'gleason_grp',
  "stromogenic_smc_loss_reactive_stroma_present",
  "non_stromogenic_smc_abundant",
  "inflammation",
  "glandular_atrophy_pin",
  "cribriform",
  "is_tumor",
  "gleason_score",
  "gleason_score_sum",
  "gleason_grp"
)

cat_name = 'gleason_grp'
order_name = 'gs_grp'
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)

var_name   <- "gleason_grp"  # x-axis categories in the heatmap
order_name <- "gs_grp"       # patient-level ordering variable

# 1) patient order (unique gs_grp per patient)
pat_order <- tmas %>%
  select(pat_id, !!sym(order_name)) %>%
  filter(!is.na(.data[[order_name]]), .data[[order_name]] != "nan") %>%
  distinct() %>%
  arrange(as.numeric(.data[[order_name]]), pat_id) %>%
  mutate(pat_id = factor(pat_id, levels = pat_id))

# 2) counts long table
pdat <- as.data.frame.matrix(table(tmas[[var_name]], tmas$pat_id)) |>
  rownames_to_column(var = var_name) |>
  pivot_longer(cols = -all_of(var_name),
               names_to = "pat_id",
               values_to = "value") |>
  semi_join(pat_order, by = "pat_id") |>
  mutate(pat_id = factor(pat_id, levels = levels(pat_order$pat_id)))

# 3) main heatmap
g_main <- ggplot(
  pdat,
  aes(y = pat_id,
      x = factor(.data[[var_name]]),
      fill = factor(value, levels = 0:4))
) +
  geom_tile(color = "black", linewidth = 0.1) +
  scale_fill_manual(
    values = c("0"="white","4"="#fde725","3"="#5ec962","2"="#21918c","1"="#440154"),
    drop = FALSE,
    name = "Count"
  ) +
  labs(y = "Patient", x = var_name) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        panel.grid = element_blank())

gs_palette <- c(
  "1" = "#f6d2d2",
  "2" = "#f1abab",
  "3" = "#f96b6b",
  "4" = "#f44336",
  "5" = "#b71c1c",
  "nan" = "#9e9e9e"
)

g_strip <- pat_order %>%
  mutate(
    gs_grp = as.character(.data[[order_name]]),
    gs_grp = factor(gs_grp, levels = names(gs_palette)),
    x = "gs_grp"
  ) %>%
  ggplot(aes(x = x, y = pat_id, fill = gs_grp)) +
  geom_tile(color = "black", linewidth = 0.1) +
  scale_fill_manual(values = gs_palette, drop = FALSE) +
  labs(x = NULL, y = NULL, fill = "GS group") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()
  )
library(patchwork)

(g_strip + g_main) + plot_layout(widths = c(1, 10))


dat = table(tmas[[var_name]],tmas$pat_id)
pdat <- as.data.frame.matrix(pdat)
pdat <- tibble::rownames_to_column(pdat, var = var_name)
pdat <- pdat %>%
  pivot_longer(
    cols = -all_of(var_name),
    names_to = "patient",
    values_to = "value"
  )

pdat = pdat[order(pdat$patient),]
g = ggplot(pdat, aes(y = patient, x = factor(.data[[var_name]]),
                     fill = factor(value, levels = 0:4))) +
  geom_tile(color = "black", linewidth = 0.1) +
  scale_fill_manual(
    values = c(
      "0" = "white",
      "4" = "#fde725",
      "3" = "#5ec962",
      "2" = "#21918c",
      "1" = "#440154"
    ),
    drop = FALSE,
    name = "Count"
  ) +
  labs(y = "Patient", x = var_name) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    # axis.ticks.x = element_blank(),
    panel.grid = element_blank()
  )
