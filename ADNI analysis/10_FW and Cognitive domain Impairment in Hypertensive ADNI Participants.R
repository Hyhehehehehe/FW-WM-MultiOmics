# =============================================================
# FW and Cognitive Impairment in Hypertensive ADNI Participants
# Analysis 1: t-test of FW across cognitive impairment groups
# Analysis 2: Boxplots of FW by cognitive impairment status
# Analysis 3: Logistic regression (FW high/low → impairment)
# Analysis 4: Forest plot of adjusted ORs
# =============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)

# =============================================================
# 1. Load data and derive cognitive impairment labels
# =============================================================
FW_data <- read.csv("FW_HTN_cognition_data_b_value.csv")

FW_data$Age_bin <- cut(
  FW_data$Age,
  breaks = c(51, 65, 75, 86, 100),
  right  = FALSE,
  labels = c("51-64", "65-74", "75-85", ">=85")
)
FW_data$Edu_bin <- cut(
  FW_data$Education,
  breaks = c(0, 13, 16, 20),
  right  = TRUE,
  labels = c("<=13", "14-16", ">=17")
)
FW_data$AgeEdu_Group <- paste0(FW_data$Edu_bin, "_", FW_data$Age_bin)

# Function: impairment defined as score < NC mean - 1 SD (per age/education stratum)
cog_impair_calc <- function(df, score_col, out_col) {
  nc_ref <- subset(df, DIAGNOSIS_Dia == "NC")
  
  ref_stats <- aggregate(
    formula(paste(score_col, "~ AgeEdu_Group")),
    data = nc_ref,
    FUN  = function(x) c(mean = mean(x, na.rm = TRUE),
                         sd   = sd(x,   na.rm = TRUE))
  )
  ref_stats$mean_cog <- ref_stats[[score_col]][, "mean"]
  ref_stats$sd_cog   <- ref_stats[[score_col]][, "sd"]
  ref_stats[[score_col]] <- NULL
  
  df_merge <- merge(df, ref_stats, by = "AgeEdu_Group", all.x = TRUE)
  
  df_merge[[out_col]] <- ifelse(
    !is.na(df_merge[[score_col]]) &
      df_merge[[score_col]] < (df_merge$mean_cog - df_merge$sd_cog),
    1, 0
  )
  df_merge[[out_col]] <- ifelse(is.na(df_merge[[score_col]]),
                                NA, df_merge[[out_col]])
  df_merge$mean_cog <- df_merge$sd_cog <- NULL
  
  df_merge
}

FW_data <- cog_impair_calc(FW_data, "PHC_MEM_PHC", "mem_imp")
FW_data <- cog_impair_calc(FW_data, "PHC_EXF_PHC", "exf_imp")
FW_data <- cog_impair_calc(FW_data, "PHC_LAN_PHC", "lan_imp")
FW_data <- cog_impair_calc(FW_data, "PHC_VSP_PHC", "vsp_imp")

fw_vars <- c("FW_total", "FW_asso", "FW_comm", "FW_proj")
cog_vars <- list(
  list(name = "MI",  var = "mem_imp"),
  list(name = "ED",  var = "exf_imp"),
  list(name = "LI",  var = "lan_imp"),
  list(name = "VSD", var = "vsp_imp")
)

# Hypertensive subgroup
htn_sub <- subset(FW_data, HTN == 1)

# =============================================================
# 2. t-test: FW across cognitive impairment groups
# =============================================================
t_test_table <- data.frame()

for (cog in cog_vars) {
  cog_name <- cog$name
  cog_var  <- cog$var
  
  for (fw in fw_vars) {
    
    sub_dat <- htn_sub[!is.na(htn_sub[[cog_var]]) & !is.na(htn_sub[[fw]]), ]
    group0  <- sub_dat[sub_dat[[cog_var]] == 0, ][[fw]]
    group1  <- sub_dat[sub_dat[[cog_var]] == 1, ][[fw]]
    
    t_res <- t.test(group0, group1)
    
    t_test_table <- rbind(t_test_table, data.frame(
      Cognitive_domain = cog_name,
      FW_index         = fw,
      Normal_N         = length(group0),
      Normal_Mean_SD   = paste0(round(mean(group0, na.rm = TRUE), 3), "±",
                                round(sd(group0,   na.rm = TRUE), 3)),
      Impair_N         = length(group1),
      Impair_Mean_SD   = paste0(round(mean(group1, na.rm = TRUE), 3), "±",
                                round(sd(group1,   na.rm = TRUE), 3)),
      t_statistic      = round(t_res$statistic, 10),
      P_value          = round(t_res$p.value,   10),
      CI_lower         = round(t_res$conf.int[1], 10),
      CI_upper         = round(t_res$conf.int[2], 10)
    ))
  }
}

print(t_test_table, row.names = FALSE)

write.csv(t_test_table, "HTN_cognition_to_FW_Results.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(FW_data, "hypertension_FW_cognition_impaired_data.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

# =============================================================
# 3. Boxplots of FW by cognitive impairment status
# =============================================================
output_dir <- "boxplots"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cog_domains <- list(
  list(name = "Memory",       var = "mem_imp"),
  list(name = "Executive",    var = "exf_imp"),
  list(name = "Language",     var = "lan_imp"),
  list(name = "Visuospatial", var = "vsp_imp")
)

custom_colors <- c(
  "1" = "#B2C4D8", "2" = "#6A8CB2",   # FW_total
  "3" = "#B8D1A3", "4" = "#7CAA56",   # FW_asso
  "5" = "#A9D3DF", "6" = "#449EB6",   # FW_comm
  "7" = "#D9A7AB", "8" = "#B65259"    # FW_proj
)

for (cog in cog_domains) {
  
  cog_name <- cog$name
  cog_var  <- cog$var
  
  temp_data <- htn_sub[!is.na(htn_sub[[cog_var]]), ]
  temp_data <- temp_data[, c(cog_var, fw_vars)]
  temp_data[[cog_var]] <- as.factor(temp_data[[cog_var]])
  
  result_clean <- data.frame(
    FW_value   = c(temp_data[[fw_vars[1]]], temp_data[[fw_vars[2]]],
                   temp_data[[fw_vars[3]]], temp_data[[fw_vars[4]]]),
    Impairment = factor(rep(temp_data[[cog_var]], 4)),
    Variable   = factor(c(
      rep(fw_vars[1], nrow(temp_data)),
      rep(fw_vars[2], nrow(temp_data)),
      rep(fw_vars[3], nrow(temp_data)),
      rep(fw_vars[4], nrow(temp_data))
    ), levels = fw_vars)
  ) %>% na.omit()
  
  result_clean <- result_clean %>%
    mutate(
      color_group = case_when(
        Variable == fw_vars[1] & Impairment == 0 ~ 1,
        Variable == fw_vars[1] & Impairment == 1 ~ 2,
        Variable == fw_vars[2] & Impairment == 0 ~ 3,
        Variable == fw_vars[2] & Impairment == 1 ~ 4,
        Variable == fw_vars[3] & Impairment == 0 ~ 5,
        Variable == fw_vars[3] & Impairment == 1 ~ 6,
        Variable == fw_vars[4] & Impairment == 0 ~ 7,
        Variable == fw_vars[4] & Impairment == 1 ~ 8,
        TRUE ~ NA_integer_
      )
    )
  
  result_clean$color_group <- as.factor(result_clean$color_group)
  
  p <- ggplot(result_clean,
              aes(x = Variable, y = FW_value, fill = color_group)) +
    geom_boxplot(alpha = 0.8, outlier.size = -1,
                 width = 0.60, color = "#525252") +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3,
                 color = "red", position = position_dodge(width = 0.6)) +
    scale_fill_manual(values = custom_colors, guide = "none") +
    scale_x_discrete(
      limits = fw_vars,
      labels = c("Total", "Association", "Commissural", "Projection")
    ) +
    labs(y = "FW value", title = paste0(cog_name, " Domain")) +
    theme_classic() +
    theme(
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      panel.background  = element_blank(),
      axis.line         = element_line(color = "black", linewidth = 0.7),
      axis.ticks.length = unit(0.15, "cm"),
      axis.ticks        = element_line(linewidth = 0.5),
      axis.text.x       = element_text(size = 10, color = "black"),
      axis.text.y       = element_text(size = 10, color = "black"),
      axis.title.x      = element_blank(),
      axis.title.y      = element_text(size = 12, color = "black", face = "bold"),
      plot.title        = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position   = "none",
      plot.margin       = margin(10, 10, 10, 10)
    )
  
  ggsave(file.path(output_dir, paste0("HTN_", cog_name, "_FW_Boxplot.pdf")),
         plot = p, width = 10, height = 12, units = "cm", dpi = 300)
}

# =============================================================
# 4. Logistic regression: FW high/low → cognitive impairment
# =============================================================
FW_data <- read.csv("hypertension_FW_cognition_impaired_data.csv")
htn_sub <- subset(FW_data, HTN == 1)

# Median split for each FW metric (outside the loop)
for (fw in fw_vars) {
  med <- median(htn_sub[[fw]], na.rm = TRUE)
  htn_sub[[paste0("HTN_", fw)]] <- ifelse(htn_sub[[fw]] < med, 0, 1)
}

exposures <- paste0("HTN_", fw_vars)
fw_labels <- fw_vars

for (exp in exposures) {
  htn_sub[[exp]] <- factor(htn_sub[[exp]], levels = c(0, 1),
                           labels = c("low", "high"))
}

results_all <- data.frame()

for (cog in cog_vars) {
  cog_name <- cog$name
  cog_var  <- cog$var
  
  for (j in seq_along(exposures)) {
    
    exp      <- exposures[j]
    fw_label <- fw_labels[j]
    med      <- median(htn_sub[[fw_label]], na.rm = TRUE)
    
    formula <- as.formula(paste(cog_var, "~", exp,
                                "+ Age + Sex + SBP + APOE4_carrier + b_value"))
    model <- glm(formula, data = htn_sub, family = binomial)
    
    model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
      filter(grepl(paste0("^", exp), term))
    
    low_dat  <- htn_sub[htn_sub[[exp]] == "low",  ]
    high_dat <- htn_sub[htn_sub[[exp]] == "high", ]
    
    results_all <- rbind(results_all, data.frame(
      FW_metric          = fw_label,
      median_cutoff      = round(med, 4),
      cognitive_outcome  = cog_name,
      total_valid_N      = nrow(htn_sub),
      low_FW_N           = nrow(low_dat),
      low_FW_impair_N    = sum(low_dat[[cog_var]]  == 1, na.rm = TRUE),
      low_FW_impair_pct  = round(mean(low_dat[[cog_var]]  == 1, na.rm = TRUE) * 100, 2),
      high_FW_N          = nrow(high_dat),
      high_FW_impair_N   = sum(high_dat[[cog_var]] == 1, na.rm = TRUE),
      high_FW_impair_pct = round(mean(high_dat[[cog_var]] == 1, na.rm = TRUE) * 100, 2),
      adjusted_OR        = round(model_summary$estimate, 3),
      OR_CI_lower        = round(model_summary$conf.low, 3),
      OR_CI_upper        = round(model_summary$conf.high, 3),
      P_value            = round(model_summary$p.value, 4)
    ))
  }
}

print(results_all)

write.csv(results_all, "HTN_FW_to_cognition_OR_Results_b_adjust.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

# =============================================================
# 5. Forest plot of adjusted ORs
# =============================================================
data <- read.csv("HTN_FW_to_cognition_OR_Results_b_adjust.csv")

data <- data %>%
  mutate(
    cognitive_outcome = factor(cognitive_outcome,
                               levels = c("VSD", "ED", "MI", "LI")),
    FW_metric         = factor(FW_metric,
                               levels = c("FW_proj", "FW_comm",
                                          "FW_asso", "FW_total"))
  )

data_clean <- data %>%
  select(FW_metric, cognitive_outcome, adjusted_OR, OR_CI_lower, OR_CI_upper) %>%
  rename(Varnames = FW_metric, Factor = cognitive_outcome,
         OR = adjusted_OR, Lower = OR_CI_lower, Upper = OR_CI_upper)

p_forest <- ggplot(data_clean,
                   aes(x = Varnames, y = OR,
                       color = Factor, shape = Factor)) +
  geom_pointrange(aes(ymin = Lower, ymax = Upper),
                  position = position_dodge(width = 0.7),
                  size = 0.5, linewidth = 0.7) +
  geom_point(position = position_dodge(width = 0.7),
             size = 2.5, stroke = 0.9) +
  scale_color_manual(values = c(
    "MI"  = "#A44A46",
    "VSD" = "#446C74",
    "ED"  = "#2A437A",
    "LI"  = "#AB779D"
  )) +
  scale_shape_manual(values = c(
    "MI"  = 16,
    "VSD" = 17,
    "ED"  = 15,
    "LI"  = 18
  )) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  labs(y = "OR (95% CI)", x = "") +
  scale_y_continuous(limits = c(0.9, 4), breaks = seq(1, 4, 1)) +
  coord_flip() +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme_classic() +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.background  = element_blank(),
    axis.line         = element_line(color = "black", linewidth = 0.7),
    axis.ticks.length = unit(0.15, "cm"),
    axis.ticks        = element_line(linewidth = 0.5),
    legend.position   = "right",
    legend.title      = element_blank(),
    plot.margin       = margin(10, 10, 10, 10)
  )

print(p_forest)

ggsave("FW_cognitive_OR.pdf",
       plot = p_forest, width = 16, height = 14, units = "cm", dpi = 300)