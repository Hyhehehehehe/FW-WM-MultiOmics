# =============================================================
# FW and Cognitive Diagnosis in Hypertensive ADNI Participants
# Analysis 1: t-test comparing FW by HTN with/without cognitive impairment
# Analysis 2: Boxplots of FW by HTN with/without cognitive impairment
# Analysis 3: Logistic regression (FW high/low → cognitive impairment)
# Analysis 4: Forest plot of ORs
# =============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)

# ---- Load data ----
FW_data <- read.csv("FW_HTN_AD_data_b_value.csv", header = TRUE)

FW_data <- FW_data[complete.cases(FW_data[, c("Age", "Sex", "Education",
                                              "APOE4_carrier")]), ]

FW_data <- FW_data %>%
  mutate(
    DIAGNOSIS_Dia = factor(DIAGNOSIS_Dia,
                           levels = c("NC", "MCI", "AD")),
    Sex           = factor(Sex, levels = c("Female", "male")),
    Ethnic        = factor(Ethnic,
                           levels = c("Hispanic-or-Latino",
                                      "Not-Hispanic-or-Latino")),
    APOE4_carrier = factor(APOE4_carrier,
                           levels = c("0_E4_alleles", "1_E4_allele",
                                      "2_E4_alleles")),
    b_value       = factor(Data_Type, levels = c("SB_first_b0", "MB"))
  )

FW_data$AD_FELLOW_UP <- FW_data$AD_FELLOW_UP / 365.25

# =============================================================
# 1. Define hypertensive subgroup and cognitive grouping
# =============================================================

htn_sub <- FW_data %>%
  filter(HTN == 1, !is.na(HTN)) %>%
  mutate(
    HTN_cog = factor(
      ifelse(DIAGNOSIS_Dia %in% c("MCI", "AD"), 1, 0),
      levels = c(0, 1),
      labels = c("HTN_non_AD", "HTN_with_AD")
    )
  )

table(htn_sub$HTN_cog)

fw_vars <- c("FW_total", "FW_asso", "FW_comm", "FW_proj")

# =============================================================
# 2. t-test: FW across HTN cognitive groups
# =============================================================

simple_results <- data.frame()

for (fw_var in fw_vars) {
  
  t_res <- t.test(as.formula(paste(fw_var, "~ HTN_cog")),
                  data = htn_sub, var.equal = TRUE)
  
  stats <- htn_sub %>%
    group_by(HTN_cog) %>%
    summarise(
      Mean = mean(!!sym(fw_var), na.rm = TRUE),
      SD   = sd(!!sym(fw_var),   na.rm = TRUE),
      N    = sum(!is.na(!!sym(fw_var))),
      .groups = "drop"
    ) %>%
    mutate(
      FW_Var   = fw_var,
      Group    = as.character(HTN_cog),
      t_stat   = t_res$statistic,
      p_value  = t_res$p.value,
      CI_lower = t_res$conf.int[1],
      CI_upper = t_res$conf.int[2]
    ) %>%
    select(Group, FW_Var, Mean, SD, N, t_stat, p_value, CI_lower, CI_upper)
  
  simple_results <- bind_rows(simple_results, stats)
}

simple_results_formatted <- simple_results %>%
  arrange(FW_Var, Group) %>%
  mutate(
    Mean     = round(Mean, 4),
    SD       = round(SD, 4),
    t_stat   = round(t_stat, 3),
    p_value  = format(p_value, scientific = TRUE, digits = 3),
    CI_lower = round(CI_lower, 4),
    CI_upper = round(CI_upper, 4)
  )

print(simple_results_formatted)

write.csv(simple_results_formatted,
          "HTN_dementia_to_FW_Results.csv",
          row.names = FALSE)

# =============================================================
# 3. Boxplots: FW by HTN cognitive status
# =============================================================

output_dir <- "boxplots"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

temp_data <- htn_sub[, c("HTN_cog", fw_vars)]
temp_data <- temp_data[complete.cases(temp_data[, fw_vars]), ]

result_clean <- data.frame(
  FW_value = c(temp_data[[fw_vars[1]]], temp_data[[fw_vars[2]]],
               temp_data[[fw_vars[3]]], temp_data[[fw_vars[4]]]),
  Group    = rep(temp_data$HTN_cog, 4),
  Variable = factor(c(
    rep(fw_vars[1], nrow(temp_data)),
    rep(fw_vars[2], nrow(temp_data)),
    rep(fw_vars[3], nrow(temp_data)),
    rep(fw_vars[4], nrow(temp_data))
  ), levels = fw_vars)
) %>% na.omit()

result_clean <- result_clean %>%
  mutate(
    color_group = case_when(
      Variable == fw_vars[1] & Group == "HTN_non_AD"  ~ 1,
      Variable == fw_vars[1] & Group == "HTN_with_AD" ~ 2,
      Variable == fw_vars[2] & Group == "HTN_non_AD"  ~ 3,
      Variable == fw_vars[2] & Group == "HTN_with_AD" ~ 4,
      Variable == fw_vars[3] & Group == "HTN_non_AD"  ~ 5,
      Variable == fw_vars[3] & Group == "HTN_with_AD" ~ 6,
      Variable == fw_vars[4] & Group == "HTN_non_AD"  ~ 7,
      Variable == fw_vars[4] & Group == "HTN_with_AD" ~ 8,
      TRUE ~ NA_integer_
    )
  )

result_clean$color_group <- as.factor(result_clean$color_group)

custom_colors <- c(
  "1" = "#B2C4D8", "2" = "#6A8CB2",   # FW_total
  "3" = "#B8D1A3", "4" = "#7CAA56",   # FW_asso
  "5" = "#A9D3DF", "6" = "#449EB6",   # FW_comm
  "7" = "#D9A7AB", "8" = "#B65259"    # FW_proj
)

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
  labs(y = "FW value",
       title = "FW Values by Hypertension with Cognitive Impairment Status") +
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

print(p)

ggsave(file.path(output_dir, "HTN_dementia_FW_Boxplot.pdf"),
       plot = p, width = 10, height = 12, units = "cm", dpi = 300)

ggsave(file.path(output_dir, "HTN_dementia_FW_Boxplot.png"),
       plot = p, width = 10, height = 12, units = "cm", dpi = 300)

# =============================================================
# 4. Logistic regression: FW high/low → cognitive impairment
# =============================================================

htn_sub <- htn_sub %>%
  mutate(
    HTN_FW  = factor(ifelse(FW_total < median(FW_total, na.rm = TRUE), 0, 1),
                     levels = c(0, 1), labels = c("low", "high")),
    HTN_AFW = factor(ifelse(FW_asso  < median(FW_asso,  na.rm = TRUE), 0, 1),
                     levels = c(0, 1), labels = c("low", "high")),
    HTN_CFW = factor(ifelse(FW_comm  < median(FW_comm,  na.rm = TRUE), 0, 1),
                     levels = c(0, 1), labels = c("low", "high")),
    HTN_PFW = factor(ifelse(FW_proj  < median(FW_proj,  na.rm = TRUE), 0, 1),
                     levels = c(0, 1), labels = c("low", "high"))
  )

write.csv(htn_sub,
          "FW_HTN_Dementia_basic_information_b_value.csv",
          row.names = FALSE)

exposures <- c("HTN_FW", "HTN_AFW", "HTN_CFW", "HTN_PFW")

results <- data.frame()

for (exp in exposures) {
  
  formula <- as.formula(paste("HTN_cog ~", exp,
                              "+ Sex + Age + SBP + APOE4_carrier + b_value"))
  
  model <- glm(formula, data = htn_sub, family = binomial)
  
  model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(grepl(paste0("^", exp), term)) %>%
    mutate(
      Exposure = exp,
      OR_95CI  = paste0(round(estimate, 3), " (",
                        round(conf.low, 3), "-",
                        round(conf.high, 3), ")")
    ) %>%
    select(Exposure, estimate, conf.low, conf.high, OR_95CI, p.value)
  
  results <- rbind(results, model_summary)
}

print(results)

write.csv(results,
          "HTN_FW_to_Dementia_OR_Results_b_adjust.csv",
          row.names = FALSE)

# =============================================================
# 5. Forest plot: FW cognitive OR (ADNI)
# =============================================================

data <- read.csv("HTN_FW_to_Dementia_OR_Results_b_adjust.csv")

custom_levels <- c("HTN_FW", "HTN_AFW", "HTN_CFW", "HTN_PFW")

data <- data %>%
  mutate(Exposure = factor(Exposure, levels = custom_levels))

data_clean <- data %>%
  select(Exposure, estimate, conf.low, conf.high) %>%
  rename(Varnames = Exposure, OR = estimate, Lower = conf.low, Upper = conf.high)

p_forest <- ggplot(data_clean,
                   aes(x = Varnames, y = OR,
                       color = Varnames, shape = Varnames)) +
  geom_pointrange(aes(ymin = Lower, ymax = Upper),
                  position = position_dodge(width = 0.7),
                  size = 0.5, linewidth = 0.7) +
  geom_point(position = position_dodge(width = 0.7),
             size = 2.5, stroke = 0.9) +
  scale_color_manual(values = c(
    "HTN_FW"  = "#A44A46",
    "HTN_AFW" = "#446C74",
    "HTN_CFW" = "#2A437A",
    "HTN_PFW" = "#AB779D"
  )) +
  scale_shape_manual(values = c(
    "HTN_FW"  = 16,
    "HTN_AFW" = 17,
    "HTN_CFW" = 15,
    "HTN_PFW" = 18
  )) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  labs(y = "OR (95% CI)", x = "") +
  scale_y_continuous(limits = c(1, 5), breaks = seq(1, 5, 1)) +
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

ggsave("FW_dementia_OR.pdf",
       plot = p_forest, width = 16, height = 14, units = "cm", dpi = 300)