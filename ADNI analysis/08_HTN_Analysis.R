# =============================================================
# 08_HTN_Analysis.R
# Hypertension-stratified analyses and logistic regression
# =============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(gtsummary)
library(broom)

# ---- 8.1 Merge hypertension status ----
FW_data <- read.csv("FW_COX_data_final_b_value.csv", header = TRUE)
PHC     <- read.csv("ADSP PHC Composite Cardiovascular Risk Scores [ADNI1,GO,2,3].csv",
                    header = TRUE)

PHC_clean <- PHC[, c(2, 5, 6, 7, 16)]

merge_FW_covariable_HTN <- FW_data %>%
  left_join(PHC_clean, by = "PTID") %>%
  mutate(HTN = case_when(
    VISCODE2_MRI == VISCODE2 ~ PHC_Hypertension,
    VISCODE2_MRI != VISCODE2 & EXAMDATE.x > EXAMDATE.y ~ PHC_Hypertension,
    VISCODE2_MRI != VISCODE2 & EXAMDATE.x < EXAMDATE.y &
      (SBP >= 140 | DBP >= 90) ~ 1,
    TRUE ~ 0
  )) %>%
  select(-c(59:62)) %>%
  rename(VISCODE = VISCODE.x, EXAMDATE = EXAMDATE.x)

write.csv(merge_FW_covariable_HTN, "FW_HTN_AD_data_b_value.csv", row.names = FALSE)

# ---- 8.2 Descriptive statistics by diagnosis ----
FW_data <- read.csv("FW_HTN_AD_data_b_value.csv", header = TRUE) %>%
  filter(complete.cases(select(., Age, Sex, Education, APOE4_carrier, b_value))) %>%
  mutate(
    DIAGNOSIS_Dia = factor(DIAGNOSIS_Dia, levels = c("NC", "MCI", "AD")),
    Sex           = factor(Sex, levels = c("Female", "male")),
    Ethnic        = factor(Ethnic,
                           levels = c("Hispanic-or-Latino",
                                      "Not-Hispanic-or-Latino")),
    history_HTN   = factor(HTN, levels = c(0, 1),
                           labels = c("non_HTN", "HTN")),
    APOE4_carrier = factor(APOE4_carrier,
                           levels = c("0_E4_alleles", "1_E4_allele",
                                      "2_E4_alleles")),
    b_value       = factor(数据信息, levels = c("SB_first_b0", "MB")),
    AD_FELLOW_UP  = AD_FELLOW_UP / 365.25
  )

descriptive_table <- FW_data %>%
  select(DIAGNOSIS_Dia, Age, Sex, BMI, SBP, DBP, Education, MRI_FELLOW_UP,
         starts_with("FW"), APOE4_carrier, history_HTN, Ethnic) %>%
  tbl_summary(
    by        = DIAGNOSIS_Dia,
    missing   = "no",
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits    = list(all_continuous() ~ 2, all_categorical() ~ 1)
  ) %>%
  add_overall() %>%
  add_p(test = list(all_continuous() ~ "aov",
                    all_categorical() ~ "chisq.test")) %>%
  modify_header(label ~ "**Variable**") %>%
  modify_caption("**Characteristics by diagnosis group**")

# ---- 8.3 Descriptive statistics by hypertension ----
FW_data <- FW_data %>% filter(!is.na(history_HTN))

descriptive_table_HTN <- FW_data %>%
  select(history_HTN, Age, Sex, BMI, DBP, Education, MRI_FELLOW_UP,
         starts_with("FW"), APOE4_carrier, Ethnic) %>%
  tbl_summary(
    by        = history_HTN,
    missing   = "no",
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits    = list(all_continuous() ~ 2, all_categorical() ~ 1)
  ) %>%
  add_overall() %>%
  add_p(test = list(all_continuous() ~ "aov",
                    all_categorical() ~ "chisq.test")) %>%
  modify_header(label ~ "**Variable**") %>%
  modify_caption("**Characteristics by hypertension status**")

# ---- 8.4 FW comparison within hypertension group ----
FW_data_HTN <- FW_data %>%
  filter(history_HTN == "HTN") %>%
  mutate(HTN_cog = factor(
    ifelse(DIAGNOSIS_Dia %in% c("MCI", "AD"), 1, 0),
    levels = c(0, 1),
    labels = c("HTN without cognitive impairment",
               "HTN with cognitive impairment")
  ))

fw_vars <- c("FW_total", "FW_asso", "FW_comm", "FW_proj")

simple_results <- lapply(fw_vars, function(fw_var) {
  t_res <- t.test(as.formula(paste(fw_var, "~ HTN_cog")),
                  data = FW_data_HTN, var.equal = TRUE)
  FW_data_HTN %>%
    group_by(HTN_cog) %>%
    summarise(Mean = mean(!!sym(fw_var), na.rm = TRUE),
              SD   = sd(!!sym(fw_var), na.rm = TRUE),
              N    = sum(!is.na(!!sym(fw_var))), .groups = "drop") %>%
    mutate(FW_Var    = fw_var,
           t_stat    = t_res$statistic,
           p_value   = t_res$p.value,
           CI_lower  = t_res$conf.int[1],
           CI_upper  = t_res$conf.int[2]) %>%
    select(Group = HTN_cog, FW_Var, Mean, SD, N, t_stat, p_value,
           CI_lower, CI_upper)
}) %>% bind_rows()

simple_results_formatted <- simple_results %>%
  mutate(Mean     = round(Mean, 4),
         SD       = round(SD, 4),
         t_stat   = round(t_stat, 3),
         p_value  = format(p_value, scientific = TRUE, digits = 3),
         CI_lower = round(CI_lower, 4),
         CI_upper = round(CI_upper, 4))

write.csv(simple_results_formatted, "HTN_dementia_to_FW_Results.csv")

# ---- 8.5 Logistic regression: HTN + FW → cognitive impairment ----
FW_data_HTN <- FW_data_HTN %>%
  mutate(
    HTN_FW  = factor(ifelse(FW_total < median(FW_total), 0, 1),
                     levels = c(0, 1), labels = c("low", "high")),
    HTN_AFW = factor(ifelse(FW_asso  < median(FW_asso),  0, 1),
                     levels = c(0, 1), labels = c("low", "high")),
    HTN_CFW = factor(ifelse(FW_comm  < median(FW_comm),  0, 1),
                     levels = c(0, 1), labels = c("low", "high")),
    HTN_PFW = factor(ifelse(FW_proj  < median(FW_proj),  0, 1),
                     levels = c(0, 1), labels = c("low", "high"))
  )

exposures <- c("HTN_FW", "HTN_AFW", "HTN_CFW", "HTN_PFW")

results <- lapply(exposures, function(exp) {
  formula <- as.formula(paste("HTN_cog ~", exp,
                              "+ Age + Sex + SBP + APOE4_carrier + b_value"))
  model <- glm(formula, data = FW_data_HTN, family = binomial)
  tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(grepl(paste0("^", exp), term)) %>%
    mutate(Exposure = exp,
           OR_95CI  = paste0(round(estimate, 3), " (",
                             round(conf.low, 3), "-",
                             round(conf.high, 3), ")")) %>%
    select(Exposure, estimate, conf.low, conf.high, OR_95CI, p.value)
}) %>% bind_rows()

write.csv(results, "HTN_FW_to_Dementia_OR_Results_b_adjust.csv", row.names = FALSE)