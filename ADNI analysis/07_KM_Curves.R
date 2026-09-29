# =============================================================
# 07_KM_Curves.R
# Kaplan-Meier curves for FW high/low groups (median split)
# Note: dichotomized FW is used ONLY for KM curves.
# =============================================================

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)

FW_data <- read.csv("FW_COX_data_final_b_value.csv", header = TRUE) %>%
  mutate(
    DIAGNOSIS_Dia = factor(DIAGNOSIS_Dia, levels = c("NC", "MCI", "AD")),
    Sex           = factor(Sex, levels = c("Female", "male")),
    Ethnic        = factor(Ethnic,
                           levels = c("Hispanic-or-Latino",
                                      "Not-Hispanic-or-Latino")),
    APOE4_carrier = factor(APOE4_carrier,
                           levels = c("0_E4_alleles", "1_E4_allele",
                                      "2_E4_alleles")),
    b_value       = factor(b_value, levels = c("SB_first_b0", "MB"))
  )

AD_FW_data <- FW_data %>%
  filter(!is.na(AD_change), DIAGNOSIS_Dia != "AD") %>%
  mutate(AD_FELLOW_UP = AD_FELLOW_UP / 365.25)

model_data <- AD_FW_data %>%
  filter(complete.cases(select(., AD_FELLOW_UP, AD_change, Age, Sex,
                               Education, APOE4_carrier, b_value))) %>%
  mutate(across(starts_with("FW"), ~ as.numeric(scale(.x))))

fw_vars <- grep("FW", names(model_data), value = TRUE)

# Median split for KM curves only
for (col in fw_vars) {
  med <- median(model_data[[col]], na.rm = TRUE)
  model_data[[paste0(col, "_group")]] <- factor(
    ifelse(model_data[[col]] > med, "High", "Low"),
    levels = c("Low", "High")
  )
}

group_cols <- grep("_group$", names(model_data), value = TRUE)

output_dir <- "KM_plot"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

for (var_name in group_cols) {
  temp_data <- model_data[, c("AD_FELLOW_UP", "AD_change", var_name)]
  names(temp_data)[3] <- "group"
  
  fit_km <- survfit(Surv(AD_FELLOW_UP, AD_change) ~ group, data = temp_data)
  
  p <- ggsurvplot(
    fit              = fit_km,
    data             = temp_data,
    pval             = TRUE,
    pval.method      = TRUE,
    conf.int         = TRUE,
    conf.int.style   = "ribbon",
    conf.int.alpha   = 0.15,
    risk.table       = TRUE,
    risk.table.col   = "strata",
    risk.table.height = 0.25,
    censor           = TRUE,
    censor.shape     = 124,
    censor.size      = 4.5,
    palette          = c("#2e9fdf", "#d37f88"),
    title            = paste0("Kaplan-Meier Curve for AD Conversion\n", var_name),
    xlab             = "Time (years)",
    ylab             = "Progression-Free Probability",
    legend.title     = var_name,
    legend.labs      = levels(temp_data$group),
    break.time.by    = 2,
    xlim             = c(0, max(temp_data$AD_FELLOW_UP, na.rm = TRUE) + 0.5),
    ylim             = c(0, 1),
    ggtheme          = theme_classic()
  )
  
  pdf(file.path(output_dir, paste0(var_name, "_NC_MCI_to_AD.pdf")),
      width = 10, height = 8)
  print(p)
  dev.off()
}