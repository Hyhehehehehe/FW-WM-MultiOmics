# =============================================================
# 05_Cox_NC_to_AD.R
# Cox regression: continuous FW → AD conversion in NC
# =============================================================

library(survival)
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
  filter(!is.na(AD_change), DIAGNOSIS_Dia == "NC") %>%
  mutate(AD_FELLOW_UP = AD_FELLOW_UP / 365.25)

model_data <- AD_FW_data %>%
  filter(complete.cases(select(., AD_FELLOW_UP, AD_change, Age, Sex,
                               Education, APOE4_carrier, b_value))) %>%
  mutate(across(starts_with("FW"), ~ as.numeric(scale(.x))))

fw_vars <- grep("FW", names(model_data), value = TRUE)

run_cox <- function(var, data, time, event) {
  formula <- as.formula(paste0("Surv(", time, ", ", event, ") ~ ", var,
                               " + Age + Sex + Education + APOE4_carrier + b_value"))
  model <- tryCatch(
    coxph(formula, data = data, control = coxph.control(iter.max = 30)),
    error = function(e) NULL, warning = function(w) NULL
  )
  if (is.null(model) || any(is.na(coef(model))) || model$iter[1] == 30) return(NULL)
  
  coef_summary <- summary(model)$coefficients
  rows <- grep(paste0("^", var, "$"), rownames(coef_summary))
  if (!length(rows)) return(NULL)
  
  row <- coef_summary[rows[1], ]
  ci  <- tryCatch(exp(confint(model)[rows[1], ]), error = function(e) c(NA, NA))
  zp  <- tryCatch(cox.zph(model), error = function(e) NULL)
  ph_p <- if (!is.null(zp)) zp$table[grep(paste0("^", var, "$"),
                                          rownames(zp$table)), "p"] else NA
  
  data.frame(
    Variable        = var,
    HR              = exp(row[1]),
    HR_CI_low       = ci[1],
    HR_CI_high      = ci[2],
    P_value         = row[5],
    PH_assumption_p = ph_p,
    Event_count     = model$nevent
  )
}

cox_results <- do.call(rbind, lapply(fw_vars, run_cox,
                                     data = model_data,
                                     time = "AD_FELLOW_UP",
                                     event = "AD_change"))

cox_results <- cox_results %>%
  mutate(HR_95CI = paste0(round(HR, 3), " (",
                          round(HR_CI_low, 3), "-",
                          round(HR_CI_high, 3), ")"),
         P_value = round(P_value, 4),
         PH_assumption_p = round(PH_assumption_p, 4))

write.csv(cox_results, "FW_continuous_Cox_NC_AD_model3_b_adjust.csv", row.names = FALSE)