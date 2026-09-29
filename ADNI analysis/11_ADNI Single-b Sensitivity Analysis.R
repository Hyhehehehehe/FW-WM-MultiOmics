# =============================================================
# ADNI Single-b Sensitivity Analysis: FW Continuous Cox Regression
# Purpose: Replace multi-b FW with single-b FW and rerun the 4 Cox models
# =============================================================

library(readxl)
library(dplyr)
library(survival)

# =============================================================
# Part 1: Replace multi-b FW with single-b FW
# =============================================================

FW_single_b <- read_excel("multi_b_to_single_b.xlsx")
FW_data     <- read.csv("FW_cox_ADNI_bvalue.csv")

fw_pairs <- list(
  list(new = "WM_mask",       old = "FW_total"),
  list(new = "3fiber_rois_1", old = "FW_asso"),
  list(new = "3fiber_rois_2", old = "FW_comm"),
  list(new = "3fiber_rois_3", old = "FW_proj")
)

for (pair in fw_pairs) {
  FW_single_sub <- FW_single_b %>% select(ID, all_of(pair$new))
  FW_data <- FW_data %>%
    left_join(FW_single_sub, by = c("ImageID" = "ID")) %>%
    mutate(!!pair$old := coalesce(.data[[pair$new]], .data[[pair$old]])) %>%
    select(-all_of(pair$new))
}


# Standardize the b_value label column
b_col <- grep("Data_Type|data_type|DataInformation", names(FW_data), value = TRUE)
if (length(b_col) == 1) {
  FW_data[[b_col]][is.na(FW_data[[b_col]])] <- "SB_first_b0"
  FW_data[[b_col]][FW_data[[b_col]] == "SB_last_b0"] <- "SB_first_b0"
  names(FW_data)[names(FW_data) == b_col] <- "b_value"
}

write.csv(FW_data,
          "FW_COX_data_final_new.csv",
          row.names = FALSE)

# =============================================================
# Part 2: Shared helper functions
# =============================================================

prepare_fw_data <- function(path) {
  df <- read.csv(path, header = TRUE)
  df$DIAGNOSIS_Dia <- factor(df$DIAGNOSIS_Dia, levels = c("NC", "MCI", "AD"))
  df$Sex           <- factor(df$Sex, levels = c("Female", "male"))
  df$APOE4_carrier <- factor(df$APOE4_carrier,
                             levels = c("0_E4_alleles", "1_E4_allele",
                                        "2_E4_alleles"))
  df$b_value       <- factor(df$b_value, levels = c("SB_first_b0", "MB"))
  df
}

run_cox_fw <- function(model_data, time_col, event_col) {
  
  fw_cols <- grep("^FW", names(model_data), value = TRUE)
  model_data <- model_data %>%
    mutate(across(all_of(fw_cols), ~ as.numeric(scale(.x))))
  
  cox_results <- data.frame()
  failed_vars <- c()
  
  for (var in fw_cols) {
    
    formula <- as.formula(paste0(
      "Surv(", time_col, ", ", event_col, ") ~ ", var,
      " + Age + Sex + Education + APOE4_carrier + b_value"
    ))
    
    model <- tryCatch(
      coxph(formula, data = model_data, control = coxph.control(iter.max = 30)),
      error   = function(e) NULL,
      warning = function(w) suppressWarnings(
        coxph(formula, data = model_data, control = coxph.control(iter.max = 30))
      )
    )
    
    if (is.null(model) || any(is.na(coef(model))) || model$iter[1] == 30) {
      failed_vars <- c(failed_vars, var)
      next
    }
    
    coef_summary <- summary(model)$coefficients
    var_rows <- grep(paste0("^", var, "$"), rownames(coef_summary))
    if (length(var_rows) == 0) {
      failed_vars <- c(failed_vars, var)
      next
    }
    
    row   <- coef_summary[var_rows[1], ]
    hr    <- exp(row[1])
    ci    <- tryCatch(exp(confint(model)[var_rows[1], ]),
                      error = function(e) c(NA, NA))
    p_val <- row[5]
    
    ph_p <- NA
    zp <- tryCatch(cox.zph(model), error = function(e) NULL)
    if (!is.null(zp)) {
      zp_row <- grep(paste0("^", var, "$"), rownames(zp$table))
      if (length(zp_row) > 0) ph_p <- zp$table[zp_row[1], "p"]
    }
    
    cox_results <- rbind(cox_results, data.frame(
      Variable        = var,
      HR              = hr,
      HR_CI_low       = ci[1],
      HR_CI_high      = ci[2],
      P_value         = p_val,
      PH_assumption_p = ph_p,
      Event_count     = model$nevent
    ))
  }
  
  cat("Successful:", nrow(cox_results),
      " | Failed:", length(failed_vars), "\n")
  if (length(failed_vars) > 0) print(failed_vars)
  
  cox_results %>%
    mutate(
      HR_95CI         = paste0(round(HR, 3), " (",
                               round(HR_CI_low, 3), "-",
                               round(HR_CI_high, 3), ")"),
      P_value         = round(P_value, 4),
      PH_assumption_p = round(PH_assumption_p, 4)
    )
}

# =============================================================
# Part 3: Four Cox models
# =============================================================

fw_path <- "FW_COX_data_final_new.csv"
out_dir <- "all_SB_continuous_with_sd"

# ---------- 3.1 NC/MCI → AD ----------
FW_data <- prepare_fw_data(fw_path)
FW_data$AD_FELLOW_UP <- FW_data$AD_FELLOW_UP / 365.25

AD_FW_data <- FW_data %>%
  filter(!is.na(AD_change), DIAGNOSIS_Dia != "AD")

model_data <- AD_FW_data[
  complete.cases(AD_FW_data[, c("AD_FELLOW_UP", "AD_change", "Age", "Sex",
                                "Education", "APOE4_carrier", "b_value")]),
]

cox_results <- run_cox_fw(model_data, "AD_FELLOW_UP", "AD_change")
print(cox_results)

write.csv(cox_results,
          file.path(out_dir, "FW_continuous_Cox.csv"),
          row.names = FALSE)

# ---------- 3.2 NC → MCI ----------
FW_data <- prepare_fw_data(fw_path)
FW_data$MCI_FELLOW_UP <- FW_data$MCI_FELLOW_UP / 365.25

MCI_FW_data <- FW_data %>%
  filter(!is.na(MCI_change), DIAGNOSIS_Dia == "NC")

model_data <- MCI_FW_data[
  complete.cases(MCI_FW_data[, c("MCI_FELLOW_UP", "MCI_change", "Age", "Sex",
                                 "Education", "APOE4_carrier", "b_value")]),
]

cox_results <- run_cox_fw(model_data, "MCI_FELLOW_UP", "MCI_change")
print(cox_results)

write.csv(cox_results,
          file.path(out_dir, "FW_continuous_Cox_MCI.csv"),
          row.names = FALSE)

# ---------- 3.3 NC → AD ----------
FW_data <- prepare_fw_data(fw_path)
FW_data$AD_FELLOW_UP <- FW_data$AD_FELLOW_UP / 365.25

AD_FW_data <- FW_data %>%
  filter(!is.na(MCI_change), DIAGNOSIS_Dia == "NC")

model_data <- AD_FW_data[
  complete.cases(AD_FW_data[, c("AD_FELLOW_UP", "AD_change", "Age", "Sex",
                                "Education", "APOE4_carrier", "b_value")]),
]

cox_results <- run_cox_fw(model_data, "AD_FELLOW_UP", "AD_change")
print(cox_results)

write.csv(cox_results,
          file.path(out_dir, "FW_continuous_Cox_NC_AD.csv"),
          row.names = FALSE)

# ---------- 3.4 MCI → AD ----------
FW_data <- prepare_fw_data(fw_path)
FW_data$AD_FELLOW_UP <- FW_data$AD_FELLOW_UP / 365.25

AD_FW_data <- FW_data %>%
  filter(!is.na(AD_change), DIAGNOSIS_Dia == "MCI")

model_data <- AD_FW_data[
  complete.cases(AD_FW_data[, c("AD_FELLOW_UP", "AD_change", "Age", "Sex",
                                "Education", "APOE4_carrier", "b_value")]),
]

cox_results <- run_cox_fw(model_data, "AD_FELLOW_UP", "AD_change")
print(cox_results)

write.csv(cox_results,
          file.path(out_dir, "FW_continuous_Cox_MCI_AD.csv"),
          row.names = FALSE)