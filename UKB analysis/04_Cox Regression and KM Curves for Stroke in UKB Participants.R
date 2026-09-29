# =============================================================
# Cox Regression and KM Curves for Stroke in UKB
# =============================================================

library(survival)
library(survminer)
library(dplyr)
library(lubridate)

# ---- Load data ----
df <- read.csv("freewater.csv")

# =============================================================
# 1. Build survival data
# =============================================================

cutoff_date <- as.Date("2021-12-31")

# Exclude participants with stroke at baseline
df_sub <- df %>% filter(I64_2 == 0)

df_surv <- df_sub %>%
  mutate(
    baseline    = Date_attending_2,
    event_date  = if_else(!is.na(Date_I64) & Date_I64 > baseline,
                          Date_I64, as.Date(NA)),
    censor_date = if_else(!is.na(Date_death_0), Date_death_0, cutoff_date),
    end_date    = pmin(event_date, censor_date, na.rm = TRUE),
    status      = if_else(!is.na(event_date) & end_date == event_date, 1, 0),
    time        = as.numeric(difftime(end_date, baseline, units = "days")) / 365.25
  ) %>%
  filter(!is.na(time) & time > 0 & !is.na(baseline))

write.csv(df_surv, "stroke_cox_data.csv", row.names = FALSE)

# =============================================================
# 2. Prepare FW variables
# =============================================================

df <- read.csv("stroke_cox_data.csv")

df_surv <- df[complete.cases(df[, c(
  "time", "status",
  "FW_2", "FW_Projection", "FW_Commissural", "FW_Association",
  "Age_attending_2", "sex", "Townsend_index_0", "BIM_PreImaging_2",
  "Smoke_status_2", "Hypertension_2", "Diabetes_2",
  "brain_volume", "WMH_2"
)]), ]

fw_vars <- c("FW_2", "FW_Projection", "FW_Commissural", "FW_Association")

# Standardize FW variables
df_surv <- df_surv %>%
  mutate(across(all_of(fw_vars), ~ as.numeric(scale(.x))))

# Median split into high/low groups
for (col in fw_vars) {
  med <- median(df_surv[[col]], na.rm = TRUE)
  df_surv[[paste0(col, "_group")]] <- ifelse(df_surv[[col]] > med, "High", "Low")
}

group_cols <- paste0(fw_vars, "_group")

df_surv <- df_surv %>%
  mutate(across(all_of(group_cols), ~ factor(.x, levels = c("Low", "High"))))

# Event counts per group
lapply(group_cols, function(g) table(df_surv[[g]], df_surv$status))

# =============================================================
# 3. Cox regression: continuous FW
# =============================================================

cox_results <- data.frame()
failed_vars <- c()

for (var in fw_vars) {
  formula <- as.formula(paste(
    "Surv(time, status) ~", var,
    "+ Age_attending_2 + sex + Townsend_index_0 + BIM_PreImaging_2",
    "+ Smoke_status_2 + Hypertension_2 + Diabetes_2 + brain_volume + WMH_2"
  ))
  
  model <- tryCatch(
    coxph(formula, data = df_surv, control = coxph.control(iter.max = 30)),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(
      coxph(formula, data = df_surv, control = coxph.control(iter.max = 30))
    )
  )
  
  if (is.null(model)) {
    failed_vars <- c(failed_vars, paste0(var, ": coxph error"))
    next
  }
  
  if (any(is.na(coef(model))) || model$iter[1] == 30) {
    failed_vars <- c(failed_vars, paste0(var, ": did not converge or NA coefficients"))
    next
  }
  
  coef_summary <- summary(model)$coefficients
  var_rows <- grep(paste0("^", var, "$"), rownames(coef_summary))
  if (length(var_rows) == 0) {
    failed_vars <- c(failed_vars, paste0(var, ": no coefficient row found"))
    next
  }
  
  row   <- coef_summary[var_rows[1], ]
  hr    <- exp(row[1])
  ci    <- tryCatch(exp(confint(model)[var_rows[1], ]),
                    error = function(e) c(NA, NA))
  p_val <- row[5]
  n_evt <- model[["nevent"]]
  
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
    n_event         = n_evt,
    PH_assumption_p = ph_p
  ))
}

cat("Successful variables:", nrow(cox_results), "\n")
cat("Skipped/failed variables:", length(failed_vars), "\n")
if (length(failed_vars) > 0) print(failed_vars[1:min(10, length(failed_vars))])

cox_results <- cox_results %>%
  mutate(
    HR_95CI         = paste0(round(HR, 3), " (",
                             round(HR_CI_low, 3), "-",
                             round(HR_CI_high, 3), ")"),
    P_value         = round(P_value, 4),
    PH_assumption_p = round(PH_assumption_p, 4)
  )

print(cox_results)

write.csv(cox_results, "FW_Cox_stroke_model.csv", row.names = FALSE)

# =============================================================
# 4. KM curves for FW high/low groups (loop)
# =============================================================

output_dir <- "KM_plots"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

for (var in group_cols) {
  
  fit_km <- survfit(
    as.formula(paste("Surv(time, status) ~", var)),
    data = df_surv
  )
  
  p <- ggsurvplot(
    fit_km,
    data             = df_surv,
    pval             = TRUE,
    pval.method      = TRUE,
    pval.coord       = c(0, 0.998),
    conf.int         = TRUE,
    risk.table       = TRUE,
    ylim             = c(0.99, 1.0),
    risk.table.col   = "strata",
    risk.table.height = 0.25,
    risk.table.y.text = FALSE,
    cumcensor        = FALSE,
    cumevents        = FALSE,
    palette          = c("#2e9fdf", "#d37f88"),
    title            = "Kaplan-Meier Curve for Stroke",
    xlab             = "Time (years)",
    ylab             = "Progression-Free Probability",
    legend.title     = var,
    legend.labs      = c("FW low", "FW High"),
    break.time.by    = 2,
    ggtheme          = theme_classic()
  )
  
  pdf(file.path(output_dir, paste0(var, "_KM.pdf")), width = 8, height = 8)
  print(p)
  dev.off()
}