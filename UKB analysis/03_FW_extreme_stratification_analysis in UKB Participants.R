# =============================================================
# FW Extreme Stratification Analysis
# Workflow:
# 1. Extreme-stratification association analysis
#    - Stroke: top 25% vs bottom 25%
#    - Dementia: top 25% vs bottom 75%
# 2. Adjusted for imaging covariates (brain volume, WMH)
# 3. Adjusted for imaging + vascular risk factors
# 4. Adjusted for imaging + vascular + APOE
# 5. Forest plots (Dementia and Stroke)
# =============================================================

library(tidyverse)
library(haven)
library(broom)
library(openxlsx)

# =============================================================
# Helper: run extreme-stratification logistic regression
# =============================================================
run_extreme_strat <- function(data, fw_vars, diseases, covariates) {
  
  results <- tibble(
    FW_Var       = character(),
    Disease      = character(),
    OR           = numeric(),
    CI_lower     = numeric(),
    CI_upper     = numeric(),
    p_value      = numeric(),
    High_N       = integer(),
    Low_N        = integer(),
    Total_N      = integer(),
    High_Disease = integer(),
    Low_Disease  = integer()
  )
  
  for (fw in fw_vars) {
    
    if (!fw %in% names(data)) next
    
    for (d in diseases) {
      
      if (!d %in% names(data)) next
      
      # ---- Stratification strategy ----
      if (d == "Stroke") {
        p25 <- quantile(data[[fw]], probs = 0.25, na.rm = TRUE)
        p75 <- quantile(data[[fw]], probs = 0.75, na.rm = TRUE)
        
        data_temp <- data %>%
          mutate(group = case_when(
            .data[[fw]] >= p75 ~ "High",
            .data[[fw]] <= p25 ~ "Low",
            TRUE ~ NA_character_
          )) %>%
          filter(!is.na(group)) %>%
          mutate(group = factor(group, levels = c("Low", "High")))
        
      } else if (d == "Dementia") {
        p75 <- quantile(data[[fw]], probs = 0.75, na.rm = TRUE)
        
        data_temp <- data %>%
          mutate(group = case_when(
            .data[[fw]] >= p75 ~ "High",
            .data[[fw]] <  p75 ~ "Low",
            TRUE ~ NA_character_
          )) %>%
          filter(!is.na(group)) %>%
          mutate(group = factor(group, levels = c("Low", "High")))
      }
      
      if (nrow(data_temp) < 10) next
      
      # ---- Logistic regression ----
      tryCatch({
        formula <- as.formula(paste(d, "~ group +", paste(covariates, collapse = " + ")))
        model   <- glm(formula, data = data_temp,
                       family = binomial(link = "logit"))
        
        model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE)
        
        if (nrow(model_summary) >= 2) {
          res <- model_summary[2, ]
          
          high_group <- data_temp %>% filter(group == "High")
          low_group  <- data_temp %>% filter(group == "Low")
          
          results <- results %>%
            add_row(
              FW_Var       = fw,
              Disease      = d,
              OR           = res$estimate,
              CI_lower     = res$conf.low,
              CI_upper     = res$conf.high,
              p_value      = res$p.value,
              High_N       = nrow(high_group),
              Low_N        = nrow(low_group),
              Total_N      = nrow(high_group) + nrow(low_group),
              High_Disease = sum(high_group[[d]] == 1, na.rm = TRUE),
              Low_Disease  = sum(low_group[[d]]  == 1, na.rm = TRUE)
            )
        }
      }, error = function(e) NULL)
    }
  }
  
  results
}

# =============================================================
# Load data
# =============================================================
data <- read_dta("FW_disease_covariate.dta")

fw_vars  <- c("FW", "FW_Association", "FW_Commissural", "FW_Projection")
diseases <- c("Stroke", "Dementia")

# =============================================================
# Model 1: Baseline (age, sex, Townsend index)
# =============================================================
results_model1 <- run_extreme_strat(
  data       = data,
  fw_vars    = fw_vars,
  diseases   = diseases,
  covariates = c("sex", "Age", "Townsend_index")
)

# =============================================================
# Model 2: + imaging covariates (brain volume, WMH)
# =============================================================
results_model2 <- run_extreme_strat(
  data       = data,
  fw_vars    = fw_vars,
  diseases   = diseases,
  covariates = c("sex", "Age", "Townsend_index", "brain_volume", "WMH")
)

# =============================================================
# Model 3: + vascular risk factors
# =============================================================
results_model3 <- run_extreme_strat(
  data       = data,
  fw_vars    = fw_vars,
  diseases   = diseases,
  covariates = c("sex", "Age", "Townsend_index",
                 "brain_volume", "WMH",
                 "Hypertension", "Diabetes", "BMI", "Smoke_status")
)

# =============================================================
# Model 4: + APOE
# =============================================================
results_model4 <- run_extreme_strat(
  data       = data,
  fw_vars    = fw_vars,
  diseases   = diseases,
  covariates = c("sex", "Age", "Townsend_index",
                 "brain_volume", "WMH",
                 "Hypertension", "Diabetes", "BMI", "Smoke_status",
                 "APOE4_group")
)

print(results_model1)
print(results_model2)
print(results_model3)
print(results_model4)

# =============================================================
# Forest plot helper
# =============================================================
plot_forest <- function(data, color_map, shape_map,
                        ylim = c(0.5, 6),
                        out_file = NULL,
                        width = 16, height = 8) {
  
  data_clean <- data %>%
    select(FW_Var, Disease, OR, CI_lower, CI_upper) %>%
    rename(Varnames = FW_Var, Factor = Disease,
           OR = OR, Lower = CI_lower, Upper = CI_upper)
  
  data_clean <- data_clean %>%
    mutate(
      lower_trunc = ifelse(Lower < ylim[1], ylim[1], Lower),
      upper_trunc = ifelse(Upper > ylim[2], ylim[2], Upper),
      left_arrow  = Lower < ylim[1],
      right_arrow = Upper > ylim[2]
    )
  
  p <- ggplot(data_clean,
              aes(x = Varnames, y = OR, color = Factor, shape = Factor)) +
    geom_pointrange(aes(ymin = lower_trunc, ymax = upper_trunc),
                    position = position_dodge(width = 0.7),
                    size = 0.5, linewidth = 0.7) +
    geom_segment(
      data = filter(data_clean, left_arrow),
      aes(y = lower_trunc, yend = Lower,
          x = Varnames, xend = Varnames),
      arrow = arrow(type = "open",
                    length = unit(0.1, "inches"), ends = "first"),
      color = "black"
    ) +
    geom_point(position = position_dodge(width = 0.7),
               size = 2.5, stroke = 0.9) +
    scale_color_manual(values = color_map) +
    scale_shape_manual(values = shape_map) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
    labs(y = "OR (95% CI)", x = "") +
    scale_y_continuous(limits = ylim, breaks = seq(1, ylim[2], 1)) +
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
  
  print(p)
  
  if (!is.null(out_file)) {
    ggsave(out_file, plot = p,
           width = width, height = height, units = "cm", dpi = 300)
  }
}

# =============================================================
# Forest plot: Dementia
# =============================================================
data_dem <- read_xlsx("FW_Disease_Analysis25%_0815.xlsx", sheet = "Sheet1")[c(1, 3, 5, 7), ]
data_dem$FW_Var <- factor(data_dem$FW_Var,
                          levels = c("P75 of FW_Projection",
                                     "P75 of FW_Commissural",
                                     "P75 of FW_Association",
                                     "P75 of FW_2"))

plot_forest(
  data      = data_dem,
  color_map = c("Dementia" = "#A44A46", "Stroke" = "#9467bd"),
  shape_map = c("Dementia" = 16,        "Stroke" = 17),
  out_file  = "FW_dementia_OR.pdf"
)

# =============================================================
# Forest plot: Stroke
# =============================================================
data_str <- read_xlsx("FW_Disease_Analysis25%_0815.xlsx", sheet = "Sheet 1")[c(1, 3, 5, 7), ]
data_str$FW_Var <- factor(data_str$FW_Var,
                          levels = c("P75 of FW_Projection",
                                     "P75 of FW_Commissural",
                                     "P75 of FW_Association",
                                     "P75 of FW_2"))

plot_forest(
  data      = data_str,
  color_map = c("Stroke" = "#ffbb78"),
  shape_map = c("Stroke" = 17),
  out_file  = "FW_stroke_OR.pdf"
)