# FW-WM and Cognition in Hypertension
# This script compares FW-WM indexes between hypertensive individuals with and
# without domain-specific cognitive impairment, and tests whether higher FW-WM
# is associated with increased odds of impairment in each cognitive domain
# Workflow:
# 1. Group comparisons of FW-WM by cognitive impairment status (t-tests)
# 2. Boxplots of FW-WM by cognitive impairment status (one per domain)
# 3. Logistic regression of cognitive impairment on dichotomized FW-WM groups
# 4. Forest plot of the adjusted odds ratios

#=============================================================================
# 1. Group comparisons of FW-WM by cognitive impairment status
#=============================================================================

library(tidyverse)
library(haven)
library(broom)

# Read the data
data <- read_dta("HBP_cognitive_FW.dta")

# Define the cognitive variables and FW indexes
cognitive_vars <- c("memory", "visuospatial", "executive", "language")
fw_vars <- c("FW", "FW_Projection", "FW_Commissural", "FW_Association")

# Create an empty result data frame
results <- tibble(
  Cognitive_Var = character(),
  Group         = character(),
  FW_Var        = character(),
  Mean          = numeric(),
  SD            = numeric(),
  N             = integer(),
  t_stat        = numeric(),
  p_value       = numeric(),
  CI_lower      = numeric(),
  CI_upper      = numeric()
)

# Loop over each cognitive variable
for (cog_var in cognitive_vars) {
  
  if (!cog_var %in% names(data)) next
  
  # Subset samples with values 0 or 1 for the current cognitive variable
  sub_data <- data %>%
    filter(!is.na(.data[[cog_var]]),
           .data[[cog_var]] %in% c(0, 1)) %>%
    mutate(group = factor(.data[[cog_var]], levels = c(0, 1),
                          labels = c("Group0", "Group1")))
  
  if (nrow(sub_data) < 10) next
  
  # Loop over each FW index
  for (fw in fw_vars) {
    
    if (!fw %in% names(sub_data)) next
    
    # Descriptive statistics by group
    group_stats <- sub_data %>%
      group_by(group) %>%
      summarise(
        Mean = mean(.data[[fw]], na.rm = TRUE),
        SD   = sd(.data[[fw]],   na.rm = TRUE),
        N    = n(),
        .groups = "drop"
      )
    
    # t-test
    t_test   <- t.test(as.formula(paste(fw, "~ group")),
                       data = sub_data, var.equal = TRUE)
    t_result <- broom::tidy(t_test)
    
    for (grp in c("Group0", "Group1")) {
      grp_stats <- group_stats %>% filter(group == grp)
      
      results <- results %>%
        add_row(
          Cognitive_Var = cog_var,
          Group         = grp,
          FW_Var        = fw,
          Mean          = grp_stats$Mean,
          SD            = grp_stats$SD,
          N             = grp_stats$N,
          t_stat        = ifelse(grp == "Group0", NA_real_, t_result$statistic),
          p_value       = ifelse(grp == "Group0", NA_real_, t_result$p.value),
          CI_lower      = ifelse(grp == "Group0", NA_real_, t_result$conf.low),
          CI_upper      = ifelse(grp == "Group0", NA_real_, t_result$conf.high)
        )
    }
  }
}

print(results, n = Inf)
write_csv(results, "HBP_cognitive_fw_t_test_results.csv")

#=============================================================================
# 2. Boxplots of FW-WM by cognitive impairment status (loop)
#=============================================================================

library(readstata13)
library(dplyr)
library(ggplot2)

data <- read.dta13("HBP_cognitive_FW.dta")
data <- as.data.frame(data)

# Domain-specific subsets
memory_fw       <- data[!is.na(data$memory_2),       c(1, 3, 13:16)]
visuospatial_fw <- data[!is.na(data$visuospatial_2), c(1, 4, 13:16)]
executive_fw    <- data[!is.na(data$executive_2),    c(1, 5, 13:16)]
language_fw     <- data[!is.na(data$language_2),     c(1, 6, 13:16)]

# Custom color scheme
custom_colors <- c(
  "1" = "#B2C4D8", "2" = "#6A8CB2",   # FW_total
  "3" = "#B8D1A3", "4" = "#7CAA56",   # FW_Association
  "5" = "#A9D3DF", "6" = "#449EB6",   # FW_Commissural
  "7" = "#D9A7AB", "8" = "#B65259"    # FW_Projection
)

# Domain configuration: name of the data frame, cognitive column, FW columns
domain_specs <- list(
  list(df = memory_fw,       cog_col = "memory_2",       out = "HBP_memory_FW.pdf"),
  list(df = executive_fw,    cog_col = "executive_2",    out = "HBP_executive_FW.pdf"),
  list(df = visuospatial_fw, cog_col = "visuospatial_2", out = "HBP_visuospatial_FW.pdf"),
  list(df = language_fw,     cog_col = "language_2",     out = "HBP_language_FW.pdf")
)

fw_columns <- c("FW_2", "FW_Association", "FW_Commissural", "FW_Projection")

for (spec in domain_specs) {
  
  df      <- spec$df
  cog_col <- spec$cog_col
  
  df[[cog_col]] <- as.factor(df[[cog_col]])
  
  result_clean <- data.frame(
    FW_value = c(df[[fw_columns[1]]], df[[fw_columns[2]]],
                 df[[fw_columns[3]]], df[[fw_columns[4]]]),
    FW       = factor(rep(df[[cog_col]], 4)),
    Variable = c(rep(fw_columns[1], nrow(df)),
                 rep(fw_columns[2], nrow(df)),
                 rep(fw_columns[3], nrow(df)),
                 rep(fw_columns[4], nrow(df)))
  ) %>%
    na.omit() %>%
    mutate(
      color = case_when(
        Variable == fw_columns[1] & FW == 0 ~ 1,
        Variable == fw_columns[1] & FW == 1 ~ 2,
        Variable == fw_columns[2] & FW == 0 ~ 3,
        Variable == fw_columns[2] & FW == 1 ~ 4,
        Variable == fw_columns[3] & FW == 0 ~ 5,
        Variable == fw_columns[3] & FW == 1 ~ 6,
        Variable == fw_columns[4] & FW == 0 ~ 7,
        Variable == fw_columns[4] & FW == 1 ~ 8,
        TRUE ~ NA_integer_
      )
    )
  
  result_clean$color <- as.factor(result_clean$color)
  
  p <- ggplot(result_clean, aes(x = Variable, y = FW_value, fill = color)) +
    geom_boxplot(alpha = 0.8, outlier.size = -1,
                 width = 0.60, color = "#525252") +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 2,
                 color = "red", position = position_dodge(width = 0.6)) +
    scale_fill_manual(values = custom_colors, name = "color", guide = "none") +
    scale_x_discrete(limits = fw_columns) +
    ylim(0.2, 0.4) +
    labs(y = "FW value") +
    theme_classic() +
    theme(
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      panel.background  = element_blank(),
      axis.line         = element_line(color = "black", linewidth = 0.7),
      axis.ticks.length = unit(0.15, "cm"),
      axis.ticks        = element_line(linewidth = 0.5),
      axis.title        = element_blank(),
      legend.position   = "none",
      legend.title      = element_blank(),
      legend.text       = element_blank(),
      plot.margin       = margin(10, 10, 10, 10)
    )
  
  print(p)
  
  ggsave(spec$out, plot = p, width = 10, height = 12, units = "cm", dpi = 300)
}

#=============================================================================
# 3. Logistic regression of cognitive impairment on dichotomized FW-WM groups
#=============================================================================

library(tidyverse)
library(broom)
library(haven)

data <- read_dta("HBP_cognitive_FW.dta")

cognitive_vars <- c("memory", "visuospatial", "executive", "language")
fw_vars        <- c("FW", "FW_Projection", "FW_Commissural", "FW_Association")
covariates     <- c("Age", "sex", "Mean_SBP")

results <- tibble()

for (fw in fw_vars) {
  
  if (!fw %in% names(data)) next
  
  fw_group_name <- paste0(fw, "_group")
  median_val    <- median(data[[fw]], na.rm = TRUE)
  
  data[[fw_group_name]] <- factor(
    ifelse(data[[fw]] <= median_val, "Low", "High"),
    levels = c("Low", "High")
  )
  
  for (cog_var in cognitive_vars) {
    
    if (!cog_var %in% names(data)) next
    
    sub_data <- data %>%
      select(all_of(c(cog_var, fw_group_name, covariates))) %>%
      na.omit() %>%
      filter(!is.na(.data[[cog_var]]),
             .data[[cog_var]] %in% c(0, 1))
    
    if (nrow(sub_data) < 50) next
    
    prop_table <- sub_data %>%
      group_by(.data[[fw_group_name]]) %>%
      summarise(
        N_impairment = sum(.data[[cog_var]]),
        N_total      = n(),
        Proportion   = mean(.data[[cog_var]]),
        .groups      = "drop"
      )
    
    low_prop  <- prop_table %>% filter(.data[[fw_group_name]] == "Low")
    high_prop <- prop_table %>% filter(.data[[fw_group_name]] == "High")
    
    formula_str <- paste(cog_var, "~", fw_group_name, "+",
                         paste(covariates, collapse = " + "))
    
    model     <- glm(as.formula(formula_str), data = sub_data,
                     family = binomial())
    model_sum <- broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
    
    high_group_term <- paste0(fw_group_name, "High")
    or_row <- model_sum %>% filter(term == high_group_term)
    
    result_row <- tibble(
      Cognitive_Var        = cog_var,
      FW_Var               = fw,
      Group_Comparison     = "High vs Low",
      N_total              = nrow(sub_data),
      N_low                = low_prop$N_total,
      N_high               = high_prop$N_total,
      Prop_low_impairment  = low_prop$Proportion,
      Prop_high_impairment = high_prop$Proportion,
      OR                   = ifelse(nrow(or_row) > 0, or_row$estimate, NA),
      OR_CI_lower          = ifelse(nrow(or_row) > 0, or_row$conf.low, NA),
      OR_CI_upper          = ifelse(nrow(or_row) > 0, or_row$conf.high, NA),
      p_value              = ifelse(nrow(or_row) > 0, or_row$p.value, NA)
    )
    
    results <- bind_rows(results, result_row)
  }
}

if (nrow(results) > 0) {
  print(results %>% select(Cognitive_Var, FW_Var, OR, p_value))
  write_csv(results, "HBP_cognitive_FW_group_comparison_results.csv")
}

#=============================================================================
# 4. Forest plot of the adjusted odds ratios
#=============================================================================

library(readxl)
library(ggplot2)
library(dplyr)

data1 <- read.csv("HBP_cognitive_FW_group_comparison_results.csv")
data  <- data1[c(1:16), ]

data <- data %>%
  mutate(
    Cognitive_Var = factor(Cognitive_Var,
                           levels = c("visuospatial_2", "executive_2",
                                      "memory_2", "language_2")),
    FW_Var        = factor(FW_Var,
                           levels = c("FW_Projection", "FW_Commissural",
                                      "FW_Association", "FW_2"))
  )

data_clean <- data %>%
  select(FW_Var, Cognitive_Var, OR, OR_CI_lower, OR_CI_upper) %>%
  rename(Varnames = FW_Var, Factor = Cognitive_Var,
         OR = OR, Lower = OR_CI_lower, Upper = OR_CI_upper)

p <- ggplot(data_clean, aes(x = Varnames, y = OR,
                            color = Factor, shape = Factor)) +
  geom_pointrange(aes(ymin = Lower, ymax = Upper),
                  position = position_dodge(width = 0.7),
                  size = 0.5, linewidth = 0.7) +
  geom_point(position = position_dodge(width = 0.7),
             size = 2.5, stroke = 0.9) +
  scale_color_manual(values = c(
    "memory_2"       = "#A44A46",
    "visuospatial_2" = "#446C74",
    "executive_2"    = "#2A437A",
    "language_2"     = "#AB779D"
  )) +
  scale_shape_manual(values = c(
    "memory_2"       = 16,
    "visuospatial_2" = 17,
    "executive_2"    = 15,
    "language_2"     = 18
  )) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  labs(y = "OR (95% CI)", x = "") +
  scale_y_continuous(limits = c(0.8, 1.6), breaks = seq(1, 1.6, 0.2)) +
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
ggsave("FW_cognitive_OR.pdf", plot = p,
       width = 16, height = 14, units = "cm", dpi = 300)