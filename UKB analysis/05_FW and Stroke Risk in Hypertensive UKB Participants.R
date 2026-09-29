# =============================================================
# FW and Stroke risk in Hypertensive UKB Participants
# 1. FW high/low groups vs stroke incidence (logistic regression + forest plot)
# 2. FW values by stroke status (boxplot + t-test)
# =============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)

# ---- Load analysis dataset ----
df <- read.csv("FW_stroke_in_hypertension_data.csv")


# =============================================================
# 1. FW high/low groups vs stroke incidence
# =============================================================

# ---- 1.1 Median split of FW metrics ----
st_cols <- grep("FW", names(df), value = TRUE)
st_cols <- st_cols[!grepl("_group$", st_cols)]

for (col in st_cols) {
  med <- median(df[[col]], na.rm = TRUE)
  df[[paste0(col, "_group")]] <- ifelse(df[[col]] > med, "High", "Low")
}

# ---- 1.2 Format variables ----
df$I64_2 <- factor(df$I64_2, levels = c(0, 1))
df$sex   <- factor(df$sex,   levels = c(0, 1))

group_cols <- c("FW_2_group", "FW_Projection_group",
                "FW_Commissural_group", "FW_Association_group")

# ---- 1.3 Helper: extract OR and 95% CI ----
extract_or_ci <- function(model, var_name) {
  coef_summary <- summary(model)$coefficients
  idx <- grep(paste0("^", var_name), rownames(coef_summary))
  if (length(idx) == 0) return(NULL)
  
  beta   <- coef_summary[idx, "Estimate"]
  se     <- coef_summary[idx, "Std. Error"]
  z_val  <- coef_summary[idx, "z value"]
  p_val  <- coef_summary[idx, "Pr(>|z|)"]
  
  data.frame(
    Variable = var_name,
    OR       = round(exp(beta), 3),
    CI_lower = round(exp(beta - 1.96 * se), 3),
    CI_upper = round(exp(beta + 1.96 * se), 3),
    P_value  = round(p_val, 4),
    Z_value  = round(z_val, 3)
  )
}

# ---- 1.4 Fit logistic models ----
results_df <- data.frame()

for (col in group_cols) {
  df[[col]] <- factor(df[[col]], levels = c("Low", "High"))
  formula <- as.formula(paste("I64_2 ~", col,
                              "+ sex + Age_attending_2 + Mean_SBP_2"))
  model <- glm(formula, data = df, family = binomial(link = "logit"))
  
  temp <- extract_or_ci(model, col)
  if (!is.null(temp)) results_df <- rbind(results_df, temp)
}

print(results_df)

write.csv(results_df,
          "FW_to_stroke_OR_results.csv",
          row.names = FALSE)


# =============================================================
# 2. Forest plot: FW → stroke OR (UKB)
# =============================================================

data <- read.csv("FW_to_stroke_OR_results.csv")

custom_levels <- c("FW_2_group",
                   "FW_Association_group",
                   "FW_Commissural_group",
                   "FW_Projection_group")

data <- data %>%
  mutate(Variable = factor(Variable, levels = custom_levels))

data_clean <- data %>%
  select(Variable, OR, CI_lower, CI_upper) %>%
  rename(Varnames = Variable, Lower = CI_lower, Upper = CI_upper)

p_forest <- ggplot(data_clean,
                   aes(x = Varnames, y = OR, color = Varnames, shape = Varnames)) +
  geom_pointrange(aes(ymin = Lower, ymax = Upper),
                  position = position_dodge(width = 0.7),
                  size = 0.5, linewidth = 0.7) +
  geom_point(position = position_dodge(width = 0.7),
             size = 2.5, stroke = 0.9) +
  scale_color_manual(values = c(
    "FW_2_group"           = "#A44A46",
    "FW_Association_group" = "#446C74",
    "FW_Commissural_group" = "#2A437A",
    "FW_Projection_group"  = "#AB779D"
  )) +
  scale_shape_manual(values = c(
    "FW_2_group"           = 16,
    "FW_Association_group" = 17,
    "FW_Commissural_group" = 15,
    "FW_Projection_group"  = 18
  )) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  labs(y = "OR (95% CI)", x = "") +
  scale_y_continuous(limits = c(1, 2.5), breaks = seq(1, 2.5, 0.5)) +
  coord_flip() +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line        = element_line(color = "black", linewidth = 0.7),
    axis.ticks.length = unit(0.15, "cm"),
    axis.ticks       = element_line(linewidth = 0.5),
    legend.position  = "right",
    legend.title     = element_blank(),
    plot.margin      = margin(10, 10, 10, 10)
  )

print(p_forest)

ggsave("FW_stroke_OR.pdf",
       plot = p_forest, width = 16, height = 14, units = "cm", dpi = 300)


# =============================================================
# 3. FW values by stroke status in hypertensive patients
# =============================================================

df_sub <- read.csv("FW_stroke_in_hypertension_data.csv")

fw_vars <- grep("FW", names(df_sub), value = TRUE)
fw_vars <- fw_vars[!grepl("_group$", fw_vars)]

output_dir <- "boxplots"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

df_sub$I64_2 <- factor(df_sub$I64_2, levels = c(0, 1),
                       labels = c("No_Stroke", "Stroke"))

print(table(df_sub$I64_2))

# ---- 3.1 Reshape to long format ----
temp_data <- df_sub[, c("I64_2", fw_vars)]
temp_data <- temp_data[complete.cases(temp_data[, fw_vars]), ]

result_clean <- data.frame(
  FW_value = c(temp_data[[fw_vars[1]]], temp_data[[fw_vars[2]]],
               temp_data[[fw_vars[3]]], temp_data[[fw_vars[4]]]),
  Group    = rep(temp_data$I64_2, 4),
  Variable = factor(c(
    rep(fw_vars[1], nrow(temp_data)),
    rep(fw_vars[2], nrow(temp_data)),
    rep(fw_vars[3], nrow(temp_data)),
    rep(fw_vars[4], nrow(temp_data))
  ), levels = fw_vars)
) %>% na.omit()

# ---- 3.2 Assign color groups ----
result_clean <- result_clean %>%
  mutate(
    color_group = case_when(
      Variable == fw_vars[1] & Group == "No_Stroke" ~ 1,
      Variable == fw_vars[1] & Group == "Stroke"    ~ 2,
      Variable == fw_vars[2] & Group == "No_Stroke" ~ 3,
      Variable == fw_vars[2] & Group == "Stroke"    ~ 4,
      Variable == fw_vars[3] & Group == "No_Stroke" ~ 5,
      Variable == fw_vars[3] & Group == "Stroke"    ~ 6,
      Variable == fw_vars[4] & Group == "No_Stroke" ~ 7,
      Variable == fw_vars[4] & Group == "Stroke"    ~ 8,
      TRUE ~ NA_integer_
    )
  )

result_clean$color_group <- as.factor(result_clean$color_group)

custom_colors <- c(
  "1" = "#B2C4D8", "2" = "#6A8CB2",   # FW_total
  "3" = "#B8D1A3", "4" = "#7CAA56",   # FW_association
  "5" = "#A9D3DF", "6" = "#449EB6",   # FW_commissural
  "7" = "#D9A7AB", "8" = "#B65259"    # FW_projection
)

# ---- 3.3 Boxplot ----
p_box <- ggplot(result_clean, aes(x = Variable, y = FW_value, fill = color_group)) +
  geom_boxplot(alpha = 0.8, outlier.size = -1, width = 0.60, color = "#525252") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3,
               color = "red", position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = custom_colors, guide = "none") +
  scale_x_discrete(
    limits = fw_vars,
    labels = c("Total", "Association", "Commissural", "Projection")
  ) +
  labs(y = "FW value",
       title = "FW Values by Stroke Status in Hypertensive Patients") +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line        = element_line(color = "black", linewidth = 0.7),
    axis.ticks.length = unit(0.15, "cm"),
    axis.ticks       = element_line(linewidth = 0.5),
    axis.text.x      = element_text(size = 10, color = "black"),
    axis.text.y      = element_text(size = 10, color = "black"),
    axis.title.x     = element_blank(),
    axis.title.y     = element_text(size = 12, color = "black", face = "bold"),
    plot.title       = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position  = "none",
    plot.margin      = margin(10, 10, 10, 10)
  )

print(p_box)

ggsave(file.path(output_dir, "Hypertension_Stroke_FW_Boxplot.pdf"),
       plot = p_box, width = 10, height = 12, units = "cm", dpi = 300)

ggsave(file.path(output_dir, "Hypertension_Stroke_FW_Boxplot.png"),
       plot = p_box, width = 10, height = 12, units = "cm", dpi = 300)

# ---- 3.4 t-test: FW by stroke status ----
t_test_results <- data.frame()

for (fw in fw_vars) {
  group0 <- result_clean[result_clean$Variable == fw &
                           result_clean$Group == "No_Stroke", "FW_value"]
  group1 <- result_clean[result_clean$Variable == fw &
                           result_clean$Group == "Stroke",    "FW_value"]
  
  t_res <- t.test(group0, group1)
  
  t_test_results <- rbind(t_test_results, data.frame(
    FW_Index       = fw,
    No_Stroke_Mean = round(mean(group0, na.rm = TRUE), 4),
    No_Stroke_SD   = round(sd(group0,   na.rm = TRUE), 4),
    No_Stroke_N    = length(group0),
    Stroke_Mean    = round(mean(group1, na.rm = TRUE), 4),
    Stroke_SD      = round(sd(group1,   na.rm = TRUE), 4),
    Stroke_N       = length(group1),
    t_statistic    = round(t_res$statistic, 4),
    df             = round(t_res$parameter, 2),
    P_value        = format(t_res$p.value, scientific = FALSE, digits = 4)
  ))
}

print(t_test_results)

write.csv(t_test_results,
          file.path(output_dir, "Hypertension_Stroke_FW_t_test_results.csv"),
          row.names = FALSE)