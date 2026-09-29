############################################################
## Mediation analysis of hypertension-related proteins linking
## hypertension to FW-WM indexes
##
## Indirect effects (Beta1 * Beta2) were tested using the
## bootstrap approach, with bootstrap samples drawn from the
## sampling distributions of the two MR effect estimates.
############################################################

library(readxl)

# Bootstrap test for the indirect (mediation) effect
bootstrap_test <- function(b1, se1, b2, se2, n_boot = 100000) {
  boots_b1 <- rnorm(n_boot, mean = b1, sd = se1)
  boots_b2 <- rnorm(n_boot, mean = b2, sd = se2)
  boots_indirect <- boots_b1 * boots_b2

  ci <- quantile(boots_indirect, probs = c(0.025, 0.975), names = FALSE)
  p_value <- 2 * min(
    mean(boots_indirect <= 0),
    mean(boots_indirect >= 0)
  )
  return(list(
    effect = b1 * b2,
    ci_lower = ci[1],
    ci_upper = ci[2],
    p = p_value
  ))
}

# Read the input data
# Required columns: X (exposure), M (mediator), Y (outcome),
# Beta1/Se1 (X -> M), Beta2/Se2 (M -> Y), Beta3 (total effect, X -> Y)
file_path <- "mediation_HBP_PRO_FW.xlsx"
data <- read_excel(file_path, sheet = "Sheet1")

# Initialize the result data frame
final_results <- data.frame()

# Analyze each candidate mediator row by row
for (i in 1:nrow(data)) {
  row <- data[i, ]

  # Extract effect estimates and standard errors
  b1 <- as.numeric(row$Beta1)
  se1 <- as.numeric(row$Se1)
  b2 <- as.numeric(row$Beta2)
  se2 <- as.numeric(row$Se2)
  total_effect <- as.numeric(row$Beta3)  # total effect

  # Bootstrap test of the indirect effect
  bootstrap <- bootstrap_test(b1, se1, b2, se2, 100000)

  # Mediation proportion (based on the point estimate)
  mediation_ratio_bootstrap <- bootstrap$effect / total_effect

  # Direct effect (total effect - indirect effect)
  direct_effect_bootstrap <- total_effect - bootstrap$effect

  # Organize the result
  result_row <- data.frame(
    RowID = i,
    X = row$X,
    M = row$M,
    Y = row$Y,
    Total_Effect = total_effect,
    Method = "Bootstrap",
    Indirect_Effect = bootstrap$effect,
    CI_Lower = bootstrap$ci_lower,
    CI_Upper = bootstrap$ci_upper,
    P_Value = bootstrap$p,
    Direct_Effect = direct_effect_bootstrap,
    Mediation_Ratio = mediation_ratio_bootstrap
  )

  final_results <- rbind(final_results, result_row)
}

# Save the results
write.csv(final_results, "mediation_HBP_PRO_FW_analysis_results.csv", row.names = FALSE)
