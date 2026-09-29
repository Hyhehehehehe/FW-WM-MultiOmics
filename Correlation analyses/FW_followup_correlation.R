# FW Follow-up Correlation
# This script evaluates the test-retest consistency between baseline and
# follow-up FW-WM measurements
# Workflow:
# 1. Calculate the Pearson correlation and ICC for each FW-WM index
# 2. Create correlation scatter plots with ICC values and panel titles
# 3. Save the combined figure in TIFF and PDF formats

# Load the required packages
library(haven)      # For reading .dta files
library(data.table) # Efficient data processing
library(ggplot2)    # Plotting
library(ggpubr)     # Plot arrangement

########## 1. Data preparation ##########
# Read the data
data_path <- "freewater.dta"
data_fw <- read_dta(data_path)
data_fw <- as.data.table(data_fw)

# Check the data
cat("Data dimensions:", dim(data_fw), "\n")

# Define the metric pairs and new label names (in the required order)
fw_metrics <- list(
  "FWtotal" = c("FW_2", "FW_3"),
  "FWassociation" = c("FW_Association", "FW_Association_3"),
  "FWcommissural" = c("FW_Commissural", "FW_Commissural_3"),
  "FWprojection" = c("FW_Projection", "FW_Projection_3")
)

# Define the display labels (in the new order)
display_labels <- list(
  "FWtotal" = list(x = "FWtotal baseline", y = "FWtotal follow-up", title = "FWtotal"),
  "FWassociation" = list(x = "FWassociation baseline", y = "FWassociation follow-up", title = "FWassociation"),
  "FWcommissural" = list(x = "FWcommissural baseline", y = "FWcommissural follow-up", title = "FWcommissural"),
  "FWprojection" = list(x = "FWprojection baseline", y = "FWprojection follow-up", title = "FWprojection")
)

########## 2. Calculate the correlation coefficients and ICC ##########
cat("\n===== Test-retest reliability analysis =====\n")

# Create the result table
results_table <- data.table(
  Metric = character(),
  n_pairs = integer(),
  Pearson_r = numeric(),
  Pearson_p = numeric(),
  ICC_estimate = numeric(),
  ICC_CI_lower = numeric(),
  ICC_CI_upper = numeric()
)

# Calculate the correlation and ICC for each metric (in the new order)
for (i in seq_along(fw_metrics)) {
  metric_name <- names(fw_metrics)[i]
  var_time1 <- fw_metrics[[i]][1]
  var_time2 <- fw_metrics[[i]][2]
  
  cat(sprintf("\nAnalyzing metric: %s (%s vs %s)\n", metric_name, var_time1, var_time2))
  
  # Extract the paired data
  paired_data <- data_fw[!is.na(get(var_time1)) & !is.na(get(var_time2)), 
                         .(time1 = get(var_time1), time2 = get(var_time2))]
  
  n_pairs <- nrow(paired_data)
  
  if (n_pairs > 1) {
    # Calculate the Pearson correlation coefficient
    cor_test <- cor.test(paired_data$time1, paired_data$time2, method = "pearson")
    pearson_r <- cor_test$estimate
    pearson_p <- cor_test$p.value
    
    # Calculate the ICC - using a simple calculation method
    # Convert the data to long format
    long_data <- rbind(
      data.table(id = 1:n_pairs, time = 1, value = paired_data$time1),
      data.table(id = 1:n_pairs, time = 2, value = paired_data$time2)
    )
    
    # Calculate the variance components
    anova_result <- aov(value ~ factor(id) + factor(time), data = long_data)
    anova_summary <- summary(anova_result)
    
    # Extract the mean squares
    MS_id <- anova_summary[[1]][1, "Mean Sq"]  # Between subjects
    MS_time <- anova_summary[[1]][2, "Mean Sq"] # Measurement occasion
    MS_error <- anova_summary[[1]][3, "Mean Sq"] # Error
    
    # Calculate ICC(3,1) - two-way mixed model, absolute agreement
    k <- 2  # Number of measurements
    n <- n_pairs  # Number of subjects
    
    # ICC(3,1) = (MS_id - MS_error) / (MS_id + (k-1)*MS_error + k*(MS_time - MS_error)/n)
    icc_value <- (MS_id - MS_error) / (MS_id + (k-1)*MS_error + k*(MS_time - MS_error)/n)
    
    # Calculate the confidence interval of the ICC (using the F distribution)
    F_value <- MS_id / MS_error
    df1 <- n - 1
    df2 <- (n - 1) * (k - 1)
    
    # Calculate the critical values of the F distribution
    alpha <- 0.05
    F_lower <- qf(alpha/2, df1, df2)
    F_upper <- qf(1 - alpha/2, df1, df2)
    
    # Calculate the confidence interval of the ICC
    icc_lower <- (F_value/F_upper - 1) / (F_value/F_upper + (k-1))
    icc_upper <- (F_value/F_lower - 1) / (F_value/F_lower + (k-1))
    
    # Ensure the ICC is within the [0, 1] range
    icc_value <- max(0, min(1, icc_value))
    icc_lower <- max(0, min(1, icc_lower))
    icc_upper <- max(0, min(1, icc_upper))
    
    # Store the results
    results_table <- rbind(results_table,
                           data.table(
                             Metric = metric_name,
                             n_pairs = n_pairs,
                             Pearson_r = round(pearson_r, 3),
                             Pearson_p = round(pearson_p, 5),
                             ICC_estimate = round(icc_value, 3),
                             ICC_CI_lower = round(icc_lower, 3),
                             ICC_CI_upper = round(icc_upper, 3)
                           ))
    
    cat(sprintf("  Sample size: %d\n", n_pairs))
    cat(sprintf("  Pearson r: %.3f (p = %.5f)\n", pearson_r, pearson_p))
    cat(sprintf("  ICC(3,1): %.3f (95%% CI: [%.3f, %.3f])\n", icc_value, icc_lower, icc_upper))
    
  } else {
    cat(sprintf("  Insufficient sample size: %d\n", n_pairs))
  }
}

########## 3. Create the correlation scatter plots (with ICC and panel titles) ##########
cat("\n===== Creating test-retest correlation scatter plots (with ICC and panel titles) =====\n")

# Scatter plot function - including the ICC value and panel title
create_scatter_plot_with_title_icc <- function(var_time1, var_time2, metric_name, results_row) {
  # Extract the paired data
  plot_data <- data_fw[!is.na(get(var_time1)) & !is.na(get(var_time2)), 
                       .(time1 = get(var_time1), time2 = get(var_time2))]
  
  # Calculate the axis ranges
  all_values <- c(plot_data$time1, plot_data$time2)
  min_val <- min(all_values, na.rm = TRUE)
  max_val <- max(all_values, na.rm = TRUE)
  range_val <- max_val - min_val
  
  # Get the display labels
  x_label <- display_labels[[metric_name]]$x
  y_label <- display_labels[[metric_name]]$y
  plot_title <- display_labels[[metric_name]]$title
  
  # Format the P value
  pearson_p_formatted <- ifelse(results_row$Pearson_p < 0.001, 
                                "p < 0.001", 
                                sprintf("p = %.3f", results_row$Pearson_p))
  
  # Format the ICC value (3 decimal places)
  icc_formatted <- sprintf("ICC = %.3f", results_row$ICC_estimate)
  
  # Create the scatter plot - add the panel title
  p <- ggplot(plot_data, aes(x = time1, y = time2)) +
    geom_point(alpha = 0.5, size = 1, color = "#4C75A0") +
    geom_smooth(method = "lm", se = TRUE, formula = y ~ x, 
                color = "#E64B35", fill = "#E64B35", alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", alpha = 0.7) +
    # Add the statistical information - including Pearson r, P value, and ICC
    annotate("text", 
             x = min_val + 0.05 * range_val,
             y = max_val - 0.05 * range_val,
             label = sprintf("r = %.3f\n%s\n%s\nn = %d", 
                             results_row$Pearson_r,
                             pearson_p_formatted,
                             icc_formatted,
                             nrow(plot_data)),
             hjust = 0, vjust = 1,
             size = 3.5,
             color = "black",
             fontface = "bold") +
    labs(
      title = plot_title,
      x = x_label,
      y = y_label
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold", margin = margin(b = 10)),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 10),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2)
    )
  
  return(p)
}

# Create a scatter plot for each metric (in the new order)
scatter_plots <- list()

for (i in seq_along(fw_metrics)) {
  metric_name <- names(fw_metrics)[i]
  var_time1 <- fw_metrics[[i]][1]
  var_time2 <- fw_metrics[[i]][2]
  
  # Get the corresponding result row
  results_row <- results_table[Metric == metric_name]
  
  if (nrow(results_row) > 0) {
    scatter_plots[[i]] <- create_scatter_plot_with_title_icc(var_time1, var_time2, metric_name, results_row)
  }
}

# Combine the scatter plots (without an overall title)
combined_scatter <- ggarrange(
  plotlist = scatter_plots,
  ncol = 2, 
  nrow = 2,
  labels = "AUTO"
)

# Note: no overall title is added here because the combined figure does not require one
print(combined_scatter)

########## 4. Save the figure ##########
output_dir <- "./followup_correlation/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Save the scatter plot
tiff(paste0(output_dir, "FW_test_retest_scatter_1.tiff"), 
     width = 12, height = 10, units = "in", res = 300, compression = "lzw")
print(combined_scatter)
dev.off()

# Save the result table (optional)
#write.csv(results_table, paste0(output_dir, "FW_reliability_results.csv"), row.names = FALSE)

