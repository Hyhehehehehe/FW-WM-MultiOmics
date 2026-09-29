############################################################
## Consistency analysis of FW values extracted from single-shell and multi-shell data
## Including:
## 1. Pearson correlation
## 2. ICC(C,1): two-way model, consistency, single measurement
## 3. Scatter plots (with regression line, identity line, and statistical annotations)
############################################################

required_packages <- c(
  "haven",
  "data.table",
  "ggplot2",
  "ggpubr",
  "irr"
)

new_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(new_packages) > 0) {
  install.packages(new_packages)
}

library(haven)
library(data.table)
library(ggplot2)
library(ggpubr)
library(irr)


########## 1. Parameter settings and data loading ##########

data_path <- "FW_ADNI_MB_SB.dta"
output_dir <- "./"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

data_fw <- read_dta(data_path)
data_fw <- as.data.table(data_fw)

cat("Data dimensions:", nrow(data_fw), "rows x", ncol(data_fw), "columns\n")


########## 2. Define the metric pairs ##########

fw_metrics <- list(
  FWtotal = c(
    single_shell = "FW_SB",
    multi_shell = "FW"
  ),
  FWassociation = c(
    single_shell = "FW_Association_SB",
    multi_shell = "FW_Association"
  ),
  FWcommissural = c(
    single_shell = "FW_Commissural_SB",
    multi_shell = "FW_Commissural"
  ),
  FWprojection = c(
    single_shell = "FW_Projection_SB",
    multi_shell = "FW_Projection"
  )
)

display_labels <- list(
  FWtotal = list(
    x = "FWTotal, Single_shell",
    y = "FWTotal, Multi_shell",
    title = "Total"
  ),
  FWassociation = list(
    x = "FWAssociation, Single_shell",
    y = "FWAssociation, Multi_shell",
    title = "Association"
  ),
  FWcommissural = list(
    x = "FWCommissural, Single_shell",
    y = "FWCommissural, Multi_shell",
    title = "Commissural"
  ),
  FWprojection = list(
    x = "FWProjection, Single_shell",
    y = "FWProjection, Multi_shell",
    title = "Projection"
  )
)


########## 3. Check whether the variables exist ##########

all_variables <- unique(unlist(fw_metrics))
missing_variables <- setdiff(all_variables, names(data_fw))

if (length(missing_variables) > 0) {
  stop(
    paste0(
      "The following variables do not exist in the data:\n",
      paste(missing_variables, collapse = ", ")
    )
  )
}

cat("All analysis variables are present.\n")


########## 4. Helper functions ##########

format_p_value <- function(p_value) {
  if (is.na(p_value)) {
    return("P = NA")
  }
  if (p_value < 0.001) {
    return("P < 0.001")
  }
  sprintf("P = %.3f", p_value)
}

calculate_relative_icc <- function(paired_data) {
  if (nrow(paired_data) < 3) {
    return(
      list(
        estimate = NA_real_,
        lower = NA_real_,
        upper = NA_real_,
        p_value = NA_real_
      )
    )
  }
  if (
    sd(paired_data[[1]], na.rm = TRUE) == 0 ||
    sd(paired_data[[2]], na.rm = TRUE) == 0
  ) {
    return(
      list(
        estimate = NA_real_,
        lower = NA_real_,
        upper = NA_real_,
        p_value = NA_real_
      )
    )
  }
  icc_result <- tryCatch(
    {
      irr::icc(
        ratings = paired_data,
        model = "twoway",
        type = "consistency",
        unit = "single",
        conf.level = 0.95
      )
    },
    error = function(e) {
      message("ICC(C,1) calculation failed: ", e$message)
      return(NULL)
    }
  )
  if (is.null(icc_result)) {
    return(
      list(
        estimate = NA_real_,
        lower = NA_real_,
        upper = NA_real_,
        p_value = NA_real_
      )
    )
  }
  list(
    estimate = as.numeric(icc_result$value),
    lower = as.numeric(icc_result$lbound),
    upper = as.numeric(icc_result$ubound),
    p_value = as.numeric(icc_result$p.value)
  )
}


########## 5. Pearson correlation and consistency ICC analysis ##########

cat("\n")
cat("====================================================\n")
cat("Single-shell vs multi-shell FW: Pearson correlation and ICC(C,1) analysis\n")
cat("====================================================\n")

results_list <- list()

for (metric_name in names(fw_metrics)) {
  
  var_single <- fw_metrics[[metric_name]]["single_shell"]
  var_multi <- fw_metrics[[metric_name]]["multi_shell"]
  
  cat("\n----------------------------------------------------\n")
  cat("Analyzing metric: ", metric_name, "\n")
  cat("Single-shell variable: ", var_single, "\n")
  cat("Multi-shell variable: ", var_multi, "\n")
  
  paired_data <- data_fw[
    !is.na(get(var_single)) &
      !is.na(get(var_multi)),
    .(
      single_shell = as.numeric(get(var_single)),
      multi_shell = as.numeric(get(var_multi))
    )
  ]
  
  paired_data <- paired_data[
    is.finite(single_shell) &
      is.finite(multi_shell)
  ]
  
  n_pairs <- nrow(paired_data)
  
  if (n_pairs < 3) {
    warning(
      paste0(
        metric_name,
        " has insufficient valid paired samples, n = ",
        n_pairs
      )
    )
    next
  }
  
  pearson_result <- cor.test(
    paired_data$single_shell,
    paired_data$multi_shell,
    method = "pearson"
  )
  
  pearson_r <- unname(pearson_result$estimate)
  pearson_p <- pearson_result$p.value
  pearson_ci_lower <- pearson_result$conf.int[1]
  pearson_ci_upper <- pearson_result$conf.int[2]
  
  icc_consistency <- calculate_relative_icc(
    paired_data = as.data.frame(paired_data)
  )
  
  results_list[[metric_name]] <- data.table(
    Metric = metric_name,
    n_pairs = n_pairs,
    Pearson_r = pearson_r,
    Pearson_CI_lower = pearson_ci_lower,
    Pearson_CI_upper = pearson_ci_upper,
    Pearson_p = pearson_p,
    ICC_C1 = icc_consistency$estimate,
    ICC_C1_CI_lower = icc_consistency$lower,
    ICC_C1_CI_upper = icc_consistency$upper,
    ICC_C1_p = icc_consistency$p_value
  )
  
  cat("Valid paired sample size: ", n_pairs, "\n")
  cat(
    sprintf(
      "Pearson r：%.4f，95%% CI [%.4f, %.4f]，%s\n",
      pearson_r,
      pearson_ci_lower,
      pearson_ci_upper,
      format_p_value(pearson_p)
    )
  )
  cat(
    sprintf(
      paste0(
        "ICC(C,1), consistency: %.4f, ",
        "95%% CI [%.4f, %.4f]，P = %.3g\n"
      ),
      icc_consistency$estimate,
      icc_consistency$lower,
      icc_consistency$upper,
      icc_consistency$p_value
    )
  )
}


########## 6. Combine and save the statistical results ##########

if (length(results_list) == 0) {
  stop("No metric was analyzed successfully.")
}

results_table <- rbindlist(
  results_list,
  use.names = TRUE,
  fill = TRUE
)

results_file <- file.path(
  output_dir,
  "FW_Single_shell_vs_Multi_shell_Pearson_ICC_C1_results.csv"
)

fwrite(
  results_table,
  results_file
)

cat("\nStatistical results have been saved to:\n")
cat(results_file, "\n")

results_display <- copy(results_table)
numeric_columns <- names(results_display)[
  vapply(results_display, is.numeric, logical(1))
]
numeric_columns <- setdiff(
  numeric_columns,
  "n_pairs"
)
results_display[
  ,
  (numeric_columns) := lapply(
    .SD,
    function(x) round(x, 5)
  ),
  .SDcols = numeric_columns
]
print(results_display)


########## 7. Calculate the global axis range ##########

all_x <- c()
all_y <- c()

for (metric_name in names(fw_metrics)) {
  var_single <- fw_metrics[[metric_name]]["single_shell"]
  var_multi <- fw_metrics[[metric_name]]["multi_shell"]
  
  plot_data <- data_fw[
    !is.na(get(var_single)) &
      !is.na(get(var_multi)),
    .(
      single_shell = as.numeric(get(var_single)),
      multi_shell = as.numeric(get(var_multi))
    )
  ]
  plot_data <- plot_data[
    is.finite(single_shell) &
      is.finite(multi_shell)
  ]
  
  all_x <- c(all_x, plot_data$single_shell)
  all_y <- c(all_y, plot_data$multi_shell)
}

global_min <- min(c(all_x, all_y), na.rm = TRUE)
global_max <- max(c(all_x, all_y), na.rm = TRUE)
value_range <- global_max - global_min

if (!is.finite(value_range) || value_range == 0) {
  value_range <- max(abs(c(global_min, global_max)), 1) * 0.1
}

padding <- value_range * 0.05
global_lower <- global_min - padding
global_upper <- global_max + padding

cat(sprintf("Unified axis range: [%.4f, %.4f]\n", global_lower, global_upper))


########## 8. Create the scatter plot function (using unified axes) ##########

create_scatter_plot <- function(
    var_single,
    var_multi,
    metric_name,
    results_row,
    global_min,
    global_max) {
  
  plot_data <- data_fw[
    !is.na(get(var_single)) &
      !is.na(get(var_multi)),
    .(
      single_shell = as.numeric(get(var_single)),
      multi_shell = as.numeric(get(var_multi))
    )
  ]
  plot_data <- plot_data[
    is.finite(single_shell) &
      is.finite(multi_shell)
  ]
  
  annotation_label <- sprintf(
    paste0(
      "r = %.3f\n",
      "ICC = %.3f\n",
      "%s"
    ),
    results_row$Pearson_r,
    results_row$ICC_C1,
    format_p_value(results_row$Pearson_p)
  )
  
  # Set the x-axis and y-axis ranges separately
  x_lower <- 0.30
  x_upper <- global_max
  
  y_lower <- global_min
  y_upper <- 0.40
  
  # The annotation is placed in the upper-left corner of the current display range
  x_pos <- x_lower + 0.03 * (x_upper - x_lower)
  y_pos <- y_upper - 0.03 * (y_upper - y_lower)
  
  ggplot(
    plot_data,
    aes(
      x = single_shell,
      y = multi_shell
    )
  ) +
    geom_point(
      alpha = 0.35,
      size = 0.8,
      color = "#4C75A0"
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      color = "#E64B35",
      fill = "#E64B35",
      alpha = 0.18,
      linewidth = 0.8
    ) +
    annotate(
      geom = "text",
      x = x_pos,
      y = y_pos,
      label = annotation_label,
      hjust = 0,
      vjust = 1,
      size = 3.3,
      fontface = "bold"
    ) +
    scale_x_continuous(
      breaks = seq(0.30, 0.55, by = 0.05),
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      breaks = seq(0.15, 0.40, by = 0.05),
      expand = expansion(mult = 0)
    ) +
    coord_fixed(
      xlim = c(x_lower, x_upper),
      ylim = c(y_lower, y_upper),
      ratio = 1,
      expand = FALSE
    ) +
    labs(
      title = display_labels[[metric_name]]$title,
      x = display_labels[[metric_name]]$x,
      y = display_labels[[metric_name]]$y
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 14,
        face = "bold",
        margin = margin(b = 10)
      ),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}


########## 9. Generate the combined scatter plot ##########

scatter_plots <- list()

for (metric_name in names(fw_metrics)) {
  
  if (!metric_name %in% results_table$Metric) {
    next
  }
  
  var_single <- fw_metrics[[metric_name]]["single_shell"]
  var_multi <- fw_metrics[[metric_name]]["multi_shell"]
  
  results_row <- results_table[
    Metric == metric_name
  ]
  
  scatter_plots[[metric_name]] <- create_scatter_plot(
    var_single = var_single,
    var_multi = var_multi,
    metric_name = metric_name,
    results_row = results_row,
    global_min = global_lower,
    global_max = global_upper
  )
}

combined_scatter <- ggarrange(
  plotlist = unname(scatter_plots),
  ncol = 2,
  nrow = 2,
  labels = "AUTO",
  align = "hv"
)

print(combined_scatter)


########## 10. Save the scatter plot in TIFF and PDF formats ##########

scatter_tiff_file <- file.path(
  output_dir,
  "FW_Single_shell_vs_Multi_shell_Pearson_ICC_C1_scatter_revised.tiff"
)

scatter_pdf_file <- file.path(
  output_dir,
  "FW_Single_shell_vs_Multi_shell_Pearson_ICC_C1_scatter_revised.pdf"
)

ggsave(
  filename = scatter_tiff_file,
  plot = combined_scatter,
  device = "tiff",
  width = 12,
  height = 10,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  filename = scatter_pdf_file,
  plot = combined_scatter,
  device = "pdf",
  width = 12,
  height = 10,
  units = "in"
)


########## 11. Analysis complete ##########

cat("\n")
cat("====================================================\n")
cat("Analysis complete\n")
cat("====================================================\n")

cat("Statistical results: ", results_file, "\n")
cat("TIFF scatter plot: ", scatter_tiff_file, "\n")
cat("PDF scatter plot: ", scatter_pdf_file, "\n")

