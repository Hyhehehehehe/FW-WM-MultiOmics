# Disease-FW Mendelian Randomization
# This script performs Mendelian randomization (MR) between neurological diseases /
# hypertension and the four FW-WM indexes (FW, FW_Association, FW_Commissural,
# FW_Projection)
# Workflow:
# 1. Run MR analysis for every exposure-outcome pair (disease as exposure,
#    FW-WM index as outcome)
# 2. Detect and remove pleiotropic outlier SNPs with MR-PRESSO
# 3. Calculate odds ratios, heterogeneity and pleiotropy tests
# 4. Generate scatter, forest, funnel, and leave-one-out plots
# 5. Draw lollipop plots of the summary MR results (one panel per FW-WM index)

#=============================================================================
# 1. Configuration of exposures and outcomes
#=============================================================================
# Exposures: neurological diseases / hypertension from local GWAS summary
# statistics (csv files in FW_related_disease/)
# Outcomes:  four FW-WM indexes from UKB GWAS summary statistics
#            (ukb_FW_Association.csv.gz, ukb_FW_Commissural.csv.gz,
#             ukb_FW_Projection.csv.gz)
# Output:    <output_base>/<exposure>_<outcome>/ for each MR pair

library(data.table)
library(TwoSampleMR)
library(tidyverse)
library(MendelianRandomization)
library(MRPRESSO)
library(ggplot2)
library(readxl)

# ---- Paths ----
setwd("G:/R/MR/GWAS")
exposure_dir <- "FW_related_disease"
outcome_dir  <- "FW_gwas"
output_base  <- "G:/R/MR/MR/Disease_FW_MR"

# ---- File lists ----
exposure_files <- list.files(exposure_dir, pattern = "\\.csv$", full.names = TRUE)
outcome_files <- c(
  "ukb_FW_Association.csv.gz",
  "ukb_FW_Commissural.csv.gz",
  "ukb_FW_Projection.csv.gz"
)

if (!dir.exists(output_base)) dir.create(output_base, recursive = TRUE)

#=============================================================================
# 2. Main loop: iterate over each outcome and each exposure
#=============================================================================
for (outcome_file in outcome_files) {
  
  outcome_path <- file.path(outcome_dir, outcome_file)
  outcome_gwas <- fread(outcome_path)
  outcome_name <- tools::file_path_sans_ext(basename(outcome_file))
  
  # ---- Inner loop: iterate over each exposure ----
  for (exposure_file in exposure_files) {
    
    exposure_name <- tools::file_path_sans_ext(basename(exposure_file))
    
    output_dir <- file.path(output_base, paste0(exposure_name, "_", outcome_name))
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    tryCatch({
      
      # ==================== Data preparation ====================
      exp <- read_exposure_data(
        filename             = exposure_file,
        sep                  = ",",
        snp_col              = "SNP",
        beta_col             = "beta",
        se_col               = "se",
        effect_allele_col    = "effect_allele",
        other_allele_col     = "other_allele",
        pval_col             = "pval.exposure"
      )
      
      valid_snps <- exp$SNP
      
      outcomeTab <- outcome_gwas[SNP %in% valid_snps] %>%
        merge(exp[, c("SNP", "effect_allele.exposure", "other_allele.exposure")],
              by = "SNP")
      
      write.csv(outcomeTab,
                file.path(output_dir, "Matched_Outcome.csv"),
                row.names = FALSE)
      
      outcome_dat <- read_outcome_data(
        snps                 = valid_snps,
        filename             = file.path(output_dir, "Matched_Outcome.csv"),
        sep                  = ",",
        snp_col              = "SNP",
        beta_col             = "beta",
        se_col               = "se",
        effect_allele_col    = "effect_allele",
        other_allele_col     = "other_allele",
        pval_col             = "pval"
      )
      
      # ==================== Core analysis ====================
      dat <- harmonise_data(exposure_dat = exp, outcome_dat = outcome_dat)
      write.csv(dat,
                file.path(output_dir, "Harmonised_Data.csv"),
                row.names = FALSE)
      
      # MR-PRESSO requires at least 5 SNPs
      num_snps <- nrow(dat)
      presso_skipped <- FALSE
      
      if (num_snps >= 5) {
        tryCatch({
          presso <- run_mr_presso(dat, NbDistribution = 1000)
          
          if (!is.null(presso[[1]]$`MR-PRESSO results`$`Outlier Test`)) {
            write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`,
                      file.path(output_dir, "PRESSO_Outliers.csv"),
                      row.names = FALSE)
            
            outlier_indices <- presso[[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
            if (!is.null(outlier_indices) && length(outlier_indices) > 0) {
              dat_clean <- dat[-outlier_indices, ]
              write.csv(dat_clean,
                        file.path(output_dir, "Cleaned_Data.csv"),
                        row.names = FALSE)
            } else {
              dat_clean <- dat
            }
          } else {
            dat_clean <- dat
          }
        }, error = function(e) {
          warning("MR-PRESSO failed for ", exposure_name, " + ", outcome_name,
                  " with error: ", e$message)
          dat_clean <<- dat
          presso_skipped <<- TRUE
        })
      } else {
        message("Skipping MR-PRESSO: ", exposure_name, " + ", outcome_name,
                " (only ", num_snps, " SNP)")
        dat_clean <- dat
        presso_skipped <- TRUE
      }
      
      if (presso_skipped) {
        skip_reason <- ifelse(num_snps < 5,
                              paste("Insufficient SNPs (n =", num_snps, ")"),
                              "MR-PRESSO failed")
        writeLines(skip_reason, file.path(output_dir, "PRESSO_SKIPPED.txt"))
      }
      
      # ==================== MR analysis ====================
      mr_results <- mr(dat_clean,
                       method_list = c("mr_ivw",
                                       "mr_egger_regression",
                                       "mr_weighted_median"))
      mr_or <- generate_odds_ratios(mr_results)
      
      write.csv(mr_results,
                file.path(output_dir, "MR_Results.csv"),
                row.names = FALSE)
      write.csv(mr_or,
                file.path(output_dir, "MR_Odds_Ratios.csv"),
                row.names = FALSE)
      
      # ==================== Quality control ====================
      heterogeneity <- mr_heterogeneity(dat_clean)
      write.csv(heterogeneity,
                file.path(output_dir, "Heterogeneity.csv"),
                row.names = FALSE)
      
      pleiotropy <- mr_pleiotropy_test(dat_clean)
      write.csv(pleiotropy,
                file.path(output_dir, "Pleiotropy.csv"),
                row.names = FALSE)
      
      # ==================== Visualization ====================
      pdf(file.path(output_dir, "MR_Scatter.pdf"), width = 8, height = 7)
      print(mr_scatter_plot(mr_results, dat_clean))
      dev.off()
      
      res_single <- mr_singlesnp(dat_clean)
      
      pdf(file.path(output_dir, "Forest_Plot.pdf"), width = 9, height = 8)
      print(mr_forest_plot(res_single))
      dev.off()
      
      pdf(file.path(output_dir, "Funnel_Plot.pdf"), width = 8, height = 7)
      print(mr_funnel_plot(res_single))
      dev.off()
      
      pdf(file.path(output_dir, "LeaveOneOut.pdf"), width = 9, height = 8)
      print(mr_leaveoneout_plot(mr_leaveoneout(dat_clean)))
      dev.off()
      
      # ==================== Progress ====================
      message(sprintf("Completed: exposure [%s] + outcome [%s] @ %s",
                      exposure_name, outcome_name,
                      format(Sys.time(), "%X")))
      
    }, error = function(e) {
      error_msg <- paste("Analysis failed:", exposure_name, "+", outcome_name,
                         ":", e$message)
      message(error_msg)
      writeLines(error_msg, file.path(output_dir, "ERROR_LOG.txt"))
    })
  }
}

message("All MR analyses completed. Results saved to: ", output_base)

#=============================================================================
# 3. Visualization of the MR results (lollipop plot)
#=============================================================================
# This section reads the aggregated MR results (Disease_FW_MR_result.xlsx)
# with one sheet per FW-WM index, and draws a lollipop plot for each sheet.
# Columns used: disease name, category (Richness), b (column 5), FDR (column 8)
# Output: Lollipop_<FW index>.pdf

setwd("G:/FW_analysis/Fig/Fig5.Disease_Pro_FW/plot_data")

sheet_names <- c("FW", "FW_Association", "FW_Projection", "FW_Commissural")

for (sheet in sheet_names) {
  
  df <- read_xlsx("Disease_FW_MR_result.xlsx", sheet = sheet)
  plot_data <- df[, c(1, 2, 5, 8)]
  colnames(plot_data) <- c("disease", "Richness", "b", "FDR")
  plot_data$b <- plot_data$b * 100
  
  p2 <- ggplot(plot_data,
               aes(x = reorder(disease, b),
                   y = b,
                   color = disease)) +
    geom_segment(aes(x = disease, xend = disease, y = 0, yend = b),
                 linewidth = 1, color = "gray60") +
    geom_point(size = 10, shape = 1) +
    geom_point(size = 8) +
    geom_text(
      aes(y = b + 0.05,
          label = case_when(
            FDR < 0.001 ~ "***",
            FDR < 0.01  ~ "**",
            FDR < 0.05  ~ "*",
            TRUE ~ ""
          )),
      size = 5, color = "black", show.legend = FALSE
    ) +
    scale_color_manual(values = c("#5ebcc2", "#46a9cb", "#5791c9", "#2D59A8",
                                  "#DC93BF", "#C77CFF", "#7a76b7", "#945893",
                                  "#9c3d62")) +
    scale_y_continuous(expand = c(0, 0), limits = c(-0.5, 1.5)) +
    theme_bw() +
    coord_flip() +
    geom_hline(aes(yintercept = 0), linewidth = 0.7, color = "gray42") +
    facet_grid(~Richness, scales = "free", space = "free_x") +
    theme(
      strip.background = element_rect(fill = "#FFF6E1"),
      strip.text       = element_text(size = 12, face = "bold", color = "#a16e53"),
      axis.text        = element_text(color = "black", size = 11),
      axis.title       = element_blank()
    )
  
  output_file <- paste0("Lollipop_", sub("FW_", "", sheet), ".pdf")
  ggsave(output_file, plot = p2, width = 9, height = 8)
  
  message("Generated: ", output_file)
}

message("All lollipop plots completed.")