# =============================================================
# Mendelian Randomization: Hypertension → Protein
# Exposure : Hypertension (ebi-a-GCST90038604)
# Outcomes : PEAR1, CATF, CR3L4, HP_Mixed_Type, HP_HPT,
#            ART4, LRIG3, MLP3A, ANGL3
# =============================================================

library(data.table)
library(TwoSampleMR)
library(tidyverse)
library(MendelianRandomization)
library(MRPRESSO)

# ---- Paths ----
setwd("G:/R/MR/GWAS")
outcome_dir <- "FW_protein_gwas"
output_base <- "G:/R/MR/MR/hypertension_protein_FW_MR"

# ---- Outcome file list ----
outcome_dir2  <- list.files(outcome_dir, pattern = "\\protein$", full.names = TRUE)
outcome_files <- list.files(outcome_dir2, full.names = TRUE)

# ---- Exposure ----
exposure_ids <- c("ebi-a-GCST90038604")   # Hypertension

if (!dir.exists(output_base)) dir.create(output_base, recursive = TRUE)

# =============================================================
# Part 1: MR analysis
# =============================================================
for (outcome_file in outcome_files) {
  
  outcome_gwas <- fread(outcome_file)
  outcome_name <- tools::file_path_sans_ext(basename(outcome_file))
  
  for (exposure_id in exposure_ids) {
    
    output_dir <- file.path(output_base, paste0(exposure_id, "_", outcome_name))
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    tryCatch({
      
      # ---- Exposure data ----
      exp_dat <- extract_instruments(
        outcomes = exposure_id,
        p1       = 5e-08,
        clump    = TRUE,
        p2       = 5e-08,
        r2       = 0.01,
        kb       = 10000
      )
      
      if (is.null(exp_dat) || nrow(exp_dat) == 0) {
        stop("Failed to retrieve exposure data from IEU: ", exposure_id)
      }
      
      write.csv(exp_dat,
                file.path(output_dir, "Exposure_Data.csv"),
                row.names = FALSE)
      
      # ---- Outcome data ----
      valid_snps <- exp_dat$SNP
      outcomeTab <- outcome_gwas[SNP %in% valid_snps]
      outcomeTab <- data.frame(outcomeTab)
      
      if (nrow(outcomeTab) == 0) {
        stop("No matched SNPs: ", exposure_id, " vs ", outcome_name)
      }
      
      write.csv(outcomeTab,
                file.path(output_dir, "Matched_Outcome.csv"),
                row.names = FALSE)
      
      outcome_dat <- format_data(
        outcomeTab,
        type                 = "outcome",
        snps                 = exp_dat$SNP,
        snp_col              = "SNP",
        beta_col             = "beta",
        se_col               = "se",
        effect_allele_col    = "effect_allele",
        other_allele_col     = "other_allele",
        pval_col             = "pval"
      )
      
      # ---- Harmonisation ----
      dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = outcome_dat)
      write.csv(dat,
                file.path(output_dir, "Harmonised_Data.csv"),
                row.names = FALSE)
      
      # ---- MR-PRESSO ----
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
          warning("MR-PRESSO failed for ", exposure_id, " + ", outcome_name,
                  " with error: ", e$message)
          dat_clean <<- dat
          presso_skipped <<- TRUE
        })
      } else {
        message("Skipping MR-PRESSO: ", exposure_id, " + ", outcome_name,
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
      
      # ---- MR analysis ----
      mr_results <- mr(dat_clean, method_list = c(
        "mr_ivw", "mr_egger_regression", "mr_weighted_median",
        "mr_simple_mode", "mr_weighted_mode"
      ))
      mr_or <- generate_odds_ratios(mr_results)
      
      write.csv(mr_results,
                file.path(output_dir, "MR_Results.csv"),
                row.names = FALSE)
      write.csv(mr_or,
                file.path(output_dir, "MR_Odds_Ratios.csv"),
                row.names = FALSE)
      
      # ---- Quality control ----
      heterogeneity <- mr_heterogeneity(dat_clean)
      write.csv(heterogeneity,
                file.path(output_dir, "Heterogeneity.csv"),
                row.names = FALSE)
      
      pleiotropy <- mr_pleiotropy_test(dat_clean)
      write.csv(pleiotropy,
                file.path(output_dir, "Pleiotropy.csv"),
                row.names = FALSE)
      
      # ---- Visualization ----
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
      
      message(sprintf("Completed: exposure [%s] + outcome [%s] (SNP = %d) @ %s",
                      exposure_id, outcome_name, num_snps,
                      format(Sys.time(), "%X")))
      
    }, error = function(e) {
      error_msg <- paste("Analysis failed:", exposure_id, "+", outcome_name,
                         ":", e$message)
      message(error_msg)
      writeLines(error_msg, file.path(output_dir, "ERROR_LOG.txt"))
    })
  }
}

message("All MR analyses completed. Results saved to: ", output_base)

# =============================================================
# Part 2: Combine MR_Odds_Ratios across all outcome folders
# =============================================================

main_folder <- output_base
top_folders <- list.files(main_folder, pattern = "\\.txt$")

# Import each MR_Odds_Ratios.csv as a named data frame
for (top_folder in top_folders) {
  
  top_folder_path <- file.path(main_folder, top_folder)
  csv_file_path   <- file.path(top_folder_path, "MR_Odds_Ratios.csv")
  
  new_variable_name <- paste0(top_folder)
  
  cat("Importing:", csv_file_path, "\n")
  cat("Variable name:", new_variable_name, "\n")
  
  assign(new_variable_name, read.csv(csv_file_path))
  
  cat("Import successful.\n\n")
}

# Tag each data frame with its dataset name
all_objects <- ls()
dataset_objects <- all_objects[sapply(all_objects, function(x) is.data.frame(get(x)))]

for (dataset_name in dataset_objects) {
  dataset <- get(dataset_name)
  dataset$DatasetName <- dataset_name
  assign(dataset_name, dataset)
}

# Combine all data frames by row
all_objects <- ls()
combined_data <- data.frame()

for (object_name in all_objects) {
  if (is.data.frame(get(object_name))) {
    combined_data <- rbind(combined_data, get(object_name))
  }
}

dir.create(file.path(output_base, "MR_results"), showWarnings = FALSE)

write.csv(combined_data,
          file.path(output_base, "MR_results", "MR_Odds_Ratios.csv"),
          row.names = FALSE)