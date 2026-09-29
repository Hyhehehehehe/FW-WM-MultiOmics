# FW-Disease Mendelian Randomization
# This script performs Mendelian randomization (MR) between the four FW-WM indexes
# (FW, FW_Association, FW_Commissural, FW_Projection) and neurological diseases / hypertension
# Workflow:
# 1. Prepare clumped exposure data for the four FW-WM indexes from the GWAS summary statistics
# 2. Merge exposure SNPs with outcome data and harmonise the alleles
# 3. Detect and remove pleiotropic outlier SNPs with MR-PRESSO
# 4. Run MR analysis and calculate odds ratios
# 5. Perform heterogeneity and pleiotropy tests
# 6. Generate scatter, forest, funnel, and leave-one-out plots
# 7. Draw the summary forest plot of the MR results (one forest plot per FW-WM index)

#=============================================================================
# 1. Preparation of exposure data for Mendelian randomization
#=============================================================================
# Genome-wide significant SNPs (P < 5e-8) for each FW-WM index are extracted
# from the GWAS summary statistics and clumped to obtain independent instruments
# Input:  GWAS summary statistics (ukb_<FW index>.csv.gz) with columns
#         SNP, CHR, BP, eaf, n, beta, se, pval, effect_allele, other_allele (in this order)
# Output: exposure_data/<FW index>_exp.csv

library(data.table)
library(TwoSampleMR)
library(MRPRESSO)

# Four FW-WM indexes analyzed in this study
fw_indexes <- c("FW", "FW_Association", "FW_Commissural", "FW_Projection")

# Output directory for the clumped exposure data
exp_dir <- "exposure_data"
if (!dir.exists(exp_dir)) dir.create(exp_dir, recursive = TRUE)

for (fw in fw_indexes) {
  message("Preparing exposure data: ", fw)
  
  # Read GWAS summary statistics and keep genome-wide significant SNPs
  gwas <- fread(paste0("ukb_", fw, ".csv.gz"))
  gwas <- subset(gwas, pval < 5e-8)
  colnames(gwas) <- c("SNP", "CHR", "BP", "eaf", "n", "beta", "se",
                      "pval.exposure", "effect_allele", "other_allele")
  
  # Clumping to obtain independent instrumental variables
  gwas_clumped <- clump_data(gwas,
                             clump_kb = 10000,
                             clump_r2 = 0.1,
                             clump_p1 = 5e-8,
                             clump_p2 = 1e-5)
  
  write.csv(gwas_clumped, file.path(exp_dir, paste0(fw, "_exp.csv")),
            row.names = FALSE)
}

#=============================================================================
# 2. Configuration of the outcome GWAS summary statistics
#=============================================================================
# Outcomes: intracerebral hemorrhage (ICH), stroke, ischemic stroke (IS),
# small vessel stroke (SVS), Alzheimer's disease (AD), dementia (DEMNAS),
# vascular dementia (VASCDEM) and multiple sclerosis (MS) from local GWAS
# summary statistics, plus hypertension (HBP; OpenGWAS ID ebi-a-GCST90038604),
# the only vascular outcome finally included in the study
#   file:      path to the summary statistics (relative path)
#   drop:      column positions removed before renaming
#   new_names: column names after removal (must contain SNP, beta,
#              se, effect_allele, other_allele and pval)
outcome_configs <- list(
  ICH = list(
    file = "GWAS_sumstats/finngen_R11_I9_ICH.gz",
    drop = c(6, 8, 12, 13),
    new_names = c("CHR", "BP", "other_allele", "effect_allele",
                  "SNP", "pval", "beta", "se", "eaf")
  ),
  Stroke = list(
    file = "GWAS_sumstats/GCST90104539_buildGRCh37_snp.tsv.gz",
    drop = c(1:4, 6:7),
    new_names = c("SNP", "CHR", "BP", "eaf", "beta", "se", "pval",
                  "OR", "ci_lower", "ci_upper", "effect_allele", "other_allele")
  ),
  IS = list(
    file = "GWAS_sumstats/GCST90104540_buildGRCh37_snp.tsv.gz",
    drop = c(1:4, 6:7),
    new_names = c("SNP", "CHR", "BP", "eaf", "beta", "se", "pval",
                  "OR", "ci_lower", "ci_upper", "effect_allele", "other_allele")
  ),
  SVS = list(
    file = "GWAS_sumstats/GCST90104543_buildGRCh37_snp.tsv.gz",
    drop = c(1:4, 6:7),
    new_names = c("SNP", "CHR", "BP", "eaf", "beta", "se", "pval",
                  "OR", "ci_lower", "ci_upper", "effect_allele", "other_allele")
  ),
  AD = list(
    file = "GWAS_sumstats/finngen_R11_G6_ALZHEIMER.gz",
    drop = c(6, 8, 12, 13),
    new_names = c("CHR", "BP", "other_allele", "effect_allele",
                  "SNP", "pval", "beta", "se", "eaf")
  ),
  DEMNAS = list(
    file = "GWAS_sumstats/finngen_R11_F5_DEMNAS.gz",
    drop = c(6, 8, 12, 13),
    new_names = c("CHR", "BP", "other_allele", "effect_allele",
                  "SNP", "pval", "beta", "se", "eaf")
  ),
  VASCDEM = list(
    file = "GWAS_sumstats/finngen_R11_F5_VASCDEM.gz",
    drop = c(6, 8, 12, 13),
    new_names = c("CHR", "BP", "other_allele", "effect_allele",
                  "SNP", "pval", "beta", "se", "eaf")
  ),
  MS = list(
    file = "GWAS_sumstats/MS1_GCST90043738_buildGRCh37.tsv.gz",
    drop = c(8:10, 14),
    new_names = c("CHR", "SNP", "BP", "effect_allele", "other_allele",
                  "N", "eaf", "beta", "se", "pval")
  )
)

#=============================================================================
# 3. Core MR pipeline once a harmonised data frame is available
#=============================================================================
run_mr_pipeline <- function(dat, file_prefix) {
  # Number of valid instrumental variables (mr_keep = TRUE)
  num_keep_snps <- sum(dat$mr_keep)
  
  if (num_keep_snps >= 4) {
    dat_keep <- dat[dat$mr_keep, ]
    
    # MR-PRESSO outlier detection (biased SNPs)
    presso_result <- tryCatch({
      run_mr_presso(dat_keep)
    }, error = function(e) {
      warning(paste("MR-PRESSO failed:", e$message))
      return(NULL)
    })
    
    if (!is.null(presso_result)) {
      write.csv(presso_result[[1]]$`MR-PRESSO results`$`Outlier Test`,
                file = paste0(file_prefix, "-PRESSO.csv"))
      
      b <- presso_result[[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
      if (!is.null(b) && is.numeric(b)) {
        dat_clean <- dat_keep[-b, ]
        message(paste("Removed", length(b), "outliers (remaining",
                      nrow(dat_clean), "SNPs)"))
      } else {
        dat_clean <- dat_keep
        message(paste("No significant outliers detected (",
                      num_keep_snps, "SNPs)"))
      }
    } else {
      dat_clean <- dat_keep
      message("Skipped outlier removal due to MR-PRESSO failure")
    }
  } else {
    # Too few valid SNPs for MR-PRESSO; use mr_keep = TRUE SNPs directly
    dat_clean <- dat[dat$mr_keep, ]
    message(paste("Skipped MR-PRESSO - only", num_keep_snps,
                  "valid instrumental variables"))
  }
  
  # Save the final data used for MR analysis
  write.csv(dat_clean, paste0(file_prefix, "-MRdat.csv"), row.names = FALSE)
  
  if (nrow(dat_clean) == 0) {
    warning("No valid SNPs remaining. Skipping MR analysis.")
    return(invisible(NULL))
  }
  
  # MR analysis
  mrResult <- mr(dat_clean)
  write.csv(mrResult, paste0(file_prefix, "-MRresult.csv"), row.names = FALSE)
  
  # Odds ratios
  mrTab <- generate_odds_ratios(mrResult)
  write.csv(mrTab, paste0(file_prefix, "-MRresultOR.csv"), row.names = FALSE)
  
  # Heterogeneity test
  heterTab <- mr_heterogeneity(dat_clean)
  write.csv(heterTab, paste0(file_prefix, "-heterogeneity.csv"), row.names = FALSE)
  
  # Pleiotropy test
  pleioTab <- mr_pleiotropy_test(dat_clean)
  write.csv(pleioTab, paste0(file_prefix, "-pleiotropy.csv"), row.names = FALSE)
  
  # Scatter plot
  pdf(paste0(file_prefix, "-scatter-plot.pdf"), width = 7.5, height = 7)
  print(mr_scatter_plot(mrResult, dat_clean))
  dev.off()
  
  # Forest plot (effect of each instrumental variable)
  res_single <- mr_singlesnp(dat_clean)
  pdf(paste0(file_prefix, "-forest.pdf"), width = 7, height = 6.5)
  print(mr_forest_plot(res_single))
  dev.off()
  
  # Funnel plot
  pdf(paste0(file_prefix, "-funnel-plot.pdf"), width = 7, height = 6.5)
  print(mr_funnel_plot(singlesnp_results = res_single))
  dev.off()
  
  # Leave-one-out sensitivity analysis
  pdf(paste0(file_prefix, "-leaveoneout.pdf"), width = 7, height = 6.5)
  print(mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat_clean)))
  dev.off()
  
  message(paste("MR analysis completed:", file_prefix,
                "with", nrow(dat_clean), "valid SNPs"))
}

#=============================================================================
# 4. MR analysis against an outcome from local GWAS summary statistics
#=============================================================================
run_mr_local_outcome <- function(exposure_a, cfg, pair_name) {
  result_dir <- file.path("Result_FW_disease", pair_name)
  if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)
  file_prefix <- file.path(result_dir, pair_name)
  
  # Import outcome summary statistics
  outcome_raw <- fread(cfg$file)
  keep_cols <- setdiff(seq_len(ncol(outcome_raw)), cfg$drop)
  outcome_raw <- outcome_raw[, keep_cols, with = FALSE]
  colnames(outcome_raw) <- cfg$new_names
  
  # Find the exposure SNPs in the outcome data
  outcomeTab <- merge(exposure_a[, "SNP", drop = FALSE],
                      outcome_raw, by = "SNP")
  write.csv(outcomeTab, paste0(file_prefix, "-merge.csv"), row.names = FALSE)
  
  outcome_dat <- read_outcome_data(snps = exposure_a$SNP,
                                   filename = paste0(file_prefix, "-merge.csv"),
                                   sep = ",",
                                   snp_col = "SNP",
                                   beta_col = "beta",
                                   se_col = "se",
                                   effect_allele_col = "effect_allele",
                                   other_allele_col = "other_allele",
                                   pval_col = "pval")
  
  # Harmonisation
  dat <- harmonise_data(exposure_dat = exposure_a, outcome_dat = outcome_dat)
  write.csv(dat, paste0(file_prefix, "-harmonise.csv"), row.names = FALSE)
  
  run_mr_pipeline(dat, file_prefix)
}

#=============================================================================
# 5. Main loop: run every outcome for each of the four FW-WM indexes
#=============================================================================

# Hypertension (HBP) was the only vascular outcome finally included
hbp_id <- "ebi-a-GCST90038604"

for (fw in fw_indexes) {
  message("==================== Exposure: ", fw, " ====================")
  
  # Import exposure data (clumped instruments prepared in Step 1)
  exposure_a <- read_exposure_data(
    filename = file.path("exposure_data", paste0(fw, "_exp.csv")),
    sep = ",",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval.exposure"
  )
  
  # Outcomes from local GWAS summary statistics
  for (outcome_name in names(outcome_configs)) {
    tryCatch({
      run_mr_local_outcome(exposure_a, outcome_configs[[outcome_name]],
                           paste0(fw, "_", outcome_name))
    }, error = function(e) {
      warning(paste("Error processing", fw, "-", outcome_name, ":", e$message))
    })
  }
  
  # Hypertension outcome from OpenGWAS
  tryCatch({
    result_dir <- file.path("Result_FW_disease", paste0(fw, "_HBP"))
    if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)
    file_prefix <- file.path(result_dir, paste0(fw, "_HBP"))
    
    out <- extract_outcome_data(snps = exposure_a$SNP, outcomes = hbp_id)
    write.csv(out, paste0(file_prefix, "-merge.csv"), row.names = FALSE)
    
    dat <- harmonise_data(exposure_dat = exposure_a, outcome_dat = out)
    write.csv(dat, paste0(file_prefix, "-harmonise.csv"), row.names = FALSE)
    
    run_mr_pipeline(dat, file_prefix)
  }, error = function(e) {
    warning(paste("Error processing", fw, "- HBP:", e$message))
  })
}

#=============================================================================
# 6. Visualization: summary forest plot of the MR results
#=============================================================================
# This section reads the aggregated MR results (FW_Disease_Converted_MR_Results_1.xlsx)
# with one sheet per FW index, and draws a forest plot for each sheet.
# The forest plot displays the odds ratio (OR) and 95% confidence interval
# for each FW-WM index and each outcome, together with the number of SNPs
# and the FDR-adjusted P value.
# Output: Forest_<FW index>_disease.pdf

library(readxl)
library(forestploter)
library(grid)
library(magrittr)
library(tidyselect)
library(dplyr)

# 定义要处理的工作表列表
sheet_names <- c("FW", "FW_Association", "FW_Projection", "FW_Commissural")

# 循环处理每个工作表
for (sheet in sheet_names) {
  
  # 加载数据
  dat <- read_excel("FW_Disease_Converted_MR_Results_1.xlsx", sheet = sheet)
  
  # 创建一个包含 20 个空格的向量（用于森林图空白列）
  dat$` ` <- paste(rep(" ", 20), collapse = "")
  
  # 格式化 OR (95% CI)
  dat$`OR(95% CI)` <- sprintf("%.3f (%.3f-%.3f)", dat$or, dat$or_lci95, dat$or_uci95)
  
  # 格式化 PFDR
  dat$PFDR <- ifelse(
    dat$PFDR < 0.0001,
    formatC(dat$PFDR, format = "e", digits = 3),
    sprintf("%.3f", dat$PFDR)
  )
  
  # 选择用于森林图显示的列
  plot_dat <- dat %>%
    select(exposure, outcome, nsnp, or, ` `, PFDR, `OR(95% CI)`,
           lo_ci, up_ci, or_lci95, or_uci95)
  
  # 森林图主题
  tm <- forest_theme(base_size       = 9,
                     ci_pch         = 15,
                     ci_col         = "#4575b4",
                     ci_lty         = 1,
                     ci_lwd         = 0.8,
                     ci_theight     = 0.2,
                     xaxis_lwd      = 0.8,
                     xaxis_cex      = 0.6,
                     refline_lwd    = 0.8,
                     refline_lty    = "dashed",
                     refline_col    = "red",
                     summary_fill   = "#4575b4",
                     summary_col    = "#4575b4",
                     footnote_cex   = 1.1,
                     footnote_fontface = "italic",
                     footnote_col   = "blue")
  
  # 根据 FW_Commissural 使用不同的 x 轴范围（与原始代码一致）
  if (sheet == "FW_Commissural") {
    xlim_use  <- c(0, 2)
    ticks_use <- c(0, 0.5, 1, 1.5, 2)
  } else {
    xlim_use  <- c(0.5, 2.5)
    ticks_use <- c(0.5, 1, 1.5, 2, 2.5)
  }
  
  # 创建森林图
  p <- forest(plot_dat[, c(1:7)],
              est        = plot_dat$or,
              lower      = plot_dat$or_lci95,
              upper      = plot_dat$or_uci95,
              sizes      = 0.25,
              ci_column  = 5,
              ref_line   = 1,
              arrow_lab  = c("risk", "protect"),
              xlim       = xlim_use,
              ticks_at   = ticks_use,
              nudge_y    = 0,
              theme      = tm)
  
  # 保存为 PDF
  output_file <- paste0("Forest_", sub("FW_", "", sheet), "_disease.pdf")
  pdf(output_file, width = 8, height = 8)
  print(p)
  dev.off()
  
  message("已生成: ", output_file)
}

message("所有森林图处理完成！")