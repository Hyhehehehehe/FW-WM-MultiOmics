# Protein-FW Mendelian Randomization
# This script performs Mendelian randomization (MR) between plasma proteins
# (cis-pQTL instruments from deCODE) and the four FW-WM indexes, and visualizes
# the results as volcano plots
# Workflow:
# 1. MR analysis: protein as exposure, FW-WM index as outcome (per array task)
#    - Instrument selection from cis-pQTL, harmonisation, F-statistic, MR-PRESSO,
#      MR estimates, heterogeneity and pleiotropy tests
# 2. Volcano plot: -log10(FDR) versus Beta for each FW-WM index
#    - FDR-adjusted significance and labeled significant proteins

#=============================================================================
# 1. MR analysis (protein as exposure, FW-WM index as outcome)
#=============================================================================
# Input:  cis-pQTL instruments from deCODE (clumped at r2 = 0.01)
#         Outcome GWAS summary statistics (one of cad, t2d, shtg)
# Output: <pheno>/decode_protein_<pheno>_mr_array<array.id>.txt

rm(list = ls())
a <- Sys.time()

library(data.table)
library(dplyr)
library(TwoSampleMR)
library(ggplot2)
library(MRPRESSO)

# ---- 1.1 Parse command-line arguments ----
array.id  <- as.numeric(commandArgs(TRUE)[1])
max_index <- as.numeric(commandArgs(TRUE)[2])
pheno     <- as.character(commandArgs(TRUE)[3])

task_count <- ceiling(max_index / 20)

if (array.id * task_count > max_index) {
  kk <- as.numeric(c(((array.id - 1) * task_count + 1) : max_index))
} else {
  kk <- as.numeric(c(((array.id - 1) * task_count + 1) : (array.id * task_count)))
}

# ---- 1.2 Load outcome GWAS ----
outgwas_file <- switch(
  pheno,
  "FW_Total"  = "G:/R/MR/GWAS/FW_gwas/ukb_FW.csv.gz",
  "FW_Association"  = "G:/R/MR/GWAS/FW_gwas/ukb_FW_Association.csv.gz",
  "FW_Commissural" = "G:/R/MR/GWAS/FW_gwas/ukb_FW_Commissural.csv.gz",
  "FW_Projection" = "G:/R/MR/GWAS/FW_gwas/ukb_FW_Projection.csv.gz"
)

out_gwas <- fread(outgwas_file)
out_gwas <- data.frame(out_gwas)
out_gwas_snp <- out_gwas$SNP

# ---- 1.3 MR analysis per protein ----
for (r2 in c(0.01)) {
  
  protein_name <- fread(paste0("clumped_decode_cis_protein_mafcut_0.01_full_name"),
                        header = FALSE)$V1[kk]
  
  cis_pqtl <- fread(paste0("0_merged_r2_", r2,
                           "_mafcut_0.01_decode_cis_pqtl.txt"))
  
  cis_pqtl <- cis_pqtl %>%
    filter(rsids %in% out_gwas_snp) %>%
    filter(!Name %in% unique(cis_pqtl$Name[which(duplicated(cis_pqtl$Name))])) %>%
    rename("outcome" = protein)
  
  protein_exp_dat <- format_data(
    as.data.frame(cis_pqtl),
    type                 = "exposure",
    phenotype_col        = "outcome",
    snp_col              = "rsids",
    beta_col             = "Beta",
    se_col               = "SE",
    effect_allele_col    = "effectAllele",
    other_allele_col     = "otherAllele",
    pval_col             = "Pval",
    eaf_col              = "effectAlleleFreq",
    samplesize_col       = "N"
  )
  
  outcome_dat <- format_data(
    out_gwas,
    type                 = "outcome",
    snps                 = protein_exp_dat$SNP,
    snp_col              = "SNP",
    beta_col             = "Beta",
    se_col               = "Se",
    effect_allele_col    = "A1",
    other_allele_col     = "A2",
    eaf_col              = "EAF",
    pval_col             = "P",
    samplesize_col       = "N"
  )
  outcome_dat$outcome <- pheno
  all_out_dat <- outcome_dat
  
  mat_mr <- NULL
  n <- 0
  
  for (i in 1:length(protein_name)) {
    
    protein <- protein_name[i]
    
    exp_dat <- protein_exp_dat %>% filter(exposure == protein)
    if (nrow(exp_dat) == 0 || is.null(exp_dat)) {
      print(paste("None iv exists in the outcome gwas dataset for exposure protein",
                  protein))
      next
    }
    
    outcome_dat <- all_out_dat %>% filter(SNP %in% exp_dat$SNP)
    
    dat <- harmonise_data(exposure_dat = exp_dat,
                          outcome_dat  = outcome_dat) %>%
      filter(mr_keep)
    
    if (nrow(dat) == 0 || is.null(dat)) {
      print(paste("None iv exists in the outcome gwas dataset after harmonising for protein",
                  protein))
      next
    }
    
    N  <- dat$samplesize.exposure
    k  <- nrow(dat)
    R2 <- sum((2 * (1 - dat$eaf.exposure) * dat$eaf.exposure * dat$beta.exposure^2) /
                (2 * (1 - dat$eaf.exposure) * dat$eaf.exposure * dat$beta.exposure^2 +
                   2 * (1 - dat$eaf.exposure) * dat$eaf.exposure * N * dat$se.exposure^2))
    F  <- (mean(N) - k - 1) / k * R2 / (1 - R2)
    
    if (F < 10) {
      print(paste("F statistic of ivs for protein", protein, "< 10"))
      next
    }
    
    print(paste("Analysis for protein", protein))
    
    presso_p <- "iv counts < 4"
    
    if (nrow(dat) > 3) {
      dat_mr_presso <- dat[, c("beta.exposure", "se.exposure", "pval.exposure",
                               "beta.outcome", "se.outcome", "pval.outcome")]
      presso <- mr_presso(BetaOutcome    = "beta.outcome",
                          BetaExposure   = "beta.exposure",
                          SdOutcome      = "se.outcome",
                          SdExposure     = "se.exposure",
                          OUTLIERtest    = TRUE,
                          DISTORTIONtest = TRUE,
                          data           = dat_mr_presso,
                          NbDistribution = 1000,
                          SignifThreshold = 0.05)
      outliers <- presso[[2]]$`Distortion Test`$`Outliers Indices`
      if (is.null(outliers) | length(outliers) == 0 |
          all(is.na(outliers)) | is.character(outliers)) {
        cat("No outlier were identified. Stopping iteration in mr_presso.\n")
        presso_p <- "No outlier were identified in mr presso"
      } else {
        presso_p <- presso[[2]]$`Distortion Test`$Pvalue
        dat <- dat[-outliers, ]
      }
    }
    
    N  <- dat$samplesize.exposure
    k  <- nrow(dat)
    R2 <- sum((2 * (1 - dat$eaf.exposure) * dat$eaf.exposure * dat$beta.exposure^2) /
                (2 * (1 - dat$eaf.exposure) * dat$eaf.exposure * dat$beta.exposure^2 +
                   2 * (1 - dat$eaf.exposure) * dat$eaf.exposure * N * dat$se.exposure^2))
    F  <- (mean(N) - k - 1) / k * R2 / (1 - R2)
    
    if (F < 10) {
      print(paste("F statistic of ivs for protein", protein,
                  "< 10 after mr presso"))
      next
    }
    
    res <- mr(dat, method_list = c("mr_wald_ratio", "mr_ivw",
                                   "mr_weighted_median", "mr_egger_regression"))
    
    mr_heterogeneity <- mr_heterogeneity(dat,
                                         method_list = c("mr_egger_regression"))
    if (nrow(mr_heterogeneity) == 0) {
      i2 <- NA
    } else {
      i2 <- (mr_heterogeneity$Q - mr_heterogeneity$Q_df) /
        mr_heterogeneity$Q * 100
    }
    
    mr_pleiotropy <- mr_pleiotropy_test(dat)
    if (is.null(mr_pleiotropy) | nrow(mr_pleiotropy) == 0 |
        ncol(mr_pleiotropy) == 0) {
      pleiotropy <- c(NA, NA, NA)
    } else {
      pleiotropy <- c(mr_pleiotropy$egger_intercept,
                      mr_pleiotropy$se, mr_pleiotropy$pval)
    }
    
    mr  <- unname(unlist(lapply(1:3, function(i) res[i, 7:9])))
    mat <- matrix(NA, 1, 18)
    mat[1, ] <- c(protein, res$outcome[1], res$nsnp[1], F, mr, i2,
                  pleiotropy, presso_p)
    mat_mr <- rbind(mat_mr, mat)
    n <- n + 1
    print(paste("Number", n, "protein"))
  }
  
  mat_mr <- data.frame(mat_mr)
  colnames(mat_mr) <- c("Exposure", "Outcome", "iv_counts", "F",
                        paste0(rep(c("wald/ivw", "weighted median", "mr egger"),
                                   each = 3),
                               c(".beta", ".se", ".pval")),
                        "I2_heterogeneity", "egger_intercept",
                        "se_pleiotropy", "P_pleiotropy", "P_mr_presso")
  
  write.table(mat_mr,
              file = paste0(pheno, "/decode_protein_", pheno,
                            "_mr_array", array.id, ".txt"),
              col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  
  message(paste("Perform MR analyses for", n, "of", task_count,
                "proteins in this array"))
}

b <- Sys.time()
cat("Elapsed Time:", as.numeric(b - a, units = "mins"), "minutes\n")

#=============================================================================
# 2. Volcano plot (Beta versus -log10(FDR))
#=============================================================================
# Input:  protein_FW_1.csv, protein_FW_Association.csv,
#         protein_FW_Projection.csv, protein_FW_Commissural.csv
# Output: pro_FW_<FW index>_volcano_<date>.pdf

rm(list = ls())

library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(tidyverse)

# ---- 2.1 File list ----
file_list <- c(
  "protein_FW_1.csv",
  "protein_FW_Association.csv",
  "protein_FW_Projection.csv",
  "protein_FW_Commissural.csv"
)

# ---- 2.2 Color gradient ----
color_gradient <- colorRampPalette(
  c("#3d78b1", "#b2daeb", "#e36d42", "#d82e1e", "darkred")
)

# ---- 2.3 Loop over each FW index ----
for (file in file_list) {
  
  # Read data and compute FDR
  data <- read.csv(file, header = TRUE)
  data$FDR <- p.adjust(data$pvalue, method = "fdr")
  
  # Compute -log10(FDR) and significance status
  data <- data %>%
    mutate(
      log_pvalue = -log10(FDR),
      status = case_when(
        FDR < 0.05 & Beta > 0 ~ "up",
        FDR < 0.05 & Beta < 0 ~ "down",
        TRUE ~ "stable"
      )
    )
  
  # Label significant proteins
  data$label <- ifelse(data$FDR < 0.05, data$Exposure, "")
  
  # Remove NA rows
  data <- na.omit(data)
  
  # Derive the FW index name for the output file
  file_type <- gsub(".csv", "", file)
  file_type <- gsub("protein_FW_", "", file_type)
  if (file_type == "protein_FW") file_type <- "total"
  
  # Volcano plot
  p <- ggplot(data, aes(x = Beta, y = log_pvalue)) +
    geom_point(aes(color = log_pvalue), alpha = 1, size = 2) +
    scale_color_gradientn(name = expression(-log[10](italic(FDR))),
                          colours = color_gradient(4),
                          limits = c(0, 3.5)) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed", color = "darkred", linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    labs(x = "Beta",
         y = expression(-log[10](italic(FDR))),
         title = paste("Volcano plot -", file_type)) +
    theme_bw() +
    ylim(c(0, 3.5)) +
    geom_text_repel(aes(label = label, color = log_pvalue), size = 3) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(hjust = 0.5)
    )
  
  # Save
  output_file <- paste0("pro_FW_", file_type, "_volcano_",
                        format(Sys.Date(), "%m%d"), ".pdf")
  ggsave(output_file, p, width = 10, height = 8)
  
  cat("Processed", file, "- saved as", output_file, "\n")
}

cat("All files processed successfully!\n")