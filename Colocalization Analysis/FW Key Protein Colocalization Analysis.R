# =============================================================
# FW Key Protein Colocalization Analysis
# Proteins: ANGPTL3, PEAR1
# GWAS outcomes: FW_total, FW_projection
# =============================================================

library(data.table)
library(biomaRt)
library(coloc)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(patchwork)
library(locuszoomr)
library(EnsDb.Hsapiens.v75)


# =============================================================
# 1. Colocalization: FW_projection + ANGPTL3
# =============================================================

# ---- 1.1 Load FW GWAS and compute MAF ----
gwas_data <- fread("G:/R/MR/GWAS/FW_gwas/ukb_FW_Projection.csv.gz", header = TRUE)
gwas_data[, MAF := pmin(eaf, 1 - eaf)]

# ---- 1.2 Retrieve gene coordinates ----
mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  GRCh    = 37
)

gene_pos <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name",
                 "start_position", "end_position"),
  filters    = "hgnc_symbol",
  values     = "ANGPTL3",
  mart       = mart
)

cis_window <- 5e5
gene_chr   <- gene_pos$chromosome_name
gene_start <- gene_pos$start_position
gene_end   <- gene_pos$end_position

# ---- 1.3 Subset FW GWAS to cis region ----
gwas_region <- gwas_data[
  CHR == as.integer(gene_chr) &
    BP >= (gene_start - cis_window) &
    BP <= (gene_end + cis_window)
]

# ---- 1.4 Load and subset pQTL ----
pQTL <- fread("G:/R/MR/GWAS/FW_protein_gwas/10391_1_ANGPTL3_ANGL3_gwas.txt.gz")
pQTL$Chrom <- as.integer(gsub("chr", "", pQTL$Chrom))

cis_pQTL <- pQTL[
  Chrom == as.integer(gene_chr) &
    Pos >= (gene_start - cis_window) &
    Pos <= (gene_end + cis_window)
]

cis_pQTL <- cis_pQTL[order(Pval), .SD[1], by = rsids]

cis_pQTL_clean <- cis_pQTL[, .(
  SNP  = rsids,
  CHR  = Chrom,
  BP   = Pos,
  beta = Beta,
  se   = SE,
  MAF  = ImpMAF,
  eaf  = effectAlleleFreq,
  N    = N,
  A1   = effectAllele,
  A2   = otherAllele,
  pval = Pval
)]

# ---- 1.5 Merge and build coloc datasets ----
merged_data <- merge(cis_pQTL_clean, gwas_region, by = "SNP")

gwas_dataset <- list(
  beta    = merged_data$beta.y,
  varbeta = merged_data$se.y^2,
  type    = "quant",
  N       = merged_data$n,
  MAF     = merged_data$MAF.y,
  snp     = merged_data$SNP
)

pQTL_dataset <- list(
  beta    = merged_data$beta.x,
  varbeta = merged_data$se.x^2,
  type    = "quant",
  N       = merged_data$N,
  MAF     = merged_data$MAF.x,
  snp     = merged_data$SNP
)

# ---- 1.6 Run coloc.abf ----
result_abf <- coloc.abf(
  dataset1 = gwas_dataset,
  dataset2 = pQTL_dataset,
  p1  = 1e-4,
  p2  = 1e-4,
  p12 = 1e-5
)

post_df <- data.frame(
  Hypothesis = c("H0 (no assoc)", "H1 (GWAS only)", "H2 (pQTL only)",
                 "H3 (both, independent)", "H4 (colocalized)"),
  Probability = c(result_abf$summary["PP.H0.abf"],
                  result_abf$summary["PP.H1.abf"],
                  result_abf$summary["PP.H2.abf"],
                  result_abf$summary["PP.H3.abf"],
                  result_abf$summary["PP.H4.abf"])
)

write.csv(post_df,
          "G:/FW_analysis/protein_gwas/coloc/FW_projection_ANGPTL3_UKB_coloc_500kb.csv",
          row.names = FALSE)


# =============================================================
# 2. Colocalization: FW_total + PEAR1
# =============================================================

# ---- 2.1 Load FW GWAS and compute MAF ----
gwas_data <- fread("G:/R/MR/GWAS/FW_gwas/ukb_FW.csv.gz", header = TRUE)
gwas_data[, MAF := pmin(eaf, 1 - eaf)]

# ---- 2.2 Retrieve gene coordinates ----
gene_pos <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name",
                 "start_position", "end_position"),
  filters    = "hgnc_symbol",
  values     = "PEAR1",
  mart       = mart
)

cis_window <- 5e5
gene_chr   <- gene_pos$chromosome_name
gene_start <- gene_pos$start_position
gene_end   <- gene_pos$end_position

# ---- 2.3 Subset FW GWAS to cis region ----
gwas_region <- gwas_data[
  CHR == as.integer(gene_chr) &
    BP >= (gene_start - cis_window) &
    BP <= (gene_end + cis_window)
]

# ---- 2.4 Load and subset pQTL ----
pQTL <- fread("G:/R/MR/GWAS/FW_protein_gwas/PEAR1.csv.gz")

cis_pQTL <- pQTL[
  CHR == as.integer(gene_chr) &
    BP >= (gene_start - cis_window) &
    BP <= (gene_end + cis_window)
]

cis_pQTL_clean <- cis_pQTL[, .(
  SNP  = SNP,
  CHR  = CHR,
  BP   = BP,
  beta = BETA,
  se   = SE,
  MAF  = IMPMAF,
  eaf  = FRQ,
  N    = N,
  A1   = A1,
  A2   = A2,
  pval = P
)]

write.csv(cis_pQTL_clean,
          "G:/FW_analysis/protein_gwas/PEAR1/PEAR1_UKB_cis_pqtl_500kb.csv",
          row.names = FALSE)

# ---- 2.5 Merge and build coloc datasets ----
merged_data <- merge(cis_pQTL_clean, gwas_region, by = "SNP")

gwas_dataset <- list(
  beta    = merged_data$beta.y,
  varbeta = merged_data$se.y^2,
  type    = "quant",
  N       = merged_data$n,
  MAF     = merged_data$MAF.y,
  snp     = merged_data$SNP
)

pQTL_dataset <- list(
  beta    = merged_data$beta.x,
  varbeta = merged_data$se.x^2,
  type    = "quant",
  N       = merged_data$N,
  MAF     = merged_data$MAF.x,
  snp     = merged_data$SNP
)

# ---- 2.6 Run coloc.abf ----
result_abf <- coloc.abf(
  dataset1 = gwas_dataset,
  dataset2 = pQTL_dataset,
  p1  = 1e-4,
  p2  = 1e-4,
  p12 = 1e-5
)

post_df <- data.frame(
  Hypothesis = c("H0 (no assoc)", "H1 (GWAS only)", "H2 (pQTL only)",
                 "H3 (both, independent)", "H4 (colocalized)"),
  Probability = c(result_abf$summary["PP.H0.abf"],
                  result_abf$summary["PP.H1.abf"],
                  result_abf$summary["PP.H2.abf"],
                  result_abf$summary["PP.H3.abf"],
                  result_abf$summary["PP.H4.abf"])
)

write.csv(post_df,
          "G:/FW_analysis/protein_gwas/coloc/FW_total_PEAR1_UKB_coloc_500kb.csv",
          row.names = FALSE)

# Identify the shared causal variant (top SNP.PP.H4)
idx     <- which.max(result_abf$results$SNP.PP.H4)
top_snp <- result_abf$results$snp[idx]
print(top_snp)


# =============================================================
# 3. LocusZoom plot: PEAR1 + FW (projection and total)
# =============================================================

# ---- 3.1 Load cis-pQTL and FW GWAS ----
pQTL       <- read.csv("G:/FW_analysis/protein_gwas/PEAR1/PEAR1_UKB_cis_pqtl_500kb.csv",
                       header = TRUE)
gwas_data1 <- fread("G:/R/MR/GWAS/FW_gwas/ukb_FW_Projection.csv.gz", header = TRUE)
gwas_data2 <- fread("G:/R/MR/GWAS/FW_gwas/ukb_FW.csv.gz",            header = TRUE)

chr_id       <- 1
region_start <- 156793656
region_end   <- 157016432

gwas_region1 <- gwas_data1 %>% filter(CHR == chr_id, BP >= region_start, BP <= region_end)
gwas_region2 <- gwas_data2 %>% filter(CHR == chr_id, BP >= region_start, BP <= region_end)
pQTL         <- pQTL       %>% filter(CHR == chr_id, BP >= region_start, BP <= region_end)

# ---- 3.2 Lead SNPs ----
gwas_lead_snp1 <- gwas_region1$SNP[which.min(gwas_region1$pval)]
gwas_lead_snp2 <- gwas_region2$SNP[which.min(gwas_region2$pval)]
pQTL_lead_snp  <- pQTL$SNP[which.min(pQTL$pval)]

# ---- 3.3 Standardize column names ----
gwas_region1$type <- "FWprojection_GWAS"
gwas_region2$type <- "FWtotal_GWAS"
pQTL$type         <- "pQTL"

colnames(pQTL)[colnames(pQTL) == "Pval"] <- "pval"

# ---- 3.4 LD reference files ----
ld_gwas <- fread("G:/Plink/FW_projection_LD_r2/rs4661077_LD_full.ld") %>%
  dplyr::select(SNP = SNP_B, ld_gwas = R2)

ld_pQTL <- fread("G:/Plink/FW_projection_LD_r2/rs4661012_LD_full.ld") %>%
  dplyr::select(SNP = SNP_B, ld_pQTL = R2)

# ---- 3.5 Combine plotting data ----
plot_dat <- bind_rows(
  gwas_region1 %>% dplyr::select(SNP, CHR, BP, beta, se, pval, type),
  gwas_region2 %>% dplyr::select(SNP, CHR, BP, beta, se, pval, type),
  pQTL         %>% dplyr::select(SNP, CHR, BP, beta, se, pval, type)
) %>%
  mutate(logP = -log10(pval))

# ---- 3.6 Attach LD per facet ----
dat_gwas1 <- plot_dat %>%
  filter(type == "FWprojection_GWAS") %>%
  left_join(ld_gwas, by = "SNP") %>%
  mutate(ld = ifelse(is.na(ld_gwas), 0, ld_gwas))

dat_gwas2 <- plot_dat %>%
  filter(type == "FWtotal_GWAS") %>%
  left_join(ld_gwas, by = "SNP") %>%
  mutate(ld = ifelse(is.na(ld_gwas), 0, ld_gwas))

dat_pQTL <- plot_dat %>%
  filter(type == "pQTL") %>%
  left_join(ld_pQTL, by = "SNP") %>%
  mutate(ld = ifelse(is.na(ld_pQTL), 0, ld_pQTL))

plot_dat_final <- bind_rows(dat_gwas1, dat_gwas2, dat_pQTL)

# ---- 3.7 Main LocusZoom plot ----
p_main <- ggplot(plot_dat_final, aes(x = BP, y = logP)) +
  geom_point(aes(fill = ld), shape = 21, size = 1.8, stroke = 0) +
  scale_fill_gradientn(
    colors = c("darkblue", "deepskyblue", "cyan", "green", "gold", "red"),
    limits = c(0, 1),
    name   = bquote(LD ~ (r^2))
  ) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "gray40") +
  geom_point(
    data = subset(plot_dat_final, SNP == gwas_lead_snp1 & type == "FWprojection_GWAS"),
    shape = 18, color = "firebrick", size = 5, stroke = 1
  ) +
  geom_point(
    data = subset(plot_dat_final, SNP == gwas_lead_snp2 & type == "FWtotal_GWAS"),
    shape = 18, color = "firebrick", size = 5, stroke = 1
  ) +
  geom_point(
    data = subset(plot_dat_final, SNP == pQTL_lead_snp & type == "pQTL"),
    shape = 18, color = "darkviolet", size = 5, stroke = 1
  ) +
  geom_text_repel(
    data = subset(plot_dat_final,
                  SNP %in% c(gwas_lead_snp1, gwas_lead_snp2, pQTL_lead_snp)),
    aes(label = SNP), size = 4, fontface = "bold", nudge_y = 2
  ) +
  facet_grid(type ~ ., scales = "free_y") +
  scale_x_continuous(labels = function(x) paste0(x / 1e6, " Mb"),
                     expand = c(0.01, 0)) +
  theme_bw() +
  theme(
    panel.grid     = element_blank(),
    strip.text     = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "gray90"),
    legend.position = "right",
    axis.text.x    = element_text(size = 6),
    axis.text.y    = element_text(size = 6)
  )

print(p_main)

ggsave(
  filename = "G:/FW_analysis/protein_gwas/coloc/PEAR1_UKB_FW_LocusZoom_plot_500kb.pdf",
  plot     = p_main,
  width    = 10, height = 10, units = "cm", dpi = 300
)

# ---- 3.8 Add gene track ----
loc <- locus(
  seqname = chr_id,
  xrange  = c(region_start, region_end),
  ens_db  = EnsDb.Hsapiens.v75
)

p_gene <- gg_genetracks(loc, filter_gene_name = "PEAR1") +
  scale_x_continuous(labels = function(x) paste0(x, " Mb"),
                     expand = c(0.01, 0)) +
  theme(axis.text.x = element_text(size = 6))

final_plot <- p_main / p_gene + plot_layout(heights = c(10, 1))
print(final_plot)

ggsave(
  filename = "G:/FW_analysis/protein_gwas/coloc/PEAR1_UKB_FW_LocusZoom_plot_gene_500kb.pdf",
  plot     = final_plot,
  width    = 10, height = 16, units = "cm", dpi = 300
)