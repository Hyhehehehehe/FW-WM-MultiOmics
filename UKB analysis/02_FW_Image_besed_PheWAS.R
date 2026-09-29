# =============================================================
# FW PheWAS: Association Analysis and Heatmap Visualization
# Workflow:
# 1. PheWAS: logistic regression of each disease on standardized FW indexes
#    (one result file per ICD-10 category per FW index)
# 2. Merge per-category results and restrict to ICD-10 chapters E, F, G, I
# 3. Draw heatmaps: full heatmap + 4 split panels by row range
# =============================================================

library(readstata13)
library(data.table)
library(dplyr)
library(readxl)
library(pheatmap)

# =============================================================
# Part 1: PheWAS regression
# =============================================================

data0 <- read.dta13("FW_disease.dta")

data0$sex <- factor(data0$sex)
for (col in 11:636) {
  data0[, col] <- factor(data0[, col])
}

exposures <- c("FW", "FW_Association", "FW_Commissural", "FW_Projection")
data_std  <- data0

for (exp in exposures) {
  mean_val <- mean(data0[[exp]], na.rm = TRUE)
  sd_val   <- sd(data0[[exp]],   na.rm = TRUE)
  data_std[[paste0(exp, "_std")]] <- (data0[[exp]] - mean_val) / sd_val
}

exposures_std <- paste0(exposures, "_std")

icd_cats <- c("A_B", "D", "E", "F", "G", "H", "I", "J",
              "K", "L", "M", "N", "O", "P", "Q")

output_dir <- "phewas_results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

for (exp in exposures_std) {
  
  fw_name <- sub("_std$", "", exp)
  
  for (icd_cat in icd_cats) {
    
    results <- data.frame()
    
    for (i in 11:636) {
      disease <- colnames(data0)[i]
      if (!grepl(paste0("^", icd_cat), disease)) next
      
      tryCatch({
        formula <- as.formula(paste(
          disease, "~", exp,
          "+ sex + Age_attending_2 + Townsend_index_0"
        ))
        
        fit <- glm(formula, data = data_std,
                   family = binomial(link = "logit"))
        
        coef_data <- coef(summary(fit))
        conf_data <- confint(fit)
        
        if (exp %in% rownames(coef_data)) {
          beta    <- coef_data[exp, "Estimate"]
          p_value <- coef_data[exp, "Pr(>|z|)"]
          or      <- exp(beta)
          or_ci   <- exp(conf_data[exp, ])
          
          results <- rbind(results, data.frame(
            Disease       = disease,
            Category      = icd_cat,
            Exposure      = exp,
            OR            = or,
            OR_95CI_Lower = or_ci[1],
            OR_95CI_Upper = or_ci[2],
            P_value       = p_value,
            stringsAsFactors = FALSE
          ))
        }
      }, error = function(e) NULL)
    }
    
    out_file <- file.path(output_dir,
                          paste0("phewas_", fw_name, "_ICD_", icd_cat, ".csv"))
    write.csv(results, file = out_file, row.names = FALSE)
    
    cat("Written:", out_file, "\n")
  }
}

# =============================================================
# Part 2: Build OR / P matrices for the four FW indexes
# =============================================================

fw_list <- c("total", "association", "commissural", "projection")
fw_file_suffix <- c("FW", "FW_Association", "FW_Commissural", "FW_Projection")

allowed_prefixes <- c("E", "F", "G", "I")

d_OR <- NULL
d_P  <- NULL

for (k in seq_along(fw_list)) {
  
  fw_label  <- fw_list[k]
  fw_suffix <- fw_file_suffix[k]
  
  dat_all <- NULL
  for (icd_cat in icd_cats) {
    f <- file.path(output_dir,
                   paste0("phewas_", fw_suffix, "_ICD_", icd_cat, ".csv"))
    if (!file.exists(f)) next
    dat_all <- rbind(dat_all, fread(f))
  }
  dat_all <- data.frame(dat_all)
  
  # Restrict to ICD-10 chapters E, F, G, I
  dat_all$FirstLetter <- toupper(substr(dat_all$Disease, 1, 1))
  dat_all <- subset(dat_all, FirstLetter %in% allowed_prefixes)
  
  dat_all$Disease <- sub("_2$", "", dat_all$Disease)
  
  # Rename Category to a readable label (Chapter)
  dat_all$Category <- dat_all$FirstLetter
  
  if (k == 1) {
    d_OR <- dat_all[, c("Disease", "Category", "OR",       "P_value")]
    d_P  <- dat_all[, c("Disease", "Category", "P_value")]
  } else {
    d_OR <- merge(d_OR, dat_all[, c("Disease", "Category", "OR",      "P_value")],
                  by = c("Disease", "Category"))
    d_P  <- merge(d_P,  dat_all[, c("Disease", "Category", "P_value")],
                  by = c("Disease", "Category"))
  }
}

# ---- Build OR matrix ----
colnames(d_OR) <- c("Disease", "Category",
                    "total", "P_total",
                    "Association", "P_Association",
                    "Commissural", "P_Commissural",
                    "Projection", "P_Projection")

dat_OR <- d_OR[, c("total", "Association", "Commissural", "Projection")]
dat1 <- as.matrix(dat_OR)
rownames(dat1) <- d_OR$Disease

# ---- Build P-value matrix ----
p_matrix <- as.matrix(d_OR[, c("P_total", "P_Association",
                               "P_Commissural", "P_Projection")])
rownames(p_matrix) <- d_OR$Disease
colnames(p_matrix) <- c("total", "Association", "Commissural", "Projection")

# ---- Annotation ----
anno <- d_OR[, c("Disease", "Category")]
rownames(anno) <- anno[, 1]
anno <- anno[, -1, drop = FALSE]
anno$Category <- gsub(" ", "_", anno$Category)

# =============================================================
# Part 3: Heatmap – full version
# =============================================================

signif_threshold <- 7.99e-05

cc1 <- c("#2266ac", "white", "#b21a2b")

or_min <- 0
or_max <- 2.5
or_white <- 1

breaks_left  <- seq(or_min, or_white, length.out = 104)
breaks_right <- seq(or_white, or_max, length.out = 103)
breaks <- unique(c(breaks_left, breaks_right))

signif_matrix <- matrix(ifelse(p_matrix <= signif_threshold, "*", ""),
                        nrow = nrow(p_matrix),
                        dimnames = dimnames(p_matrix))

ann_colors <- list(
  Category = c(E = "#931635",
               F = "#9467BD",
               G = "#8C564B",
               I = "#FFBB78")
)

pdf("heatmap_full.pdf", width = 6, height = 80)
pheatmap(dat1,
         color = colorRampPalette(cc1)(length(breaks) - 1),
         breaks = breaks,
         border_color = "black",
         cluster_rows = FALSE, cluster_cols = FALSE,
         legend_breaks = c(0, 1, 2.5),
         annotation_row = anno,
         cellwidth = 10, cellheight = 10,
         annotation_colors = ann_colors,
         display_numbers = signif_matrix,
         fontsize_number = 15, number_color = "black",
         angle_col = "90",
         legend_labels = c("0", "1", "2.5"))
dev.off()

# =============================================================
# Part 4: Heatmap – four split panels
# =============================================================

# Load disease abbreviation table
df2 <- read_xlsx("Psychiatric.xlsx", sheet = "Sheet2")[, -3]
df3 <- read_xlsx("Endocrine.xlsx",   sheet = "Sheet2")[, -3]
df4 <- read_xlsx("Nervous.xlsx",     sheet = "Sheet2")[, -3]
df6 <- read_xlsx("Cardiovascular.xlsx", sheet = "Sheet2")[, -3]

colnames(df2) <- c("ICD-10", "Abbreviated_Name", "disease")
colnames(df3) <- c("ICD-10", "Abbreviated_Name", "disease")
colnames(df4) <- c("ICD-10", "Abbreviated_Name", "disease")
colnames(df6) <- c("ICD-10", "Abbreviated_Name", "disease")

df <- rbind(df2, df3, df4, df6)
df <- df[, -2]

full_color_map <- c(
  E = "#279127",
  F = "#931635",
  G = "#9467BD",
  I = "#FFBB78"
)

row_index <- c("1:52", "53:104", "105:156", "157:207")

for (rr in 1:4) {
  
  dat1_1 <- dat1[eval(parse(text = row_index[rr])), , drop = FALSE]
  
  current_icds <- rownames(dat1_1)
  disease_names <- df$disease[match(current_icds, df$`ICD-10`)]
  disease_names[is.na(disease_names)] <- current_icds[is.na(disease_names)]
  rownames(dat1_1) <- disease_names
  
  p_matrix_1 <- p_matrix[eval(parse(text = row_index[rr])), , drop = FALSE]
  rownames(p_matrix_1) <- disease_names
  
  anno_1 <- anno[eval(parse(text = row_index[rr])), , drop = FALSE]
  rownames(anno_1) <- disease_names
  anno_1$Category <- gsub(" ", "_", anno_1$Category)
  
  current_categories <- unique(anno_1$Category)
  current_colors <- full_color_map[names(full_color_map) %in% current_categories]
  ann_colors_rr <- list(Category = current_colors)
  
  signif_matrix_1 <- matrix(ifelse(p_matrix_1 <= signif_threshold, "*", ""),
                            nrow = nrow(p_matrix_1))
  rownames(signif_matrix_1) <- disease_names
  
  pdf(paste0("heatmap", rr, ".pdf"), width = 6, height = 25)
  pheatmap(dat1_1,
           color = colorRampPalette(cc1)(length(breaks) - 1),
           breaks = breaks,
           border_color = "black",
           cluster_rows = FALSE, cluster_cols = FALSE,
           legend_breaks = c(0, 1, 2.5),
           cellwidth = 20, cellheight = 20,
           annotation_row = anno_1,
           annotation_colors = ann_colors_rr,
           display_numbers = signif_matrix_1,
           fontsize_number = 15, number_color = "black",
           angle_col = "90",
           legend_labels = c("0", "1", "2.5"))
  dev.off()
  
  cat("Panel", rr, "saved.\n")
}

cat("All heatmaps completed.\n")