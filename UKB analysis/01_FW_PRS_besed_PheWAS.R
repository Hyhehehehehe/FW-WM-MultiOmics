# FW PheWAS: Association Analysis, Manhattan Plot, and Volcano Plot
# This script performs a phenome-wide association study (PheWAS) of the four
# FW-WM indexes (FW, FW_Association, FW_Commissural, FW_Projection) against
# ICD-10 diseases in UKB, and visualizes the results as Manhattan and volcano plots
# Workflow:
# 1. PheWAS: logistic regression of each ICD-10 disease on standardized FW indexes
#    (one result file per FW index and ICD-10 chapter)
# 2. Manhattan plots (one panel per FW index)
# 3. Volcano plot based on OR and P value across all ICD-10 chapters

#=============================================================================
# 1. PheWAS regression
#=============================================================================

library(readstata13)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(tibble)

# ---- 1.1 Load data ----
data0 <- read.dta13("FW_disease.dta")

# Define categorical variables
data0$sex <- factor(data0$sex)
for (col in 11:636) {
  data0[, col] <- factor(data0[, col])
}

# ---- 1.2 Standardize FW indexes ----
exposures <- c("FW", "FW_Association", "FW_Commissural", "FW_Projection")
data_std  <- data0

for (exp in exposures) {
  mean_val <- mean(data0[[exp]], na.rm = TRUE)
  sd_val   <- sd(data0[[exp]],   na.rm = TRUE)
  data_std[[paste0(exp, "_std")]] <- (data0[[exp]] - mean_val) / sd_val
}

exposures_std <- paste0(exposures, "_std")

# ---- 1.3 ICD-10 chapter list ----
icd_cats <- c("A_B", "D", "E", "F", "G", "H", "I", "J",
              "K", "L", "M", "N", "O", "P", "Q")

# ---- 1.4 Loop: write one result file per FW index and ICD-10 chapter ----
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
          se      <- coef_data[exp, "Std. Error"]
          p_value <- coef_data[exp, "Pr(>|z|)"]
          or      <- exp(beta)
          or_ci   <- exp(conf_data[exp, ])
          
          results <- rbind(results, data.frame(
            Disease       = disease,
            Exposure      = exp,
            Beta          = beta,
            OR            = or,
            OR_95CI_Lower = or_ci[1],
            OR_95CI_Upper = or_ci[2],
            P_value       = p_value,
            stringsAsFactors = FALSE
          ))
        }
      }, error = function(e) NULL)
    }
    
    out_file <- file.path(
      output_dir,
      paste0("ukb_phewas_", fw_name, "_ICD_", icd_cat, ".csv")
    )
    write.csv(results, file = out_file, row.names = FALSE)
    
    cat("Written:", out_file, "\n")
  }
}

#=============================================================================
# 2. Manhattan plots
#=============================================================================

# ---- 2.1 Color palette and ICD-10 chapter labels ----
color_palette <- c(
  "#1F77B4", "#FF7F0E", "#C49C94", "#D62728",
  "#9467BD", "#8C564B", "#E377C2", "#FFBB78",
  "#BCBD22", "#17BECF", "#FFCC33", "#AEC7E8",
  "#98DF8A", "#FF9896", "#C5B0D5", "#2CA02C"
)

cat_labels <- c(
  "A" = "ICD-10 A: Certain infectious and parasitic diseases",
  "B" = "ICD-10 B: Certain infectious and parasitic diseases",
  "D" = "ICD-10 D: Blood, blood-forming organs and certain immune disorders",
  "E" = "ICD-10 E: Endocrine, nutritional and metabolic diseases",
  "F" = "ICD-10 F: Mental and behavioural disorders",
  "G" = "ICD-10 G: Nervous system disorders",
  "H" = "ICD-10 H: Eye and adnexa disorders; Ear and mastoid process disorders",
  "I" = "ICD-10 I: Circulatory system disorders",
  "J" = "ICD-10 J: Respiratory system disorders",
  "K" = "ICD-10 K: Digestive system disorders",
  "L" = "ICD-10 L: Skin and subcutaneous tissue disorders",
  "M" = "ICD-10 M: Musculoskeletal system and connective tissue disorders",
  "N" = "ICD-10 N: Genitourinary system disorders",
  "O" = "ICD-10 O: Pregnancy, childbirth and the puerperium",
  "P" = "ICD-10 P: Certain conditions originating in the perinatal period",
  "Q" = "ICD-10 Q: Congenital disruptions and chromosomal abnormalities"
)

# ---- 2.2 Helper: draw Manhattan plot for one FW index ----
plot_phewas_manhattan <- function(fw_name, icd_cats, input_dir, out_dir) {
  
  dat <- NULL
  for (icd_cat in icd_cats) {
    f <- file.path(
      input_dir,
      paste0("ukb_phewas_", fw_name, "_ICD_", icd_cat, ".csv")
    )
    if (!file.exists(f)) next
    dat <- rbind(dat, fread(f), fill = TRUE)
  }
  dat <- data.frame(dat)
  
  if (nrow(dat) == 0) {
    cat("No valid rows for", fw_name, "- skipping Manhattan plot\n")
    return(invisible(NULL))
  }
  
  dat$cat <- gsub("ICD-10|[0-9]", "", dat$Disease)
  dat <- na.omit(dat)
  
  dat_man <- dat[, c("cat", "Disease", "P_value")]
  colnames(dat_man) <- c("Cat", "Disease", "P")
  dat_man <- dat_man %>%
    group_by(Cat) %>%
    arrange(P, .by_group = TRUE) %>%
    ungroup() %>%
    mutate(Pos = row_number())
  
  P_b <- 0.05 / nrow(dat_man)
  
  dat_mark <- dat_man %>%
    filter(P < P_b) %>%
    arrange(P) %>%
    slice_head(n = 10)
  
  X_axis <- dat_man %>%
    group_by(Cat) %>%
    summarize(center = (max(Pos) + min(Pos)) / 2, .groups = "drop")
  
  dat_man$Cat <- trimws(dat_man$Cat)
  unique_cats <- unique(dat_man$Cat)
  color_mapping <- setNames(color_palette[seq_along(unique_cats)], unique_cats)
  
  p <- ggplot(dat_man, aes(x = Pos, y = -log10(P),
                           color = Cat, size = -log10(P))) +
    geom_point() +
    scale_size_continuous(range = c(1, 4)) +
    scale_x_continuous(label = X_axis$Cat, breaks = X_axis$center,
                       expand = c(0, 0)) +
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(0, 150),
      breaks = seq(0, 150, by = 10)
    ) +
    geom_hline(yintercept = -log10(P_b),
               color = "black", linewidth = 0.4,
               linetype = 2, alpha = 0.75) +
    labs(x = "ICD code", y = expression(-log[10](italic(P)))) +
    scale_color_manual(name = NULL, labels = cat_labels,
                       values = color_mapping) +
    theme(
      panel.background  = element_rect(fill = "white", color = "white"),
      panel.border      = element_rect(color = "white", fill = NA, linewidth = 0.5),
      axis.line.y       = element_line(),
      axis.line.x       = element_line(),
      axis.text.y.right = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.line.y.right  = element_blank(),
      panel.grid        = element_blank(),
      legend.position   = "top",
      axis.text         = element_text(size = 8)
    ) +
    guides(color = guide_legend(nrow = 8, keywidth = 2)) +
    geom_text_repel(
      data = dat_mark,
      aes(x = Pos, y = -log10(P), label = Disease),
      color = "black",
      size = 3,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "grey50",
      max.overlaps = 20
    )
  
  ggsave(file.path(out_dir,
                   paste0("phewas_", fw_name, "_manhattan_plot.pdf")),
         plot = p, width = 25, height = 20, units = "cm")
}

# ---- 2.3 Draw Manhattan plot for each FW index ----
out_dir_plot <- "phewas_manhattan"
if (!dir.exists(out_dir_plot)) dir.create(out_dir_plot, recursive = TRUE)

for (fw in c("fw", "fw_Association", "fw_Commissural", "fw_Projection")) {
  plot_phewas_manhattan(fw, icd_cats, output_dir, out_dir_plot)
}

#=============================================================================
# 3. Volcano plot
#=============================================================================
# This section reads all per-chapter PheWAS result files for FW_total and
# combines them into a single table for the volcano plot.
# Output: Fw_volcano_OR.pdf

# ---- 3.1 Read and combine all per-chapter results ----
dat <- NULL
for (icd_cat in icd_cats) {
  f <- file.path(output_dir,
                 paste0("ukb_phewas_fw_ICD_", icd_cat, ".csv"))
  if (!file.exists(f)) next
  dat <- rbind(dat, fread(f))
}

dat <- as.data.frame(na.omit(dat))
dat$OR   <- exp(dat$Beta)
dat$cat  <- gsub("ICD-10|[0-9]", "", dat$Disease) %>% str_trim(side = "both")
dat <- dat %>%
  mutate(
    sig   = ifelse(P_value < 5.95e-05, "sig", "non_sig"),
    group = cat
  )

# ---- 3.2 Strip and background ranges ----
strip_ymin <- 0.98
strip_ymax <- 1.02

lower_dat <- dat %>% filter(OR < strip_ymin)
upper_dat <- dat %>% filter(OR > strip_ymax)

range_dat <- dat %>%
  group_by(group) %>%
  summarise(
    min_or = min(OR, na.rm = TRUE),
    max_or = max(OR, na.rm = TRUE),
    .groups = "drop"
  )

cat_order <- c("A", "B", "D", "E", "F", "G", "H", "I",
               "J", "K", "L", "M", "N", "O", "P", "Q")

custom_colors <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728",
  "#9467BD", "#8C564B", "#E377C2", "#FFBB78",
  "#BCBD22", "#17BECF", "#FFCC33", "#AEC7E8",
  "#98DF8A", "#FF9896", "#C5B0D5", "#C49C94"
)

bg_lower <- range_dat %>%
  mutate(
    x    = as.integer(factor(group, levels = cat_order)),
    xmin = x - 0.4,
    xmax = x + 0.4,
    ymin = min_or,
    ymax = strip_ymin,
    fill = "lightgray",
    alpha = 0.6
  )

bg_upper <- range_dat %>%
  mutate(
    x    = as.integer(factor(group, levels = cat_order)),
    xmin = x - 0.4,
    xmax = x + 0.4,
    ymin = strip_ymax,
    ymax = max_or,
    fill = "lightgray",
    alpha = 0.6
  )

strip_data <- tibble(
  group = cat_order,
  x     = as.integer(factor(group, levels = cat_order)),
  xmin  = x - 0.5,
  xmax  = x + 0.5,
  ymin  = strip_ymin,
  ymax  = strip_ymax,
  fill  = custom_colors
)

# ---- 3.3 Volcano plot ----
p <- ggplot() +
  geom_rect(data = bg_lower,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = bg_lower$fill, alpha = bg_lower$alpha,
            inherit.aes = FALSE) +
  geom_rect(data = bg_upper,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = bg_upper$fill, alpha = bg_upper$alpha,
            inherit.aes = FALSE) +
  geom_rect(data = strip_data,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = strip_data$fill, alpha = 0.8,
            inherit.aes = FALSE) +
  geom_jitter(data = lower_dat,
              aes(x = group, y = OR, color = sig, size = OR),
              width = 0.3, alpha = 1) +
  geom_jitter(data = upper_dat,
              aes(x = group, y = OR, color = sig, size = OR),
              width = 0.3, alpha = 1) +
  geom_text(data = strip_data,
            aes(x = x, y = (strip_ymin + strip_ymax) / 2, label = group),
            color = "white", size = 3.5, fontface = "bold") +
  scale_x_discrete(limits = cat_order) +
  scale_color_manual(values = c("sig" = "#e2492e", "non_sig" = "#618fc0")) +
  scale_size_continuous(range = c(1, 4)) +
  theme_minimal() +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "right",
    plot.title       = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x     = element_blank(),
    axis.title.y     = element_text(size = 13, color = "black", face = "bold")
  ) +
  labs(
    y     = "Odds Ratio (OR)",
    title = "OR of Free Water PRS Associated with Diseases by System",
    color = "P-value Significance",
    size  = "OR Magnitude"
  )

print(p)

ggsave(p,
       file = "Fw_volcano_OR.pdf",
       width = 25, height = 20, units = "cm")