# FW PheWAS OR Correlation (Motion Adjustment)
# This script evaluates the robustness of the imaging-derived PheWAS results to
# head motion by comparing the OR estimates before and after additional
# adjustment for a diffusion MRI head-motion metric
# Workflow:
# 1. Read and standardize the original and motion-adjusted PheWAS result files
# 2. Match the two result sets and calculate the Pearson correlation of the ORs
# 3. Generate four aligned scatter panels and export TIFF and PDF figures

packages <- c(
  "readr", "dplyr", "stringr", "ggplot2", "gridExtra", "grid", "gtable"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(grid)
library(gtable)


## ============================================================
## 1. Set the paths
## ============================================================

base_dir <- "."

file_original <- file.path(base_dir, "FW_disease.csv")
file_adjust_motion <- file.path(base_dir, "FW_disease_adjust_motion.csv")

out_dir <- file.path(
  base_dir,
  "fw_disease_adjust_motion_comparison_outputs",
  "FW_vs_adjust_motion"
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


## ============================================================
## 2. Helper functions
## ============================================================

pick_col <- function(dat, candidates, required = TRUE) {

  hit <- intersect(candidates, colnames(dat))

  if (length(hit) == 0) {
    if (required) {
      stop("Column not found: ", paste(candidates, collapse = " / "))
    } else {
      return(NA_character_)
    }
  }

  hit[1]
}

to_numeric <- function(x) {
  suppressWarnings(as.numeric(gsub("^<", "", as.character(x))))
}

normalize_exposure <- function(x) {

  x %>%
    as.character() %>%
    str_remove("_std$") %>%
    str_remove("_e$") %>%
    str_remove("_v$") %>%
    str_remove("_adjust_motion$") %>%
    str_remove("_motion$") %>%
    str_replace("^FW_2$", "FW")
}

format_p_value <- function(p) {

  if (is.na(p)) {
    return("P = NA")
  }

  if (p < 0.001) {
    return("P < 0.001")
  }

  paste0("P = ", formatC(p, format = "f", digits = 3))
}


## ============================================================
## 3. Read and standardize the PheWAS files
## ============================================================

read_phewas <- function(file) {

  if (!file.exists(file)) {
    stop("File does not exist: ", file)
  }

  dat_raw <- read_csv(
    file,
    show_col_types = FALSE,
    guess_max = 100000
  )

  colnames(dat_raw) <- colnames(dat_raw) %>%
    str_replace("^\ufeff", "") %>%
    str_trim()

  disease_col <- pick_col(
    dat_raw,
    c("Disease", "Phenotype", "ICD.code", "ICD_code")
  )

  exposure_col <- pick_col(
    dat_raw,
    c("Exposure")
  )

  p_col <- pick_col(
    dat_raw,
    c("P_value", "P", "p")
  )

  or_col <- pick_col(
    dat_raw,
    c("OR", "or")
  )

  dat_raw %>%
    transmute(
      Disease = as.character(.data[[disease_col]]),
      Exposure = as.character(.data[[exposure_col]]),
      Group = normalize_exposure(.data[[exposure_col]]),
      OR = to_numeric(.data[[or_col]]),
      P = to_numeric(.data[[p_col]])
    ) %>%
    filter(
      !is.na(Disease),
      !is.na(Group),
      !is.na(OR),
      !is.na(P)
    ) %>%
    group_by(Group, Disease) %>%
    arrange(P, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()
}


## ============================================================
## 4. Read and match the two files
## ============================================================

dat_original <- read_phewas(file_original) %>%
  select(
    Group,
    Disease,
    Original_OR = OR
  )

dat_adjust_motion <- read_phewas(file_adjust_motion) %>%
  select(
    Group,
    Disease,
    Adjust_motion_OR = OR
  )

plot_df <- inner_join(
  dat_original,
  dat_adjust_motion,
  by = c("Group", "Disease")
) %>%
  filter(
    is.finite(Original_OR),
    is.finite(Adjust_motion_OR)
  )

if (nrow(plot_df) < 3) {
  stop("Fewer than 3 matched rows are available for the correlation analysis.")
}


## ============================================================
## 5. Calculate the Pearson correlation of the OR values across groups
## ============================================================

cor_results <- plot_df %>%
  group_by(Group) %>%
  group_modify(
    ~ {
      dat <- .x %>%
        filter(
          is.finite(Original_OR),
          is.finite(Adjust_motion_OR)
        )

      if (nrow(dat) < 3 ||
          sd(dat$Original_OR) == 0 ||
          sd(dat$Adjust_motion_OR) == 0) {
        return(
          tibble(
            r = NA_real_,
            P = NA_real_
          )
        )
      }

      test_result <- cor.test(
        dat$Original_OR,
        dat$Adjust_motion_OR,
        method = "pearson"
      )

      tibble(
        r = unname(test_result$estimate),
        P = test_result$p.value
      )
    }
  ) %>%
  ungroup() %>%
  mutate(
    label = paste0(
      "r = ",
      ifelse(is.na(r), "NA", sprintf("%.3f", r)),
      "\n",
      vapply(P, format_p_value, character(1))
    )
  )


## ============================================================
## 6. Unify the axis ranges and ticks
## ============================================================

common_range <- range(
  c(plot_df$Original_OR, plot_df$Adjust_motion_OR),
  na.rm = TRUE
)

range_width <- diff(common_range)

if (!is.finite(range_width) || range_width == 0) {
  range_width <- 0.01
}

range_padding <- range_width * 0.05

common_limits <- c(
  common_range[1] - range_padding,
  common_range[2] + range_padding
)

common_breaks <- pretty(common_limits, n = 5)
common_breaks <- common_breaks[
  common_breaks >= common_limits[1] &
    common_breaks <= common_limits[2]
]

preferred_order <- c(
  "FW",
  "FW_Association",
  "FW_Commissural",
  "FW_Projection"
)

group_order <- c(
  preferred_order[preferred_order %in% unique(plot_df$Group)],
  setdiff(unique(plot_df$Group), preferred_order)
)


## ============================================================
## 7. Function for plotting a single panel
## ============================================================

make_panel <- function(group_name) {

  panel_data <- plot_df %>%
    filter(Group == group_name)

  panel_cor <- cor_results %>%
    filter(Group == group_name)

  annotation_x <- common_limits[1] + 0.02 * diff(common_limits)
  annotation_y <- common_limits[2] - 0.02 * diff(common_limits)

  ggplot(
    panel_data,
    aes(
      x = Original_OR,
      y = Adjust_motion_OR
    )
  ) +
    geom_point(
      alpha = 0.45,
      size = 1.4
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.7
    ) +
    annotate(
      "text",
      x = annotation_x,
      y = annotation_y,
      label = panel_cor$label,
      hjust = 0,
      vjust = 1,
      size = 5
    ) +
    scale_x_continuous(
      limits = common_limits,
      breaks = common_breaks,
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      limits = common_limits,
      breaks = common_breaks,
      expand = expansion(mult = 0)
    ) +
    coord_fixed(ratio = 1) +
    labs(
      title = group_name,
      x = "OR in original model",
      y = "OR in head motion-adjusted model"
    ) +
    theme_classic(base_size = 15) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 18,
        margin = margin(b = 10)
      ),
      axis.title.x = element_text(
        size = 14,
        margin = margin(t = 8)
      ),
      axis.title.y = element_text(
        size = 14,
        margin = margin(r = 8)
      ),
      axis.text = element_text(
        size = 12,
        color = "black"
      ),
      axis.line = element_line(linewidth = 0.7),
      axis.ticks = element_line(linewidth = 0.7),
      plot.margin = margin(
        t = 12,
        r = 12,
        b = 12,
        l = 12
      )
    )
}


## ============================================================
## 8. Generate the 4 panels and force alignment
## ============================================================

panel_list <- lapply(group_order, make_panel)

# Convert to grob
grob_list <- lapply(panel_list, ggplotGrob)

# Unify the column widths and row heights of all panels to ensure strict alignment of the four panels
max_width <- do.call(unit.pmax, lapply(grob_list, function(g) g$widths))
max_height <- do.call(unit.pmax, lapply(grob_list, function(g) g$heights))

grob_list_aligned <- lapply(
  grob_list,
  function(g) {
    g$widths <- max_width
    g$heights <- max_height
    g
  }
)

p_or <- arrangeGrob(
  grobs = grob_list_aligned,
  ncol = 2,
  nrow = 2
)


## ============================================================
## 9. Export TIFF and PDF
## ============================================================

or_tiff <- file.path(
  out_dir,
  "FW_vs_adjust_motion_OR_Pearson_correlation_aligned.tiff"
)

or_pdf <- file.path(
  out_dir,
  "FW_vs_adjust_motion_OR_Pearson_correlation_aligned.pdf"
)

tiff(
  filename = or_tiff,
  width = 12,
  height = 12,
  units = "in",
  res = 600,
  compression = "lzw"
)
grid.newpage()
grid.draw(p_or)
dev.off()

pdf(
  file = or_pdf,
  width = 12,
  height = 12,
  useDingbats = FALSE
)
grid.newpage()
grid.draw(p_or)
dev.off()

message("Aligned TIFF saved to: ", or_tiff)
message("Aligned PDF saved to: ", or_pdf)
