# =============================================================
# 01_Covariates_Matching.R
# Match blood pressure, BMI, demographics, and APOE genotype
# =============================================================

library(readr)
library(readxl)
library(readstata13)
library(dplyr)
library(tidyr)

order_cols <- c("PHASE_MRI", "PTID", "RID_MRI", "VISCODE", "VISCODE2_MRI", "EXAMDATE",
                "Group", "Age", "Sex", "FW_total", "FW_asso", "FW_comm", "FW_proj",
                "FIELD_STRENGTH_MRI", "MMCONDCT_MRI", "MMREASON_MRI", "SITEID_MRI")

# ---- 1.1 Blood pressure and BMI ----
BP_BMI  <- read.dta13("Vital Signs_ADNI1,GO,2,3,4.dta")
FW_data <- read.csv("FW_information.csv", header = TRUE)

# Fill in unmatched records
miss_FW_data <- FW_data %>%
  filter(is.na(VISCODE2_MRI)) %>%
  select(PTID, VISCODE, EXAMDATE, Group, Age, Sex,
         FW_total, FW_asso, FW_comm, FW_proj, FIELD_STRENGTH_MRI,
         MMCONDCT_MRI, MMREASON_MRI, SITEID_MRI, PHASE_MRI, RID_MRI, VISCODE2_MRI) %>%
  left_join(BP_BMI, by = c("PTID", "VISCODE")) %>%
  select(all_of(order_cols))

FW_data <- bind_rows(FW_data %>% filter(!is.na(VISCODE2_MRI)), miss_FW_data)
write.csv(FW_data, "FW_information2.csv", row.names = FALSE)

# ---- 1.2 Recompute BMI (single height per subject) ----
FW_data <- read.csv("FW_information2.csv", header = TRUE)
BP_BMI  <- read.dta13("Vital Signs_ADNI1,GO,2,3,4.dta") %>%
  select(1:7, 9, 13, 11, 12) %>%
  group_by(PTID) %>%
  mutate(height = if (all(is.na(height))) NA else first(height[!is.na(height)])) %>%
  ungroup() %>%
  mutate(BMI    = round(weight / height^2, 2),
         weight = round(weight, 2),
         height = round(height, 2))

names(BP_BMI) <- paste0(names(BP_BMI), "_signs")

# ---- 1.3 Merge FW with BP/BMI ----
FW_data$VISCODE_ST <- ifelse(FW_data$VISCODE == "scmri", "sc", FW_data$VISCODE)

merged_FW_data <- left_join(
  FW_data, BP_BMI,
  by = c("PTID" = "PTID_signs", "PHASE_MRI" = "PHASE_signs",
         "VISCODE_ST" = "VISCODE_signs")
)

# Fill unmatched via VISCODE2
miss_merged_FW_data <- merged_FW_data %>%
  filter(is.na(RID_signs)) %>%
  select(-c(19:26)) %>%
  mutate(VISCODE2_MRI = ifelse(VISCODE2_MRI == "scmri", "sc", VISCODE2_MRI))

merged_FW_data_miss <- left_join(
  miss_merged_FW_data, BP_BMI,
  by = c("PTID" = "PTID_signs", "PHASE_MRI" = "PHASE_signs",
         "VISCODE2_MRI" = "VISCODE2_signs")
) %>% mutate(VISCODE2_signs = VISCODE2_MRI) %>% select(all_of(order_cols))

merged_FW_data2 <- bind_rows(
  merged_FW_data %>% filter(!is.na(RID_signs)) %>% select(all_of(order_cols)),
  merged_FW_data_miss
)

write.csv(merged_FW_data2, "FW_BP_BMI.csv", row.names = FALSE)

# ---- 1.4 Demographics (education and ethnicity) ----
merged_FW_data <- read.csv("FW_BP_BMI.csv", header = TRUE)
Demographics <- read.csv("Subject Demographics [ADNI1,GO,2,3,4].csv", header = TRUE) %>%
  select(1:6, 13, 24, 25, 78) %>%
  group_by(PTID) %>%
  mutate(PTEDUCAT = max(PTEDUCAT, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(PTEDUCAT = ifelse(PTEDUCAT == "-Inf", NA, PTEDUCAT),
         VISDATE   = as.Date(VISDATE)) %>%
  group_by(PTID) %>%
  arrange(VISDATE, .by_group = TRUE) %>%
  mutate(PTETHCAT = {
    valid <- PTETHCAT[!PTETHCAT %in% c(-4) & !is.na(PTETHCAT)]
    if (length(valid) > 0) first(valid) else first(PTETHCAT)
  }) %>%
  ungroup() %>%
  select(PTID, PTEDUCAT, PTETHCAT, PTRACCAT)

Demographics_unique <- Demographics %>%
  group_by(PTID) %>% slice(1) %>% ungroup()
names(Demographics_unique) <- paste0(names(Demographics_unique), "_demo")

merged_FW_data1 <- left_join(
  merged_FW_data, Demographics_unique, by = c("PTID" = "PTID_demo")
) %>% select(all_of(order_cols))

write.csv(merged_FW_data1, "FW_covariable.csv", row.names = FALSE)

# ---- 1.5 APOE genotype ----
FW_covariable <- read.csv("FW_covariable.csv", header = TRUE)
APOE <- read.csv("ApoE Genotyping Results [ADNI1,GO,2,3,4].csv", header = TRUE) %>%
  select(PTID, GENOTYPE)

merge_FW_covariable <- left_join(FW_covariable, APOE, by = "PTID")
write.csv(merge_FW_covariable, "FW_covariable_APOE.csv", row.names = FALSE)