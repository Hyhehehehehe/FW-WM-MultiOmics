# =============================================================
# 02_Imaging_and_Diagnosis_Merging.R
# Merge WMH imaging metrics and derive diagnosis conversion dates
# =============================================================

library(readr)
library(readxl)
library(readstata13)
library(dplyr)
library(tidyr)

# ---- 2.1 Merge WMH metrics ----
merge_FW_covariable <- read.csv("FW_covariable_APOE.csv")
WMH <- read.csv("UCD - White Matter Hyperintensity Volumes_ADNIGO,2,3,4.csv")
names(WMH) <- paste0(names(WMH), "_WMH")

merged_FW_data <- left_join(
  merge_FW_covariable, WMH,
  by = c("PTID" = "PTID_WMH", "PHASE_MRI" = "PHASE_WMH",
         "VISCODE" = "VISCODE_WMH")
)

miss_merged_FW_data <- merged_FW_data %>%
  filter(is.na(RID_WMH)) %>%
  select(-c(27:50))

merged_FW_data_miss <- left_join(
  miss_merged_FW_data, WMH,
  by = c("PTID" = "PTID_WMH", "PHASE_MRI" = "PHASE_WMH",
         "VISCODE2_MRI" = "VISCODE2_WMH")
)

cols_keep <- c("PHASE_MRI", "PTID", "RID_MRI", "VISCODE", "VISCODE2_MRI",
               "Group", "EXAMDATE", "Age", "Sex", "PTEDUCAT_demo",
               "PTETHCAT_demo", "PTRACCAT_demo", "weight_signs", "height_signs",
               "BMI_signs", "VSBPSYS_signs", "VSBPDIA_signs",
               "FW_total", "FW_asso", "FW_comm", "FW_proj",
               "FIELD_STRENGTH_MRI", "MMCONDCT_MRI", "MMREASON_MRI", "SITEID_MRI",
               "GENOTYPE", "MANUFACTURER_WMH", "MANUFACTURERSMODELNAME_WMH",
               "CEREBRUM_TCV_WMH", "CEREBRUM_TCB_WMH", "CEREBRUM_TCC_WMH",
               "CEREBRUM_GRAY_WMH", "CEREBRUM_WHITE_WMH",
               "LEFT_HIPPO_WMH", "RIGHT_HIPPO_WMH", "TOTAL_HIPPO_WMH",
               "TOTAL_CSF_WMH", "TOTAL_GRAY_WMH", "TOTAL_WHITE_WMH",
               "TOTAL_WMH_WMH", "TOTAL_BRAIN_WMH")

merged_FW_data2 <- bind_rows(
  merged_FW_data %>% filter(!is.na(RID_WMH)) %>% select(all_of(cols_keep)),
  merged_FW_data_miss %>% mutate(VISCODE2_signs = VISCODE2_MRI) %>%
    select(all_of(cols_keep))
)

write.csv(merged_FW_data2, "FW_covariable_APOE_WMH.csv", row.names = FALSE)

# ---- 2.2 Diagnosis conversion dates ----
Diagnosis <- read.csv("Diagnostic Summary [ADNI1,GO,2,3,4].csv", header = TRUE) %>%
  select(1:7)
names(Diagnosis) <- paste0(names(Diagnosis), "_Dia")

Diagnosis <- Diagnosis %>%
  mutate(EXAMDATE_Dia = as.Date(EXAMDATE_Dia, format = "%Y-%m-%d"),
         DIAGNOSIS_Dia = as.numeric(DIAGNOSIS_Dia))

# Fill missing EXAMDATE using Vital Signs VISDATE
BP_BMI <- read.dta13("Vital Signs_ADNI1,GO,2,3,4.dta") %>% select(1:6)
merged_Diagnosis <- Diagnosis %>%
  left_join(BP_BMI,
            by = c("PTID_Dia" = "PTID", "PHASE_Dia" = "PHASE",
                   "RID_Dia" = "RID", "VISCODE_Dia" = "VISCODE")) %>%
  mutate(EXAMDATE_Dia = as.Date(EXAMDATE_Dia),
         VISDATE      = as.Date(VISDATE))

Diagnosis$EXAMDATE_Dia <- ifelse(is.na(Diagnosis$EXAMDATE_Dia),
                                 merged_Diagnosis$VISDATE,
                                 Diagnosis$EXAMDATE_Dia)
Diagnosis <- Diagnosis %>%
  filter(!is.na(EXAMDATE_Dia), !is.na(DIAGNOSIS_Dia))

write.csv(Diagnosis, "Diagnosis_processed.csv", row.names = FALSE)

# ---- 2.3 Derive AD_DATE, MCI_DATE and change labels ----
Diagnosis <- read.csv("Diagnosis_processed.csv", header = TRUE) %>%
  mutate(EXAMDATE_Dia = as.Date(EXAMDATE_Dia))

Diagnosis <- Diagnosis %>%
  group_by(PTID_Dia) %>%
  arrange(EXAMDATE_Dia, .by_group = TRUE) %>%
  mutate(
    first_diag = first(DIAGNOSIS_Dia[!is.na(DIAGNOSIS_Dia)]),
    ad_date_candidate = if (any(DIAGNOSIS_Dia == 3, na.rm = TRUE)) {
      first(EXAMDATE_Dia[DIAGNOSIS_Dia == 3 & !is.na(DIAGNOSIS_Dia)])
    } else {
      if (all(is.na(EXAMDATE_Dia))) NA else max(EXAMDATE_Dia, na.rm = TRUE)
    },
    mci_date_candidate = if (any(DIAGNOSIS_Dia == 2, na.rm = TRUE)) {
      first(EXAMDATE_Dia[DIAGNOSIS_Dia == 2 & !is.na(DIAGNOSIS_Dia)])
    } else {
      if (all(is.na(EXAMDATE_Dia))) NA else max(EXAMDATE_Dia, na.rm = TRUE)
    }
  ) %>%
  mutate(
    AD_DATE  = ifelse(first_diag != 3, ad_date_candidate, NA) %>% as.Date(),
    MCI_DATE = ifelse(first_diag == 1, mci_date_candidate, NA) %>% as.Date(),
    AD_change  = ifelse(first_diag == 3, NA,
                        ifelse(any(DIAGNOSIS_Dia == 3, na.rm = TRUE), 1, 0)),
    MCI_change = ifelse(first_diag != 1, NA,
                        ifelse(any(DIAGNOSIS_Dia == 2, na.rm = TRUE), 1, 0))
  ) %>%
  ungroup()

write.csv(Diagnosis, "Diagnosis_changed_label.csv", row.names = FALSE)

# ---- 2.4 Merge FW with diagnosis ----
Diagnosis <- read.csv("Diagnosis_changed_label.csv")
FW_data   <- read.csv("FW_covariable_APOE_WMH.csv") %>%
  mutate(EXAMDATE = as.Date(EXAMDATE))

merged_FW_data_dia <- FW_data %>%
  left_join(Diagnosis,
            by = c("PTID" = "PTID_Dia", "RID_MRI" = "RID_Dia",
                   "VISCODE" = "VISCODE_Dia", "PHASE_MRI" = "PHASE_Dia"))

# Fill unmatched via VISCODE2
miss_merge_FW_data_dia <- merged_FW_data_dia %>%
  filter(is.na(VISCODE2_Dia)) %>%
  select(-c(42:51)) %>%
  mutate(VISCODE2_MRI = ifelse(VISCODE2_MRI == "scmri", "bl", VISCODE2_MRI))

merged_FW_data_miss_dia <- left_join(
  miss_merge_FW_data_dia, Diagnosis,
  by = c("PTID" = "PTID_Dia", "RID_MRI" = "RID_Dia",
         "VISCODE2_MRI" = "VISCODE2_Dia", "PHASE_MRI" = "PHASE_Dia")
) %>% mutate(VISCODE2_Dia = VISCODE2_MRI)

order_cols <- c("PHASE_MRI", "PTID", "RID_MRI", "VISCODE", "VISCODE2_MRI",
                "EXAMDATE", "EXAMDATE_Dia", "Group", "DIAGNOSIS_Dia",
                "first_diag", "ad_date_candidate", "AD_DATE", "AD_change",
                "mci_date_candidate", "MCI_DATE", "MCI_change",
                "Age", "Sex", "PTEDUCAT_demo", "PTETHCAT_demo", "PTRACCAT_demo",
                "weight_signs", "height_signs", "BMI_signs",
                "VSBPSYS_signs", "VSBPDIA_signs",
                "FW_total", "FW_asso", "FW_comm", "FW_proj",
                "FIELD_STRENGTH_MRI", "MMCONDCT_MRI", "MMREASON_MRI", "SITEID_MRI",
                "GENOTYPE", "MANUFACTURER_WMH", "MANUFACTURERSMODELNAME_WMH",
                "CEREBRUM_TCV_WMH", "CEREBRUM_TCB_WMH", "CEREBRUM_TCC_WMH",
                "CEREBRUM_GRAY_WMH", "CEREBRUM_WHITE_WMH",
                "LEFT_HIPPO_WMH", "RIGHT_HIPPO_WMH", "TOTAL_HIPPO_WMH",
                "TOTAL_CSF_WMH", "TOTAL_GRAY_WMH", "TOTAL_WHITE_WMH",
                "TOTAL_WMH_WMH", "TOTAL_BRAIN_WMH")

merged_FW_data_dia2 <- bind_rows(
  merged_FW_data_dia %>% filter(!is.na(VISCODE2_Dia)) %>% select(all_of(order_cols)),
  merged_FW_data_miss_dia %>% select(all_of(order_cols))
)

# Compute follow-up time from baseline MRI
merged_FW_data_dia2 <- merged_FW_data_dia2 %>%
  mutate(EXAMDATE = as.Date(EXAMDATE)) %>%
  group_by(PTID) %>%
  mutate(min_MRI_DATE = min(EXAMDATE, na.rm = TRUE),
         max_MRI_DATE = max(EXAMDATE, na.rm = TRUE)) %>%
  ungroup()

write.csv(merged_FW_data_dia2, "FW_diagnosis_COX_data.csv", row.names = FALSE)

# ---- 2.5 Derive follow-up time and clean negative values ----
FW_data <- read.csv("FW_diagnosis_COX_data.csv", header = TRUE) %>%
  mutate(across(c(EXAMDATE, EXAMDATE_Dia, min_MRI_DATE, max_MRI_DATE,
                  AD_DATE, MCI_DATE, ad_date_candidate, mci_date_candidate),
                as.Date))

FW_data <- FW_data %>%
  mutate(AD_FELLOW_UP  = AD_DATE  - min_MRI_DATE,
         MCI_FELLOW_UP = MCI_DATE - min_MRI_DATE)

ptid_counts  <- FW_data %>% count(PTID)
single_ptids <- ptid_counts %>% filter(n == 1) %>% pull(PTID)

FW_data <- FW_data %>%
  mutate(
    AD_FELLOW_UP = case_when(
      !is.na(AD_FELLOW_UP) & !is.na(DIAGNOSIS_Dia) &
        AD_FELLOW_UP <= 0 & DIAGNOSIS_Dia == 3 ~ NA_real_,
      PTID %in% single_ptids & AD_FELLOW_UP <= 0 &
        DIAGNOSIS_Dia %in% c(1, 2) ~ NA_real_,
      AD_FELLOW_UP <= 0 ~ NA_real_,
      TRUE ~ as.numeric(AD_FELLOW_UP)
    ),
    MCI_FELLOW_UP = case_when(
      !is.na(MCI_FELLOW_UP) & !is.na(DIAGNOSIS_Dia) &
        MCI_FELLOW_UP <= 0 & DIAGNOSIS_Dia != 1 ~ NA_real_,
      PTID %in% single_ptids & MCI_FELLOW_UP <= 0 &
        DIAGNOSIS_Dia == 1 ~ NA_real_,
      TRUE ~ as.numeric(MCI_FELLOW_UP)
    ),
    AD_change  = ifelse(is.na(AD_FELLOW_UP),  NA, AD_change),
    MCI_change = ifelse(is.na(MCI_FELLOW_UP), NA, MCI_change)
  )

order_cols <- c("PHASE_MRI", "PTID", "RID_MRI", "VISCODE", "VISCODE2_MRI",
                "EXAMDATE", "min_MRI_DATE", "max_MRI_DATE", "EXAMDATE_Dia",
                "Group", "DIAGNOSIS_Dia", "first_diag", "ad_date_candidate",
                "AD_DATE", "AD_FELLOW_UP", "AD_change",
                "mci_date_candidate", "MCI_DATE", "MCI_change", "MCI_FELLOW_UP",
                "Age", "Sex", "PTEDUCAT_demo", "PTETHCAT_demo", "PTRACCAT_demo",
                "weight_signs", "height_signs", "BMI_signs",
                "VSBPSYS_signs", "VSBPDIA_signs",
                "FW_total", "FW_asso", "FW_comm", "FW_proj",
                "FIELD_STRENGTH_MRI", "MMCONDCT_MRI", "MMREASON_MRI", "SITEID_MRI",
                "GENOTYPE", "MANUFACTURER_WMH", "MANUFACTURERSMODELNAME_WMH",
                "CEREBRUM_TCV_WMH", "CEREBRUM_TCB_WMH", "CEREBRUM_TCC_WMH",
                "CEREBRUM_GRAY_WMH", "CEREBRUM_WHITE_WMH",
                "LEFT_HIPPO_WMH", "RIGHT_HIPPO_WMH", "TOTAL_HIPPO_WMH",
                "TOTAL_CSF_WMH", "TOTAL_GRAY_WMH", "TOTAL_WHITE_WMH",
                "TOTAL_WMH_WMH", "TOTAL_BRAIN_WMH")

FW_data <- FW_data[, order_cols]
write.csv(FW_data, "FW_COX_data_processed.csv", row.names = FALSE)

# ---- 2.6 Final covariate formatting (baseline row per subject) ----
FW_data <- read.csv("FW_COX_data_processed.csv", header = TRUE) %>%
  mutate(across(c(EXAMDATE, EXAMDATE_Dia, min_MRI_DATE, max_MRI_DATE,
                  AD_DATE, MCI_DATE, ad_date_candidate, mci_date_candidate),
                as.Date))

FW_data_min_date <- FW_data %>%
  group_by(PTID) %>%
  arrange(EXAMDATE, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  filter(!is.na(DIAGNOSIS_Dia)) %>%
  mutate(
    APOE4_carrier = factor(
      case_when(GENOTYPE == "4/4" ~ "2_E4_alleles",
                GENOTYPE %in% c("2/4", "3/4") ~ "1_E4_allele",
                TRUE ~ "0_E4_alleles"),
      levels = c("0_E4_alleles", "1_E4_allele", "2_E4_alleles")),
    DIAGNOSIS_Dia = factor(DIAGNOSIS_Dia, levels = c(1, 2, 3),
                           labels = c("NC", "MCI", "AD")),
    Sex = factor(Sex, levels = c("F", "M"), labels = c("Female", "male")),
    PTETHCAT_demo = factor(PTETHCAT_demo, levels = c("1", "2"),
                           labels = c("Hispanic-or-Latino",
                                      "Not-Hispanic-or-Latino")),
    MRI_FELLOW_UP = as.numeric(max_MRI_DATE - min_MRI_DATE) / 365.25
  )

# Rename to publication-friendly variable names
colnames(FW_data_min_date) <- c(
  "PHASE_MRI", "PTID", "RID_MRI", "VISCODE", "VISCODE2_MRI", "EXAMDATE",
  "min_MRI_DATE", "max_MRI_DATE", "EXAMDATE_Dia", "Group", "DIAGNOSIS_Dia",
  "first_diag", "ad_date_candidate", "AD_DATE", "AD_FELLOW_UP", "AD_change",
  "mci_date_candidate", "MCI_DATE", "MCI_change", "MCI_FELLOW_UP",
  "Age", "Sex", "Education", "Ethnic", "PTRACCAT_demo",
  "weight_signs", "height_signs", "BMI", "SBP", "DBP",
  "FW_total", "FW_asso", "FW_comm", "FW_proj",
  "FIELD_STRENGTH", "QC", "QC_reason", "SITEID", "GENOTYPE",
  "MANUFACTURER", "MANUFACTURERSMODELNAME",
  "CEREBRUM_TCV", "CEREBRUM_TCB", "CEREBRUM_TCC", "CEREBRUM_GRAY",
  "CEREBRUM_WHITE", "LEFT_HIPPO", "RIGHT_HIPPO", "TOTAL_HIPPO", "TOTAL_CSF",
  "TOTAL_GRAY", "TOTAL_WHITE", "TOTAL_WMH", "TOTAL_BRAIN",
  "APOE4_carrier", "MRI_FELLOW_UP"
)

write.csv(FW_data_min_date, "FW_COX_data_final.csv", row.names = FALSE)