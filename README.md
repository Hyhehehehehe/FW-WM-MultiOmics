# FW-WM Multi-Omics
This repository contains the source code for the imaging, genetic, phenotypic, longitudinal, and proteomic analyses presented in the study:

**“Multi-omics characterization of white matter free water reveals genetic and vascular determinants of neurological disorders”**

The analytical workflow integrates data from the UK Biobank and Alzheimer’s Disease Neuroimaging Initiative (ADNI) to characterize white matter free water (FW-WM) across multiple white matter fiber systems and to investigate its genetic architecture, neurological disease associations, vascular determinants, and candidate molecular pathways.

The repository includes code for FW-WM estimation, genome-wide association studies (GWAS), polygenic risk score (PRS) analyses, phenome-wide association studies (PheWAS), longitudinal survival analyses, Mendelian randomization (MR), proteomic analyses, colocalization, mediation analyses, and sensitivity analyses.

## Repository Structure
```text
.
├── FW estimation/
│   ├── FW_multishell.py
│   └── FW_singleshell.py
│
├── GWAS and PRS/
│   ├── 1_run_gwas_fw.sh
│   ├── 2_merge_fw.sh
│   └── 3_PRScs.sh
│
├── UKB analysis/
│   ├── 01_FW_PRS_besed_PheWAS.R
│   ├── 02_FW_Image_besed_PheWAS.R
│   ├── 03_FW_extreme_stratification_analysis in UKB Participants.R
│   ├── 04_Cox Regression and KM Curves for Stroke in UKB Participants.R
│   ├── 05_FW and Stroke Risk in Hypertensive UKB Participants.R
│   └── 06_FW and Cognitive domain Impairment in Hypertensive UKB Participants.R
│
├── ADNI analysis/
│   ├── 01_Covariates_Matching.R
│   ├── 02_Imaging_and_Diagnosis_Merging.R
│   ├── 03_Cox_AD_Conversion.R
│   ├── 04_Cox_MCI_Conversion.R
│   ├── 05_Cox_NC_to_AD.R
│   ├── 06_Cox_MCI_to_AD.R
│   ├── 07_KM_Curves.R
│   ├── 08_HTN_Analysis.R
│   ├── 09_FW and Cognitive Diagnosis in Hypertensive ADNI Participants.R
│   ├── 10_FW and Cognitive domain Impairment in Hypertensive ADNI Participants.R
│   └── 11_ADNI Single-b Sensitivity Analysis.R
│
├── MR analysis/
│   ├── 01_Disease_FW_MR.R
│   ├── 02_FW_Disease_MR.R
│   ├── 03_decode_protein_FW_MR.R
│   └── 04_hypertension_protein_FW_MR.R
│
├── Colocalization Analysis/
│   └── FW Key Protein Colocalization Analysis.R
│
├── Mediation analysis/
│   └── FW_mediation_analysis.R
│
├── Correlation analyses/
│   ├── FW_eroded_correlation.R
│   ├── FW_followup_correlation.R
│   ├── FW_PheWAS_OR_correlation_adjust_motion.R
│   ├── FW_ventricle_correlation.R
│   └── MB_SB_FW_correlation.R
│
└── README.md
```
