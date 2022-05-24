# ++++++++++++++++++ #
# Survival Dataframe #
# ++++++++++++++++++ #

# Build survival dataframe with following info: 
# - survival time
# - survival status


# input variables ----

# tumor type
# it can be either 
# "COAD": colon adenocarcinoma 
# "LUAD": lung adenocarcinoma
tumor_type <- "COAD"

# type of survival analysis
# "OS" stands for overall survival
# "DFS" stands for disease-free survival
survival_analysis <- "DFS"

# type of sample (we only select primary tumor samples)
sample_type <- "Primary Tumor"

# selected tumor stages
stages <- c("Stage IA", 
            "Stage II", 
            "Stage IIA", 
            "Stage IIB", 
            "Stage IIC", 
            "Stage III", 
            "Stage IIIA", 
            "Stage IIIB", 
            "Stage IIIC",
            "Stage IVA", 
            "Stage IVB")

# tumor stages taken into account
stages_info <- "I-IV"

# clinical dataset
clinical_path <- file.path(
  "./data/clinical", 
  tumor_type,
  paste0(
    "nationwidechildrens.org_clinical_patient_", 
    tolower(tumor_type),
    ".txt")
  )

# biospecimen dataset
samples_path <- file.path(
  "./data/sample", 
  tumor_type, 
  paste0("nationwidechildrens.org_biospecimen_sample_", 
  tolower(tumor_type), 
  ".txt")
  )


# import functions ----

# samples filtering 
source("./R/functions/samples_filtering.R")


# import dataset ----

# survival table
surv_df <- read.table(
  file.path(
    "./data/cBioPortal_survival-tables",
    paste0("TCGA-", 
           tumor_type), 
    paste0(
      survival_analysis, 
      "_", 
      tumor_type, 
      ".txt")
    ),
  sep = "\t",
  stringsAsFactors = F,
  col.names = c(
    "tcga_study", 
    "bcr_patient_barcode", 
    paste(survival_analysis, 
          "status", 
          sep = "_"), 
    "months")
  )

# use "filter_samples" function from "./R/functions/samples_filtering.R"
# to select desired patients
patientsIds <- filter_samples(
  clinical_input_path = clinical_path, 
  samples_input_path = samples_path, 
  sample_type = sample_type, 
  stages = stages)


# Data wrangling ----

## survival table ##

# split status column at ":"
surv_df_split <- cbind(
  surv_df[, c(1,2)], 
  do.call(
    'rbind', 
    strsplit(surv_df[[paste(survival_analysis, "status", sep = "_")]], 
             ':', 
             fixed=TRUE)),
  surv_df[, 4])

# rename columns 3 and 4
colnames(surv_df_split) <- c(
  colnames(surv_df)[1:2], 
  paste(survival_analysis, "status", "int", sep = "_"), 
  paste(survival_analysis, "status", "chr", sep = "_"), 
  "months")

# make column 3 an integers column
class(surv_df_split[[
  paste(survival_analysis, "status", "int", sep = "_")]]) <- "integer"

# add column with years
surv_df_split$years <- surv_df_split$months/12

# subset patients 
surv_df_split_sub <- surv_df_split[
  surv_df_split$bcr_patient_barcode %in% patientsIds, ]

# save surv_df_split
saveRDS(
  surv_df_split_sub,
  file = file.path(
    "./output",
    tumor_type,
    paste0(
      "surv-",
      tolower(survival_analysis),
      "_",
      tolower(tumor_type),
      ".rds")
  )
)
