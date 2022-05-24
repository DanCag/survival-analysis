# +++++++++++++++++++++++++ #  
# Immune cells infiltration #
# +++++++++++++++++++++++++ #   

# Integrate immune cells infiltration into survival table 

# input parameters ----

# tumor type
# it can be either 
# "COAD": colon adenocarcinoma 
# "LUAD": lung adenocarcinoma
tumor_type <- "LUAD"

# type of survival analysis
# "OS" stands for overall survival
# "DFS" stands for disease-free survival
survival_analysis <- "OS"

# tumor stages taken into account
stages_info <- "I-IV"

# signature used for inferring cells fractions
# it can be either "CD8-no-split" or "CD8-split"
signature_type <- "CD8-split"

# CIBERSORTx parameters
cibersort_info <- "rm-batch-S-mode_1000-permutations"

# cells fractions path
cells_fractions_path <- file.path(
    "./data/CIBERSORTx",
    tumor_type,
    paste0(
      "CIBERSORTx_",
      tumor_type,
      "_",
      signature_type,
      "_", 
      cibersort_info,
      ".txt")
  )


# import functions ----

# functions for integrating features
source("./R/functions/process_cells-fractions.R")


# import datasets ----

# survival table 
surv_df <- readRDS(
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

# cells fractions
cells_fractions <- process_cells_fractions(
  cells_fractions_path = cells_fractions_path)


# Build survival dataframe + immune cells infiltrations ----

# merge surv_df with cells_fractions
# by common column "bcr_patient_barcode"
surv_features_df <- merge(surv_df, 
                          cells_fractions, 
                          by = "bcr_patient_barcode")


# save surv_features_df
saveRDS(
  surv_features_df,
  file = file.path(
    "./output",
    tumor_type,
    paste0(
      "surv-",
      tolower(survival_analysis),
      "_",
      tolower(tumor_type),
      "_",
      paste("cells-fractions",
            collapse = '-'),
      "_",
      signature_type,
      "_",
      stages_info,
      ".rds")
    )
  )
