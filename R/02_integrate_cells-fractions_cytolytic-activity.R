# ++++++++++++++++++++++++++++++ #  
# Immune cells infiltration and  #
# cytolytic activity information #
# ++++++++++++++++++++++++++++++ #  

# Integrate immune cells infiltration and cytolytic activity information
# into survival table 


# input parameters ----

# signature used for inferring cells fractions
# it is "LM22"
signature <- "LM22"

# feature we will integrate
feature <- "cells-fractions_cytolytic-activity"

# save parameters that I will need for script "03_km.R"
save(signature, 
     feature, 
     file = "./output/parameters/integrate_parameters.RData")

# load build parameters obtained from script "01_build_survival-dataframe.R"
load("./output/parameters/build_parameters.RData")

# CIBERSORTx parameters
cibersort_info <- "1000-permutations"


# import packages ----
library(survminer)
source("./R/functions/process_cells-fractions.R")


# import datasets ----

# survival table
surv_df <- readRDS(
  file = file.path(
    "./output", 
    tumor_type,
    paste0(
      "surv-dfs_", 
      tolower(tumor_type), 
      ".rds")
  )
)

## cells fractions ##
cells_fractions <- process_cells_fractions(
  cells_fractions_path = file.path(
    "./data/CIBERSORTx",
    tumor_type,
    paste0(
      "CIBERSORTx_",
      tumor_type,
      "_",
      signature,
      "_", 
      cibersort_info,
      ".txt")
    )
  )

## cytolytic activity ## 

# geometric mean of 
# - GZMA (ENSG00000145649.7) 
# - Perforin (ENSG00000180644.6) 
# expression level (TPM)

cyt <- read.delim(
  file.path(
    "./data/cytolyitic-activity",
    tumor_type,
    paste0(
      "cyt_tcga-", 
      tolower(tumor_type), 
      ".tsv")),        
  stringsAsFactors = F
)


# Data wrangling ----

## cytolytic activity ##

# add new column with trimmed barcode
cyt$bcr_patient_barcode <- gsub(
  "\\.", "-",
  substr(cyt$Mixture, 1, 12)
)

# find positions of surv_df patients in cyt dataframe
match_patients <- match(surv_df$bcr_patient_barcode, 
                        cyt$bcr_patient_barcode)

# subset cyt dataframe by picking matching samples
cyt_subset <- cyt[match_patients, ] 


# Build dataframe with survival + cells fraction + cytolytic activity ----

# survival + cells-fractions + cytolytic activity
surv_cells_cyt <- Reduce(
  merge, 
  list(surv_df, cells_fractions, cyt_subset)
  )

# add column with Neutrophils presence/absence
surv_cells_cyt$Neutrophils_presence <- factor(
  ifelse(
    surv_cells_cyt$Neutrophils == 0, "No", "Yes")
  )

# define low and high cytolytic activity
# based on maximally selected rank statistics
cutpoint_cyt <- surv_cutpoint(
  surv_cells_cyt, 
  time = "years",
  event = paste(
    survival_analysis, 
    "status", 
    "int", 
    sep = "_"),
  variables = "Cyt")

# add column where we split cytolytic activity into "High" or "Low"
surv_cells_cyt$Cyt_categorical <- ifelse(
  surv_cells_cyt$Cyt > cutpoint_cyt$cutpoint$cutpoint, "High", "Low")

# add col with info about neutrophils presence and cytolytic activity 
surv_cells_cyt$Neutrophils_Cyt <- factor(
  paste0(
    surv_cells_cyt$Neutrophils_presence, 
    "-", 
    surv_cells_cyt$Cyt_categorical))

# save surv_exp 
saveRDS(
  object = surv_cells_cyt, 
  file = file.path(
    "./output", 
    tumor_type,
    paste0(
      "surv-", 
      tolower(survival_analysis), 
      "_",
      tolower(tumor_type), 
      "_",
      feature,
      "_",
      signature, 
      "_",
      stages_info, 
      ".rds")
  )
)