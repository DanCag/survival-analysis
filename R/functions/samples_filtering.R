# ++++++++++++++++++++ #
# Filter tumor samples #
# ++++++++++++++++++++ #

filter_samples <- function(clinical_input_path, 
                           samples_input_path,
                           sample_type,
                           stages) {
  
  ## Filter TCGA samples ##
  
  ## Input:
  # - clinical_input_path: chr, path of "clinical" file
  # - samples_input_path: chr, path of "samples" file
  # - sample_type: chr, type of sample you want to select
  #                ("Primary Tumor", "Metastatic", ...)
  # - stages: chr, stages you want to select (I, II, ...)
  
  # import datasets ----
  
  # clinical info 
  clinical <- read.delim(clinical_input_path) [-c(1,2),]
  
  # samples info 
  samples <- read.delim(samples_input_path,
                        stringsAsFactors = F)[-1, ]
  
  # Data wrangling ----
  
  # select primary tumor sample ids
  # exclude metastases and normal tissues (blood derived or solid)
  primary_sampleIds <- samples$bcr_sample_barcode[
    samples$sample_type == sample_type]
  
  # extract bcr_patient_barcode from bcr_sample_barcode 
  primary_patientIds <- substr(primary_sampleIds, 1, 12)
  
  # get duplicated primary_patientIDs
  duplicatedIds <- primary_patientIds[duplicated(primary_patientIds)]
  
  # remove duplicated IDs
  primary_patientIds_unique <- unique(primary_patientIds)

  # samples with wanted stages  
  patientIds_desidered_stages <- clinical$bcr_patient_barcode[
    clinical$ajcc_pathologic_tumor_stage %in% stages]
  
  # intersect primary and desidred stages patients
  patientIds_intersection <- intersect(
    primary_patientIds_unique, 
    patientIds_desidered_stages)

}
