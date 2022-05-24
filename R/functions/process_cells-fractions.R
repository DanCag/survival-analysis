# ++++++++++++++++++++++++++++ #
# Process cells-fractions file #
# ++++++++++++++++++++++++++++ #

process_cells_fractions <- function(
    cells_fractions_path, 
    sample_type_code = NULL) {
  
  ## Wrangling and cleaning of CIBERSORT cells-fractions input ##
  
  ## Input:
  
  # - cells_fractions_path: chr, path of CIBERSORT cells-fractions file
  # - sample_type_code: chr, (default is NULL). 
  #                     Specify if you want to pick primary, metastases, ..
  
  
  # import dataset ---- 
  
  cells_fractions <- read.delim(cells_fractions_path)
  
  
  # Data wrangling ----
  
  # replace "." with "-" in barcode
  cells_fractions$Mixture <- gsub("\\.", "-", cells_fractions$Mixture)
  
  # trim tcga barcode to get
  
  # - bcr_sample_barcode
  cells_fractions$bcr_sample_barcode <- substr(
    cells_fractions$Mixture, 1, 16)
  
  # - bcr_patient_barcode
  cells_fractions$bcr_patient_barcode <- substr(
    cells_fractions$Mixture, 1, 12)
  
  
  # Data cleaning ----
  
  # remove duplicated patient ids
  cells_fractions_no_dup <- cells_fractions[
    !(duplicated(cells_fractions$bcr_patient_barcode) |
        duplicated(cells_fractions$bcr_patient_barcode, fromLast = TRUE)), ]
  
  # order it to have different barcodes
  # (Mixture, bcr_sample_barcode, bcr_patient_barcode) close each other
  cells_fractions_no_dup <- cells_fractions_no_dup[
    , c(1, 
        (ncol(cells_fractions_no_dup)) - 1, 
        ncol(cells_fractions_no_dup), 
        2:(ncol(cells_fractions_no_dup) -2))]
  
  # pick samples matching sample_type_code
  if (!is.null(sample_type_code)) {
    
    # pick only samples matching sample_type_code
    ids <- grep(paste0(sample_type_code, "."),
                cells_fractions_no_dup$bcr_sample_barcode, 
                value = T)
    
    # subset cells_fractions
    cells_fractions_no_dup_ids <- cells_fractions_no_dup[
      cells_fractions_no_dup$bcr_sample_barcode %in% ids, ]
    
    # processed cells_fractions
    cells_fractions_final <- cells_fractions_no_dup_ids
    
  } else {
    
    # processed cells_fractions
    cells_fractions_final <- cells_fractions_no_dup
  }
  
  return(cells_fractions_final)
}