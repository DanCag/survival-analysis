# ++++++++++++ #
# Kaplan-Meier #
# ++++++++++++ #

# Kaplan-Meier analysis includes
# - maximally selected ranks statistics cutpoint
#   (see commented lines 127-130 of the script to plot cutpoint)
# - KM plots


# clear the environment
rm(list=ls())

# input parameters ----

# variable we use to stratify patients
strat <- "CD8_Tem_GZMK_high"

# load build parameters obtained from script "01_build_survival-dataframe.R"
load(file = "./output/parameters/build_parameters.RData")

# load integrate parameters obtained from script "02_"
load(file = "./output/parameters/integrate_parameters.RData")

# flag to decide how to stratify patients:
# - if FALSE: patients are stratified based on cutpoint computed with 
#             maximally selected ranks statistics
# - if TRUE: patients are stratified based on cutpoint that keeps 
#            the same high/low patients proportion
#            as found for CD8_Tem population
#            when using the "CD8-no-split" signature.
#            54% of patients in "high" group (>= cutpoint),
#            46% of patients in the "low" group (< cutpoint)
keep_CD8_Tem_proportion <- FALSE


# import packages ----

library(survival)
library(survminer)
source("./R/functions/km_fit.R")


# import datasets ----

# survival dataframe (survival table + molecular features)
surv_mf <- readRDS(
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


# Kaplan Meier ----

if (is.numeric(surv_mf[[strat]]) ) {
  
  # if stratification variable is numeric
  
  if (keep_CD8_Tem_proportion) {
    
    # if keep_CD8_Tem_proportion is TRUE
    
    cutpoint <- quantile(surv_mf[[strat]],
                         probs = c(0.46))
    
  } else {
      
    # define cutpoint with maximally selected ranks statistics
    cutpoint_l <- surv_cutpoint(
      surv_mf, 
      time = "years",
      event = paste(survival_analysis, 
                    "status", 
                    "int", 
                    sep = "_"),
      variables = strat)
    
    # extract cutpoint from cutpoint list
    cutpoint <- cutpoint_l$cutpoint$cutpoint
    
    }
  
  # km fit object 
  km_fit_object <- km_fit(
    df = surv_mf,
    time = "years", 
    censor = paste(survival_analysis,
                   "status",
                   "int",
                   sep = "_"),
    strat = strat, 
    cutpoint = cutpoint
  )
  
  # km plot
  ggsurvplot(km_fit_object, 
             data = surv_mf,
             conf.int = F,
             pval = T, 
             risk.table = T,
             legend.title = "Infiltration",
             legend.labs = c("high", "low"),
             xlim = c(0, 5), 
             break.x.by = 1,
             risk.table.y.text.col = T,
             risk.table.y.text = FALSE,
             title = strat, 
             subtitle = survival_analysis,
             xlab = "Years",
             font.x = 14,
             font.y = 14)
  
  # if you want to check how cutpoint is defined and plot it, 
  # uncomment the following 2 lines
  # summary(cutpoint_l)   # print cutpoint value
  # plot(cutpoint_l, strat, palette = "npg") # Plot cutpoint for strat feature
  
  
} else {
  
  # The only case where the stratification variable is not numeric
  # is when we consider the neutrophils infiltration status
  # coupled with cytolytic activity (strat = "Neutrophils_Cyt")
  
  # km plot
  ggsurvplot(
    fit = survfit(
      Surv(time = years, 
           event = DFS_status_int) ~ Neutrophils_Cyt, 
      data = surv_mf), 
    xlab = "years", 
    ylab = paste(survival_analysis, 
                 "probability"), 
    pval = T,
    risk.table = T,
    legend.title = "Neutrophils Infiltration - Cytolytic activity",
    legend.labs = c("No Neutr-High Cyt", 
                    "No Neutr-Low Cyt",
                    "Yes Neutr-High Cyt", 
                    "Yes Neutr-Low Cyt"), 
    xlim = c(0, 5), 
    break.x.by = 1,
    risk.table.y.text.col = T, 
    risk.table.y.text = FALSE ) + 
    ggtitle("Neutrophils infiltration + Cytolytic activity", 
            subtitle = survival_analysis) +
    xlab("Years")
  
}
