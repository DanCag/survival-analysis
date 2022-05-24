# ++++++++++++ #
# Kaplan-Meier #
# ++++++++++++ #

# Kaplan-Meier fit object
km_fit <- function(df, time, censor, stratification, cutpoint) {
  
  ## It takes in input a df with time, censor and stratification info ##
  ## and returns a Kaplan-Meier fit object                            ##
  
  # Input:
  # - df: dataframe, dataframe with survival info + features
  # - time: chr, name of column with time info
  # - censor: chr, name of column with censor info
  # - stratification: chr, name of column you want to use for stratification
  # - cutpoint: num, value to split patients into "low" and "high" group
  
  # stratify based on cutpoint #  
  strat <- factor(ifelse(df[[stratification]] <= cutpoint, "low", "high"))
  
  # fit #
  
  # formula
  form <- as.formula("Surv(time = df[[time]], event = df[[censor]]) ~ strat")

  # fit the survival object
  fit <- survfit(form, data = df)
  
  # manually put the formula into fit$call$formula
  fit$call$formula <- form
  
  return(fit)
}