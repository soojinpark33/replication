#### Clear workspace ####
rm(list = ls())

# Load necessary libraries and functions
library(dplyr)
setwd("/Users/soojinp/Library/CloudStorage/OneDrive-UniversityofCalifornia,Riverside/replication/Code")
source("replication_functions.R")

load("data_final.rda")

#::::::::::::::::::::::::::::::::::::::::::::::::
#::::: 1. Prepare and Clean the Initial Data :::::
#::::::::::::::::::::::::::::::::::::::::::::::::
# Define variable lists for easier management
vars_to_scale <- c(
  "age", "soc_stat", "know_psy", "know_sopsy", "know_conhyp",
  "cont_hmls", "ext_cont_hmls", "intgr_anxty_hmls", "empathy",
  "RWA", "SDO", "pol_orient", "viv_imag", "motiv_part",
  "intr_stud_incent", "intr_stud_topic", "noise_lvl", "distr",
  "temptr", "awar_hmls", "PO14_01"
)

data_prepared <- d %>%
  # Filter for studies 1 and 3
  filter(BV02 %in% c(1, 3)) %>%
  # Use mutate(across(...)) for all transformations
  mutate(
    # Scale all specified metric variables
    across(all_of(vars_to_scale), ~ as.numeric(scale(.))),
    
    # Center baseline covariates to their means
    gender_I = as.numeric(gender) - mean(as.numeric(gender), na.rm = TRUE),
    ger_skil_I = as.numeric(ger_skil) - mean(as.numeric(ger_skil), na.rm = TRUE),
    prfsn_I = as.numeric(prfsn) - mean(as.numeric(prfsn), na.rm = TRUE),
    
    # Recode setting variables (0 = reference)
    locat_I = as.numeric(locat) - 1,
    device_I = as.numeric(device) - 1,
    pres_othr_I = as.numeric(pres_othr_I) - 1,
    
    # Create the combined study/treatment group indicator 'XR'
    XR = case_when(
      BV03 == 2 & BV02 == 1 ~ "00",
      BV03 == 2 & BV02 == 3 ~ "01",
      BV03 == 1 & BV02 == 1 ~ "10",
      TRUE ~ "11" # All other cases
    )
  ) %>%
  # Order the data by the new group indicator
  arrange(XR)


#::::::::::::::::::::::::::::::::::::::::::::::::
#::::: 2. Create the Final Analysis Dataset  :::::
#::::::::::::::::::::::::::::::::::::::::::::::::

# Define final variable sets
Zvars <- c(
  "gender_I", "ger_skil_I", "age", "prfsn_I", "soc_stat", "know_psy",
  "know_sopsy", "know_conhyp", "cont_hmls", "ext_cont_hmls",
  "intgr_anxty_hmls", "empathy", "RWA", "SDO", "pol_orient",
  "viv_imag", "motiv_part", "intr_stud_incent", "intr_stud_topic"
)
Wvars <- c("locat_I", "pres_othr_I", "noise_lvl", "temptr")

# Create the final analysis-ready data frame 'df'
df <- data_prepared %>%
  # Impute missing numeric values with the column median
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
  # Create binary indicators, outcome, and Wprime
  mutate(
    R = ifelse(XR %in% c("00", "10"), 0, 1),
    X = ifelse(XR %in% c("00", "01"), 0, 1),
    Y = PO14_01,
    Wprime = device_I
  ) %>%
  # Select only the necessary variables for the analysis
  dplyr::select(Y, R, X, Wprime, all_of(Zvars), all_of(Wvars))

# Clean up the environment by removing intermediate objects
rm(d, data_prepared, vars_to_scale)


#:::::::::::::::::::::::::::::::::
#:::::::POC:::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::
poc_results <- bootstrap_poc(
  B = 1000,
  df = df,
  Y = "Y",
  R = "R",
  X = "X",
  Zvars = Zvars,
  Wvars = Wvars,
  Wprime = "Wprime",
  w_fix = 0
)

#--- View the Final Results ---
print(poc_results)

#:::::::::::::::::::::::::::::::::
#:::::::Imputation::::::::::::::::
#:::::::::::::::::::::::::::::::::
imputation_results <- bootstrap_imputation(
  B = 1000, 
  df = df,
  Y = "Y",
  R = "R",
  X = "X",
  Zvars = Zvars,
  Wvars = Wvars,
  Wprime = "Wprime",
  w_fix = 0
)

# View the final results
print(imputation_results)

#:::::::::::::::::::::::::::::::::
#:::::::LASSO::::::::::::::::
#:::::::::::::::::::::::::::::::::
# Hybrid LASSO (hybrid is now the default in replication_functions.R).
lasso_hybrid_results_min <- bootstrap_lasso(
  B = 1000,
  df = df,
  Y = "Y",
  R = "R",
  X = "X",
  Zvars = Zvars,
  Wvars = Wvars,
  Wprime = "Wprime",
  w_fix = 0,
  lambda_rule = "min",
  lasso_type = "double"
)
print(lasso_hybrid_results_min)

lasso_hybrid_results_1se <- bootstrap_lasso(
  B = 1000,
  df = df,
  Y = "Y",
  R = "R",
  X = "X",
  Zvars = Zvars,
  Wvars = Wvars,
  Wprime = "Wprime",
  w_fix = 0,
  lambda_rule = "1se",
  lasso_type = "double"
)
print(lasso_hybrid_results_1se)

# Selection diagnostics for each LASSO run (double x {min, 1se}).
lasso_diagnostics_double <- collect_lasso_diagnostics(
  df = df,
  Y = "Y",
  R = "R",
  X = "X",
  Zvars = Zvars,
  Wvars = Wvars,
  Wprime = "Wprime",
  lasso_types = "double"
)
print(as.data.frame(lasso_diagnostics_double), row.names = FALSE)
