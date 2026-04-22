# replication
The codes are to replicate analyses for "Identifying Causes of Treatment Effect Heterogeneity Across Replication Studies: Definitions, Identification, and Estimation of Remaining Heterogeneity."

Abstract: Replication studies are crucial for evaluating treatment efficacy and exploring the conditions that moderate treatment effects. Understanding this treatment effect heterogeneity across studies is essential for building robust theories. However, identifying the causes of this heterogeneity is challenging, as multiple study characteristics often vary simultaneously, even in replication studies in which many aspects are deliberately held constant across studies. Traditional meta-regression typically lacks a clear causal framework, obscuring the causal interpretation of moderating effects. While a more recent causal method provides a way to identify the replication effect after adjusting for unintended differences, it relies on the strong assumption that all effect moderators that cause the treatment effect heterogeneity across studies must be measured and uses a parametric estimation model. We take a different view on that by introducing a causal framework that decomposes sources of treatment effect heterogeneity. Our approach uses a Directed Acyclic Graph (DAG) to identify these sources and to develop strategies to single out their contributions. By focusing on  different (causal and non-causal) quantities, the proposed strategies require weaker assumptions than the earlier causal approach. We also propose robust estimation strategies using a variable selection method designed to model these higher-order interactions with minimal bias. We illustrate our method using a real-world example of replication studies from social psychology. Through a simulation and case study, we provide useful recommendations for applied researchers on designing future replication studies and identifying moderators that cause effect heterogeneity.

Here, we provide `R` codes to reproduce our simulation study and replicate our data analysis. 

## Case Study

* `case_main_code.R` 
 
   This `R` file replicates Table 2 of our study.

* `replication_functions.R`  
 
   This `R` file includes source functions required to run full analyses for our case study.
  
## Simulation Study

* `Simulation_main_code.R`  

   This `R` file contains the simulation codes for our proposed analysis. This code replicates Figure 2 of our paper.

* `Simulation_functions.R`  
 
   This `R` file includes source functions required to run our simulation codes.
  
# Replication Analysis Guide

To use the source code, load the functions into your R session:

```r
# Load the core functions
source("replication_functions.R")

# 1. Coefficient-based approach (relies on Eq. (12) and (13))
poc_results <- bootstrap_poc(
  B = 1000, 
  df = df, 
  Y = "Y", 
  R = "R", 
  X = "X", 
  Zvars = c("age", "gender"),  
  Wvars = c("noise", "loc"), 
  Wprime = "device", 
  w_fix = 0
)
# 2. Imputation estimator (relies on Eq. (12) and (13))
imp_results <- bootstrap_imputation(
  B = 1000, 
  df = df, 
  Y = "Y", 
  R = "R", 
  X = "X", 
  Zvars = c("age", "gender"), 
  Wvars = c("noise", "loc"), 
  Wprime = "device", 
  w_fix = 0
)
# 3. Double Lasso Imputation
lasso_results <- bootstrap_lasso(
  B = 1000, 
  df = df, 
  Y = "Y", 
  R = "R", 
  X = "X", 
  Zvars = c("age", "gender"), 
  Wvars = c("noise", "loc"), 
  Wprime = "device", 
  w_fix = 0, 
  lambda_rule = "min"
)
