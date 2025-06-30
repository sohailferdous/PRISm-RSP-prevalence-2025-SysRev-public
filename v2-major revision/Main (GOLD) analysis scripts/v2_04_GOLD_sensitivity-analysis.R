# Global prevalence of preserved ratio impaired spirometry and restrictive spirometry pattern in the 
# general population: a systematic review and meta-analysis 
#
# SCRIPT INTRODUCTION
# 
# This script performs the same function as v2_gold_sens_prevalence script, only with studies with 
# QA score >= 7.
# 
# This script is divided into the following broad sections:
# 1. Data import and preparation
# 2. Model building
# 3. Model summary checking
# 4. Output display
# 5. Publication bias checking
# 
# IMPORTANT - since some variable names may be the same as v2_GOLD_prevalence script, it is highly
# recommended to restart your R session (or clear the environment) if compiling this script 
# immediately after v2_GOLD_prevalence.
# 
# END OF INTRODUCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Please ensure working directory contains v2_data_prevalence.csv file.
#
# Install and import packages
#
# install.packages("tidyverse")
# install.packages("meta")
# install.packages("metafor")
# install.packages("remotes")
# remotes::install_github("MathiasHarrer/dmetar")
#
{                       # Compile here to finish package loading, data import and preparation, and model building in one step
#
library(tidyverse)
library(meta)
library(metafor)
library(dmetar)
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. DATA IMPORT AND PREPARATION ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
{                       # Compile here to finish the entire import and preparation in one step
#
#* 1.1 Data import and filtering by gold_sens definition -------------------------------------------
#* 
overall <- read.csv("v2_data_prevalence.csv")

overall <- overall %>%
  mutate(geo_location = factor(geo_location),
         id = factor(id),
         country = factor(country),
         income_class = factor(income_class),
         diag = factor(diag),
         diag_defn = factor(diag_defn),
         bronchod = factor(bronchod),
         overall.prev = factor(overall.prev),
         gender.prev = factor(gender.prev),
         comb.prism.rsp.prev = factor(comb.prism.rsp.prev),
         smoking.prev = factor(smoking.prev))

# Create gold_sens dataset - for main analysis
gold_sens <- filter(overall, diag_defn == "gold" & qa_score >=7) %>% # only take gold_sens studies with QA score >=7
  filter(!is.na(case), !is.na(sample_size), !is.na(geo_location))

#* 1.2 Creating subsets for meta-analysis ----------------------------------------------------------
#
#** 1.2.1 Overall prevalence
#
# prism
prism_gold_sens <-filter(gold_sens, diag == "prism" & overall.prev == "y")

# rsp
rsp_gold_sens <-filter(gold_sens, diag == "rsp" & overall.prev == "y")

# Combined prism and rsp 
comb_gold_sens <-filter(gold_sens, comb.prism.rsp.prev == "y")

#** 1.2.2 Subgroups - sex
#
# Males - combined, prism and rsp groups
male_gold_sens_comb <- filter(gold_sens, gender.prev == "y") %>% 
  drop_na(male_total, male_cases)

male_gold_sens_prism <- filter(gold_sens, diag == "prism" & gender.prev == "y") %>% 
  drop_na(male_total, male_cases)

male_gold_sens_rsp <- filter(gold_sens, diag == "rsp" & gender.prev == "y") %>% 
  drop_na(male_total, male_cases)

# females - combined, prism and rsp groups
female_gold_sens_comb <- filter(gold_sens, gender.prev == "y") %>% 
  drop_na(female_total, female_cases)

female_gold_sens_prism <- filter(gold_sens, diag == "prism" & gender.prev == "y") %>% 
  drop_na(female_total, female_cases)

female_gold_sens_rsp <- filter(gold_sens, diag == "rsp" & gender.prev == "y") %>% 
  drop_na(female_total, female_cases)

#** 1.2.3 Subgroups - smoking status
#
# Non-smokers - combined, prism and rsp groups
non_smoker_gold_sens_comb <- filter(gold_sens, smoking.prev == "y") %>% 
  drop_na(tot_non,case_non)

non_smoker_gold_sens_prism <- filter(gold_sens, diag == "prism" & smoking.prev == "y") %>% 
  drop_na(tot_non,case_non)

non_smoker_gold_sens_rsp <- filter(gold_sens, diag == "rsp" & smoking.prev == "y") %>% 
  drop_na(tot_non,case_non)

# Ex-smokers - combined, prism and rsp groups
ex_smoker_gold_sens_comb <- filter(gold_sens, smoking.prev == "y") %>% 
  drop_na(tot_ex,case_ex)

ex_smoker_gold_sens_prism <- filter(gold_sens, diag == "prism" & smoking.prev == "y") %>% 
  drop_na(tot_ex,case_ex)

ex_smoker_gold_sens_rsp <- filter(gold_sens, diag == "rsp" & smoking.prev == "y") %>% 
  drop_na(tot_ex,case_ex)

# Current smokers - combined, prism and rsp groups
cur_smoker_gold_sens_comb <- filter(gold_sens, smoking.prev == "y") %>% 
  drop_na(tot_cur,case_cur)

cur_smoker_gold_sens_prism <- filter(gold_sens, diag == "prism" & smoking.prev == "y") %>% 
  drop_na(tot_cur,case_cur)

cur_smoker_gold_sens_rsp <- filter(gold_sens, diag == "rsp" & smoking.prev == "y") %>% 
  drop_na(tot_cur,case_cur)

#* 1.3 Calculating effect sizes and creating variance-covariance matrix for MLMA analysis ----------
#
# prism
prism_es <- escalc(
  measure = "PLO",                 # Logit transformation for proportions
  xi = case,                       # Number of cases (events)
  ni = sample_size,                # Total sample size
  data = prism_gold_sens)

V_matrix_prism <- diag(prism_es$vi)

# rsp
rsp_es <- escalc(
  measure = "PLO",
  xi = case,
  ni = sample_size,
  data = rsp_gold_sens)

V_matrix_rsp <- diag(rsp_es$vi)

# comb
comb_es <- escalc(
  measure = "PLO",
  xi = case,
  ni = sample_size,
  data = comb_gold_sens)

V_matrix_comb <- diag(comb_es$vi)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. MODEL BUILDING ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
{                                # Compile here to finish model building in one step
#
#* 2.1 MLMA modelling ------------------------------------------------------------------------------
#
#** 2.1.1 Cluster variable - WHO geographical region ###############################################
# 
#*** 2.1.1.1 prism ####
prism_region_gold_sens_mlma_model <- rma.mv(
  yi = yi,                         # Effect size (logit transformed prevalence)
  V = V_matrix_prism,              # Variance of the effect size
  random = ~ 1 | geo_location/id,  # Random intercepts: region and study level
  data = prism_es,
  method = "REML"                  # Estimation method
)

# Generate prediction intervals
pred_prism_region <- predict(
  prism_region_gold_sens_mlma_model,
  pi.type = "normal"
)

# Back-transform logit estimates to normal prevalence
pred_prev_prism_region <- data.frame(
  Predicted_Prevalence = plogis(pred_prism_region$pred) * 100,
  CI_Lower = plogis(pred_prism_region$ci.lb) * 100,
  CI_Upper = plogis(pred_prism_region$ci.ub) * 100,
  PI_Lower = plogis(pred_prism_region$pi.lb) * 100,
  PI_Upper = plogis(pred_prism_region$pi.ub) * 100
)

#*** 2.1.1.2 rsp ####
rsp_region_gold_sens_mlma_model <- rma.mv(
  yi = yi,
  V = V_matrix_rsp,
  random = ~ 1 | geo_location/id,
  data = rsp_es,
  method = "REML"
)

# Generate prediction intervals
pred_rsp_region <- predict(
  rsp_region_gold_sens_mlma_model,
  pi.type = "normal"
)

# Back-transform logit estimates to normal prevalence
pred_prev_rsp_region <- data.frame(
  Predicted_Prevalence = plogis(pred_rsp_region$pred) * 100,
  CI_Lower = plogis(pred_rsp_region$ci.lb) * 100,
  CI_Upper = plogis(pred_rsp_region$ci.ub) * 100,
  PI_Lower = plogis(pred_rsp_region$pi.lb) * 100,
  PI_Upper = plogis(pred_rsp_region$pi.ub) * 100
)

#*** 2.1.1.3 combined ####
comb_region_gold_sens_mlma_model <- rma.mv(
  yi = yi,
  V = V_matrix_comb,
  random = ~ 1 | geo_location/id,
  data = comb_es,
  method = "REML"
)

# Generate prediction intervals
pred_comb_region <- predict(
  comb_region_gold_sens_mlma_model,
  pi.type = "normal"
)

# Back-transform logit estimates to normal prevalence
pred_prev_comb_region <- data.frame(
  Predicted_Prevalence = plogis(pred_comb_region$pred) * 100,
  CI_Lower = plogis(pred_comb_region$ci.lb) * 100,
  CI_Upper = plogis(pred_comb_region$ci.ub) * 100,
  PI_Lower = plogis(pred_comb_region$pi.lb) * 100,
  PI_Upper = plogis(pred_comb_region$pi.ub) * 100
)

#*** 2.1.1.4 Extracting WHO region subgroup prevalence from MLMA model =============================
#
# 2.1.1.4.1 prism regional prevalence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Extract random effects at geo_location level and fixed intercept
geo_location_prism_effects <- ranef(prism_region_gold_sens_mlma_model)$geo_location
fixed_intercept_prism_region <- prism_region_gold_sens_mlma_model$b[1, 1]

# Calculate predicted logit values
predicted_logit_prism_region <- fixed_intercept_prism_region + geo_location_prism_effects$intrcpt

# 95% CI and PI calculation
z_critical <- qnorm(0.975)  # value of 1.96 assigned for 95% CI

se_random_effects_prism_region <- geo_location_prism_effects$se
ci_lower_logit_prism_region <- predicted_logit_prism_region - z_critical * se_random_effects_prism_region
ci_upper_logit_prism_region <- predicted_logit_prism_region + z_critical * se_random_effects_prism_region

# For PI
residual_variance <- prism_region_gold_sens_mlma_model$sigma2[1]
se_pred_interval <- sqrt(se_random_effects_prism_region^2 + residual_variance)

pi_lower_logit_prism_region <- predicted_logit_prism_region - z_critical * se_pred_interval
pi_upper_logit_prism_region <- predicted_logit_prism_region + z_critical * se_pred_interval

# Back-transform prevalence predictions and CIs
predicted_prevalence_prism_region <- plogis(predicted_logit_prism_region)
ci_lower_prevalence_prism_region <- plogis(ci_lower_logit_prism_region)
ci_upper_prevalence_prism_region <- plogis(ci_upper_logit_prism_region)
pi_lower_prevalence_prism_region <- plogis(pi_lower_logit_prism_region)
pi_upper_prevalence_prism_region <- plogis(pi_upper_logit_prism_region)

# Data frame creation for results
predictions_df_prism_region <- data.frame(
  geo_location = rownames(geo_location_prism_effects),
  predicted_prevalence_prism_region = predicted_prevalence_prism_region,
  ci_lower_prevalence_prism_region = ci_lower_prevalence_prism_region,
  ci_upper_prevalence_prism_region = ci_upper_prevalence_prism_region,
  pi_lower_prevalence_prism_region = pi_lower_prevalence_prism_region,
  pi_upper_prevalence_prism_region = pi_upper_prevalence_prism_region
)

# Store prevalence estimates descending order of prevalence
predictions_df_prism_region <- predictions_df_prism_region[order(-predictions_df_prism_region$predicted_prevalence_prism_region), ]

# 2.1.1.4.2 rsp regional prevalence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Extract random effects at geo_location level and fixed intercept
geo_location_rsp_effects <- ranef(rsp_region_gold_sens_mlma_model)$geo_location
fixed_intercept_rsp_region <- rsp_region_gold_sens_mlma_model$b[1, 1]

# Calculate predicted logit values
predicted_logit_rsp_region <- fixed_intercept_rsp_region + geo_location_rsp_effects$intrcpt

# 95% CI and PI calculation
z_critical <- qnorm(0.975)  # value of 1.96 assigned for 95% CI

se_random_effects_rsp_region <- geo_location_rsp_effects$se
ci_lower_logit_rsp_region <- predicted_logit_rsp_region - z_critical * se_random_effects_rsp_region
ci_upper_logit_rsp_region <- predicted_logit_rsp_region + z_critical * se_random_effects_rsp_region

# For PI
residual_variance <- rsp_region_gold_sens_mlma_model$sigma2[1]
se_pred_interval <- sqrt(se_random_effects_rsp_region^2 + residual_variance)

pi_lower_logit_rsp_region <- predicted_logit_rsp_region - z_critical * se_pred_interval
pi_upper_logit_rsp_region <- predicted_logit_rsp_region + z_critical * se_pred_interval

# Back-transform prevalence predictions and CIs
predicted_prevalence_rsp_region <- plogis(predicted_logit_rsp_region)
ci_lower_prevalence_rsp_region <- plogis(ci_lower_logit_rsp_region)
ci_upper_prevalence_rsp_region <- plogis(ci_upper_logit_rsp_region)
pi_lower_prevalence_rsp_region <- plogis(pi_lower_logit_rsp_region)
pi_upper_prevalence_rsp_region <- plogis(pi_upper_logit_rsp_region)

# Data frame creation for results
predictions_df_rsp_region <- data.frame(
  geo_location = rownames(geo_location_rsp_effects),
  predicted_prevalence_rsp_region = predicted_prevalence_rsp_region,
  ci_lower_prevalence_rsp_region = ci_lower_prevalence_rsp_region,
  ci_upper_prevalence_rsp_region = ci_upper_prevalence_rsp_region,
  pi_lower_prevalence_rsp_region = pi_lower_prevalence_rsp_region,
  pi_upper_prevalence_rsp_region = pi_upper_prevalence_rsp_region
)

# Store prevalence estimates descending order of prevalence
predictions_df_rsp_region <- predictions_df_rsp_region[order(-predictions_df_rsp_region$predicted_prevalence_rsp_region), ]

# 2.1.1.4.3 combined regional prevalence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#
# Extract random effects at geo_location level and fixed intercept
geo_location_comb_effects <- ranef(comb_region_gold_sens_mlma_model)$geo_location
fixed_intercept_comb_region <- comb_region_gold_sens_mlma_model$b[1, 1]

# Calculate predicted logit values
predicted_logit_comb_region <- fixed_intercept_comb_region + geo_location_comb_effects$intrcpt

# 95% CI and PI calculation
z_critical <- qnorm(0.975)  # value of 1.96 assigned for 95% CI

se_random_effects_comb_region <- geo_location_comb_effects$se
ci_lower_logit_comb_region <- predicted_logit_comb_region - z_critical * se_random_effects_comb_region
ci_upper_logit_comb_region <- predicted_logit_comb_region + z_critical * se_random_effects_comb_region

# For PI
residual_variance <- comb_region_gold_sens_mlma_model$sigma2[1]
se_pred_interval <- sqrt(se_random_effects_comb_region^2 + residual_variance)

pi_lower_logit_comb_region <- predicted_logit_comb_region - z_critical * se_pred_interval
pi_upper_logit_comb_region <- predicted_logit_comb_region + z_critical * se_pred_interval

# Back-transform prevalence predictions and CIs
predicted_prevalence_comb_region <- plogis(predicted_logit_comb_region)
ci_lower_prevalence_comb_region <- plogis(ci_lower_logit_comb_region)
ci_upper_prevalence_comb_region <- plogis(ci_upper_logit_comb_region)
pi_lower_prevalence_comb_region <- plogis(pi_lower_logit_comb_region)
pi_upper_prevalence_comb_region <- plogis(pi_upper_logit_comb_region)

# Data frame creation for results
predictions_df_comb_region <- data.frame(
  geo_location = rownames(geo_location_comb_effects),
  predicted_prevalence_comb_region = predicted_prevalence_comb_region,
  ci_lower_prevalence_comb_region = ci_lower_prevalence_comb_region,
  ci_upper_prevalence_comb_region = ci_upper_prevalence_comb_region,
  pi_lower_prevalence_comb_region = pi_lower_prevalence_comb_region,
  pi_upper_prevalence_comb_region = pi_upper_prevalence_comb_region
)

# Store prevalence estimates descending order of prevalence
predictions_df_comb_region <- predictions_df_comb_region[order(-predictions_df_comb_region$predicted_prevalence_comb_region), ]

#** 2.1.2 Cluster variable - World Bank income-level ###############################################
#
#*** 2.1.2.1 prism ====
prism_income_gold_sens_mlma_model <- rma.mv(
  yi = yi,
  V = V_matrix_prism,
  random = ~ 1 | income_class/id,
  data = prism_es,
  method = "REML"
)

# Generate prediction intervals
pred_prism_income <- predict(
  prism_income_gold_sens_mlma_model,
  pi.type = "normal"
)

# Back-transform logit estimates to normal prevalence
pred_prev_prism_income <- data.frame(
  Predicted_Prevalence = plogis(pred_prism_income$pred) * 100,
  CI_Lower = plogis(pred_prism_income$ci.lb) * 100,
  CI_Upper = plogis(pred_prism_income$ci.ub) * 100,
  PI_Lower = plogis(pred_prism_income$pi.lb) * 100,
  PI_Upper = plogis(pred_prism_income$pi.ub) * 100
)

#*** 2.1.2.2 rsp ====
rsp_income_gold_sens_mlma_model <- rma.mv(
  yi = yi,
  V = V_matrix_rsp,
  random = ~ 1 | income_class/id,
  data = rsp_es,
  method = "REML"
)

# Generate prediction intervals
pred_rsp_income <- predict(
  rsp_income_gold_sens_mlma_model,
  pi.type = "normal"
)

# Back-transform logit estimates to normal prevalence
pred_prev_rsp_income <- data.frame(
  Predicted_Prevalence = plogis(pred_rsp_income$pred) * 100,
  CI_Lower = plogis(pred_rsp_income$ci.lb) * 100,
  CI_Upper = plogis(pred_rsp_income$ci.ub) * 100,
  PI_Lower = plogis(pred_rsp_income$pi.lb) * 100,
  PI_Upper = plogis(pred_rsp_income$pi.ub) * 100
)

#*** 2.1.2.3 combined ====
comb_income_gold_sens_mlma_model <- rma.mv(
  yi = yi,
  V = V_matrix_comb,
  random = ~ 1 | income_class/id,
  data = comb_es,
  method = "REML"
)

# Generate prediction intervals
pred_comb_income <- predict(
  comb_income_gold_sens_mlma_model,
  pi.type = "normal"
)

# Back-transform logit estimates to normal prevalence
pred_prev_comb_income <- data.frame(
  Predicted_Prevalence = plogis(pred_comb_income$pred) * 100,
  CI_Lower = plogis(pred_comb_income$ci.lb) * 100,
  CI_Upper = plogis(pred_comb_income$ci.ub) * 100,
  PI_Lower = plogis(pred_comb_income$pi.lb) * 100,
  PI_Upper = plogis(pred_comb_income$pi.ub) * 100
)

#*** 2.1.2.4 Extracting World Bank income level subgroup prevalence from MLMA model ================
#
# 2.1.2.4.1 prism income prevalence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Extract random effects at income_class level and fixed intercept
income_class_prism_effects <- ranef(prism_income_gold_sens_mlma_model)$income_class
fixed_intercept_prism_income <- prism_income_gold_sens_mlma_model$b[1, 1]

# Calculate predicted logit values
predicted_logit_prism_income <- fixed_intercept_prism_income + income_class_prism_effects$intrcpt

# 95% CI and PI calculation
z_critical <- qnorm(0.975)  # value of 1.96 assigned for 95% CI

se_random_effects_prism_income <- income_class_prism_effects$se
ci_lower_logit_prism_income <- predicted_logit_prism_income - z_critical * se_random_effects_prism_income
ci_upper_logit_prism_income <- predicted_logit_prism_income + z_critical * se_random_effects_prism_income

# For PI
residual_variance <- prism_income_gold_sens_mlma_model$sigma2[1]
se_pred_interval <- sqrt(se_random_effects_prism_income^2 + residual_variance)

pi_lower_logit_prism_income <- predicted_logit_prism_income - z_critical * se_pred_interval
pi_upper_logit_prism_income <- predicted_logit_prism_income + z_critical * se_pred_interval

# Back-transform prevalence predictions and CIs
predicted_prevalence_prism_income <- plogis(predicted_logit_prism_income)
ci_lower_prevalence_prism_income <- plogis(ci_lower_logit_prism_income)
ci_upper_prevalence_prism_income <- plogis(ci_upper_logit_prism_income)
pi_lower_prevalence_prism_income <- plogis(pi_lower_logit_prism_income)
pi_upper_prevalence_prism_income <- plogis(pi_upper_logit_prism_income)

# Data frame creation for results
predictions_df_prism_income <- data.frame(
  income_class = rownames(income_class_prism_effects),
  predicted_prevalence_prism_income = predicted_prevalence_prism_income,
  ci_lower_prevalence_prism_income = ci_lower_prevalence_prism_income,
  ci_upper_prevalence_prism_income = ci_upper_prevalence_prism_income,
  pi_lower_prevalence_prism_income = pi_lower_prevalence_prism_income,
  pi_upper_prevalence_prism_income = pi_upper_prevalence_prism_income
)

# Store prevalence estimates descending order of prevalence
predictions_df_prism_income <- predictions_df_prism_income[order(-predictions_df_prism_income$predicted_prevalence_prism_income), ]

# 2.1.2.4.1 rsp income prevalence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Extract random effects at income_class level and fixed intercept
income_class_rsp_effects <- ranef(rsp_income_gold_sens_mlma_model)$income_class
fixed_intercept_rsp_income <- rsp_income_gold_sens_mlma_model$b[1, 1]

# Calculate predicted logit values
predicted_logit_rsp_income <- fixed_intercept_rsp_income + income_class_rsp_effects$intrcpt

# 95% CI and PI calculation
z_critical <- qnorm(0.975)  # value of 1.96 assigned for 95% CI

se_random_effects_rsp_income <- income_class_rsp_effects$se
ci_lower_logit_rsp_income <- predicted_logit_rsp_income - z_critical * se_random_effects_rsp_income
ci_upper_logit_rsp_income <- predicted_logit_rsp_income + z_critical * se_random_effects_rsp_income

# For PI
residual_variance <- rsp_income_gold_sens_mlma_model$sigma2[1]
se_pred_interval <- sqrt(se_random_effects_rsp_income^2 + residual_variance)

pi_lower_logit_rsp_income <- predicted_logit_rsp_income - z_critical * se_pred_interval
pi_upper_logit_rsp_income <- predicted_logit_rsp_income + z_critical * se_pred_interval

# Back-transform prevalence predictions and CIs
predicted_prevalence_rsp_income <- plogis(predicted_logit_rsp_income)
ci_lower_prevalence_rsp_income <- plogis(ci_lower_logit_rsp_income)
ci_upper_prevalence_rsp_income <- plogis(ci_upper_logit_rsp_income)
pi_lower_prevalence_rsp_income <- plogis(pi_lower_logit_rsp_income)
pi_upper_prevalence_rsp_income <- plogis(pi_upper_logit_rsp_income)

# Data frame creation for results
predictions_df_rsp_income <- data.frame(
  income_class = rownames(income_class_rsp_effects),
  predicted_prevalence_rsp_income = predicted_prevalence_rsp_income,
  ci_lower_prevalence_rsp_income = ci_lower_prevalence_rsp_income,
  ci_upper_prevalence_rsp_income = ci_upper_prevalence_rsp_income,
  pi_lower_prevalence_rsp_income = pi_lower_prevalence_rsp_income,
  pi_upper_prevalence_rsp_income = pi_upper_prevalence_rsp_income
)

# Store prevalence estimates descending order of prevalence
predictions_df_rsp_income <- predictions_df_rsp_income[order(-predictions_df_rsp_income$predicted_prevalence_rsp_income), ]

# 2.1.2.4.1 combined income prevalence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Extract random effects at income_class level and fixed intercept
income_class_comb_effects <- ranef(comb_income_gold_sens_mlma_model)$income_class
fixed_intercept_comb_income <- comb_income_gold_sens_mlma_model$b[1, 1]

# Calculate predicted logit values
predicted_logit_comb_income <- fixed_intercept_comb_income + income_class_comb_effects$intrcpt

# 95% CI and PI calculation
z_critical <- qnorm(0.975)  # value of 1.96 assigned for 95% CI

se_random_effects_comb_income <- income_class_comb_effects$se
ci_lower_logit_comb_income <- predicted_logit_comb_income - z_critical * se_random_effects_comb_income
ci_upper_logit_comb_income <- predicted_logit_comb_income + z_critical * se_random_effects_comb_income

# For PI
residual_variance <- comb_income_gold_sens_mlma_model$sigma2[1]
se_pred_interval <- sqrt(se_random_effects_comb_income^2 + residual_variance)

pi_lower_logit_comb_income <- predicted_logit_comb_income - z_critical * se_pred_interval
pi_upper_logit_comb_income <- predicted_logit_comb_income + z_critical * se_pred_interval

# Back-transform prevalence predictions and CIs
predicted_prevalence_comb_income <- plogis(predicted_logit_comb_income)
ci_lower_prevalence_comb_income <- plogis(ci_lower_logit_comb_income)
ci_upper_prevalence_comb_income <- plogis(ci_upper_logit_comb_income)
pi_lower_prevalence_comb_income <- plogis(pi_lower_logit_comb_income)
pi_upper_prevalence_comb_income <- plogis(pi_upper_logit_comb_income)

# Data frame creation for results
predictions_df_comb_income <- data.frame(
  income_class = rownames(income_class_comb_effects),
  predicted_prevalence_comb_income = predicted_prevalence_comb_income,
  ci_lower_prevalence_comb_income = ci_lower_prevalence_comb_income,
  ci_upper_prevalence_comb_income = ci_upper_prevalence_comb_income,
  pi_lower_prevalence_comb_income = pi_lower_prevalence_comb_income,
  pi_upper_prevalence_comb_income = pi_upper_prevalence_comb_income
)

# Store prevalence estimates descending order of prevalence
predictions_df_comb_income <- predictions_df_comb_income[order(-predictions_df_comb_income$predicted_prevalence_comb_income), ]

#* 2.2 Random-effects modelling --------------------------------------------------------------------
#
#** 2.2.1 Overall prevalence - standard random-effects model (for supplementary) ###################
#
#*** 2.2.1.1 prism =================================================================================
meta_prism_gold_sens <- prism_gold_sens %>%
  metaprop(
    event = case,                       # Column for number of cases             
    n = sample_size,                    # Column for total sample size           
    studlab = paste0(author),           # Label studies by the author names 
    method = "Inverse",                 # Inverse variance method       
    sm = "PLOGIT",                      # Use logit-transformed proportions            
    random = TRUE,                      # Use random effects model            
    common = FALSE,                     # Do not report common (fixed) effects model           
    warn = TRUE,                        # Display warnings              
    prediction = TRUE,                  # Report prediction intervals        
    pscale = 100)                       # Increase proportion by 100x (converting prevalence to percentage)

#*** 2.2.1.2 rsp ===================================================================================
meta_rsp_gold_sens <- rsp_gold_sens %>%
  metaprop(
    event = case,            
    n = sample_size,          
    studlab = paste0(author), 
    method = "Inverse",       
    sm = "PLOGIT",           
    random = TRUE,   
    common = FALSE,
    warn = TRUE,
    prediction = TRUE,
    pscale=100)

#*** 2.2.1.3 combined ==============================================================================
meta_comb_gold_sens <- comb_gold_sens %>%
  metaprop(
    event = case,
    n = sample_size, 
    studlab = paste0(author),
    method = "Inverse",
    sm = "PLOGIT",
    random = TRUE,
    common = FALSE,
    warn = TRUE,
    prediction = TRUE,
    pscale=100)

#** 2.2.2 Subgroup prevalence - sex ################################################################
#
#*** 2.2.2.1 Males =================================================================================
#
# prism
meta_male_gold_sens_prism <-male_gold_sens_prism %>%
  metaprop(
    event = male_cases,            
    n = male_total,         
    studlab = paste0(author),
    method = "Inverse",    
    sm = "PLOGIT",        
    random = TRUE,
    common = FALSE,
    warn = TRUE,
    prediction = TRUE,
    pscale=100)

# rsp
meta_male_gold_sens_rsp <-male_gold_sens_rsp %>%
  metaprop(
    event = male_cases,            
    n = male_total,         
    studlab = paste0(author),
    method = "Inverse",    
    sm = "PLOGIT",        
    random = TRUE,
    common = FALSE,
    warn = TRUE,
    prediction = TRUE,
    pscale=100)

# combined
meta_male_gold_sens_comb <-male_gold_sens_comb %>%
  metaprop(
    event = male_cases,            
    n = male_total,         
    studlab = paste0(author),
    method = "Inverse",    
    sm = "PLOGIT",        
    random = TRUE,
    common = FALSE,
    warn = TRUE,
    prediction = TRUE,
    pscale=100)

#*** 2.2.2.2 Females ===============================================================================
#
# prism
meta_female_gold_sens_prism <-female_gold_sens_prism %>%
  metaprop(
    event = female_cases,            
    n = female_total,         
    studlab = paste0(author),
    method = "Inverse",    
    sm = "PLOGIT",        
    random = TRUE,   
    common = FALSE,
    warn = TRUE,
    prediction = TRUE,
    pscale=100)

# rsp
meta_female_gold_sens_rsp <-female_gold_sens_rsp %>%
  metaprop(
    event = female_cases,            
    n = female_total,         
    studlab = paste0(author),
    method = "Inverse",    
    sm = "PLOGIT",        
    random = TRUE,   
    common = FALSE,
    warn = TRUE,
    prediction = TRUE,
    pscale=100)

# combined
meta_female_gold_sens_comb <-female_gold_sens_comb %>%
  metaprop(
    event = female_cases,            
    n = female_total,         
    studlab = paste0(author),
    method = "Inverse",    
    sm = "PLOGIT",        
    random = TRUE,   
    common = FALSE,
    warn = TRUE,
    prediction = TRUE,
    pscale=100)

#** 2.2.3 Subgroup prevalence - smoking habits #####################################################
#
#*** 2.2.3.1 Non-smokers ===========================================================================
#
# prism
meta_non_smoker_gold_sens_prism <- non_smoker_gold_sens_prism %>% 
  metaprop(
    event = case_non,            
    n = tot_non,         
    studlab = paste0(author),
    method = "Inverse",    
    sm = "PLOGIT",        
    random = TRUE,   
    common = FALSE,
    warn = TRUE,
    prediction = TRUE,
    pscale=100)

# rsp
meta_non_smoker_gold_sens_rsp <- non_smoker_gold_sens_rsp %>% 
  metaprop(
    event = case_non,            
    n = tot_non,         
    studlab = paste0(author),
    method = "Inverse",    
    sm = "PLOGIT",        
    random = TRUE,   
    common = FALSE,
    warn = TRUE,
    prediction = TRUE,
    pscale=100)

# combined
meta_non_smoker_gold_sens_comb <- non_smoker_gold_sens_comb %>% 
  metaprop(
    event = case_non,            
    n = tot_non,         
    studlab = paste0(author),
    method = "Inverse",    
    sm = "PLOGIT",        
    random = TRUE,   
    common = FALSE,
    warn = TRUE,
    prediction = TRUE,
    pscale=100)

#*** 2.2.3.2 Ex-smokers ============================================================================
#
# prism
meta_ex_smoker_gold_sens_prism <- ex_smoker_gold_sens_prism %>% 
  metaprop(
    event = case_ex,            
    n = tot_ex,         
    studlab = paste0(author),
    method = "Inverse",    
    sm = "PLOGIT",        
    random = TRUE,   
    common = FALSE,
    warn = TRUE,
    prediction = TRUE,
    pscale=100)

# rsp
meta_ex_smoker_gold_sens_rsp <- ex_smoker_gold_sens_rsp %>% 
  metaprop(
    event = case_ex,            
    n = tot_ex,         
    studlab = paste0(author),
    method = "Inverse",    
    sm = "PLOGIT",        
    random = TRUE,   
    common = FALSE,
    warn = TRUE,
    prediction = TRUE,
    pscale=100)

# combined
meta_ex_smoker_gold_sens_comb <- ex_smoker_gold_sens_comb %>% 
  metaprop(
    event = case_ex,            
    n = tot_ex,         
    studlab = paste0(author),
    method = "Inverse",    
    sm = "PLOGIT",        
    random = TRUE,   
    common = FALSE,
    warn = TRUE,
    prediction = TRUE,
    pscale=100)

#*** 2.2.3.3 Current smokers =======================================================================
#
# prism
meta_cur_smoker_gold_sens_prism <- cur_smoker_gold_sens_prism %>% 
  metaprop(
    event = case_cur,            
    n = tot_cur,         
    studlab = paste0(author),
    method = "Inverse",    
    sm = "PLOGIT",        
    random = TRUE,  
    common = FALSE,
    warn = TRUE,
    prediction = TRUE,
    pscale=100)

# rsp
meta_cur_smoker_gold_sens_rsp <- cur_smoker_gold_sens_rsp %>% 
  metaprop(
    event = case_cur,            
    n = tot_cur,         
    studlab = paste0(author),
    method = "Inverse",    
    sm = "PLOGIT",        
    random = TRUE,  
    common = FALSE,
    warn = TRUE,
    prediction = TRUE,
    pscale=100)

# combined
meta_cur_smoker_gold_sens_comb <- cur_smoker_gold_sens_comb %>% 
  metaprop(
    event = case_cur,            
    n = tot_cur,         
    studlab = paste0(author),
    method = "Inverse",    
    sm = "PLOGIT",        
    random = TRUE,  
    common = FALSE,
    warn = TRUE,
    prediction = TRUE,
    pscale=100)
}
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. MODEL SUMMARY CHECKING ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# prism region
summary(prism_region_gold_sens_mlma_model)
# Extract total and level-specific I2 values
i2_prism_region_gold_sens_mlma_model<-var.comp(prism_region_gold_sens_mlma_model)
summary(i2_prism_region_gold_sens_mlma_model)

# rsp region
summary(rsp_region_gold_sens_mlma_model)
i2_rsp_region_gold_sens_mlma_model<-var.comp(rsp_region_gold_sens_mlma_model)
summary(i2_rsp_region_gold_sens_mlma_model)

# comb region
summary(comb_region_gold_sens_mlma_model)
i2_comb_region_gold_sens_mlma_model<-var.comp(comb_region_gold_sens_mlma_model)
summary(i2_comb_region_gold_sens_mlma_model)

# prism income
summary(prism_income_gold_sens_mlma_model)
i2_prism_income_gold_sens_mlma_model<-var.comp(prism_income_gold_sens_mlma_model)
summary(i2_prism_income_gold_sens_mlma_model)

# rsp income
summary(rsp_income_gold_sens_mlma_model)
i2_rsp_income_gold_sens_mlma_model<-var.comp(rsp_income_gold_sens_mlma_model)
summary(i2_rsp_income_gold_sens_mlma_model)

# comb income
summary(comb_income_gold_sens_mlma_model)
i2_comb_income_gold_sens_mlma_model<-var.comp(comb_income_gold_sens_mlma_model)
summary(i2_comb_income_gold_sens_mlma_model)

# overall standard random-effects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# prism
summary(meta_prism_gold_sens)

forest(meta_prism_gold_sens,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"),
       leftlabs = c("Study","Cases","Total"))

# rsp
summary(meta_rsp_gold_sens)

forest(meta_rsp_gold_sens,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"),
       leftlabs = c("Study","Cases","Total"))

# combined
summary(meta_comb_gold_sens)

forest(meta_comb_gold_sens,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"),
       leftlabs = c("Study","Cases","Total"))

# Gender - males ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# males - prism
summary(meta_male_gold_sens_prism)

forest(meta_male_gold_sens_prism,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"),
       leftlabs = c("Study","Cases","Total"))

# males - rsp
summary(meta_male_gold_sens_rsp)

forest(meta_male_gold_sens_rsp,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"),
       leftlabs = c("Study","Cases","Total"))

# men - combined
summary(meta_male_gold_sens_comb)

forest(meta_male_gold_sens_comb,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"),
       leftlabs = c("Study","Cases","Total"))

# Gender - females
# females - prism
summary(meta_female_gold_sens_prism)

forest(meta_female_gold_sens_prism,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"),
       leftlabs = c("Study","Cases","Total"))

# females - rsp
summary(meta_female_gold_sens_rsp)

forest(meta_female_gold_sens_rsp,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"),
       leftlabs = c("Study","Cases","Total"))

# females - combined
summary(meta_female_gold_sens_comb)

forest(meta_female_gold_sens_comb,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"),
       leftlabs = c("Study","Cases","Total"))

# Smoking habit ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Non-smokers
summary(meta_non_smoker_gold_sens_prism)

forest(meta_non_smoker_gold_sens_prism,
        rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"),
        leftlabs = c("Study","Cases","Total"))

summary(meta_non_smoker_gold_sens_rsp)

forest(meta_non_smoker_gold_sens_rsp,
        rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"),
        leftlabs = c("Study","Cases","Total"))

summary(meta_non_smoker_gold_sens_comb)

forest(meta_non_smoker_gold_sens_comb,
        rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"),
        leftlabs = c("Study","Cases","Total"))
 
# ex-smokers
summary(meta_ex_smoker_gold_sens_prism)

forest(meta_ex_smoker_gold_sens_prism,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"),
       leftlabs = c("Study","Cases","Total"))

summary(meta_ex_smoker_gold_sens_rsp)

forest(meta_ex_smoker_gold_sens_rsp,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"),
       leftlabs = c("Study","Cases","Total"))
summary(meta_ex_smoker_gold_sens_comb)

forest(meta_ex_smoker_gold_sens_comb,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"),
       leftlabs = c("Study","Cases","Total"))

# current smokers
summary(meta_cur_smoker_gold_sens_prism)

forest(meta_cur_smoker_gold_sens_prism,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"),
       leftlabs = c("Study","Cases","Total"))

summary(meta_cur_smoker_gold_sens_rsp)

forest(meta_cur_smoker_gold_sens_rsp,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"),
       leftlabs = c("Study","Cases","Total"))

summary(meta_cur_smoker_gold_sens_comb)

forest(meta_cur_smoker_gold_sens_comb,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"),
       leftlabs = c("Study","Cases","Total"))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. OUTPUT DISPLAY ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#* 4.1 Overall prevalence from MLMA models ---------------------------------------------------------
#
{                        # Compile here to generate table in one step
#
overall_prevalence_table <- bind_rows(
  data.frame(Type = "prism region", pred_prev_prism_region),
  data.frame(Type = "rsp region", pred_prev_rsp_region),
  data.frame(Type = "combined region", pred_prev_comb_region),
  data.frame(Type = "prism income", pred_prev_prism_income),
  data.frame(Type = "rsp income", pred_prev_rsp_income),
  data.frame(Type = "combined income", pred_prev_comb_income)
)

overall_prevalence_table_rounded <- overall_prevalence_table %>%
  mutate(across(c(Predicted_Prevalence, CI_Lower, CI_Upper, PI_Lower, PI_Upper), ~round(.x, 2)))

# Create formatted columns
overall_prevalence_table_formatted <- overall_prevalence_table %>%
  mutate(
    Prevalence_CI = sprintf("%.2f (%.2f-%.2f)", Predicted_Prevalence, CI_Lower, CI_Upper),
    Prediction_Interval = sprintf("%.2f-%.2f", PI_Lower, PI_Upper)
  ) %>%
  select(Type, Prevalence_CI, Prediction_Interval)

# Print the formatted result
print(overall_prevalence_table_formatted)
}

# who prism
# Clean the column names and create a table for prevalence by multiplying values *100 and rounding off to two decimal places
predictions_df_prism_region <- predictions_df_prism_region %>%
  rename(
    Region = geo_location,
    Prevalence = predicted_prevalence_prism_region,
    CI_Lower = ci_lower_prevalence_prism_region,
    CI_Upper = ci_upper_prevalence_prism_region,
    PI_Lower = pi_lower_prevalence_prism_region,
    PI_Upper = pi_upper_prevalence_prism_region
  ) %>%
  mutate(across(c(Prevalence, CI_Lower, CI_Upper, PI_Lower, PI_Upper), ~ round(.x * 100, 2))) %>%
  mutate(
    Prevalence_CI = paste0(Prevalence, " (", CI_Lower, "-", CI_Upper, ") ",
                           PI_Lower, "-", PI_Upper)
  ) %>%
  select(Region, Prevalence_CI)

print(predictions_df_prism_region)

# who rsp
# Clean the column names and create a table for prevalence by multiplying values *100 and rounding off to two decimal places
predictions_df_rsp_region <- predictions_df_rsp_region %>%
  rename(
    Region = geo_location,
    Prevalence = predicted_prevalence_rsp_region,
    CI_Lower = ci_lower_prevalence_rsp_region,
    CI_Upper = ci_upper_prevalence_rsp_region,
    PI_Lower = pi_lower_prevalence_rsp_region,
    PI_Upper = pi_upper_prevalence_rsp_region
  ) %>%
  mutate(across(c(Prevalence, CI_Lower, CI_Upper, PI_Lower, PI_Upper), ~ round(.x * 100, 2))) %>%
  mutate(
    Prevalence_CI = paste0(Prevalence, " (", CI_Lower, "-", CI_Upper, ") ",
                              PI_Lower, "-", PI_Upper)
  ) %>%
  select(Region, Prevalence_CI)

print(predictions_df_rsp_region)

# who combined
# Clean the column names and create a table for prevalence by multiplying values *100 and rounding off to two decimal places
predictions_df_comb_region <- predictions_df_comb_region %>%
  rename(
    Region = geo_location,
    Prevalence = predicted_prevalence_comb_region,
    CI_Lower = ci_lower_prevalence_comb_region,
    CI_Upper = ci_upper_prevalence_comb_region,
    PI_Lower = pi_lower_prevalence_comb_region,
    PI_Upper = pi_upper_prevalence_comb_region
  ) %>%
  mutate(across(c(Prevalence, CI_Lower, CI_Upper, PI_Lower, PI_Upper), ~ round(.x * 100, 2))) %>%
  mutate(
    Prevalence_CI = paste0(Prevalence, " (", CI_Lower, "-", CI_Upper, ") ",
                              PI_Lower, "-", PI_Upper)
  ) %>%
  select(Region, Prevalence_CI)

print(predictions_df_comb_region)

# income prism
# Clean the column names and create a table for prevalence by multiplying values *100 and rounding off to two decimal places
predictions_df_prism_income <- predictions_df_prism_income %>%
  rename(
    Income = income_class,
    Prevalence = predicted_prevalence_prism_income,
    CI_Lower = ci_lower_prevalence_prism_income,
    CI_Upper = ci_upper_prevalence_prism_income,
    PI_Lower = pi_lower_prevalence_prism_income,
    PI_Upper = pi_upper_prevalence_prism_income
  ) %>%
  mutate(across(c(Prevalence, CI_Lower, CI_Upper, PI_Lower, PI_Upper), ~ round(.x * 100, 2))) %>%
  mutate(
    Prevalence_CI = paste0(Prevalence, " (", CI_Lower, "-", CI_Upper, ") ",
                              PI_Lower, "-", PI_Upper)
  ) %>%
  select(Income, Prevalence_CI)

print(predictions_df_prism_income)

# income rsp
# Clean the column names and create a table for prevalence by multiplying values *100 and rounding off to two decimal places
predictions_df_rsp_income <- predictions_df_rsp_income %>%
  rename(
    Income = income_class,
    Prevalence = predicted_prevalence_rsp_income,
    CI_Lower = ci_lower_prevalence_rsp_income,
    CI_Upper = ci_upper_prevalence_rsp_income,
    PI_Lower = pi_lower_prevalence_rsp_income,
    PI_Upper = pi_upper_prevalence_rsp_income
  ) %>%
  mutate(across(c(Prevalence, CI_Lower, CI_Upper, PI_Lower, PI_Upper), ~ round(.x * 100, 2))) %>%
  mutate(
    Prevalence_CI = paste0(Prevalence, " (", CI_Lower, "-", CI_Upper, ") ",
                              PI_Lower, "-", PI_Upper)
  ) %>%
  select(Income, Prevalence_CI)

print(predictions_df_rsp_income)

# income combined
# Clean the column names and create a table for prevalence by multiplying values *100 and rounding off to two decimal places
predictions_df_comb_income <- predictions_df_comb_income %>%
  rename(
    Income = income_class,
    Prevalence = predicted_prevalence_comb_income,
    CI_Lower = ci_lower_prevalence_comb_income,
    CI_Upper = ci_upper_prevalence_comb_income,
    PI_Lower = pi_lower_prevalence_comb_income,
    PI_Upper = pi_upper_prevalence_comb_income
  ) %>%
  mutate(across(c(Prevalence, CI_Lower, CI_Upper, PI_Lower, PI_Upper), ~ round(.x * 100, 2))) %>%
  mutate(
    Prevalence_CI = paste0(Prevalence, " (", CI_Lower, "-", CI_Upper, ") ",
                              PI_Lower, "-", PI_Upper)
  ) %>%
  select(Income, Prevalence_CI)

print(predictions_df_comb_income)

#* 4.2 All prevalences from random-effects models --------------------------------------------------
#
{              # Compile here to generate table in one step
#
# Inverse logit transformation for proportions
inv_logit <- function(x) {
 plogis(x) * 100
  }
  
# Create a function to extract data from models
extract_prevalence_results <- function(model, data, model_name) {
k <- model$k
    
# Back-transform from logit to percentage
est <- round(inv_logit(model$TE.random), 2)
ci_lb <- round(inv_logit(model$lower.random), 2)
ci_ub <- round(inv_logit(model$upper.random), 2)
    
pi_lb <- round(inv_logit(model$lower.predict), 2)
pi_ub <- round(inv_logit(model$upper.predict), 2)
    
estimate_ci <- paste0(est, " (", ci_lb, "-", ci_ub, ")")
prediction_interval <- paste0(pi_lb, "-", pi_ub)

i2 <- round(model$I2*100, 1)
    
results <- data.frame(
  Model = model_name,
  k = k,
  I2 = i2,
  Estimate_CI = estimate_ci,
  Prediction_Interval = prediction_interval,
  row.names = NULL
)
    
return(results)
}
  
# Apply to each model
res01 <- extract_prevalence_results(meta_prism_gold_sens, prism_gold_sens, "prism gold_sens")
res02 <- extract_prevalence_results(meta_rsp_gold_sens, rsp_gold_sens, "rsp gold_sens")
res03 <- extract_prevalence_results(meta_comb_gold_sens, comb_gold_sens, "combined gold_sens")
res04 <- extract_prevalence_results(meta_male_gold_sens_prism, male_gold_sens_prism, "male prism")  
res05 <- extract_prevalence_results(meta_male_gold_sens_rsp, male_gold_sens_rsp, "male rsp")  
res06 <- extract_prevalence_results(meta_male_gold_sens_comb, male_gold_sens_comb, "male comb")  
res07 <- extract_prevalence_results(meta_female_gold_sens_prism, female_gold_sens_prism, "female prism")  
res08 <- extract_prevalence_results(meta_female_gold_sens_rsp, female_gold_sens_rsp, "female rsp")
res09 <- extract_prevalence_results(meta_female_gold_sens_comb, female_gold_sens_comb, "female comb")
res10 <- extract_prevalence_results(meta_non_smoker_gold_sens_prism, non_smoker_gold_sens_prism, "non-smoker prism")
res11 <- extract_prevalence_results(meta_non_smoker_gold_sens_rsp, non_smoker_gold_sens_rsp, "non-smoker rsp")
res12 <- extract_prevalence_results(meta_non_smoker_gold_sens_comb, non_smoker_gold_sens_comb, "non-smoker comb")
res13 <- extract_prevalence_results(meta_ex_smoker_gold_sens_prism, ex_smoker_gold_sens_prism, "ex-smoker prism")
res14 <- extract_prevalence_results(meta_ex_smoker_gold_sens_rsp, ex_smoker_gold_sens_rsp, "ex-smoker rsp")
res15 <- extract_prevalence_results(meta_ex_smoker_gold_sens_comb, ex_smoker_gold_sens_comb, "ex-smoker comb")
res16 <- extract_prevalence_results(meta_cur_smoker_gold_sens_prism, cur_smoker_gold_sens_prism, "current smoker prism")
res17 <- extract_prevalence_results(meta_cur_smoker_gold_sens_rsp, cur_smoker_gold_sens_rsp, "current smoker rsp")
res18 <- extract_prevalence_results(meta_cur_smoker_gold_sens_comb, cur_smoker_gold_sens_comb, "current smoker comb")

# Combine into final table
random_effects_results_table <- bind_rows(res01, res02, res03, res04, res05, res06, res07, res08, res09, res10, res11, 
                                res12, res13, res14, res15, res16, res17, res18)
  
# View final table
print(random_effects_results_table)
  
# Clean up
rm(res01, res02, res03, res04, res05, res06, res07, res08, res09, res10, res11, res12, res13, res14, res15, 
     res16, res17, res18)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. PUBLICATION BIAS CHECKING ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get raw effect sizes (yi) and standard errors
prism_es$sei <- sqrt(diag(V_matrix_prism))
prism_es$precision <- 1 / prism_es$sei

rsp_es$sei <- sqrt(diag(V_matrix_rsp))
rsp_es$precision <- 1 / rsp_es$sei

comb_es$sei <- sqrt(diag(V_matrix_comb))
comb_es$precision <- 1 / comb_es$sei

#* 5.1 prism region MLMA region ====================================================================
#
# Funnel plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
# Get the pooled effect estimate from your MLMA model
mu_hat_prism_region <- as.numeric(coef(prism_region_gold_sens_mlma_model))

# Define SE range for drawing funnel lines
se_seq_prism_region <- seq(min(prism_es$sei), max(prism_es$sei), length.out = 100)
upper_prism_region <- mu_hat_prism_region + 1.96 * se_seq_prism_region
lower_prism_region <- mu_hat_prism_region - 1.96 * se_seq_prism_region

# 4. Create the funnel plot
plot(prism_es$yi, prism_es$sei,
     xlab = "Prevalence (logit-transformed)",
     ylab = "Standard Error",
     main = "Funnel plot for publication bias",
     pch = 19,
     ylim = rev(range(prism_es$sei)))  # Reverse y-axis

# Add funnel
lines(upper_prism_region, se_seq_prism_region, lty = 2, col = "black")
lines(lower_prism_region, se_seq_prism_region, lty = 2, col = "black")
abline(v = mu_hat_prism_region, col = "red", lty = 2)  # Pooled estimate line

# Add study labels (customize this column as needed)
text(prism_es$yi, prism_es$sei,
     labels = paste(prism_es$author, prism_es$year, sep = " "),
     cex = 0.7, pos = 4, offset = 0.3)

rm(mu_hat_prism_region, se_seq_prism_region, upper_prism_region, lower_prism_region)
}

# Egger's regression test ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
egger_model_prism_region <- rma.mv(
  yi = yi,
  V = V_matrix_prism,
  mods = ~ precision,                  # Precision as moderator
  random = ~ 1 | geo_location/id,
  data = prism_es,
  method = "REML"
)

summary(egger_model_prism_region)

#* 5.2 rsp region MLMA region ======================================================================
#
# Funnel plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
# Get the pooled effect estimate from your MLMA model
mu_hat_rsp_region <- as.numeric(coef(rsp_region_gold_sens_mlma_model))

# Define SE range for drawing funnel lines
se_seq_rsp_region <- seq(min(rsp_es$sei), max(rsp_es$sei), length.out = 100)
upper_rsp_region <- mu_hat_rsp_region + 1.96 * se_seq_rsp_region
lower_rsp_region <- mu_hat_rsp_region - 1.96 * se_seq_rsp_region

# 4. Create the funnel plot
plot(rsp_es$yi, rsp_es$sei,
     xlab = "Prevalence (logit-transformed)",
     ylab = "Standard Error",
     main = "Funnel plot for publication bias",
     pch = 19,
     ylim = rev(range(rsp_es$sei)))  # Reverse y-axis

# Add funnel
lines(upper_rsp_region, se_seq_rsp_region, lty = 2, col = "black")
lines(lower_rsp_region, se_seq_rsp_region, lty = 2, col = "black")
abline(v = mu_hat_rsp_region, col = "red", lty = 2)  # Pooled estimate line

# Add study labels (customize this column as needed)
text(rsp_es$yi, rsp_es$sei,
     labels = paste(rsp_es$author, rsp_es$year, sep = " "),
     cex = 0.7, pos = 4, offset = 0.3)

rm(mu_hat_rsp_region, se_seq_rsp_region, upper_rsp_region, lower_rsp_region)
}

# Egger's regression test ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
egger_model_rsp_region <- rma.mv(
  yi = yi,
  V = V_matrix_rsp,
  mods = ~ precision,                  # Precision as moderator
  random = ~ 1 | geo_location/id,
  data = rsp_es,
  method = "REML"
)

summary(egger_model_rsp_region)

#* 5.1 comb region MLMA region =====================================================================
#
# Funnel plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
# Get the pooled effect estimate from your MLMA model
mu_hat_comb_region <- as.numeric(coef(comb_region_gold_sens_mlma_model))

# Define SE range for drawing funnel lines
se_seq_comb_region <- seq(min(comb_es$sei), max(comb_es$sei), length.out = 100)
upper_comb_region <- mu_hat_comb_region + 1.96 * se_seq_comb_region
lower_comb_region <- mu_hat_comb_region - 1.96 * se_seq_comb_region

# 4. Create the funnel plot
plot(comb_es$yi, comb_es$sei,
     xlab = "Prevalence (logit-transformed)",
     ylab = "Standard Error",
     main = "Funnel plot for publication bias",
     pch = 19,
     ylim = rev(range(comb_es$sei)))  # Reverse y-axis

# Add funnel
lines(upper_comb_region, se_seq_comb_region, lty = 2, col = "black")
lines(lower_comb_region, se_seq_comb_region, lty = 2, col = "black")
abline(v = mu_hat_comb_region, col = "red", lty = 2)  # Pooled estimate line

# Add study labels (customize this column as needed)
text(comb_es$yi, comb_es$sei,
     labels = paste(comb_es$author, comb_es$year, sep = " "),
     cex = 0.7, pos = 4, offset = 0.3)

rm(mu_hat_comb_region, se_seq_comb_region, upper_comb_region, lower_comb_region)
}

# Egger's regression test ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
egger_model_comb_region <- rma.mv(
  yi = yi,
  V = V_matrix_comb,
  mods = ~ precision,                  # Precision as moderator
  random = ~ 1 | geo_location/id,
  data = comb_es,
  method = "REML"
)

summary(egger_model_comb_region)