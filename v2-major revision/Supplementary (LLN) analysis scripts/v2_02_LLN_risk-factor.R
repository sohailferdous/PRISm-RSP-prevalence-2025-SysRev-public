# Global prevalence of preserved ratio impaired spirometry and restrictive spirometry pattern in the
# general population: a systematic review and meta-analysis
#
#
# This script performs meta-analyses for risk factors of combined LLN-PRISm and LLN-RSP only. This 
# has been done due to insufficient studies when looking at individual categories. Effect estimate 
# is in the form of odds ratios. Meta-analysis of only unadjusted (univariable) ORs have been 
# performed due to insufficient studies adjusting for similar confounders.
# #
# This script is divided into the following broad sections:
# 1. Data import and preparation
# 2. Model building
#
# Please ensure working directory contains v2_data_LLN_risk-factor.csv file.
#
# Install and import packages
# install.packages("tidyverse")
# install.packages("meta")
#
#
library(tidyverse)
library(meta)
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. DATA IMPORT AND PREPARATION ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data <- read.csv ("v2_data_LLN_risk-factor.csv")

data <- data %>% 
  mutate(Study.ID = factor(Study.ID),
         diag = factor(diag),
         diag_defn = factor(diag_defn),
         risk_fac = factor(risk_fac))

#* 1.1 Creating subsets ----------------------------------------------------------------------------
#
#** 1.1.1 Combined LLN-PRISm+RSP subsets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#*** 1.1.1.1 Univariable odds ratios (OR) subset
comb_uni <- data %>%
  drop_na(uni_odds_ratio)

#*** Calculate columns for log odds, log CIs, and standard errors 
comb_uni$logOR <- log(comb_uni$uni_odds_ratio)
comb_uni$logLowerCI <- log(comb_uni$uni_lower_ci)
comb_uni$logUpperCI <- log(comb_uni$uni_upper_ci)
comb_uni$SE <- (comb_uni$logUpperCI - comb_uni$logLowerCI) / (2 * 1.96)

#*** Checking number of studies per risk factor
# table(comb_uni$risk_fac)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. MODEL BUILDING ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# I have not made easily copyable tables for these models - results have been manually entered
# from forest plots.
#
#* 2.1 Asthma (vs non-asthma) ----------------------------------------------------------------------
asthma_comb_uni <- filter(comb_uni, risk_fac=="asthma")

meta_asthma_comb_uni <- metagen(TE = logOR,
                                 seTE = SE,
                                 studlab = author,
                                 data = asthma_comb_uni,
                                 sm = "OR",
                                 method.tau = "REML",
                                 random = TRUE,
                                 common = FALSE,
                                 prediction = TRUE,
                                 warn = TRUE)

forest(meta_asthma_comb_uni)

#* 2.2 Current smoking (vs never/non-smoking) ------------------------------------------------------
current_comb_uni <- filter(comb_uni, risk_fac=="current")

meta_current_comb_uni <- metagen(TE = logOR,
                                  seTE = SE,
                                  studlab = author,
                                  data = current_comb_uni,
                                  sm = "OR",
                                  method.tau = "REML",
                                  random = TRUE,
                                  common = FALSE,
                                  prediction = TRUE,
                                  warn = TRUE)

forest(meta_current_comb_uni)

#* 2.3 Ex smoking (vs never/non-smoking) -----------------------------------------------------------
ex_comb_uni <- filter(comb_uni, risk_fac=="ex")

meta_ex_comb_uni <- metagen(TE = logOR,
                             seTE = SE,
                             studlab = author,
                             data = ex_comb_uni,
                             sm = "OR",
                             method.tau = "REML",
                             random = TRUE,
                             common = FALSE,
                             prediction = TRUE,
                             warn = TRUE)

forest(meta_ex_comb_uni)

#* 2.4 Diabetes (vs non-diabetes) ------------------------------------------------------------------
diabetes_comb_uni <- filter(comb_uni, risk_fac=="diabetes")

meta_diabetes_comb_uni <- metagen(TE = logOR,
                                   seTE = SE,
                                   studlab = author,
                                   data = diabetes_comb_uni,
                                   sm = "OR",
                                   method.tau = "REML",
                                   random = TRUE,
                                   common = FALSE,
                                   prediction = TRUE,
                                   warn = TRUE)

forest(meta_diabetes_comb_uni)

#* 2.5 Male gender (vs females) --------------------------------------------------------------------
male_comb_uni <- filter(comb_uni, risk_fac=="male")

meta_male_comb_uni <- metagen(TE = logOR,
                               seTE = SE,
                               studlab = author,
                               data = male_comb_uni,
                               sm = "OR",
                               method.tau = "REML",
                               random = TRUE,
                               common = FALSE,
                               prediction = TRUE,
                               warn = TRUE)

forest(meta_male_comb_uni)

#* 2.6 Hypertension (vs non-hypertension) ----------------------------------------------------------
hypertension_comb_uni <- filter(comb_uni, risk_fac=="hypertension")

meta_hypertension_comb_uni <- metagen(TE = logOR,
                                       seTE = SE,
                                       studlab = author,
                                       data = hypertension_comb_uni,
                                       sm = "OR",
                                       method.tau = "REML",
                                       random = TRUE,
                                       common = FALSE,
                                       prediction = TRUE,
                                       warn = TRUE)

forest(meta_hypertension_comb_uni)

#* 2.7 Obesity (>= 30 kg/m2 vs normal BMI - 18.5 to 24.9 kg/m2) ------------------------------------
obese_comb_uni <- filter(comb_uni, risk_fac=="obese")

meta_obese_comb_uni <- metagen(TE = logOR,
                                seTE = SE,
                                studlab = author,
                                data = obese_comb_uni,
                                sm = "OR",
                                method.tau = "REML",
                                random = TRUE,
                                common = FALSE,
                                prediction = TRUE,
                                warn = TRUE)

forest(meta_obese_comb_uni)

#* 2.8 Overweight (25-29.9 kg/m2 vs normal BMI - 18.5 to 24.9 kg/m2) -------------------------------
overweight_comb_uni <- filter(comb_uni, risk_fac=="overweight")

meta_overweight_comb_uni <- metagen(TE = logOR,
                                     seTE = SE,
                                     studlab = author,
                                     data = overweight_comb_uni,
                                     sm = "OR",
                                     method.tau = "REML",
                                     random = TRUE,
                                     common = FALSE,
                                     prediction = TRUE,
                                     warn = TRUE)

forest(meta_overweight_comb_uni)
