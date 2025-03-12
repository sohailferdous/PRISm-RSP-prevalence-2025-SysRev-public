# Global prevalence of preserved ratio impaired spirometry and restrictive spirometry pattern in the general 
# population: a systematic review and meta-analysis
#
# Supplementary analysis
#
# This script performs meta-analyses for risk factors of combined LLN-PRISm and LLN-RSP only. This has been 
# done due to insufficient studies when looking at individual categories. Effect estimate is in the form of 
# odds ratios. Meta-analysis of unadjusted (univariable) and adjusted (multivariatble) ORs have been 
# performed. For adjusted ORs - they were only recorded if they were adjusted for similar confounders
#
# Set working directory to directory with data_LLN_risk-factor.csv file
#
# Install and import packages
install.packages("tidyverse")
install.packages("meta")
#
#
library(tidyverse)
library(meta)
#
#
# Data import and preparation ################################################################################
data <- read.csv ("data_LLN_risk-factor.csv")

data <- data %>% 
  mutate(Study.ID = factor(Study.ID),
         diag = factor(diag),
         diag_defn = factor(diag_defn),
         risk_fac = factor(risk_fac))

# Creating subsets ###########################################################################################
#
# Only univariate combined meta-analysis has been performed for LLN risk factors because of insufficient 
# studies
#
# 1. Combined PRISM+RSP subsets ####
#
# Univariable odds ratios (OR) subset
comb_uni <- data %>%
  drop_na(uni_odds_ratio)

# Calculate columns for log odds, log CIs, and standard errors 
comb_uni$logOR <- log(comb_uni$uni_odds_ratio)
comb_uni$logLowerCI <- log(comb_uni$uni_lower_ci)
comb_uni$logUpperCI <- log(comb_uni$uni_upper_ci)
comb_uni$SE <- (comb_uni$logUpperCI - comb_uni$logLowerCI) / (2 * 1.96)

table(comb_uni$risk_fac)

# Meta-analysis for odds ratio for risk factors ##############################################################

# Asthma (vs non-asthma) ####
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

# Current smoking (vs never/non-smoking) ####
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

# Ex smoking (vs never/non-smoking) ####
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

# Diabetes (vs non-diabetes) ####
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

# Male gender (vs females) ####
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

# Hypertension (vs non-hypertension) ####
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

# Obesity (>= 30 kg/m2 vs normal BMI - 18.5 to 24.9 kg/m2) ####
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

# Overweight (25-29.9 kg/m2 vs normal BMI - 18.5 to 24.9 kg/m2) ####
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
