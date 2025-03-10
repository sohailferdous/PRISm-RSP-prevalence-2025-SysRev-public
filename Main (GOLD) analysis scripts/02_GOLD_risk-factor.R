# Global prevalence of preserved ratio impaired spirometry and restrictive spirometry pattern in the general 
# population: a systematic review and meta-analysis 
#
# This script performs meta-analyses for risk factors of GOLD-PRISm and GOLD-RSP. Effect estimate is in the
# form of odds ratios. Meta-analysis of unadjusted (univariable) and adjusted (multivariatble) ORs have been
# performed.For adjusted ORs - they were only recorded if they were adjusted for similar confounders
#
# Set working directory to directory with data_GOLD_risk-factor.csv file
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
data <- read.csv ("data_GOLD_risk-factor.csv")

data <- data %>% 
  mutate(Study.ID = factor(Study.ID),
         diag = factor(diag),
         diag_defn = factor(diag_defn),
         risk_fac = factor(risk_fac))

# Creating subsets ###########################################################################################
#
# 1. PRISm subsets ####
#
# Univariable odds ratios (OR) subset
prism_uni <- data %>% 
  filter(diag == "prism") %>% 
  drop_na(uni_odds_ratio)

# Calculating columns for log odds, log CIs, and standard errors 
prism_uni$logOR <- log(prism_uni$uni_odds_ratio)
prism_uni$logLowerCI <- log(prism_uni$uni_lower_ci)
prism_uni$logUpperCI <- log(prism_uni$uni_upper_ci)
prism_uni$SE <- (prism_uni$logUpperCI - prism_uni$logLowerCI) / (2 * 1.96)

# Checking number of studies per risk factor
table(prism_uni$risk_fac)

# Multivariable OR subset
prism_multi <- data %>% 
  filter(diag == "prism") %>% 
  drop_na(multi_odds_ratio)

# Calculating columns for log odds, log CIs, and standard errors 
prism_multi$logOR <- log(prism_multi$multi_odds_ratio)
prism_multi$logLowerCI <- log(prism_multi$multi_lower_ci)
prism_multi$logUpperCI <- log(prism_multi$multi_upper_ci)
prism_multi$SE <- (prism_multi$logUpperCI - prism_multi$logLowerCI) / (2 * 1.96)

# Checking number of studies per risk factor
table(prism_multi$risk_fac)

# 2. RSP subsets ####
#
# Univariable OR subset
rsp_uni <- data %>% 
  filter(diag == "rsp") %>% 
  drop_na(uni_odds_ratio)

# Calculate columns for log odds, log CIs, and standard errors 
rsp_uni$logOR <- log(rsp_uni$uni_odds_ratio)
rsp_uni$logLowerCI <- log(rsp_uni$uni_lower_ci)
rsp_uni$logUpperCI <- log(rsp_uni$uni_upper_ci)
rsp_uni$SE <- (rsp_uni$logUpperCI - rsp_uni$logLowerCI) / (2 * 1.96)

# Checking number of studies per risk factor
table(rsp_uni$risk_fac)

# Multivariable OR subset
rsp_multi <- data %>% 
  filter(diag == "rsp") %>% 
  drop_na(multi_odds_ratio)

# Calculate columns for log odds, log CIs, and standard errors 
rsp_multi$logOR <- log(rsp_multi$multi_odds_ratio)
rsp_multi$logLowerCI <- log(rsp_multi$multi_lower_ci)
rsp_multi$logUpperCI <- log(rsp_multi$multi_upper_ci)
rsp_multi$SE <- (rsp_multi$logUpperCI - rsp_multi$logLowerCI) / (2 * 1.96)

# Checking number of studies per risk factor
table(rsp_multi$risk_fac)

# Meta-analysis for odds ratio for risk factors ##############################################################
# A. PRISm univariable risk factors ####
#
# 1. Age - 50-59 (vs 40-49 years)
age50.59_prism_uni <- filter(prism_uni, risk_fac=="50-59")

meta_age50.59_prism_uni <- metagen(TE = logOR,             # Treatment effect (logOR)
                                seTE = SE,                 # Standard error
                                studlab = author,          # Study labels
                                data = age50.59_prism_uni, # Dataset
                                sm = "OR",                 # Summary measure: Odds Ratio
                                method.tau = "REML",       # Restricted maximum likelihood method to calculate variance
                                random = TRUE,             # Random effects model
                                common = FALSE,            # Fixed effects model
                                prediction = TRUE,         # Generate prediction interval
                                warn = TRUE)               # Report warnings

forest(meta_age50.59_prism_uni)

# 2. Age 60-69 (vs 40-49 years)
age60.69_prism_uni <- filter(prism_uni, risk_fac=="60-69")

meta_age60.69_prism_uni <- metagen(TE = logOR,
                                   seTE = SE,
                                   studlab = author,
                                   data = age60.69_prism_uni,
                                   sm = "OR",
                                   method.tau = "REML",
                                   random = TRUE,
                                   common = FALSE,
                                   prediction = TRUE,
                                   warn = TRUE)

forest(meta_age60.69_prism_uni)

# 3. Asthma (vs non-asthma)
asthma_prism_uni <- filter(prism_uni, risk_fac=="asthma")

meta_asthma_prism_uni <- metagen(TE = logOR,
                                 seTE = SE,
                                 studlab = author,
                                 data = asthma_prism_uni,
                                 sm = "OR",
                                 method.tau = "REML",
                                 random = TRUE,
                                 common = FALSE,
                                 prediction = TRUE,
                                 warn = TRUE)

forest(meta_asthma_prism_uni)

# 4. Current smoking (vs never/non-smoking)
current_prism_uni <- filter(prism_uni, risk_fac=="current")

meta_current_prism_uni <- metagen(TE = logOR,
                                 seTE = SE,
                                 studlab = author,
                                 data = current_prism_uni,
                                 sm = "OR",
                                 method.tau = "REML",
                                 random = TRUE,
                                 common = FALSE,
                                 prediction = TRUE,
                                 warn = TRUE)

forest(meta_current_prism_uni)

# 5. Ex smoking (vs never/non-smoking)
ex_prism_uni <- filter(prism_uni, risk_fac=="ex")

meta_ex_prism_uni <- metagen(TE = logOR,
                                  seTE = SE,
                                  studlab = author,
                                  data = ex_prism_uni,
                                  sm = "OR",
                                  method.tau = "REML",
                                  random = TRUE,
                                  common = FALSE,
                                  prediction = TRUE,
                                  warn = TRUE)

forest(meta_ex_prism_uni)

# 6. Diabetes (vs non-diabetes)
diabetes_prism_uni <- filter(prism_uni, risk_fac=="diabetes")

meta_diabetes_prism_uni <- metagen(TE = logOR,
                                   seTE = SE,
                                   studlab = author,
                                   data = diabetes_prism_uni,
                                   sm = "OR",
                                   method.tau = "REML",
                                   random = TRUE,
                                   common = FALSE,
                                   prediction = TRUE,
                                   warn = TRUE)

forest(meta_diabetes_prism_uni)

# 7. Female gender (vs males)              # Not reported
female_prism_uni <- filter(prism_uni, risk_fac=="female")

meta_female_prism_uni <- metagen(TE = logOR,
                                   seTE = SE,
                                   studlab = author,
                                   data = female_prism_uni,
                                   sm = "OR",
                                   method.tau = "REML",
                                   random = TRUE,
                                   common = FALSE,
                                   prediction = TRUE,
                                   warn = TRUE)

forest(meta_female_prism_uni)

# 8. Male gender (vs females)
male_prism_uni <- filter(prism_uni, risk_fac=="male")

meta_male_prism_uni <- metagen(TE = logOR,
                                   seTE = SE,
                                   studlab = author,
                                   data = male_prism_uni,
                                   sm = "OR",
                                   method.tau = "REML",
                                   random = TRUE,
                                   common = FALSE,
                                   prediction = TRUE,
                                   warn = TRUE)

forest(meta_male_prism_uni)

# 9. Hypertension (vs non-hypertension)
hypertension_prism_uni <- filter(prism_uni, risk_fac=="hypertension")

meta_hypertension_prism_uni <- metagen(TE = logOR,
                               seTE = SE,
                               studlab = author,
                               data = hypertension_prism_uni,
                               sm = "OR",
                               method.tau = "REML",
                               random = TRUE,
                               common = FALSE,
                               prediction = TRUE,
                               warn = TRUE)

forest(meta_hypertension_prism_uni)

# 10. Obesity (>= 30 kg/m2 vs normal BMI - 18.5 to 24.9 kg/m2)
obese_prism_uni <- filter(prism_uni, risk_fac=="obese")

meta_obese_prism_uni <- metagen(TE = logOR,
                                       seTE = SE,
                                       studlab = author,
                                       data = obese_prism_uni,
                                       sm = "OR",
                                       method.tau = "REML",
                                       random = TRUE,
                                       common = FALSE,
                                       prediction = TRUE,
                                       warn = TRUE)

forest(meta_obese_prism_uni)

# 11. History of stroke (no history)
stroke_prism_uni <- filter(prism_uni, risk_fac=="stroke")

meta_stroke_prism_uni <- metagen(TE = logOR,
                                seTE = SE,
                                studlab = author,
                                data = stroke_prism_uni,
                                sm = "OR",
                                method.tau = "REML",
                                random = TRUE,
                                common = FALSE,
                                prediction = TRUE,
                                warn = TRUE)

forest(meta_stroke_prism_uni)

# B. PRISm multivariable risk factors ####
#
# Not enough studies for any risk factor (minimum 3)
#
# C. RSP univariable risk factors ####
#
# 1. Current smoking (vs never/non-smoking)
current_rsp_uni <- filter(rsp_uni, risk_fac=="current")

meta_current_rsp_uni <- metagen(TE = logOR,
                                  seTE = SE,
                                  studlab = author,
                                  data = current_rsp_uni,
                                  sm = "OR",
                                  method.tau = "REML",
                                  random = TRUE,
                                  common = FALSE,
                                  prediction = TRUE,
                                  warn = TRUE)

forest(meta_current_rsp_uni)

# 2. Ex smoking (vs never/non-smoking)
ex_rsp_uni <- filter(rsp_uni, risk_fac=="ex")

meta_ex_rsp_uni <- metagen(TE = logOR,
                             seTE = SE,
                             studlab = author,
                             data = ex_rsp_uni,
                             sm = "OR",
                             method.tau = "REML",
                             random = TRUE,
                             common = FALSE,
                             prediction = TRUE,
                             warn = TRUE)

forest(meta_ex_rsp_uni)

# 3. Diabetes (vs non-diabetes)
diabetes_rsp_uni <- filter(rsp_uni, risk_fac=="diabetes")

meta_diabetes_rsp_uni <- metagen(TE = logOR,
                                   seTE = SE,
                                   studlab = author,
                                   data = diabetes_rsp_uni,
                                   sm = "OR",
                                   method.tau = "REML",
                                   random = TRUE,
                                   common = FALSE,
                                   prediction = TRUE,
                                   warn = TRUE)

forest(meta_diabetes_rsp_uni)

# 4. Female gender (vs males)          # Not reported
female_rsp_uni <- filter(rsp_uni, risk_fac=="female")

meta_female_rsp_uni <- metagen(TE = logOR,
                                 seTE = SE,
                                 studlab = author,
                                 data = female_rsp_uni,
                                 sm = "OR",
                                 method.tau = "REML",
                                 random = TRUE,
                                 common = FALSE,
                                 prediction = TRUE,
                                 warn = TRUE)

forest(meta_female_rsp_uni)

# 5. Male gender (vs females)
male_rsp_uni <- filter(rsp_uni, risk_fac=="male")

meta_male_rsp_uni <- metagen(TE = logOR,
                               seTE = SE,
                               studlab = author,
                               data = male_rsp_uni,
                               sm = "OR",
                               method.tau = "REML",
                               random = TRUE,
                               common = FALSE,
                               prediction = TRUE,
                               warn = TRUE)

forest(meta_male_rsp_uni)

# 6. Underweight (<18.5 kg/m2 vs normal BMI - 18.5 to 24.9 kg/m2)
underweight_rsp_uni <- filter(rsp_uni, risk_fac=="underweight")

meta_underweight_rsp_uni <- metagen(TE = logOR,
                                seTE = SE,
                                studlab = author,
                                data = underweight_rsp_uni,
                                sm = "OR",
                                method.tau = "REML",
                                random = TRUE,
                                common = FALSE,
                                prediction = TRUE,
                                warn = TRUE)

forest(meta_underweight_rsp_uni)

# 7. Overweight (25-29.9 kg/m2 vs normal BMI - 18.5 to 24.9 kg/m2)
overweight_rsp_uni <- filter(rsp_uni, risk_fac=="overweight")

meta_overweight_rsp_uni <- metagen(TE = logOR,
                                seTE = SE,
                                studlab = author,
                                data = overweight_rsp_uni,
                                sm = "OR",
                                method.tau = "REML",
                                random = TRUE,
                                common = FALSE,
                                prediction = TRUE,
                                warn = TRUE)

forest(meta_overweight_rsp_uni)

# 8. Obesity (>= 30 kg/m2 vs normal BMI - 18.5 to 24.9 kg/m2)
obese_rsp_uni <- filter(rsp_uni, risk_fac=="obese")

meta_obese_rsp_uni <- metagen(TE = logOR,
                                seTE = SE,
                                studlab = author,
                                data = obese_rsp_uni,
                                sm = "OR",
                                method.tau = "REML",
                                random = TRUE,
                                common = FALSE,
                                prediction = TRUE,
                                warn = TRUE)

forest(meta_obese_rsp_uni)

# D. RSP multivariable risk factors ####
#
# Only risk factors which adjusted for similar factors were included in the analysis. Meta-analysis was
# conducted only if there were a minimum 3 studies
#
# 1. Current smoking (vs never/non-smoking)
current_rsp_multi <- filter(rsp_multi, risk_fac=="current")

meta_current_rsp_multi <- metagen(TE = logOR,
                                seTE = SE,
                                studlab = author,
                                data = current_rsp_multi,
                                sm = "OR",
                                method.tau = "REML",
                                random = TRUE,
                                common = FALSE,
                                prediction = TRUE,
                                warn = TRUE)

forest(meta_current_rsp_multi)

# 2. Male gender (vs females)
male_rsp_multi <- filter(rsp_multi, risk_fac=="male")

meta_male_rsp_multi <- metagen(TE = logOR,
                             seTE = SE,
                             studlab = author,
                             data = male_rsp_multi,
                             sm = "OR",
                             method.tau = "REML",
                             random = TRUE,
                             common = FALSE,
                             prediction = TRUE,
                             warn = TRUE)

forest(meta_male_rsp_multi)



# 3. Obesity (>= 30 kg/m2 vs normal BMI - 18.5 to 24.9 kg/m2)
obese_rsp_multi <- filter(rsp_multi, risk_fac=="obese")

meta_obese_rsp_multi <- metagen(TE = logOR,
                              seTE = SE,
                              studlab = author,
                              data = obese_rsp_multi,
                              sm = "OR",
                              method.tau = "REML",
                              random = TRUE,
                              common = FALSE,
                              prediction = TRUE,
                              warn = TRUE)

forest(meta_obese_rsp_multi)
