# Global prevalence of preserved ratio impaired spirometry and restrictive spirometry pattern in the 
# general population: a systematic review and meta-analysis 
#
# This script performs meta-analyses for risk factors of GOLD-PRISm and GOLD-RSP. Effect estimate is 
# in the form of odds ratios. Meta-analysis of unadjusted (univariable) and adjusted (multivarirable) 
# ORs have been performed.For adjusted ORs - they were only performed if they were adjusted for 
# similar confounders.
#
# This script is divided into the following broad sections:
# 1. Data import and preparation
# 2. Model building
# 3. Output display
# 
# Please ensure working directory contains v2_data_GOLD_risk-factor.csv file.
#
# Install and import packages
# install.packages("tidyverse")
# install.packages("meta")
#
#
{                  # Compile here to complete data preparation and model building steps in one go
library(tidyverse)
library(meta)
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. DATA IMPORT AND PREPARATION ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
{                  # Compile here to complete data import and preparation (step 1) in one go
#
data <- read.csv ("v2_data_GOLD_risk-factor.csv")

data <- data %>% 
  mutate(Study.ID = factor(Study.ID),
         diag = factor(diag),
         diag_defn = factor(diag_defn),
         risk_fac = factor(risk_fac))

#* 1.1 Creating subsets ----------------------------------------------------------------------------
#
#** 1.1.1 PRISm subsets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#*** 1.1.1.1 Univariable odds ratios (OR) subset
prism_uni <- data %>% 
  filter(diag == "prism") %>% 
  drop_na(uni_odds_ratio)

#*** Calculating columns for log odds, log CIs, and standard errors 
prism_uni$logOR <- log(prism_uni$uni_odds_ratio)
prism_uni$logLowerCI <- log(prism_uni$uni_lower_ci)
prism_uni$logUpperCI <- log(prism_uni$uni_upper_ci)
prism_uni$SE <- (prism_uni$logUpperCI - prism_uni$logLowerCI) / (2 * 1.96)

#*** Checking number of studies per risk factor
table(prism_uni$risk_fac)

#*** 1.1.1.2 Multivariable OR subset
prism_multi <- data %>% 
  filter(diag == "prism") %>% 
  drop_na(multi_odds_ratio)

#*** Calculating columns for log odds, log CIs, and standard errors 
prism_multi$logOR <- log(prism_multi$multi_odds_ratio)
prism_multi$logLowerCI <- log(prism_multi$multi_lower_ci)
prism_multi$logUpperCI <- log(prism_multi$multi_upper_ci)
prism_multi$SE <- (prism_multi$logUpperCI - prism_multi$logLowerCI) / (2 * 1.96)

#*** Checking number of studies per risk factor
# table(prism_multi$risk_fac)

#** 1.1.2 RSP subsets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#*** 1.1.2.1 Univariable odds ratios (OR) subset
rsp_uni <- data %>% 
  filter(diag == "rsp") %>% 
  drop_na(uni_odds_ratio)

#*** Calculate columns for log odds, log CIs, and standard errors 
rsp_uni$logOR <- log(rsp_uni$uni_odds_ratio)
rsp_uni$logLowerCI <- log(rsp_uni$uni_lower_ci)
rsp_uni$logUpperCI <- log(rsp_uni$uni_upper_ci)
rsp_uni$SE <- (rsp_uni$logUpperCI - rsp_uni$logLowerCI) / (2 * 1.96)

#*** Checking number of studies per risk factor
table(rsp_uni$risk_fac)

#*** 1.1.2.2 Multivariable OR subset
rsp_multi <- data %>% 
  filter(diag == "rsp") %>% 
  drop_na(multi_odds_ratio)

#*** Calculate columns for log odds, log CIs, and standard errors 
rsp_multi$logOR <- log(rsp_multi$multi_odds_ratio)
rsp_multi$logLowerCI <- log(rsp_multi$multi_lower_ci)
rsp_multi$logUpperCI <- log(rsp_multi$multi_upper_ci)
rsp_multi$SE <- (rsp_multi$logUpperCI - rsp_multi$logLowerCI) / (2 * 1.96)

#*** Checking number of studies per risk factor
# table(rsp_multi$risk_fac)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. MODEL BUILDING ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#* 2.1 PRISm univariable risk factors --------------------------------------------------------------
#
#** 2.1.1 Age - 50-59 (vs 40-49 years)
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

# summary(meta_age50.59_prism_uni)
# 
# forest(meta_age50.59_prism_uni)

#** 2.1.2 Age 60-69 (vs 40-49 years)
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

# summary(meta_age60.69_prism_uni)
# 
# forest(meta_age60.69_prism_uni)

#** 2.1.3 Age 70-79 (vs 40-49 years)
age70.79_prism_uni <- filter(prism_uni, risk_fac=="70-79")

meta_age70.79_prism_uni <- metagen(TE = logOR,
                                   seTE = SE,
                                   studlab = author,
                                   data = age70.79_prism_uni,
                                   sm = "OR",
                                   method.tau = "REML",
                                   random = TRUE,
                                   common = FALSE,
                                   prediction = TRUE,
                                   warn = TRUE)

# summary(meta_age70.79_prism_uni)
# 
# forest(meta_age70.79_prism_uni)

#** 2.1.4 Asthma (vs non-asthma)
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

# summary(meta_asthma_prism_uni)
# 
# forest(meta_asthma_prism_uni)

#** 2.1.5 Current smoking (vs never/non-smoking)
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

# summary(meta_current_prism_uni)
# 
# forest(meta_current_prism_uni)

#** 2.1.6 Ex smoking (vs never/non-smoking)
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

# summary(meta_ex_prism_uni)
# 
# forest(meta_ex_prism_uni)

#** 2.1.7 Diabetes (vs non-diabetes)
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

# summary(meta_diabetes_prism_uni)
# 
# forest(meta_diabetes_prism_uni)

#** 2.1.8 Male gender (vs females)
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

# summary(meta_male_prism_uni)
# 
# forest(meta_male_prism_uni)

#** 2.1.9 Hypertension (vs non-hypertension)
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

# summary(meta_hypertension_prism_uni)
# 
# forest(meta_hypertension_prism_uni)

#** 2.1.10 Obesity (>= 30 kg/m2 vs normal BMI - 18.5 to 24.9 kg/m2)
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

# summary(meta_obese_prism_uni)
# 
# forest(meta_obese_prism_uni)

#** 2.1.11 Overweight (25-29.9 kg/m2)
overweight_prism_uni <- filter(prism_uni, risk_fac=="overweight")

meta_overweight_prism_uni <- metagen(TE = logOR,
                                seTE = SE,
                                studlab = author,
                                data = overweight_prism_uni,
                                sm = "OR",
                                method.tau = "REML",
                                random = TRUE,
                                common = FALSE,
                                prediction = TRUE,
                                warn = TRUE)

# summary(meta_overweight_prism_uni)
# 
# forest(meta_overweight_prism_uni)

#** 2.1.12 Underweight (<18.5 kg/m2)
underweight_prism_uni <- filter(prism_uni, risk_fac=="underweight")

meta_underweight_prism_uni <- metagen(TE = logOR,
                                seTE = SE,
                                studlab = author,
                                data = underweight_prism_uni,
                                sm = "OR",
                                method.tau = "REML",
                                random = TRUE,
                                common = FALSE,
                                prediction = TRUE,
                                warn = TRUE)

# summary(meta_underweight_prism_uni)
# 
# forest(meta_underweight_prism_uni)

#** 2.1.13 History of stroke (no history)
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

# summary(meta_stroke_prism_uni)
# 
# forest(meta_stroke_prism_uni)

#** 2.1.14 History of cvd (no history)
cvd_prism_uni <- filter(prism_uni, risk_fac=="cvd")

meta_cvd_prism_uni <- metagen(TE = logOR,
                                 seTE = SE,
                                 studlab = author,
                                 data = cvd_prism_uni,
                                 sm = "OR",
                                 method.tau = "REML",
                                 random = TRUE,
                                 common = FALSE,
                                 prediction = TRUE,
                                 warn = TRUE)

# summary(meta_cvd_prism_uni)
# 
# forest(meta_cvd_prism_uni)

#** 2.1.15 History of tuberculosis (no history)
tuberculosis_prism_uni <- filter(prism_uni, risk_fac=="tuberculosis")

meta_tuberculosis_prism_uni <- metagen(TE = logOR,
                                 seTE = SE,
                                 studlab = author,
                                 data = tuberculosis_prism_uni,
                                 sm = "OR",
                                 method.tau = "REML",
                                 random = TRUE,
                                 common = FALSE,
                                 prediction = TRUE,
                                 warn = TRUE)

# summary(meta_tuberculosis_prism_uni)
# 
# forest(meta_tuberculosis_prism_uni)

#* 2.2 PRISm multivariable risk factors ------------------------------------------------------------
#
#** 2.2.1 Current smoking (vs never/non-smoking)
current_prism_multi <- filter(prism_multi, risk_fac=="current")

meta_current_prism_multi <- metagen(TE = logOR,
                                  seTE = SE,
                                  studlab = author,
                                  data = current_prism_multi,
                                  sm = "OR",
                                  method.tau = "REML",
                                  random = TRUE,
                                  common = FALSE,
                                  prediction = TRUE,
                                  warn = TRUE)

# summary(meta_current_prism_multi)
# 
# forest(meta_current_prism_multi)

#** 2.2.2 Ex smoking (vs never/non-smoking)
ex_prism_multi <- filter(prism_multi, risk_fac=="ex")

meta_ex_prism_multi <- metagen(TE = logOR,
                             seTE = SE,
                             studlab = author,
                             data = ex_prism_multi,
                             sm = "OR",
                             method.tau = "REML",
                             random = TRUE,
                             common = FALSE,
                             prediction = TRUE,
                             warn = TRUE)

# summary(meta_ex_prism_multi)
# 
# forest(meta_ex_prism_multi)


#** 2.2.3 Diabetes (vs non-diabetes)
diabetes_prism_multi <- filter(prism_multi, risk_fac=="diabetes")

meta_diabetes_prism_multi <- metagen(TE = logOR,
                                   seTE = SE,
                                   studlab = author,
                                   data = diabetes_prism_multi,
                                   sm = "OR",
                                   method.tau = "REML",
                                   random = TRUE,
                                   common = FALSE,
                                   prediction = TRUE,
                                   warn = TRUE)

# summary(meta_diabetes_prism_multi)
# 
# forest(meta_diabetes_prism_multi)

#** 2.2.4 Hypertension (vs non-hypertension)
hypertension_prism_multi <- filter(prism_multi, risk_fac=="hypertension")

meta_hypertension_prism_multi <- metagen(TE = logOR,
                                       seTE = SE,
                                       studlab = author,
                                       data = hypertension_prism_multi,
                                       sm = "OR",
                                       method.tau = "REML",
                                       random = TRUE,
                                       common = FALSE,
                                       prediction = TRUE,
                                       warn = TRUE)

# summary(meta_hypertension_prism_multi)
# 
# forest(meta_hypertension_prism_multi)

#* 2.3 RSP univariable risk factors ----------------------------------------------------------------
#
#** 2.3.1 Current smoking (vs never/non-smoking)
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

# summary(meta_current_rsp_uni)
# 
# forest(meta_current_rsp_uni)

#** 2.3.2 Ex smoking (vs never/non-smoking)
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

# summary(meta_ex_rsp_uni)
# 
# forest(meta_ex_rsp_uni)

#** 2.3.3 Diabetes (vs non-diabetes)
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

# summary(meta_diabetes_rsp_uni)
# 
# forest(meta_diabetes_rsp_uni)

#** 2.3.4 Male gender (vs females)
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

# summary(meta_male_rsp_uni)
# 
# forest(meta_male_rsp_uni)

#** 2.3.5 Underweight (<18.5 kg/m2 vs normal BMI - 18.5 to 24.9 kg/m2)
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

# summary(meta_underweight_rsp_uni)
# 
# forest(meta_underweight_rsp_uni)

#** 2.3.6 Overweight (25-29.9 kg/m2 vs normal BMI - 18.5 to 24.9 kg/m2)
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

# summary(meta_overweight_rsp_uni)
# 
# forest(meta_overweight_rsp_uni)

#** 2.3.7 Obesity (>= 30 kg/m2 vs normal BMI - 18.5 to 24.9 kg/m2)
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

# summary(meta_obese_rsp_uni)
# 
# forest(meta_obese_rsp_uni)

#** 2.3.8 Hypertension (vs non-hypertension)
hypertension_rsp_uni <- filter(rsp_uni, risk_fac=="hypertension")

meta_hypertension_rsp_uni <- metagen(TE = logOR,
                                       seTE = SE,
                                       studlab = author,
                                       data = hypertension_rsp_uni,
                                       sm = "OR",
                                       method.tau = "REML",
                                       random = TRUE,
                                       common = FALSE,
                                       prediction = TRUE,
                                       warn = TRUE)

# summary(meta_hypertension_rsp_uni)
# 
# forest(meta_hypertension_rsp_uni)

#* 2.4 RSP univariable risk factors ----------------------------------------------------------------
#
#** 2.4.1 Current smoking (vs never/non-smoking)
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

# summary(meta_current_rsp_multi)
# 
# forest(meta_current_rsp_multi)

#** 2.4.2 Male gender (vs females)
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

# summary(meta_male_rsp_multi)
# 
# forest(meta_male_rsp_multi)

#** 2.4.3 Obesity (>= 30 kg/m2 vs normal BMI - 18.5 to 24.9 kg/m2)
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

# summary(meta_obese_rsp_multi)
# 
# forest(meta_obese_rsp_multi)

#** 2.4.4 Diabetes (vs non-diabetes)
diabetes_rsp_multi <- filter(rsp_multi, risk_fac=="diabetes")

meta_diabetes_rsp_multi <- metagen(TE = logOR,
                                     seTE = SE,
                                     studlab = author,
                                     data = diabetes_rsp_multi,
                                     sm = "OR",
                                     method.tau = "REML",
                                     random = TRUE,
                                     common = FALSE,
                                     prediction = TRUE,
                                     warn = TRUE)

# summary(meta_diabetes_rsp_multi)
# 
# forest(meta_diabetes_rsp_multi)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. OUTPUT DISPLAY ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
# Create a function to extract model results
  extract_risk_results <- function(model, model_name) {
    k <- model$k
    
    # Extract odds ratio and CI
    est <- round(exp(model$TE.random), 2)
    ci_lb <- round(exp(model$lower.random), 2)
    ci_ub <- round(exp(model$upper.random), 2)
    
    # Prediction interval
    pi_lb <- round(exp(model$lower.predict), 2)
    pi_ub <- round(exp(model$upper.predict), 2)
    
    # I2
    i2 <- round(model$I2*100, 1)
    
    # Format output
    estimate_ci <- paste0(est, " (", ci_lb, "-", ci_ub, ")")
    prediction_interval <- paste0(pi_lb, "-", pi_ub)
    
    results <- data.frame(
      Model = model_name,
      k = k,
      I2 = i2,
      OR_CI = estimate_ci,
      Prediction_Interval = prediction_interval,
      row.names = NULL
    )
    
    return(results)
  }

res01 <- extract_risk_results(meta_age50.59_prism_uni, "50-59 prism uni")
res02 <- extract_risk_results(meta_age60.69_prism_uni, "60-69 prism uni")
res03 <- extract_risk_results(meta_age70.79_prism_uni, "70-79 prism uni")
res04 <- extract_risk_results(meta_asthma_prism_uni, "asthma prism uni")
res05 <- extract_risk_results(meta_current_prism_uni, "current smoker prism uni")
res06 <- extract_risk_results(meta_ex_prism_uni, "ex smoker prism uni")
res07 <- extract_risk_results(meta_diabetes_prism_uni, "diabetes prism uni")
res08 <- extract_risk_results(meta_male_prism_uni, "male prism uni")
res09 <- extract_risk_results(meta_hypertension_prism_uni, "hypertension prism uni")
res10 <- extract_risk_results(meta_obese_prism_uni, "obese prism uni")
res11 <- extract_risk_results(meta_overweight_prism_uni, "overweight prism uni")
res12 <- extract_risk_results(meta_underweight_prism_uni, "underweight prism uni")
res13 <- extract_risk_results(meta_stroke_prism_uni, "stroke prism uni")
res14 <- extract_risk_results(meta_cvd_prism_uni, "cvd prism uni")
res15 <- extract_risk_results(meta_tuberculosis_prism_uni, "tuberculosis prism uni")
res16 <- extract_risk_results(meta_current_prism_multi, "current smoker prism multi")
res17 <- extract_risk_results(meta_ex_prism_multi, "ex smoker prism multi")
res18 <- extract_risk_results(meta_diabetes_prism_multi, "diabetes prism multi")
res19 <- extract_risk_results(meta_hypertension_prism_multi, "hypertension prism multi")
res20 <- extract_risk_results(meta_current_rsp_uni, "current smoker rsp uni")
res21 <- extract_risk_results(meta_ex_rsp_uni, "ex smoker rsp uni")
res22 <- extract_risk_results(meta_diabetes_rsp_uni, "diabetes rsp uni")
res23 <- extract_risk_results(meta_male_rsp_uni, "males rsp uni")
res24 <- extract_risk_results(meta_obese_rsp_uni, "obese rsp uni")
res25 <- extract_risk_results(meta_overweight_rsp_uni, "overweight rsp uni")
res26 <- extract_risk_results(meta_underweight_rsp_uni, "underweight rsp uni")
res27 <- extract_risk_results(meta_hypertension_rsp_uni, "hypertension rsp uni")
res28 <- extract_risk_results(meta_current_rsp_multi, "current smoking rsp multi")
res29 <- extract_risk_results(meta_male_rsp_multi, "males rsp multi")
res30 <- extract_risk_results(meta_obese_rsp_multi, "obese rsp multi")
res31 <- extract_risk_results(meta_diabetes_rsp_multi, "diabetes rsp multi")

# Combine results
or_summary_table <- bind_rows(res01, res02, res03, res04, res05, res06, res07, res08, res09, res10, res11, res12, res13, res14, res15, 
                              res16, res17, res18, res19, res20, res21, res22, res23, res24, res25, res26, res27,
                              res28, res29, res30, res31)
print(or_summary_table)

rm(res01, res02, res03, res04, res05, res06, res07, res08, res09, res10, res11, res12, res13, res14, res15, 
   res16, res17, res18, res19, res20, res21, res22, res23, res24, res25, res26, res27,
   res28, res29, res30, res31)
}
