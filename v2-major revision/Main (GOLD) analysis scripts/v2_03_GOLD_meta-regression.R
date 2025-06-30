# Global prevalence of preserved ratio impaired spirometry and restrictive spirometry pattern in the 
# general population: a systematic review and meta-analysis 
#
# This script performs meta-regression analysis to explore heterogeneity during calculation of 
# overall prevalence of GOLD-PRISm and GOLD-RSP.
#
# Variables were for multivariate regressions were selected by manual backward elimination to build 
# a final model which best explained our data (quantified by AIC and R2 values). As per reveiwer
# recommendation, multivariate models with all covariates have been reported in supplementary
# materials.
#
# This script is divided into the following broad sections:
# 1. Data import and preparation
# 2. PRISm meta-regressions
# 3. RSP meta-regressions
#
# Please ensure working directory contains v2_data_prevalence.csv file.
#
# Install and import packages
# install.packages("tidyverse")
# install.packages("metafor")
#
{               # Compile here to finish data import and preparation (step 1) in one go
library(tidyverse)
library(metafor)
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. DATA IMPORT AND PREPARATION ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
overall <- read.csv("v2_data_prevalence.csv")

overall <- overall %>%
  drop_na(id) %>% 
  mutate(geo_location = relevel(factor(geo_location), ref= "WPR"),
         id = factor(id),
         country = factor(country),
         income_class = relevel(factor(income_class),ref="high"),
         diag = factor(diag),
         diag_defn = factor(diag_defn),
         bronchod = factor(bronchod),
         overall.prev = factor(overall.prev),
         gender.prev = factor(gender.prev),
         comb.prism.rsp.prev = factor(comb.prism.rsp.prev),
         smoking.prev = factor(smoking.prev))

summary(overall)

# Create GOLD dataset
gold <- filter(overall, diag_defn == "gold")

# Creating sub-datasets
# PRISm
prism_gold <-filter(gold, diag == "prism" & overall.prev == "y")

# RSP
rsp_gold <-filter(gold, diag == "rsp" & overall.prev == "y")

# Combined PRISm and RSP 
comb_gold <-filter(gold, comb.prism.rsp.prev == "y")

# Logit transforming prevalence and calculating variance
prism_gold$logit_prev <- log(prism_gold$overall_reg / (1 - prism_gold$overall_reg))
prism_gold$var_logit <- 1 / (prism_gold$sample_size * prism_gold$overall_reg * (1 - prism_gold$overall_reg))

rsp_gold$logit_prev <- log(rsp_gold$overall_reg / (1 - rsp_gold$overall_reg))
rsp_gold$var_logit <- 1 / (rsp_gold$sample_size * rsp_gold$overall_reg * (1 - rsp_gold$overall_reg))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. PRISm META-REGRESSIONS ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{                 # Compile here to finish generating univariate regression models in one go
#* 2.1 Univariate regressions ----------------------------------------------------------------------
#
#** 2.1.1 publication year
reg_prism_gold_pub_year <- rma(yi = logit_prev, 
                               vi = var_logit, 
                               mods = ~ pub_year, 
                               data = prism_gold)

#** 2.1.2 geo_location
reg_prism_gold_geo_location <- rma(yi = logit_prev, 
                               vi = var_logit, 
                               mods = ~ geo_location, 
                               data = prism_gold)

#** 2.1.3 bronchod
reg_prism_gold_bronchod <- rma(yi = logit_prev, 
                                   vi = var_logit, 
                                   mods = ~ bronchod, 
                                   data = prism_gold)

#** 2.1.4 mean_age
reg_prism_gold_mean_age <- rma(yi = logit_prev, 
                               vi = var_logit, 
                               mods = ~ mean_age, 
                               data = prism_gold)

#** 2.1.5 Proportion of females in study
reg_prism_gold_females_percent <- rma(yi = logit_prev, 
                               vi = var_logit, 
                               mods = ~ females_percent, 
                               data = prism_gold)

#** 2.1.6 Proportion of current smokers in study
reg_prism_gold_current_smoker_percent <- rma(yi = logit_prev, 
                                      vi = var_logit, 
                                      mods = ~ current_smoker_percent, 
                                      data = prism_gold)

#** 2.1.7 income_class
reg_prism_gold_income_class <- rma(yi = logit_prev, 
                                   vi = var_logit, 
                                   mods = ~ income_class, 
                                   data = prism_gold)

}

#** 2.1.8 Create a table to display results from each model
{
#*** 2.1.8.1 Function to extract and format model results
extract_model_results <- function(model, mod_name) {
  sum_model <- summary(model)
  
  coefs <- as.data.frame(sum_model$beta)
  pvals <- as.data.frame(sum_model$pval)
  
#*** 2.1.8.2 Creating table
  results <- data.frame(
    Moderator = rownames(coefs),
    Estimate = round(coefs[, 1], 2),
    CI_LB = round(sum_model$ci.lb, 2),
    CI_UB = round(sum_model$ci.ub, 2),
    Pval = round(pvals[, 1], 2),
    R2 = round(sum_model$R2, 2),
    Model = mod_name,
    row.names = NULL
  )
  
#*** 2.1.8.3 Removing intrcpt
  results <- results %>% filter(Moderator != "intrcpt")
  
#*** 2.1.8.4 Changing format to make it easily copyable (no this is not AI lol)
  results$Formatted <- paste0(
    results$Estimate, " (", results$CI_LB, " to ", results$CI_UB, ")"
  )
  
  results <- results %>%
    select(Model, Moderator, Formatted, Pval, R2)
  
  return(results)
}

#*** 2.1.8.5 Applying function to each model
{
res1 <- extract_model_results(reg_prism_gold_pub_year, "publication Year")
res2 <- extract_model_results(reg_prism_gold_geo_location, "geo Location")
res3 <- extract_model_results(reg_prism_gold_bronchod, "broncho use")
res4 <- extract_model_results(reg_prism_gold_mean_age, "mean age")
res5 <- extract_model_results(reg_prism_gold_females_percent, "female percent")
res6 <- extract_model_results(reg_prism_gold_current_smoker_percent, "current smoker percent")
res7 <- extract_model_results(reg_prism_gold_income_class, "income-level")

#*** 2.1.8.6 Combining and viewing results, followed by clean-up
univariate_prism_meta_reg <- bind_rows(res1, res2, res3, res4, res5, res6, res7)
print(univariate_prism_meta_reg)

rm(res1, res2, res3, res4, res5, res6, res7)
}
}

#* 2.2 Multivariate regressions --------------------------------------------------------------------
#** 2.2.1 Full model for supplementary
reg_prism_gold_multi_full <- rma(yi = logit_prev, 
                               vi = var_logit, 
                               mods = ~ pub_year+geo_location+bronchod+mean_age+females_percent+
                              current_smoker_percent+income_class, 
                               data = prism_gold)

#** Applying function to model and displaying results
{
  multivariate_prism_full_reg <- extract_model_results(reg_prism_gold_multi_full, "multivariate prism full model")
  print(multivariate_prism_full_reg)
}

#* 2.2.2 Best performing model for main text
reg_prism_gold_multi_best <- rma(yi = logit_prev,
                            vi = var_logit, 
                            mods = ~ mean_age+income_class+pub_year, 
                            data = prism_gold)

#** Applying function to model and displaying results
{
  multivariate_prism_best_reg <- extract_model_results(reg_prism_gold_multi_best, "multivariate prism full model")
  print(multivariate_prism_best_reg)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. RSP META-REGRESSIONS ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{                 # Compile here to finish generating univariate regression models in one go
#* 3.1 Univariate regressions ----------------------------------------------------------------------
#
#** 3.1.1 publication year
reg_rsp_gold_pub_year <- rma(yi = logit_prev, 
                               vi = var_logit, 
                               mods = ~ pub_year, 
                               data = rsp_gold)

#** 3.1.2 geo_location
reg_rsp_gold_geo_location <- rma(yi = logit_prev, 
                                   vi = var_logit, 
                                   mods = ~ geo_location, 
                                   data = rsp_gold)

#** 3.1.3 bronchod
reg_rsp_gold_bronchod <- rma(yi = logit_prev, 
                               vi = var_logit, 
                               mods = ~ bronchod, 
                               data = rsp_gold)

#** 3.1.4 mean_age
reg_rsp_gold_mean_age <- rma(yi = logit_prev, 
                               vi = var_logit, 
                               mods = ~ mean_age, 
                               data = rsp_gold)

#** 3.1.5 Proportion of females in study
reg_rsp_gold_females_percent <- rma(yi = logit_prev, 
                                      vi = var_logit, 
                                      mods = ~ females_percent, 
                                      data = rsp_gold)

#** 3.1.6 Proportion of current smokers in study
reg_rsp_gold_current_smoker_percent <- rma(yi = logit_prev, 
                                             vi = var_logit, 
                                             mods = ~ current_smoker_percent, 
                                             data = rsp_gold)

#** 3.1.7 income_class
reg_rsp_gold_income_class <- rma(yi = logit_prev, 
                                   vi = var_logit, 
                                   mods = ~ income_class, 
                                   data = rsp_gold)
}

#** 3.1.8 Create table with all results using previously created function
{
  res1 <- extract_model_results(reg_rsp_gold_pub_year, "publication Year")
  res2 <- extract_model_results(reg_rsp_gold_geo_location, "geo Location")
  res3 <- extract_model_results(reg_rsp_gold_bronchod, "broncho use")
  res4 <- extract_model_results(reg_rsp_gold_mean_age, "mean age")
  res5 <- extract_model_results(reg_rsp_gold_females_percent, "female percent")
  res6 <- extract_model_results(reg_rsp_gold_current_smoker_percent, "current smoker percent")
  res7 <- extract_model_results(reg_rsp_gold_income_class, "income-level")
  
#** 3.1.8.1 Combine and view results
  univariate_rsp_meta_reg <- bind_rows(res1, res2, res3, res4, res5, res6, res7)
  print(univariate_rsp_meta_reg)
  
  rm(res1, res2, res3, res4, res5, res6, res7)
}

#* 3.2 Multivariate regressions --------------------------------------------------------------------
#** 3.2.1 Full model for supplementary
reg_rsp_gold_multi_full <- rma(yi = logit_prev, 
                            vi = var_logit, 
                            mods = ~ pub_year+geo_location+bronchod+mean_age+females_percent+
                              current_smoker_percent+income_class, 
                            data = rsp_gold)

#** Applying function to model and displaying results
{
  multivariate_rsp_full_reg <- extract_model_results(reg_rsp_gold_multi_full, "multivariate rsp full model")
  print(multivariate_rsp_full_reg)
}

#** 3.2.2 Best performing model for main-text
reg_rsp_gold_multi_best <- rma(yi = logit_prev, 
                          vi = var_logit, 
                          mods = ~ geo_location+bronchod+mean_age+current_smoker_percent,
                          data = rsp_gold)

#** Applying function to model and displaying results
{
  multivariate_rsp_best_reg <- extract_model_results(reg_rsp_gold_multi_best, "multivariate rsp full model")
  print(multivariate_rsp_best_reg)
}