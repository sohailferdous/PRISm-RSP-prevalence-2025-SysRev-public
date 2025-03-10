# Global prevalence of preserved ratio impaired spirometry and restrictive spirometry pattern in the general 
# population: a systematic review and meta-analysis
#
# Supplementary analysis
#
# This script performs meta-regression analysis to explore heterogeneity during calculation of overall 
# prevalence of LLN-PRISm and LLN-RSP
#
# Variables were included in the multivariable meta-regression only if they had a significance level of p<0.20
# in their univariable meta-regression
#
# Set working directory to directory with data_prevalence.csv file
#
# Install and import packages
install.packages("tidyverse")
install.packages("meta")
install.packages("metafor")
#
#
library(tidyverse)
library(meta)
library(metafor)
#
# Data import and preparation ################################################################################
overall <- read.csv("data_prevalence.csv")

overall <- overall %>%
  drop_na(id) %>% 
  mutate(pub_year = publication.year,
         geo_location = relevel(factor(geo_location), ref= "WPR"),
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

summary(overall)

# Create LLN dataset
lln <- filter(overall, diag_defn == "lln")

# CREATING SUB-DATASETS ######################################################################################
#
# Only heterogeneity for combined prevalence of LLN-PRISm and LLN-RSP was explored via meta-regression due to
# insufficient studies
#
# Combined PRISm and RSP subset
comb_lln <-filter(lln, comb.prism.rsp.prev == "y")

# OVERALL Combined prevalence ################################################################################
# Logit transforming prevalence and calculating variance
comb_lln$logit_prev <- log(comb_lln$overall_reg / (1 - comb_lln$overall_reg))
comb_lln$var_logit <- 1 / (comb_lln$sample_size * comb_lln$overall_reg * (1 - comb_lln$overall_reg))

# Univariate regressions ####
# 1. publication year
reg_comb_lln_pub_year <- rma(yi = logit_prev, 
                               vi = var_logit, 
                               mods = ~ pub_year, 
                               data = comb_lln)
summary(reg_comb_lln_pub_year)

# 2. geo_location
# Having to make a new subset because of study 132
geo_location_comb_lln <-filter(comb_lln, id != 132) # Study ID 132 has been removed as it is conducted in multiple WHO locations

# Logit transforming prevalence and calculating variance
geo_location_comb_lln$logit_prev <- log(geo_location_comb_lln$overall_reg / (1 - geo_location_comb_lln$overall_reg))
geo_location_comb_lln$var_logit <- 1 / (geo_location_comb_lln$sample_size * geo_location_comb_lln$overall_reg * (1 - geo_location_comb_lln$overall_reg))

reg_comb_lln_geo_location <- rma(yi = logit_prev, 
                                 vi = var_logit, 
                                 mods = ~ geo_location, 
                                 data = geo_location_comb_lln)
summary(reg_comb_lln_geo_location)

# 3. bronchod
reg_comb_lln_bronchod <- rma(yi = logit_prev, 
                               vi = var_logit, 
                               mods = ~ bronchod, 
                               data = comb_lln)
summary(reg_comb_lln_bronchod)

# 4. mean_age
reg_comb_lln_mean_age <- rma(yi = logit_prev, 
                               vi = var_logit, 
                               mods = ~ mean_age, 
                               data = comb_lln)
summary(reg_comb_lln_mean_age)

# 5. Proportion of females in study
reg_comb_lln_females_percent <- rma(yi = logit_prev, 
                                      vi = var_logit, 
                                      mods = ~ females_percent, 
                                      data = comb_lln)
summary(reg_comb_lln_females_percent)

# 6. Proportion of current smokers in study
reg_comb_lln_current_smoker_percent <- rma(yi = logit_prev, 
                                             vi = var_logit, 
                                             mods = ~ current_smoker_percent, 
                                             data = comb_lln)
summary(reg_comb_lln_current_smoker_percent)