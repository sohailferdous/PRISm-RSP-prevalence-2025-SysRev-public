# Global prevalence of preserved ratio impaired spirometry and restrictive spirometry pattern in the general 
# population: a systematic review and meta-analysis 
#
# This script performs meta-regression analysis to explore heterogeneity during calculation of overall 
# prevalence of GOLD-PRISm and GOLD-RSP
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

# Create GOLD dataset
gold <- filter(overall, diag_defn == "gold")

# Creating sub-datasets ######################################################################################
# PRISm
prism_gold <-filter(gold, diag == "prism" & overall.prev == "y")

# RSP
rsp_gold <-filter(gold, diag == "rsp" & overall.prev == "y")

# Combined PRISm and RSP 
comb_gold <-filter(gold, comb.prism.rsp.prev == "y")

# Overall PRISm prevalence ###################################################################################
# Logit transforming prevalence and calculating variance
prism_gold$logit_prev <- log(prism_gold$overall_reg / (1 - prism_gold$overall_reg))
prism_gold$var_logit <- 1 / (prism_gold$sample_size * prism_gold$overall_reg * (1 - prism_gold$overall_reg))

# Univariate regressions ####
# 1. publication year
reg_prism_gold_pub_year <- rma(yi = logit_prev, 
                               vi = var_logit, 
                               mods = ~ pub_year, 
                               data = prism_gold)
summary(reg_prism_gold_pub_year)

# 2. geo_location
reg_prism_gold_geo_location <- rma(yi = logit_prev, 
                               vi = var_logit, 
                               mods = ~ geo_location, 
                               data = prism_gold)
summary(reg_prism_gold_geo_location)

# 3. bronchod
reg_prism_gold_bronchod <- rma(yi = logit_prev, 
                                   vi = var_logit, 
                                   mods = ~ bronchod, 
                                   data = prism_gold)
summary(reg_prism_gold_bronchod)

# 4. mean_age
reg_prism_gold_mean_age <- rma(yi = logit_prev, 
                               vi = var_logit, 
                               mods = ~ mean_age, 
                               data = prism_gold)
summary(reg_prism_gold_mean_age)

# 5. Proportion of females in study
reg_prism_gold_females_percent <- rma(yi = logit_prev, 
                               vi = var_logit, 
                               mods = ~ females_percent, 
                               data = prism_gold)
summary(reg_prism_gold_females_percent)

# 6. Proportion of current smokers in study
reg_prism_gold_current_smoker_percent <- rma(yi = logit_prev, 
                                      vi = var_logit, 
                                      mods = ~ current_smoker_percent, 
                                      data = prism_gold)
summary(reg_prism_gold_current_smoker_percent)

# Multivariate regression
reg_prism_gold_multi <- rma(yi = logit_prev, 
                               vi = var_logit, 
                               mods = ~ mean_age+bronchod+females_percent, 
                               data = prism_gold)
summary(reg_prism_gold_multi)


# Overall RSP prevalence #####################################################################################
# Logit transforming prevalence and calculating variance
rsp_gold$logit_prev <- log(rsp_gold$overall_reg / (1 - rsp_gold$overall_reg))
rsp_gold$var_logit <- 1 / (rsp_gold$sample_size * rsp_gold$overall_reg * (1 - rsp_gold$overall_reg))

# Univariate regressions ####
# 1. publication year
reg_rsp_gold_pub_year <- rma(yi = logit_prev, 
                               vi = var_logit, 
                               mods = ~ pub_year, 
                               data = rsp_gold)
summary(reg_rsp_gold_pub_year)

# 2. geo_location
# Having to make a new subset because of study 603
geo_location_rsp_gold <-filter(rsp_gold, id != 603) # Study ID 603 has been removed as it is conducted in multiple WHO locations

summary(geo_location_rsp_gold$geo_location)

# Logit transforming prevalence and calculating variance
geo_location_rsp_gold$logit_prev <- log(geo_location_rsp_gold$overall_reg / (1 - geo_location_rsp_gold$overall_reg))
geo_location_rsp_gold$var_logit <- 1 / (geo_location_rsp_gold$sample_size * geo_location_rsp_gold$overall_reg * (1 - geo_location_rsp_gold$overall_reg))

reg_rsp_gold_geo_location <- rma(yi = logit_prev, 
                                   vi = var_logit, 
                                   mods = ~ geo_location, 
                                   data = geo_location_rsp_gold)
summary(reg_rsp_gold_geo_location)

# 3. bronchod
reg_rsp_gold_bronchod <- rma(yi = logit_prev, 
                               vi = var_logit, 
                               mods = ~ bronchod, 
                               data = rsp_gold)
summary(reg_rsp_gold_bronchod)

# 4. mean_age
reg_rsp_gold_mean_age <- rma(yi = logit_prev, 
                               vi = var_logit, 
                               mods = ~ mean_age, 
                               data = rsp_gold)
summary(reg_rsp_gold_mean_age)

# 5. Proportion of females in study
#  Making a new subset to remove studies with no data on female proportion
females_percent_rsp_gold <- rsp_gold %>% 
  drop_na(females_percent)

# Logit transforming prevalence and calculating variance
females_percent_rsp_gold$logit_prev <- log(females_percent_rsp_gold$overall_reg / (1 - females_percent_rsp_gold$overall_reg))
females_percent_rsp_gold$var_logit <- 1 / (females_percent_rsp_gold$sample_size * females_percent_rsp_gold$overall_reg * (1 - females_percent_rsp_gold$overall_reg))

reg_rsp_gold_females_percent <- rma(yi = logit_prev, 
                               vi = var_logit, 
                               mods = ~ females_percent, 
                               data = females_percent_rsp_gold)
summary(reg_rsp_gold_females_percent)

# 6. Proportion of current smokers in study
#  Making a new subset to remove studies with no data on current smoker proportion
current_smoker_percent_rsp_gold <- rsp_gold %>% 
  drop_na(current_smoker_percent)

current_smoker_percent_rsp_gold$logit_prev <- log(current_smoker_percent_rsp_gold$overall_reg / (1 - current_smoker_percent_rsp_gold$overall_reg))
current_smoker_percent_rsp_gold$var_logit <- 1 / (current_smoker_percent_rsp_gold$sample_size * current_smoker_percent_rsp_gold$overall_reg * (1 - current_smoker_percent_rsp_gold$overall_reg))

reg_rsp_gold_current_smoker_percent <- rma(yi = logit_prev, 
                                             vi = var_logit, 
                                             mods = ~ current_smoker_percent, 
                                             data = current_smoker_percent_rsp_gold)
summary(reg_rsp_gold_current_smoker_percent)

# Multivariate regression
reg_rsp_gold_multi <- rma(yi = logit_prev, 
                            vi = var_logit, 
                            mods = ~ mean_age+bronchod, 
                            data = rsp_gold)
summary(reg_rsp_gold_multi)
