# Global prevalence of preserved ratio impaired spirometry and restrictive spirometry pattern in the general 
# population: a systematic review and meta-analysis 
#
# This script performs meta-analyses for overall and subgroup prevalence of GOLD-PRISm, GOLD-RSP and combined 
# GOLD-PRISm and GOLD-RSP
#
# Set working directory to directory with v1_data_prevalence.csv file
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
overall <- read.csv("v1_data_prevalence.csv")

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

summary(overall)

# Create GOLD dataset
gold <- filter(overall, diag_defn == "gold") # the main study analysis will only use GOLD dataset

# Creating subsets for meta-analysis #########################################################################
#
# 1. Overall prevalence ####
#
# prism
prism_gold <-filter(gold, diag == "prism" & overall.prev == "y")

# rsp
rsp_gold <-filter(gold, diag == "rsp" & overall.prev == "y")

# Combined prism and rsp 
comb_gold <-filter(gold, comb.prism.rsp.prev == "y")

# 2. For prevalence by sex ####
#
# Males - combined, prism and rsp groups
male_gold_comb <- filter(gold, gender.prev == "y") %>% 
  drop_na(male_total, male_cases)

male_gold_prism <- filter(gold, diag == "prism" & gender.prev == "y") %>% 
  drop_na(male_total, male_cases)

male_gold_rsp <- filter(gold, diag == "rsp" & gender.prev == "y") %>% 
  drop_na(male_total, male_cases)

# females - combined, prism and rsp groups
female_gold_comb <- filter(gold, gender.prev == "y") %>% 
  drop_na(female_total, female_cases)

female_gold_prism <- filter(gold, diag == "prism" & gender.prev == "y") %>% 
  drop_na(female_total, female_cases)

female_gold_rsp <- filter(gold, diag == "rsp" & gender.prev == "y") %>% 
  drop_na(female_total, female_cases)

# 3. For prevalence by smoking status ####
#
# Non-smokers - combined, prism and rsp groups
non_smoker_gold_comb <- filter(gold, smoking.prev == "y") %>% 
  drop_na(tot_non,case_non)

non_smoker_gold_prism <- filter(gold, diag == "prism" & smoking.prev == "y") %>% 
  drop_na(tot_non,case_non)

non_smoker_gold_rsp <- filter(gold, diag == "rsp" & smoking.prev == "y") %>% 
  drop_na(tot_non,case_non)

# Ex-smokers - combined, prism and rsp groups
ex_smoker_gold_comb <- filter(gold, smoking.prev == "y") %>% 
  drop_na(tot_ex,case_ex)

ex_smoker_gold_prism <- filter(gold, diag == "prism" & smoking.prev == "y") %>% 
  drop_na(tot_ex,case_ex)

ex_smoker_gold_rsp <- filter(gold, diag == "rsp" & smoking.prev == "y") %>% 
  drop_na(tot_ex,case_ex)

# Current smokers - combined, prism and rsp groups
cur_smoker_gold_comb <- filter(gold, smoking.prev == "y") %>% 
  drop_na(tot_cur,case_cur)

cur_smoker_gold_prism <- filter(gold, diag == "prism" & smoking.prev == "y") %>% 
  drop_na(tot_cur,case_cur)

cur_smoker_gold_rsp <- filter(gold, diag == "rsp" & smoking.prev == "y") %>% 
  drop_na(tot_cur,case_cur)

# 4. WHO location ####
# Groups were only created for regions with sufficient studies to calculate GOLD prevalence (n>=3)
#
# Region of the Americas (AMR) - combined, prism and rsp groups
amro_gold_comb <- filter(gold, geo_location == "AMR")

amro_gold_prism <- filter(gold, diag == "prism" & geo_location == "AMR"& overall.prev == "y")

amro_gold_rsp <- filter(gold, diag == "rsp" & geo_location == "AMR"& overall.prev == "y")

# European region (EUR) - combined, prism and rsp groups
euro_gold_comb <- filter(gold, geo_location == "EUR" & id != 439)

euro_gold_prism <- filter(gold, diag == "prism" & geo_location == "EUR" & overall.prev == "y")

euro_gold_rsp <- filter(gold, diag == "rsp" & geo_location == "EUR" & overall.prev == "y")

# Western Pacific region (WPR) - combined, prism and rsp groups
wpro_gold_comb <- filter(gold, geo_location == "WPR")

wpro_gold_prism <- filter(gold, diag == "prism" & geo_location == "WPR" & overall.prev == "y")

wpro_gold_rsp <- filter(gold, diag == "rsp" & geo_location == "WPR" & overall.prev == "y")

# African region (AFR)) - combined, prism and rsp groups
afro_gold_comb <- filter(gold, geo_location == "AFR")

afro_gold_prism <- filter(gold, diag == "prism" & geo_location == "AFR" & overall.prev == "y")

afro_gold_rsp <- filter(gold, diag == "rsp" & geo_location == "AFR" & overall.prev == "y")

# ANALYSIS CODE ##############################################################################################
#
# 1. Overall prevalence ####
#
# Combined PRISm and RSP
#
meta_comb_gold <- comb_gold %>%
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
    pscale=100)                         # Increase proportion by 100x (converting prevalence to percentage)

forest(meta_comb_gold, 
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# Funnel plot
funnel(meta_comb_gold, 
       studlab = TRUE, 
       xlab = "Prevalence (logit-transformed)",
       ylab = "Standard Error",
       main = "Funnel Plot for Publication Bias",
       xlim = c(-4, -0.5))

# Egger's test for asymmetry
egger_test_meta_comb_gold <- metabias(meta_comb_gold, method = "linreg")
print(egger_test_meta_comb_gold)

# PRISm
#
meta_prism_gold <- prism_gold %>%
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
    pscale = 100)             

forest(meta_prism_gold, 
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# Funnel plot
funnel(meta_prism_gold, 
       studlab = TRUE, 
       xlab = "Prevalence (logit-transformed)",
       ylab = "Standard Error",
       main = "Funnel Plot for Publication Bias",
       xlim = c(-3.5, -0.5))

# Egger's test for asymmetry
egger_test_meta_prism_gold <- metabias(meta_prism_gold, method = "linreg")
print(egger_test_meta_prism_gold)

# RSP
#
meta_rsp_gold <- rsp_gold %>%
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

forest(meta_rsp_gold, 
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# Funnel plot
funnel(meta_rsp_gold, 
       studlab = TRUE, 
       xlab = "Prevalence (logit-transformed)",
       ylab = "Standard Error",
       main = "Funnel Plot for Publication Bias",
       xlim = c(-4, -0.5))

# Egger's test for asymmetry
egger_test_meta_rsp_gold <- metabias(meta_rsp_gold, method = "linreg")
print(egger_test_meta_rsp_gold)

# 2. Prevalence by sex ####
#
# Males
#
# Combined
meta_male_gold_comb <-male_gold_comb %>%
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

forest(meta_male_gold_comb,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# PRISm
meta_male_gold_prism <-male_gold_prism %>%
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

forest(meta_male_gold_prism,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# RSP
meta_male_gold_rsp <-male_gold_rsp %>%
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

forest(meta_male_gold_rsp,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# Females
#
# Combined
meta_female_gold_comb <-female_gold_comb %>%
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

forest(meta_female_gold_comb,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# PRISm
meta_female_gold_prism <-female_gold_prism %>%
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

forest(meta_female_gold_prism,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# RSP
meta_female_gold_rsp <-female_gold_rsp %>%
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

forest(meta_female_gold_rsp,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# 3. Prevalence by smoking status ####
#
# Non-smokers
#
# Combined
meta_non_smoker_gold_comb <- non_smoker_gold_comb %>% 
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

forest(meta_non_smoker_gold_comb,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# PRISm
meta_non_smoker_gold_prism <- non_smoker_gold_prism %>% 
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

forest(meta_non_smoker_gold_prism,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# RSP
meta_non_smoker_gold_rsp <- non_smoker_gold_rsp %>% 
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

forest(meta_non_smoker_gold_rsp,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# Ex-smokers
#
# Combined
meta_ex_smoker_gold_comb <- ex_smoker_gold_comb %>% 
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

forest(meta_ex_smoker_gold_comb,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# PRISm
meta_ex_smoker_gold_prism <- ex_smoker_gold_prism %>% 
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

forest(meta_ex_smoker_gold_prism,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# RSP
meta_ex_smoker_gold_rsp <- ex_smoker_gold_rsp %>% 
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

forest(meta_ex_smoker_gold_rsp,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# Current smokers
#
# Combined
meta_cur_smoker_gold_comb <- cur_smoker_gold_comb %>% 
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

forest(meta_cur_smoker_gold_comb,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# PRISm
meta_cur_smoker_gold_prism <- cur_smoker_gold_prism %>% 
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

forest(meta_cur_smoker_gold_prism,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# RSP
meta_cur_smoker_gold_rsp <- cur_smoker_gold_rsp %>% 
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

forest(meta_cur_smoker_gold_rsp,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# 4. Prevalence by WHO geographical location ####
#
# AMRO
#
# Combined
meta_amro_gold_comb <- amro_gold_comb %>% 
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

forest(meta_amro_gold_comb,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# PRISm
meta_amro_gold_prism <- amro_gold_prism %>% 
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

forest(meta_amro_gold_prism,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# RSP
meta_amro_gold_rsp <- amro_gold_rsp %>% 
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

forest(meta_amro_gold_rsp,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# EURO
#
# Combined
meta_euro_gold_comb <- euro_gold_comb %>% 
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

forest(meta_euro_gold_comb,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# PRISm
meta_euro_gold_prism <- euro_gold_prism %>% 
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

forest(meta_euro_gold_prism,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# RSP
meta_euro_gold_rsp <- euro_gold_rsp %>% 
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

forest(meta_euro_gold_rsp,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# WPRO
#
# Combined
meta_wpro_gold_comb <- wpro_gold_comb %>% 
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

forest(meta_wpro_gold_comb,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# PRISm
meta_wpro_gold_prism <- wpro_gold_prism %>% 
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

forest(meta_wpro_gold_prism,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# RSP
meta_wpro_gold_rsp <- wpro_gold_rsp %>% 
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

forest(meta_wpro_gold_rsp,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))
