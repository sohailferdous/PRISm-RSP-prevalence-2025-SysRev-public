# Global prevalence of preserved ratio impaired spirometry and restrictive spirometry pattern in the general 
# population: a systematic review and meta-analysis
#
# Supplementary analysis
#
# This script performs meta-analyses for overall and subgroup prevalence of LLN-PRISm, LLN-RSP and combined 
# LLN-PRISm and LLN-RSP
#
# Set working directory to directory with data_prevalence.csv file
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
overall <- read.csv("data_prevalence.csv")

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

# Create LLN dataset
lln <- filter(overall, diag_defn == "lln")

# 1. Overall prevalence ####
#
# prism
prism_lln <-filter(lln, diag == "prism" & overall.prev == "y")

# rsp
rsp_lln <-filter(lln, diag == "rsp" & overall.prev == "y")

# Combined prism and rsp 
comb_lln <-filter(lln, comb.prism.rsp.prev == "y")

# 2. For prevalence by sex ####
#
# Males - combined, prism and rsp groups
male_lln_comb <- filter(lln, gender.prev == "y") %>% 
  drop_na(male_total, male_cases)

male_lln_prism <- filter(lln, diag == "prism" & gender.prev == "y") %>% 
  drop_na(male_total, male_cases)

male_lln_rsp <- filter(lln, diag == "rsp" & gender.prev == "y") %>% 
  drop_na(male_total, male_cases)

# females - combined, prism and rsp groups
female_lln_comb <- filter(lln, gender.prev == "y") %>% 
  drop_na(female_total, female_cases)

female_lln_prism <- filter(lln, diag == "prism" & gender.prev == "y") %>% 
  drop_na(female_total, female_cases)

female_lln_rsp <- filter(lln, diag == "rsp" & gender.prev == "y") %>% 
  drop_na(female_total, female_cases)

# 3. For prevalence by smoking status ####
#
# Non-smokers - combined, prism and rsp groups
non_smoker_lln_comb <- filter(lln, smoking.prev == "y") %>% 
  drop_na(tot_non,case_non)

non_smoker_lln_prism <- filter(lln, diag == "prism" & smoking.prev == "y") %>% 
  drop_na(tot_non,case_non)

non_smoker_lln_rsp <- filter(lln, diag == "rsp" & smoking.prev == "y") %>% 
  drop_na(tot_non,case_non)

# Ex-smokers - combined, prism and rsp groups
ex_smoker_lln_comb <- filter(lln, smoking.prev == "y") %>% 
  drop_na(tot_ex,case_ex)

ex_smoker_lln_prism <- filter(lln, diag == "prism" & smoking.prev == "y") %>% 
  drop_na(tot_ex,case_ex)

ex_smoker_lln_rsp <- filter(lln, diag == "rsp" & smoking.prev == "y") %>% 
  drop_na(tot_ex,case_ex)

# Current smokers - combined, prism and rsp groups
cur_smoker_lln_comb <- filter(lln, smoking.prev == "y") %>% 
  drop_na(tot_cur,case_cur)

cur_smoker_lln_prism <- filter(lln, diag == "prism" & smoking.prev == "y") %>% 
  drop_na(tot_cur,case_cur)

cur_smoker_lln_rsp <- filter(lln, diag == "rsp" & smoking.prev == "y") %>% 
  drop_na(tot_cur,case_cur)

# 4. WHO location ####
#
# Region of the Americas (AMR) - combined, prism and rsp groups
amro_lln_comb <- filter(lln, geo_location == "AMR" & comb.prism.rsp.prev == "y")

amro_lln_prism <- filter(lln, diag == "prism" & geo_location == "AMR"& overall.prev == "y")

amro_lln_rsp <- filter(lln, diag == "rsp" & geo_location == "AMR"& overall.prev == "y")

# European region (EUR) - combined, prism and rsp groups
euro_lln_comb <- filter(lln, geo_location == "EUR" & id != 158 & comb.prism.rsp.prev == "y")

euro_lln_prism <- filter(lln, diag == "prism" & geo_location == "EUR" & overall.prev == "y")

euro_lln_rsp <- filter(lln, diag == "rsp" & geo_location == "EUR" & overall.prev == "y")

# Western Pacific region (WPR) - combined, prism and rsp groups
wpro_lln_comb <- filter(lln, geo_location == "WPR" & comb.prism.rsp.prev == "y")

wpro_lln_prism <- filter(lln, diag == "prism" & geo_location == "WPR" & overall.prev == "y")

wpro_lln_rsp <- filter(lln, diag == "rsp" & geo_location == "WPR" & overall.prev == "y")

# African region (AFR)) - combined, prism and rsp groups
afro_lln_comb <- filter(lln, geo_location == "AFR" & comb.prism.rsp.prev == "y")

afro_lln_prism <- filter(lln, diag == "prism" & geo_location == "AFR" & overall.prev == "y")

afro_lln_rsp <- filter(lln, diag == "rsp" & geo_location == "AFR" & overall.prev == "y")

# Analysis code ##############################################################################################
#
# 1. Overall prevalence ####
#
# Combined PRISm and RSP
#
meta_comb_lln <- comb_lln %>%
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

forest(meta_comb_lln, 
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# Funnel plot
funnel(meta_comb_lln, 
       studlab = TRUE, 
       xlab = "Prevalence (logit-transformed)",
       ylab = "Standard Error",
       main = "Funnel Plot for Publication Bias",
       xlim = c(-4, -0.5))

# Egger's test for asymmetry
egger_test_meta_comb_lln <- metabias(meta_comb_lln, method = "linreg")
print(egger_test_meta_comb_lln)

# PRISm
#
meta_prism_lln <- prism_lln %>%
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

forest(meta_prism_lln, 
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# Funnel plot
funnel(meta_prism_lln, 
       studlab = TRUE, 
       xlab = "Prevalence (logit-transformed)",
       ylab = "Standard Error",
       main = "Funnel Plot for Publication Bias",
       xlim = c(-3.5, -0.5))

# Egger's test for asymmetry
egger_test_meta_prism_lln <- metabias(meta_prism_lln, method = "linreg")
print(egger_test_meta_prism_lln)

# RSP
#
meta_rsp_lln <- rsp_lln %>%
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

forest(meta_rsp_lln, 
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# Funnel plot
funnel(meta_rsp_lln, 
       studlab = TRUE, 
       xlab = "Prevalence (logit-transformed)",
       ylab = "Standard Error",
       main = "Funnel Plot for Publication Bias",
       xlim = c(-4, -0.5))

# Egger's test for asymmetry
egger_test_meta_rsp_lln <- metabias(meta_rsp_lln, method = "linreg")
print(egger_test_meta_rsp_lln)

# 2. Prevalence by sex ####
#
# Males
#
# Combined
meta_male_lln_comb <-male_lln_comb %>%
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

forest(meta_male_lln_comb,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# PRISm
meta_male_lln_prism <-male_lln_prism %>%
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

forest(meta_male_lln_prism,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# RSP
meta_male_lln_rsp <-male_lln_rsp %>%
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

forest(meta_male_lln_rsp,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# Females
#
# Combined
meta_female_lln_comb <-female_lln_comb %>%
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

forest(meta_female_lln_comb,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# PRISm
meta_female_lln_prism <-female_lln_prism %>%
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

forest(meta_female_lln_prism,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# RSP
meta_female_lln_rsp <-female_lln_rsp %>%
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

forest(meta_female_lln_rsp,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# 3. Prevalence by smoking status ####
#
# Non-smokers
#
# Combined
meta_non_smoker_lln_comb <- non_smoker_lln_comb %>% 
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

forest(meta_non_smoker_lln_comb,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# PRISm
meta_non_smoker_lln_prism <- non_smoker_lln_prism %>% 
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

forest(meta_non_smoker_lln_prism,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# RSP
meta_non_smoker_lln_rsp <- non_smoker_lln_rsp %>% 
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

forest(meta_non_smoker_lln_rsp,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# Ex-smokers
#
# Combined
meta_ex_smoker_lln_comb <- ex_smoker_lln_comb %>% 
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

forest(meta_ex_smoker_lln_comb,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# PRISm
meta_ex_smoker_lln_prism <- ex_smoker_lln_prism %>% 
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

forest(meta_ex_smoker_lln_prism,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# RSP
meta_ex_smoker_lln_rsp <- ex_smoker_lln_rsp %>% 
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

forest(meta_ex_smoker_lln_rsp,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# Current smokers
#
# Combined
meta_cur_smoker_lln_comb <- cur_smoker_lln_comb %>% 
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

forest(meta_cur_smoker_lln_comb,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# PRISm
meta_cur_smoker_lln_prism <- cur_smoker_lln_prism %>% 
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

forest(meta_cur_smoker_lln_prism,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# RSP
meta_cur_smoker_lln_rsp <- cur_smoker_lln_rsp %>% 
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

forest(meta_cur_smoker_lln_rsp,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# 4. Prevalence by WHO geographical location ####
#
# AMRO
#
# Combined
meta_amro_lln_comb <- amro_lln_comb %>% 
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

forest(meta_amro_lln_comb,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# RSP
meta_amro_lln_rsp <- amro_lln_rsp %>% 
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

forest(meta_amro_lln_rsp,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# EURO
#
# Combined
meta_euro_lln_comb <- euro_lln_comb %>% 
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

forest(meta_euro_lln_comb,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# RSP
meta_euro_lln_rsp <- euro_lln_rsp %>% 
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

forest(meta_euro_lln_rsp,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# WPRO
#
# Combined
meta_wpro_lln_comb <- wpro_lln_comb %>% 
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

forest(meta_wpro_lln_comb,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# RSP
meta_wpro_lln_rsp <- wpro_lln_rsp %>% 
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

forest(meta_wpro_lln_rsp,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# AFRO
#
# Combined
meta_afro_lln_comb <- afro_lln_comb %>% 
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

forest(meta_afro_lln_comb,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))

# RSP
meta_afro_lln_rsp <- afro_lln_rsp %>% 
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

forest(meta_afro_lln_rsp,
       rightlabs = c("Prevalence (%)", "95%-CI", "Weight (%)"), 
       leftlabs = c("Study","Cases","Total"))
