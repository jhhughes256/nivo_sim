# Nivolumab PKPD Model - Population
# ------------------------------------------------------------------------------
# Create population of representative patients for simulations
# Updated to correlate ECOG and TUM_0
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare workspace
# Clear workspace
  # rm(list=ls(all=TRUE))
  # graphics.off()

# Set working directory
# if not working with RStudio project place working directory here
# setwd("C:/.../nivo_sim/")

# # Load package libraries
  # library(dplyr)	# dplyr required for mrgsolve
  # library(mrgsolve)	  # Metrum differential equation solver for pharmacometrics
  library(MASS)  # mvrnorm in trunc_mvrnorm
  library(MBESS)  # cor2cov in trunc_mvrnorm

# Source external scripts
  # source("scripts/functions_utility.R")  # functions utility
  source("models/NivoPKTS_Final.R")  # PopPK model script 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Define demographic data
  set.seed(123456) 
  
# Age - 61.12 years (11.12) [23 - 87]
# Weight - 79.09 kg (19.28) [34 - 168]
# Albumin - 3.92 g/dl (0.45) [2.5 - 5.1]
# Sex - 66.7% male, 33.3% female
# ECOG - 38.73% 0, 61.26% >0
# Number of individuals
  nid <- 1000  # Number of individuals
  ID <- 1:nid  # Sequence of individual ID's

# Continuous
  mean_AGE <- 61.12
  sd_AGE <- 11.12
  range_AGE <- c(23, 87)
  
  mean_BWT <- 79.09
  sd_BWT <- 19.28
  range_BWT <- c(34, 168)
  
  mean_ALB <- 3.92
  sd_ALB <- 0.45
  range_ALB <- c(2.5, 5.1)
  
# Categorical
  sex_prob <- 0.667  # males
  ecog_prob <- 0.6126  # >0
  
# Non-literature Values
# Correlate SeCr with SEX (20% correlation)
# SeCr (males) - 90.6 umol/L
# SeCr (females) - 71.5 umol/L
# give mean eGFR when using mean covariate values
  mean_SECR_M <- 90.6*1.04
  mean_SECR_F <- 71.5*1.04
  sd_SECR <- 28
  range_SECR <- c(0, 300)
  
# Correlate baseline tumour size with ECOG (20% correlation)
# Final distribution designed to have
# Average tumour size - 54.6 mm (80% CV) [10 - 300]
  mean_BTS_0 <- 43.4
  sd_BTS_0 <- 35
  mean_BTS_1 <- 62.2
  sd_BTS_1 <- 52.6
  range_BTS <- c(10, 300)
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Sample covariate values from distributions
# Age, Albumin and initial tumour size - normal distribution
  AGE <- trunc_rnorm(n = nid, mean = mean_AGE, sd = sd_AGE, range = range_AGE)
  ALB <- trunc_rnorm(n = nid, mean = mean_ALB, sd = sd_ALB, range = range_ALB)
  TUM_0 <- rnorm(n = nid, mean = 54.6, sd = 10.92)
  
# Sex & ECOG - binomial distribution
  SEX <- rbinom(nid, 1, sex_prob)
  ECOG <- rbinom(nid, 1, ecog_prob)
 
# Weight, Serum Creatinine and Base Tumour Size - log-normal distribution
  BWT <- trunc_rnorm(n = nid, mean = mean_BWT, sd = sd_BWT, 
    range = range_BWT, log = T)
  
  SECR_M <- trunc_rnorm(n = length(SEX[SEX == 1]), mean = mean_SECR_M, 
    sd = sd_SECR, range = range_SECR, log = T)
  SECR_F <- trunc_rnorm(n = length(SEX[SEX == 0]), mean = mean_SECR_F, 
    sd = sd_SECR, range = range_SECR, log = T)
  
  BTS_0 <- trunc_rnorm(n = length(ECOG[ECOG == 0]), mean = mean_BTS_0, 
    sd = sd_BTS_0, range = range_BTS, log = T)
  BTS_1 <- trunc_rnorm(n = length(ECOG[ECOG == 1]), mean = mean_BTS_1, 
    sd = sd_BTS_1, range = range_BTS, log = T)
  
# Bind covariate data.frame
  cov_df <- data.frame(AGE, ALB, BWT, SEX, ECOG)
  cov_df$SECR <- 0
  cov_df$SECR[cov_df$SEX == 1] <- SECR_M
  cov_df$SECR[cov_df$SEX == 0] <- SECR_F
  cov_df$TUM_0 <- 0
  cov_df$TUM_0[cov_df$ECOG == 1] <- BTS_1
  cov_df$TUM_0[cov_df$ECOG == 0] <- BTS_0
 
# Calculated eGFR
  cov_df$GFR <- dplyr::select(cov_df, SECR, AGE, SEX) %>%
    purrr::pmap(function(SECR, AGE, SEX) ckdepi_fn(SECR, AGE, SEX, 0)) %>%
    unlist()
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Population parameter variability
# Extract omega block values from model and allocate individual distributions
  ETABSV <- mrgsolve::omat(mod) %>%  # Omega block values from model
    as.matrix() %>%  # convert to matrix
    diag()
  ETA_df <- ETABSV %>%  # extract diagonal elements
    purrr::map_dfc(function(Z) rnorm(nid, mean = 0, sd = sqrt(Z)))  # allocate
  names(ETA_df) <- paste0("ETA", 1:9)
  ETA_df$UEVENT <- runif(nid, 0, 1)
  ETA_df$UCENSOR <- runif(nid, 0, 1)

# Create data frame of individuals with varying demographics and ETA values
  pop_df <- dplyr::bind_cols(ID = ID, cov_df, ETA_df) %>%
    dplyr::select(-SECR)
  