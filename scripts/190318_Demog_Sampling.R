# Determine covariate distributions for Nivo Population
# -----------------------------------------------------------------------------
# Create a population that is representative of the data used to develop the
#   Nivolumab models shown in:

# Bajaj G, Wang X, Agrawal S, Gupta M, Roy A, Feng Y. Model-Based Population
#   Pharmacokinetic Analysis of Nivolumab in Patients With Solid Tumors. CPT
#   Pharmacometrics Syst Pharmacol. 2017;6(1):58-66.
# Liu C, Yu J, Li H, Liu J, Xu Y, Song P, et al. Association of time-varying
#   clearance of nivolumab with disease dynamics and its implications on 
#   exposure response analysis. Clin Pharmacol Ther. 2017;101(5):657-66.

# This population is required to meet the following criteria:
# - Must include weight, height, age, sex, ECOG, serum creatinine and albumin
# - Albumin variability representative of metastatic melanoma patients
# - Weight and height are appropriately correlated
# - eGFR as determined by CKD-EPI representative of metastatic melanoma patients
# - Cancer type as metastatic melanoma
# - Unknown ADA assay results
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare work environment
# Clear workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Load libraries
  library(MASS)
  library(MBESS)

# Load functions utility
  source("C:/Users/hugjh001/Documents/nivo_sim/scripts/functions_utility.R")

# Define known demographic data from Bajaj et al.
# Age - 61.12 years (11.12) [23 - 87]
# Weight - 79.09 kg (19.28) [34 - 168]
# Sex - 66.7% male, 33.3% female
# ECOG - 38.73% 0, 61.26% >0
# Number of individuals
  nid <- 100000

# Continuous
  mean.AGE <- 61.12
  sd.AGE <- 11.12
  range.AGE <- c(23, 87)
  
  mean.WT <- 79.09
  sd.WT <- 19.28
  range.WT <- c(34, 168)
  
# Categorical
  sex.prob <- 0.667  # males
  ecog.prob <- 0.6126  # >0

# To determine eGFR height and SeCr must be determined
# eGFR - 78.49 mL/min/1.73m^2 (21.63)
  
# To determine albumin, values from Datta M et al. used
# Stage III Melanoma - 4.15 g/dl (0.33) [2.5 - 5.1]
# Stage IV Melanoma - 3.92 g/dl (0.45) [2.5 - 5.1]
  mean.ALB <- 3.92
  sd.ALB <- 0.45
  range.ALB <- c(2.5, 5.1)
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Determine mean height
# Based off of weight and BMI from Cortellini et al.
# Weight - 71 kg [35 - 139]
# BMI - 24.9 [13.5 - 46.6]
# BMI Category - Underweight 4.1%, Normal 46.3%, Overweight 38.6%, Obese 11%
# As BMI = wt/ht^2; ht = sqrt(wt/bmi)
  mean.HT <- sqrt(71/24.9)*100
  
# Determine standard deviation of height distribution
# Create vectors for mean, sd and range (upper and lower)
# Range required as a truncated multivariate distribution is used
  mean.WT.HT <- c(71, mean.HT)
  sd.WT.HT <- c(16, 20)
  lower.WT.HT <- c(range.WT[1], 150)
  upper.WT.HT <- c(range.WT[2], 195)
  
# Set up correlation matrix for height and weight
  corr.WT.HT <- matrix(c(
    1.0, 0.7,
    0.7, 1.0
  ), nrow = 2, ncol = 2)  # symmetrical correlation matrix
  
# Use trunc.mvrnorm function
  corrsim <- trunc.mvrnorm(n = nid, mean = mean.WT.HT, sd = sd.WT.HT, 
    corr_mat = corr.WT.HT, lower = lower.WT.HT, upper = upper.WT.HT, log = T)
  
# Determine BMI for samples
  samp.bmi <- apply(corrsim, 1, bmi.fn)
  
# Look at results
  hist(corrsim[,2])
  c(
    mean.WT = 71,
    sd.WT = 16,
    mean.HT = mean.HT,
    sd.HT = sd.WT.HT[2]
  )
  c(
    "UW" = length(samp.bmi[samp.bmi <= 18.5])/1000,
    "N" = length(samp.bmi[samp.bmi > 18.5 & samp.bmi <= 25])/1000,
    "OW" = length(samp.bmi[samp.bmi > 25 & samp.bmi < 30])/1000,
    "O" = length(samp.bmi[samp.bmi >= 30])/1000
  )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Determine mean serum creatinine
# CKD-EPI equation used to determine eGFR
  
# eGFR = 141 * min(\frac{[Cr]}{\kappa}, 1)^{\alpha} *
#        max(\frac{[Cr]}{\kappa},1)^{-1.209} * 0.993^{Age} * 
#        1.018^{female} * 1.159^{black} * 1.73/BSA
  
# When substituting [Cr] for 0.9 and 0.7 for males and females respectively
#   the mean eGFR is underpredicted allowing for the following simplication

# eGFR = 141 * \bigg(\frac{[Cr]}{\kappa}\bigg)^{-1.209} *
#        0.993^{age} * 1.018^{female} * 1.159^{black}
  
# And therefore you can solve for [Cr]
  
# [Cr] = \kappa * \bigg(\frac{eGFR}
#        {141*0.993^{age}*1.018^{female}*1.159^{black}}\bigg)^{-\frac{1}{1.209}}
  
# The weighted mean can be determined using the mean serum creatinine determined
#   for both males and females along with the proportion of males to females
#   from Bajaj et al.
  
# While determination of the mean for serum creatinine can be analytical, the
#   standard deviation and range are determined empirically, by altering values
#   for SECR and HT until eGFR matches the distributions from Bajaj et al.

# Define demographic data
# Some objects have been defined above
  mean.SECR <- 78.49
  sd.SECR <- 10
  range.SECR <- c(0, 150)
  
  mean.HT <- sqrt(71/24.9)*100
  sd.HT <- 10
  
# Weight and Height - multivariate log-normal distribution
# Some objects have been defined above
  mean.WT.HT <- c(mean.WT, mean.HT)
  sd.WT.HT <- c(sd.WT, sd.HT)
  
  corrsim <- trunc.mvrnorm(n = nid, mean = mean.WT.HT, sd = sd.WT.HT, 
    corr_mat = corr.WT.HT, lower = lower.WT.HT, upper = upper.WT.HT, log = T)
  
  WT <- corrsim[,1]
  HT <- corrsim[,2]
  
# Age - normal distribution
  AGE <- trunc.rnorm(n = nid, mean = mean.AGE, sd = sd.AGE, range = range.AGE)
  
# Sex - binomial distribution
  SEX <- rbinom(nid, 1, sex.prob)
  
# Serum Creatinine - log-normal distribution
  SECR <- trunc.rnorm(n = nid, mean = mean.SECR, sd = sd.SECR, 
    range = range.SECR, log = T)
  
# eGFR
  cov.df <- data.frame(WT, HT, AGE, SEX, SECR)
  eGFR <- apply(cov.df, 1, function(df) {
    ckdepi.fn(df["SECR"], df["AGE"], df["SEX"], 0)*1.73/bsa.fn(df["WT"], df["HT"])
  })
  
  mean(eGFR)
  sd(eGFR)
  length(eGFR[eGFR < 39])/nid
  length(eGFR[eGFR >= 30 & eGFR < 60])/nid
  length(eGFR[eGFR >= 60 & eGFR < 90])/nid
  length(eGFR[eGFR >= 90])/nid
  