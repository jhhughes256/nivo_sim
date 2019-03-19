# Nivolumab PKPD Model - Population
# ------------------------------------------------------------------------------
# Create population of representative patients for simulations
# Cleaned up for use within git directory
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare workspace
# Clear workspace
  # rm(list=ls(all=TRUE))
  # graphics.off()

# Set working directory
# if not working with RStudio project place working directory here
# setwd("C:/.../nivo_sim/")

# # Load package libraries
  # library(dplyr)	    # New plyr - required for mrgsolve
  # library(mrgsolve)	  # Metrum differential equation solver for pharmacometrics
  library(MASS)	      # mvrnorm function
  library(MBESS)	    # cor2cov function

# Source PopPK model script
  source("models/NivoPKTS_wPD.R")
  
# Source functions utility
  source("scripts/functions_utility.R")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Set up input data frame
# Define values for individual to simulate
  set.seed(123456) 
  n <- 10  # Number of individuals
  ID <- 1:n  # Sequence of individual ID's

  SEX <- rbinom(n, size = 1, prob = 0.66)  # Sex (1 = male, 0 = female)
  ALB <- rnorm(n, mean = 0.4, sd = 0)  # Albumin 
  RCC <- rbinom(n, size = 1, prob = 0.5)  # Renal cell carninoma
  OTHERRCC <- -RCC + 1  # those not with RCC have other cancer
  ADApos <- rbinom(n, size = 1, prob = 0.5) # Anti-Drug Antibodies
  ADAunk <- -ADApos + 1  # those not ADA positive are unknown
  TS <- 60  # tumour size
  SQNSQ <- 1  # squamous+non-squamous 
  ECOG <- 0  # ECOG Performance Status
  UEVENT  <- runif(n, 0, 1)  # Chance of event
  UCENSOR <- runif(n, 0, 1)  # Chance of dropout
  
# Correlation matrix for WT, AGE, GFR
# Define mean and standard deviation for correlated values
  mean.AGE <- 60  # Age (years)
  sd.AGE <- 15
  mean.BWT <- 80  # Weight (kg)
  sd.BWT <- 20
  mean.GFR <- 80  # Glomerular Filtration Rate (mL/min)
  sd.GFR <- 15
  
# Define correlation between values and define correlation matrix
  x11 <- 1  # Correlation AGE vs AGE
  x12 <- 0.1  # Correlation AGE vs BWT
  x13 <- 0.5  # Correlation AGE vs GFR
  x22 <- 1  # Correlation BWT vs BWT
  x23 <- 0.1  # Correlation BWT vs GFR
  x33 <- 1  # Correlation GFR vs GFR
  corr1 <- matrix(c(x11, x12, x13, x12, x22, x23, x13, x23, x33), 3, 3)
  
# Create truncated multivariate normal distribution
# Function converts correlation to covariance matrix
# Removes individuals with values outside of clinical relevance
# BWT 45-120; AGE 35-90; GFR 50-130
  lower.lim <- c(45, 35, 50)
  upper.lim <- c(120, 90, 130)
  cor.values <- trunc_mvrnorm(n, mean = c(mean.AGE, mean.BWT, mean.GFR), 
    sd = c(sd.AGE, sd.BWT, sd.GFR), corr_mat = corr1, lower.lim, upper.lim
  )

# Populaton parameter variability
  omega.block <- as.matrix(omat(mod))  # Omega block values from model
  ZCL <- omega.block[1, 1]  # ZCL (variance) from model
  ZVC <- omega.block[2, 2]  # ZVC (variance) from model
  ZVP <- omega.block[3, 3]  # ZVP (variance) from model
  ZEMAX <- omega.block[5, 5]  # ZEMAX (variance) from model
  ZEC50 <- omega.block[6, 6]  # ZEC50 (variance) from model
  ZTG <- omega.block[6, 6]  # ZTG (variance) from model
  ZR <- omega.block[7, 7]  # ZR (variance) from model
  ZHZ <- omega.block[8, 8]  # ZHZ (variance) from model
  ZTMAX <- omega.block[8, 8]  # ZTMAX (variance) from model

# Allocate individuals distribution (standard deviation)
  ETA1 <- rnorm(n, mean = 0, sd = sqrt(ZCL))
  ETA2 <- rnorm(n, mean = 0, sd = sqrt(ZVC))
  ETA3 <- rnorm(n, mean = 0, sd = sqrt(ZVP))
  ETA4 <- rnorm(n, mean = 0, sd = sqrt(ZEMAX))
  ETA5 <- rnorm(n, mean = 0, sd = sqrt(ZEC50))
  ETA6 <- rnorm(n, mean = 0, sd = sqrt(ZTG)) 
  ETA7 <- rnorm(n, mean = 0, sd = sqrt(ZR)) 
  ETA8 <- rnorm(n, mean = 0, sd = sqrt(ZHZ)) 
  ETA9 <- rnorm(n, mean = 0, sd = sqrt(ZTMAX))


# Create data frame of individuals with varying demographics and ETA values
  
  
# Check population data
  dim(pop.data)
  plot(density(pop.data$UEVENT))
  plot(density(pop.data$UCENSOR))

# Subset population data for regimen
  pop.data <- subset(pop.data, 
    select = c(
      ID, ECOG, UEVENT, UCENSOR, 
      ETA1, ETA2, ETA3, ETA4, ETA5, ETA6, ETA7, ETA8, ETA9
    )
  )