# Nivolumab PK Model - Population
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
  source("models/NivoPK.R")
  
# Source functions utility
  source("scripts/functions_utility.R")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Set up input data frame
# Define values for individual to simulate
  set.seed(123456) 
  n <- 12  # Number of individuals
  ID <- 1:n  # Sequence of individual ID's
  
  SEX <- rbinom(n, size = 1, prob = 0.66)  # Sex (1 = male, 0 = female)
  PS <- rbinom(n, size = 1, prob = 0.64)  # PS (0 or 1)
  ALB <- rnorm(n, mean = 0.4, sd = 0)  # Albumin 
  RCC <- rbinom(n, size = 1, prob = 0.5)  # Renal cell carninoma
  OTHERRCC <- -RCC + 1  # those not with RCC have other cancer
  ADApos <- rbinom(n, size = 1, prob = 0.5) # Anti-Drug Antibodies
  ADAunk <- -ADApos + 1  # those not ADA positive are unknown
  TS <- 60  # tumour size
  SQNSQ <- 1  # squamous+non-squamous 
  
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
  ZTMAX <- omega.block[4, 4]  # ZTMAX (variance) from model
  
# Allocate individuals distribution (standard deviation)
  ETA1 <- rnorm(n, mean = 0, sd = sqrt(ZCL))
  ETA2 <- rnorm(n, mean = 0, sd = sqrt(ZVC))
  ETA3 <- rnorm(n, mean = 0, sd = sqrt(ZVP))
  ETA4 <- rnorm(n, mean = 0, sd = sqrt(ZTMAX))

# Create data frame of individuals with varying demographics and ETA values
  pop.data <- data.frame(ID)
  pop.data$SEX <- SEX
  pop.data$BWT <- cor.values[, 2]
  pop.data$AGE <- cor.values[, 1]
  pop.data$ALB <- ALB
  pop.data$TS <- TS
  pop.data$GFR <- cor.values[, 3] 
  pop.data$GFR[pop.data$SEX == 0] <- cor.values[pop.data$SEX == 0, 3]*0.85
  pop.data$PS <- PS
  pop.data$RCC <- RCC
  pop.data$OTHERRCC <- OTHERRCC
  pop.data$ADApos <- ADApos
  pop.data$ADAunk <- ADAunk
  pop.data$SQNSQ <- SQNSQ
  
  pop.data$ETA1 <- ETA1
  pop.data$ETA2 <- ETA2

# Check population data
  dim(pop.data)
  hist(pop.data$SEX)
  plot(density(pop.data$BWT))
  plot(density(pop.data$AGE))
  plot(density(pop.data$GFR))

# Subset population data for regimen
  pop.data <- subset(pop.data, 
    select = c(ID, SEX, BWT, AGE, GFR, ALB, ETA1, ETA2)
  )