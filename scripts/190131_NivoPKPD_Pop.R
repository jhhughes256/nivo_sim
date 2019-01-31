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
  source("models/NivoPKPD_Hours.R")
  
# Source functions utility
  source("scripts/functions_utility.R")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Set up input data frame
# Define values for individual to simulate
  set.seed(123456) 
  n <- 12  # Number of individuals
  ID <- 1:n  # Sequence of individual ID's

  ECOG <- rbinom(n, size = 1, prob = 0.64)  # ECOG Performance Status
  UEVENT  <- runif(n, 0, 1)  # Chance of event
  UCENSOR <- runif(n, 0, 1)  # Chance of dropout

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

# Allocate individuals distribution (standard deviation)
  ETA1 <- rnorm(n, mean = 0, sd = sqrt(ZCL))
  ETA2 <- rnorm(n, mean = 0, sd = sqrt(ZVC))
  ETA3 <- rnorm(n, mean = 0, sd = sqrt(ZVP))
  ETA4 <- rnorm(n, mean = 0, sd = sqrt(ZEMAX))
  ETA5 <- rnorm(n, mean = 0, sd = sqrt(ZEC50))
  ETA6 <- rnorm(n, mean = 0, sd = sqrt(ZTG)) 
  ETA7 <- rnorm(n, mean = 0, sd = sqrt(ZR)) 
  ETA8 <- rnorm(n, mean = 0, sd = sqrt(ZHZ)) 

# Create data frame of individuals with varying demographics and ETA values
  pop.data <- data.frame(ID)
  pop.data$ECOG <- ECOG
  pop.data$UEVENT <- UEVENT
  pop.data$UCENSOR <- UCENSOR
  
  pop.data$ETA1 <- ETA1
  pop.data$ETA2 <- ETA2
  pop.data$ETA3 <- ETA3
  pop.data$ETA4 <- ETA4
  pop.data$ETA5 <- ETA5
  pop.data$ETA6 <- ETA6
  pop.data$ETA7 <- ETA7
  pop.data$ETA8 <- ETA8

# Check population data
  dim(pop.data)
  hist(pop.data$ECOG)
  plot(density(pop.data$UEVENT))
  plot(density(pop.data$UCENSOR))

# Subset population data for regimen
  pop.data <- subset(pop.data, 
    select = c(
      ID, ECOG, UEVENT, UCENSOR, 
      ETA1, ETA2, ETA3, ETA4, ETA5, ETA6, ETA7, ETA8
    )
  )