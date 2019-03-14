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
  source("models/NivoPKPD.R")
  
# Source functions utility
  source("scripts/functions_utility.R")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Set up input data frame
# Define values for individual to simulate
  set.seed(123456) 
  n <- 1000  # Number of individuals
  ID <- 1:n  # Sequence of individual ID's

  ECOG <- 0  # ECOG Performance Status
  UEVENT  <- runif(n, 0, 1)  # Chance of event
  UCENSOR <- runif(n, 0, 1)  # Chance of dropout

# Populaton parameter variability
# Source etas from NONMEM output
  nmdata <- read.csv("output/nm_placebo_haz_etas.nm7/nm_placebo_haz_etas.fit.csv",
      stringsAsFactors = F, na.strings = ".")
  etadata <- nmdata[nmdata$TIME == 0,]

# Create data frame of individuals with varying demographics and ETA values
  pop.data <- data.frame(ID)
  pop.data$ECOG <- ECOG
  pop.data$UEVENT <- etadata$UEVENT
  pop.data$UCENSOR <- etadata$UCENSOR
  
  pop.data$ETA1 <- etadata$ETA1
  pop.data$ETA2 <- etadata$ETA2
  pop.data$ETA3 <- etadata$ETA3
  pop.data$ETA4 <- etadata$ETA4
  pop.data$ETA5 <- etadata$ETA5
  pop.data$ETA6 <- etadata$ETA6
  pop.data$ETA7 <- etadata$ETA7
  pop.data$ETA8 <- etadata$ETA8

# Check population data
  dim(pop.data)
  plot(density(pop.data$UEVENT))
  plot(density(pop.data$UCENSOR))

# Subset population data for regimen
  pop.data <- subset(pop.data, 
    select = c(
      ID, ECOG, UEVENT, UCENSOR, 
      ETA1, ETA2, ETA3, ETA4, ETA5, ETA6, ETA7, ETA8
    )
  )