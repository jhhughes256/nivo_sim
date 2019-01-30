#Nivolumab - Create population of representative patients for simulations

# ------------------------------------------------------------------------------
# Clear workspace
rm(list=ls(all=TRUE))
graphics.off()

# Set working directory
setwd("C:/Users/Ugo/Desktop/TestPK/")

# # Load pakage libraries
library(dplyr)	    # New plyr - required for mrgsolve
library(mrgsolve)	  # Metrum differential equation solver for pharmacometrics
library(MASS)	      # mvrnorm function
library(MBESS)	    # cor2cov function

# Source PopPK model script
source("NivolumabModelSurvivalHours.R")

# ------------------------------------------------------------------------------
# Set up an input data frame
# Define values for individual to simulate
set.seed(123456) 
n           <-  12                                                             # Number of individuals
ID          <-  1:n                                                               # Sequence of individual ID's


ECOG  <- rbinom(n, size = 1, prob = 0.64)
CMT3 <- 60
UEVENT  <- runif(n,0,1)
UCENSOR <- runif(n,0,1)

# Populaton parameter variability
omega.block <- as.matrix(omat(mod))                                              # Omega block values from model
ZCL       <- omega.block[1,1]                                                    # ZCL (variance) from model
ZVC       <- omega.block[2,2]                                                    # ZVC (variance) from model
ZVP       <- omega.block[3,3]                                                    # ZVP (variance) from model
ZEMAX     <- omega.block[4,4]                                                    # ZTMAX (variance) from model
ZEC50     <- omega.block[5,5]
ZTG       <- omega.block[6,6]
ZR        <- omega.block[7,7]
ZHZ       <- omega.block[8,8]

ETA1        <- rnorm(n,mean = 0,sd = sqrt(ZCL))                                # Allocate individuals distribution of ZCL (SD)
ETA2        <- rnorm(n,mean = 0,sd = sqrt(ZVC))                                # Allocate individuals distribution of ZVC (SD)
ETA3        <- rnorm(n,mean = 0,sd = sqrt(ZVP))                                # Allocate individuals distribution of ZVP (SD)
ETA4        <- rnorm(n,mean = 0,sd = sqrt(ZEMAX))                              # Allocate individuals distribution of ZTMAX (SD)
ETA5        <- rnorm(n,mean = 0,sd = sqrt(ZEC50)) 
ETA6        <- rnorm(n,mean = 0,sd = sqrt(ZTG)) 
ETA7        <- rnorm(n,mean = 0,sd = sqrt(ZR)) 
ETA8        <- rnorm(n,mean = 0,sd = sqrt(ZHZ)) 

# Create data frame of individuals with varying demographics and ETA values
population.data                                     <- data.frame(ID)
population.data$ECOG                                <- ECOG
population.data$CMT3                                <- CMT3
population.data$UEVENT                              <- UEVENT
population.data$UCENSOR                             <- UCENSOR
population.data$ETA1                                <- ETA1
population.data$ETA2                                <- ETA2
population.data$ETA3                                <- ETA3
population.data$ETA4                                <- ETA4
population.data$ETA5                                <- ETA5
population.data$ETA6                                <- ETA6
population.data$ETA7                                <- ETA7
population.data$ETA8                                <- ETA8


population.data <- subset(population.data, select=c(ID,ECOG,CMT3,UEVENT,UCENSOR,ETA1, ETA2, ETA3, ETA4,ETA5,ETA6,ETA7,ETA8))