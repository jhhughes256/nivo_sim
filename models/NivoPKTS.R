# Carboplatin - Create population of representative patients for simulations

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
source("NivolumabwithTS.R")

# ------------------------------------------------------------------------------
# Set up an input data frame
# Define values for individual to simulate
set.seed(123456) 
n           <-  12                                                              # Number of individuals
ID          <-  1:n                                                               # Sequence of individual ID's

BWT         <- rnorm(n,mean = 80,sd = 20)                                        # Weight (kg)
SEX         <- rbinom(n, size = 1, prob = 0.66)                                  # Sex (1 = male, 0 = female)
AGE         <- rnorm(n, mean = 60, sd = 15)                                      # Age (years)
PS          <- rbinom(n, size = 1, prob = 0.64)                                  # PS (0 or 1)
ALB         <- rnorm(n, mean = 0.4, sd = 0)                                      # Albumin 
GFR         <- rnorm(n, mean = 80, sd = 15)                                      # GFR
RCC         <- rbinom(n, size =1, prob =0.5)
if (RCC == 1) {
  OTHERRCC <- 0} else {
    OTHERRCC <- 1}
ADApos        <- rbinom(n, size =1, prob =0.5)
if (ADApos == 1) {
  ADAunk <- 0} else {
    ADAunk <- 1}

SQNSQ <-1
# Correlation matrix for WT, AGE, GFR
mean.AGE    <-  60                                                               # Age (years)
sd.AGE      <-  15
mean.BWT    <-  80                                                              # Weight (kg)
sd.BWT      <-  20
mean.GFR    <-  80
sd.GFR      <-  15

x11         <-  1                                                                 # Correlation AGE vs AGE
x12         <-  0.1                                                               # Correlation AGE vs BWT
x13         <-  0.5                                                               # Correlation AGE vs GFR
x22         <-  1                                                                 # Correlation BWT vs BWT
x23         <-  0.1
x33         <-  1

corr1       <-  matrix(c(x11,x12,x13,x12,x22,x23,x13,x23,x33),3,3)
cov.mat1    <-  cor2cov(cor.mat = corr1,sd = c(sd.AGE, sd.BWT, sd.GFR))
cor.values <-  mvrnorm(n, mu=c(mean.AGE, mean.BWT, mean.GFR), Sigma = cov.mat1)


# Populaton parameter variability
omega.block <- as.matrix(omat(mod))                                              # Omega block values from model
ZCL       <- omega.block[1,1]                                                    # ZCL (variance) from model
ZVC       <- omega.block[2,2]                                                    # ZVC (variance) from model
ZVP       <- omega.block[3,3]                                                    # ZVP (variance) from model
ZTMAX     <- omega.block[4,4]                                                    # ZTMAX (variance) from model
ZEMAX     <- omega.block[5,5]                                                    # ZEMAX (variance) from model
ZEC50     <- omega.block[6,6]                                                    # ZEC50 (variance) from model

ETA1        <- rnorm(n,mean = 0,sd = sqrt(ZCL))                                # Allocate individuals distribution of ZCL (SD)
ETA2        <- rnorm(n,mean = 0,sd = sqrt(ZVC))                                # Allocate individuals distribution of ZVC (SD)
ETA3        <- rnorm(n,mean = 0,sd = sqrt(ZVP))                                # Allocate individuals distribution of ZVP (SD)
ETA4        <- rnorm(n,mean = 0,sd = sqrt(ZTMAX))                              # Allocate individuals distribution of ZTMAX (SD)
ETA5        <- rnorm(n,mean = 0,sd = sqrt(ZEMAX))                              # Allocate individuals distribution of ZEMAX (SD)
ETA6        <- rnorm(n,mean = 0,sd = sqrt(ZEC50))                              # Allocate individuals distribution of ZEC50 (SD)

# Create data frame of individuals with varying demographics and ETA values
population.data                                     <- data.frame(ID)
population.data$SEX                                 <- SEX
population.data$BWT                                 <- cor.values[,2]
population.data$AGE                                 <- cor.values[,1]
population.data$ALB                                 <- ALB
population.data$TS                                  <- TS
population.data$GFR                                 <- cor.values[,3]
population.data$GFR[population.data$SEX == 0]       <- cor.values[population.data$SEX == 0 ,3]*0.85    
population.data$PS                                  <- PS
population.data$RCC                                 <- RCC
population.data$OTHERRCC                            <- OTHERRCC
population.data$ADApos                              <- ADApos
population.data$ADAunk                              <- ADAunk
population.data$SQNSQ                               <- SQNSQ
# 
population.data$ETA1                                <- ETA1
population.data$ETA2                                <- ETA2
population.data$ETA3                                <- ETA3
population.data$ETA4                                <- ETA4
population.data$ETA5                                <- ETA5
population.data$ETA6                                <- ETA6

# Remove individuals with values outside of clinical relevance
population.data                                     <- population.data[population.data$BWT  > 45  &  population.data$BWT <120,]
population.data                                     <- population.data[population.data$AGE  > 35  &  population.data$AGE < 80,]
population.data                                     <- population.data[population.data$GFR  > 50  &  population.data$GFR <130,]



dim(population.data)


hist(population.data$SEX)
plot(density(population.data$BWT))
plot(density(population.data$AGE))
plot(density(population.data$GFR))

population.data <- subset(population.data, select=c(ID, SEX, BWT, AGE, GFR,ALB, ETA1, ETA2, ETA3, ETA4, ETA5,ETA6))