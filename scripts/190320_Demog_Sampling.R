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
  library(ggplot2)
  
# Define colourblind palette and custom palette
  cbPalette <- data.frame(
    grey = "#999999",
    orange = "#E69F00",
    skyblue = "#56B4E9",
    green = "#009E73",
    yellow = "#F0E442",
    blue = "#0072B2",
    red = "#D55E00",
    pink = "#CC79A7",
    stringsAsFactors = F
  )

# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))

# Load functions utility
  source("scripts/functions_utility.R")

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

# To determine eGFR SeCr must be determined
# Units of eGFR as determined by CKD-EPI is mL/min/1.73m^2
# No need to know BSA
# eGFR - 78.49 mL/min/1.73m^2 (21.63)
  
# To determine albumin, values from Datta M et al. used
# Stage III Melanoma - 4.15 g/dl (0.33) [2.5 - 5.1]
# Stage IV Melanoma - 3.92 g/dl (0.45) [2.5 - 5.1]
  mean.ALB <- 3.92
  sd.ALB <- 0.45
  range.ALB <- c(2.5, 5.1)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Determine mean serum creatinine
# CKD-EPI equation used to determine eGFR
  
# eGFR = 141 * min(\frac{[Cr]}{\kappa}, 1)^{\alpha} *
#        max(\frac{[Cr]}{\kappa},1)^{-1.209} * 0.993^{Age} * 
#        1.018^{female} * 1.159^{black} * 1.73/BSA
  
# where $[Cr]$ is the serum creatinine concentration in mg/dL, $\kappa$ is $0.9$ 
# for males and $0.7$ for females, $female$ and $black$ are binary values.
# 
# To simplify the equation the eGFR that correlates with a $[Cr]$ of 0.9 and 0.7 
#   mg/dL was determined for the average individual from the model demographics 
#   (Age = 61.12, black = 0) for both males and females respectively.
# 
# The values were 91.8 mL/min for males and 93.4 mL/min respectively. This is 
#   above the average eGFR from Bajaj et al. of 78.49 mL/min/1.73m2, therefore 
#   $[Cr]$ is more than 0.9 and 0.7 mg/dL for males and females. This simplifies
#   the equation to:

# eGFR = 141 * \bigg(\frac{[Cr]}{\kappa}\bigg)^{-1.209} *
#        0.993^{age} * 1.018^{female} * 1.159^{black}
  
# And therefore you can solve for [Cr]
  
# [Cr] = \kappa * \bigg(\frac{eGFR}
#        {141*0.993^{age}*1.018^{female}*1.159^{black}}\bigg)^{-\frac{1}{1.209}}
  
# This means the average value for $[Cr]$ is 1.024 mg/dL for males and 
#   0.809 mg/dL for females. 
  
# Two distributions of SeCr could be made
  
# While determination of the mean for serum creatinine can be analytical, the
#   standard deviation and range are determined empirically, by altering values
#   for SECR and HT until eGFR matches the distributions from Bajaj et al.

# Define demographic data
# Some objects have been defined above
  mean.SECR.M <- 90.6*0.8
  mean.SECR.F <- 71.5*0.8
  sd.SECR <- 40
  range.SECR <- c(25, 500)
  
  mean.eGFR <- 78.49
  sd.eGFR <- 21.63
  range.eGFR <- c(0, 120)
  
# Age - normal distribution
  AGE <- trunc.rnorm(n = nid, mean = mean.AGE, sd = sd.AGE, range = range.AGE)
  
# Sex - binomial distribution
  SEX <- rbinom(nid, 1, sex.prob)
  
# Serum Creatinine - log-normal distribution
  SECR.M <- trunc.rnorm(n = length(SEX[SEX == 1]), mean = mean.SECR.M, 
    sd = sd.SECR, range = range.SECR, log = T)
  SECR.F <- trunc.rnorm(n = length(SEX[SEX == 0]), mean = mean.SECR.F, 
    sd = sd.SECR, range = range.SECR, log = T)

# Bind into covariate table
  cov.df <- data.frame(AGE, SEX, SECR = 0)
  cov.df$SECR[cov.df$SEX == 1] <- SECR.M
  cov.df$SECR[cov.df$SEX == 0] <- SECR.F
  
# Calculated eGFR
  cov.df$cGFR <- apply(cov.df, 1, function(df) {
    ckdepi.fn(df["SECR"], df["AGE"], df["SEX"], 0)
  })
  
# Randomly Sampled GFR
  cov.df$rGFR <- trunc.rnorm(n = nid, mean = mean.eGFR, sd = sd.eGFR, 
    range = range.eGFR)
  
# Randomly Sampled GFR correlated to AGE
  mean.AGE.GFR <- c(mean.AGE, mean.eGFR)
  sd.AGE.GFR <- c(sd.AGE, sd.eGFR)
  
  corr.AGE.GFR <- matrix(c(
     1.00, -0.34,
    -0.34,  1.00
  ), nrow = 2, ncol = 2)  # symmetrical correlation matrix
  
  lower.AGE.GFR <- c(range.AGE[1], range.eGFR[1])
  upper.AGE.GFR <- c(range.AGE[2], range.eGFR[2])
  
  corrsim <- trunc.mvrnorm(n = nid, mean = mean.AGE.GFR, sd = sd.AGE.GFR, 
    corr_mat = corr.AGE.GFR, lower = lower.AGE.GFR, upper = upper.AGE.GFR, 
    log = F)
  cov.df$mvAGE <- corrsim[, 1]
  cov.df$mvGFR <- corrsim[, 2]
  
  p <- NULL
  p <- ggplot(data = cov.df)
  p <- p + geom_histogram(aes(cGFR, fill = cbPalette$red), 
    bins = 30, alpha = 0.3)
  p <- p + geom_histogram(aes(rGFR, fill = cbPalette$blue), 
    bins = 30, alpha = 0.3)
  p <- p + geom_histogram(aes(mvGFR, fill = cbPalette$green), 
    bins = 30, alpha = 0.3)
  p <- p + scale_fill_identity(name = "Dist", guide = "legend", 
    labels = c("Rand. GFR", "Corr. GFR", "Calc. GFR"))
  p
  
# So what does the relationship between serum creatinine and eGFR look like
  p <- NULL
  p <- ggplot(data = cov.df[cov.df$SEX == 0,])
  p <- p + geom_point(aes(x = SECR, y = cGFR), colour = cbPalette$red, 
    shape = 1, alpha = 0.25)
  p <- p + scale_x_continuous("Serum Creatinine (umol/L)")
  p
  
# Maybe the correlation between required covariates and eGFR can be obtained
#   from this work and used to sample from a large multivariate distribution
  with(cov.df, cor.test(SEX, cGFR))
  with(cov.df, cor.test(mvAGE, mvGFR))
  
  
  