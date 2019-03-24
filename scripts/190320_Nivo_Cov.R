# Sample Covariate Values from Random Distributions for Nivo Population
# -----------------------------------------------------------------------------
# Based on the work done in ******_Demog_Sampling.R and nivo_population.Rmd
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare work environment
# Clear workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Load libraries
  library(MASS)  # mvrnorm in trunc_mvrnorm
  library(MBESS)  # cor2cov in trunc_mvrnorm
  library(magrittr)  # %>%
  library(ggplot2)  # standard plotting functions
  library(cowplot)  # plot_grid
  
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

# Load functions utility - trunc_rnorm, trunc_mvrnorm, ckdepi_fn
  source("scripts/functions_utility.R")
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Define demographic data
  set.seed(123456) 
  
# Age - 61.12 years (11.12) [23 - 87]
# Weight - 79.09 kg (19.28) [34 - 168]
# Albumin - 3.92 g/dl (0.45) [2.5 - 5.1]
# Sex - 66.7% male, 33.3% female
# ECOG - 38.73% 0, 61.26% >0
# Number of individuals
  nid <- 10000

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
  
# Experimental Values (give mean eGFR when using mean covariate values)
# SeCr (males) - 90.6 umol/L
# SeCr (females) - 71.5 umol/L
  mean_SECR_M <- 90.6*1.04
  mean_SECR_F <- 71.5*1.04
  sd_SECR <- 28
  range_SECR_M <- c(0, 300)
  range_SECR_F <- c(0, 300)
  range_SECR <- c(range_SECR_F[1], range_SECR_M[2])
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Sample covariate values from distributions
# Age, Albumin and initial tumour size - normal distribution
  AGE <- trunc_rnorm(n = nid, mean = mean_AGE, sd = sd_AGE, range = range_AGE)
  ALB <- trunc_rnorm(n = nid, mean = mean_ALB, sd = sd_ALB, range = range_ALB)
  TUM_0 <- rnorm(n = nid, mean = 54.6, sd = 10.92)
  
# Sex & ECOG - binomial distribution
  SEX <- rbinom(nid, 1, sex_prob)
  ECOG <- rbinom(nid, 1, ecog_prob)
 
# Weight & Serum Creatinine - log-normal distribution
  BWT <- trunc_rnorm(n = nid, mean = mean_BWT, sd = sd_BWT, 
    range = range_BWT, log = T)
  SECR_M <- trunc_rnorm(n = length(SEX[SEX == 1]), mean = mean_SECR_M, 
    sd = sd_SECR, range = range_SECR_M, log = T)
  SECR_F <- trunc_rnorm(n = length(SEX[SEX == 0]), mean = mean_SECR_F, 
    sd = sd_SECR, range = range_SECR_F, log = T)
  
# Bind covariate data.frame
  cov_df <- data.frame(AGE, ALB, BWT, SEX, ECOG, TUM_0)
  cov_df$SECR <- 0
  cov_df$SECR[cov_df$SEX == 1] <- SECR_M
  cov_df$SECR[cov_df$SEX == 0] <- SECR_F
 
# Calculated eGFR
  cov_df$GFR <- dplyr::select(cov_df, SECR, AGE, SEX) %>%
    purrr::pmap(function(SECR, AGE, SEX) ckdepi_fn(SECR, AGE, SEX, 0)) %>%
    unlist()
  
# Save cov_df as output
  readr::write_csv(x = cov_df, path = "output/Nivo_Cov.csv")
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Covariate Population Diagnostics
# Check mean and standard deviation of GFR
# Bajaj et al. eGFR - 78.49 (21.63)
  with(cov_df, c(mean(GFR), sd(GFR)))
  
  with(cov_df, c(
    "Normal" = length(GFR[GFR > 90])*100/nid,
    "Mild" = length(GFR[GFR > 60 & GFR <= 90])*100/nid,
    "Moderate" = length(GFR[GFR > 30 & GFR <= 60])*100/nid,
    "Severe" = length(GFR[GFR <= 30])*100/nid
  ))
  
# Plot distributions
# Add SEXf column
  cov_df$SEXf <- factor(cov_df$SEX)
  levels(cov_df$SEXf) <- c("Female", "Male")
  
  p <- NULL
  p <- ggplot(data = cov_df)
  plot_AGE <- p + geom_histogram(aes(AGE), bins = 50, fill = cbPalette$green)
  plot_AGE <- plot_AGE + geom_vline(xintercept = range_AGE, linetype = "dashed")
  
  plot_ALB <- p + geom_histogram(aes(ALB), bins = 50, fill = cbPalette$orange)
  plot_ALB <- plot_ALB + geom_vline(xintercept = range_ALB, linetype = "dashed")
  
  p_M <- NULL
  p_M <- ggplot(data = cov_df[cov_df$SEX == 1, ])
  plot_SECR_M <- p_M + geom_histogram(aes(SECR), bins = 50, 
    fill = cbPalette$blue)
  plot_SECR_M <- plot_SECR_M + geom_vline(xintercept = range_SECR_M, 
    linetype = "dashed")
  plot_SECR_M <- plot_SECR_M + xlab("SECR (male)")
  
  p_F <- NULL
  p_F <- ggplot(data = cov_df[cov_df$SEX == 0, ])
  plot_SECR_F <- p_F + geom_histogram(aes(SECR), bins = 50, 
    fill = cbPalette$pink)
  plot_SECR_F <- plot_SECR_F + geom_vline(xintercept = range_SECR_F, 
    linetype = "dashed")
  plot_SECR_F <- plot_SECR_F + xlab("SECR (female)")
  
  plot_BWT <- p + geom_histogram(aes(BWT), bins = 50, fill = cbPalette$skyblue)
  plot_BWT <- plot_BWT + geom_vline(xintercept = range_BWT, linetype = "dashed")
  
  plot_GFR <- p + geom_histogram(aes(GFR), fill = cbPalette$red, bins = 50)
  
  plot_hist <- plot_grid(ncol = 2,
    plot_AGE, plot_ALB, plot_SECR_M, plot_SECR_F,
    plot_BWT, plot_GFR)
  plot_hist
  
# Plot relationship between age/secr and GFR
  plot_AGE_GFR <- p + geom_point(aes(x = AGE, y = GFR), 
    colour = cbPalette$red, shape = 1, alpha = 0.3)
  plot_AGE_GFR <- plot_AGE_GFR + xlab("Age (years)")
  plot_AGE_GFR <- plot_AGE_GFR + ylab("eGFR (mL/min/1.73m2)")
  plot_AGE_GFR
  
  plot_SECR_GFR <- p + geom_point(aes(x = SECR, y = GFR, colour = SEXf), 
    shape = 1, alpha = 0.2)
  plot_SECR_GFR <- plot_SECR_GFR + xlab("Serum Creatinine (umol/L)")
  plot_SECR_GFR <- plot_SECR_GFR + ylab("eGFR (mL/min/1.73m2)")
  plot_SECR_GFR <- plot_SECR_GFR + guides(
    colour = guide_legend(override.aes = list(alpha=1))
  )
  plot_SECR_GFR