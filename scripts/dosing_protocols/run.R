# Nivolumab PKPD Model - Base Run Script
# ------------------------------------------------------------------------------
# Script that runs other scripts for use with server. The objective of this
#   simulation study is to determine whether there is better dosing protocols
#   that could be used for nivolumab. In doing so we believe that money can be
#   saved!
#
# Via simulation, this script compares the following protocols
#   1. Flat dosing
#   2. Model-based dosing (Bayes)
#   3. Therapeutic drug monitoring
#   4. ???
#   5. Profit!
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Prepare workspace
# Clear workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set working directory
# if not working with RStudio project place working directory here
  setwd("E:/Hughes/Git/nivo_sim/scripts/dosing_protocols")

# Load package libraries
# Commented packages not required for scripts, written to document used packages
  library(dplyr)	# mutate, select, also required for mrgsolve
  library(mrgsolve)	  # Metrum differential equation solver for pharmacometrics
  # library(tidyr)  # gather, spread
  # library(readr)  # speedy reading and writing of files
  # library(purrr)  # functional programming
  # library(furrr)  # parallel processing with purrr
  # library(tibble)  # alternative data.frame data objects
  # library(ggplot2)  # graphical
  # library(MASS)  # mvrnorm in trunc_mvrnorm from functions_utility
  # library(MBESS)  # cor2cov in trunc_mvrnorm from functions_utility

# Source external scripts
  source("functions_utility.R")  # functions utility
  source("model.R")  # PopPK model script

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Begin simulation study
# Set seed for random number generator
  set.seed(123456)

# Number of individuals
  nid <- 10000  # Number of individuals
  ID <- 1:nid  # Sequence of individual ID's

# Generate Population - new object: pop_df
# Script contains information on population as well as number of simulated???
#   individuals. Seed should be set before simulating altered population.
  flag <- 0
  source("population.R")

# Run flat dosing simulation - new objects: output_flat_df & trough_flat_df
# The initial doses are used as induction doses for the model based dosing script
  source("flat_dosing.R")

# Run Bayes algorithm for model based dosing - new objects: output_bayes_df
# Script is set up to split the individuals into groups to allow parallel
#   processing of patients.
  # source("model_based.R")
