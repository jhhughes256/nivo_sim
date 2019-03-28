# Nivolumab Model and Bayes Validation - Run Validation
# ------------------------------------------------------------------------------
# Control script for validation of mrgsolve model and Bayes estimation 
#   function against NONMEM. Structure of validation is as follows:
#
#  1. Source model
#  2. Create study population
#  3. Create input data set for mrgsolve and NONMEM simulation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare workspace
# Clear workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set working directory
# if not working with RStudio project place working directory here
# setwd("C:/.../nivo_sim/")

# Load package libraries
  library(dplyr)	# dplyr required for mrgsolve
  library(mrgsolve)	  # Metrum differential equation solver for pharmacometrics
  library(ggplot2)  # graphical package
  # library(MASS)  # mvrnorm in trunc_mvrnorm
  # library(MBESS)  # cor2cov in trunc_mvrnorm
  
# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))
  
# Set colourblind palette
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
  myPalette <- with(cbPalette, c(blue, red))
  
# Source external scripts
  script_path <- "scripts/model_bayes_validation/"
  source("scripts/functions_utility.R")  # functions utility

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Validation
# 1. Source model
# New objects: mod
  source("models/NivoPKTS_Final.R")  # PopPK model script 
  
# 2. Create study population
# New objects: pop_df
  source(paste0(script_path, "create_population.R"))

# 3. Create input data set for mrgsolve and NONMEM simulation
# New objects: 
  source(paste0(script_path))