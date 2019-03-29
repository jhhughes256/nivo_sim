# Nivolumab Model and Bayes Validation - Run Validation
# ------------------------------------------------------------------------------
# Control script for validation of mrgsolve model and Bayes estimation 
#   function against NONMEM. Structure of validation is as follows:
#
# Validation of mrgsolve model
#  1. Source model
#  2. Create study population
#  3. Create input data set for mrgsolve and NONMEM simulation
#  4. Run NONMEM simulation externally
#  5. Process mrgsolve and NONMEM output and validate mrgsolve model
#
# Validation of Bayes function
#  1. Use model, population and outupt from mrgsolve validation
#  2. Create input data set for mrgsolve and NONMEM Bayes estimation and run R Bayes
#  3. Run R Bayes function
#  4. Run NONMEM Bayes estimation externally
#  5. Process mrgsolve and NONMEM output and validate Bayes function

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare workspace
# Clear workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set working directory
# if not working with RStudio project then define working directory
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
  
# Source external scripts
  script_path <- "scripts/model_bayes_validation/"
  source("scripts/functions_utility.R")  # functions utility

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Validation of mrgsolve model
# 1. Source model
# New objects: mod
  source("models/NivoPKTS_Final.R")  # PopPK model script 
  
# 2. Create study population
# New objects: pop_df
  source(paste0(script_path, "create_population.R"))

# 3. Create input data set for mrgsolve and NONMEM simulation
# New objects: input_mrgsim_df input_nonmem_df output_mrgsim_df
  source(paste0(script_path, "create_sim_input.R"))
  
# 4. Run NONMEM simulation externally
  
# 5. Process mrgsolve and NONMEM output and validate mrgsolve model
# New objects: output_nonmem_df p_CvT p_TvT
  source(paste0(script_path, "process_sim_output.R"))
  
# View plots
# Plot of concentration versus time for mrgsolve (orange) & NONMEM (blue)
  p_CvT  
  
# Plot of tumour size versus time for mrgsolve (orange) & NONMEM (blue)
  p_TvT  

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Validation of Bayes function
# 1. Use model, population and outupt from mrgsolve validation
# Reused objects: mod pop_df output_mrgsim_df
  
# 2. Create input data set for R and NONMEM Bayes estimation and run R Bayes
# New objects: input_rbayes_df input_nmbayes_df
  source(paste0(script_path, "create_bayes_input.R"))
  
# 4. Run NONMEM Bayes estimation externally
  
# 5. Process mrgsolve and NONMEM output and validate Bayes function
# New objects:
  source(paste0(script_path, "process_bayes_output.R"))
  