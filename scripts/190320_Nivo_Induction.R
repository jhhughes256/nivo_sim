# Nivolumab PK Model with Tumour Growth - Dosing Regimen 
# -----------------------------------------------------------------------------
# Simulation regimen for comparison to NONMEM model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Clear workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set working directory
# if not working with RStudio project place working directory here
# setwd("C:/.../nivo_sim/")

# Load package libraries
  library(dplyr)  # Split and rearrange data - required for mrgsolve
  library(mrgsolve)  # Metrum differential equation solver for pharmacometrics
  library(ggplot2)  # Graphical package

# Source PopPK model script
  source("scripts/functions_utility.R")  # functions utility
  source("scripts/190320_Nivo_Population.R")  # population data (sources model too)
  regimen.name <- "190320"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Replicate test population for concentration dataset
  
# Duplicate population data for number of times and order data
# purrr::map used to run rep.int on each column
  conc_df <- data.frame(
    ID = rep(ID, each = length(conc_times)),
    time = conc_times,
    amt = 0,  
    cmt = 1, 
    evid = 0, 
    rate = 0
  )

# Replicate test population for dose dataset
# Define dose time, dose and infusion duration
  dose_times <- c(0, 14, 28)
  dose_amt <- 240  # mg
  
# Duplicated population data for number of dose times and order data
  dose_df <- data.frame(
    ID = rep(ID, each = length(dose_times)),
    time = dose_times,
    amt = 240, 
    cmt = 1, 
    evid = 1, 
    rate = -2
  )

# Combine into simulation dataset
  input_df <- dplyr::bind_rows(dose_df, conc_df) %>%
    dplyr::arrange(ID, time, evid) %>%
    dplyr::bind_cols(EPS_df)
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Simulate standard dosing regimen
# Pipe dataset to model
  output_df <- mod %>%
    mrgsolve::data_set(input_df) %>%  # set input data (observations/events)
    mrgsolve::idata_set(pop_df) %>%  # set individual data (sets tumour size)
    mrgsolve::carry.out(amt, evid) %>%  # copy to simulated output
    mrgsolve::mrgsim() %>%  # simulate using mrgsolve
    tibble::as_tibble()
  
# Extract trough data from output
  trough_df <- dplyr::filter(output_df, evid == 0) %>%  # remove dose rows
    dplyr::filter(time %% 14 == 0) %>%  # filter to trough times
    dplyr::select(ID, time, DV, AUC, TUM) %>%  # select important columns
    dplyr::group_by(ID) %>%  # group by ID (initiate ddply)
    dplyr::mutate(Cavg = c(0, diff(AUC))) %>%  # calculate delta AUC (ddply .fun)
    dplyr::ungroup()  # ungroup (collect ddply output)
  
    
    
