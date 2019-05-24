# Nivolumab PK Model with Tumour Growth - Flat Dosing Regimen/Induction
# -----------------------------------------------------------------------------
# Simulation of dosing 240 mg every 2 weeks
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

# Read in data
  pop_df <- readr::read_rds("pop_df.rds")
  flat_df <- readr::read_rds("flat_dosing.rds")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create simulation input dataset
# Define individuals
  nid <- 10000  # Number of individuals
  ID <- 1:nid  # Sequence of individual ID's

# Define time points
  conc_times <- seq(from = 0, 364, by = 14)  # 1 year of half daily data
  dose_times <- 0:25*14  # 26 doses separated by 14 days

# Duplicate population data for number of times and alter dose times to include
#   dosage events for mrgsolve
  input_flat_df <- data.frame(
    ID = rep(ID, each = length(conc_times)),
    time = conc_times) %>%
    dplyr::mutate(
      amt = dplyr::if_else(time %in% dose_times, 120, 0),
      cmt = 1,
      evid = dplyr::if_else(amt != 0, 1, 0),
      rate = dplyr::if_else(amt != 0, -2, 0)
    )  %>%
    dplyr::bind_cols(flat_df["EPS1"])

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulate standard dosing regimen
# Pipe dataset to model
  output_flat_df <- mod %>%
    mrgsolve::data_set(input_flat_df) %>%  # set input data (observations/events)
    mrgsolve::idata_set(pop_df) %>%  # set individual data (sets tumour size)
    mrgsolve::carry_out(amt, evid, rate, cmt) %>%  # copy to simulated output
    mrgsolve::mrgsim() %>%  # simulate using mrgsolve
    tibble::as_tibble()

# Extract trough data from output
  trough_flat_df <- output_flat_df %>%
    dplyr::filter(time %in% c(dose_times, 364)) %>%  # filter to trough times
    dplyr::select(ID, time, amt, evid, rate, cmt, DV, AUC, TUM,
      AGE, ALB, BWT, GFR, SEX, ECOG, EPS1) %>%  # select important columns
    dplyr::group_by(ID) %>%  # group by ID (initiate ddply)
    dplyr::mutate(Cavg = c(0, diff(AUC))) %>%  # calculate delta AUC (ddply .fun)
    dplyr::ungroup()  # ungroup (collect ddply output)

# Save data as .RDS for trough_simulation.rmd
  # readr::write_rds(pop_df, "pop_df.rds")
  readr::write_rds(trough_flat_df, "flat_dosing_120.rds")
