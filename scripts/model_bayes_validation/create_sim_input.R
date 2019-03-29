# Nivolumab Model and Bayes Validation - Simulation Input Data
# ------------------------------------------------------------------------------
# Create input for both mrgsolve and NONMEM models for the validation of 
#   mrgsolve model and create simulation data from mrgsolve for the validation
#   of the Bayes estimation function against NONMEM
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare workspace
# Clear workspace
  # rm(list=ls(all=TRUE))
  # graphics.off()

# Set working directory
# if not working with RStudio project place working directory here
# setwd("C:/.../nivo_sim/")

# # Load package libraries
  # library(dplyr)	# dplyr required for mrgsolve
  # library(mrgsolve)  # Metrum differential equation solver for pharmacometrics
  # library(MASS)  # mvrnorm in trunc_mvrnorm
  # library(MBESS)  # cor2cov in trunc_mvrnorm

# Source external scripts
  # script_path <- "scripts/model_bayes_validation/"
  # source("scripts/functions_utility.R")  # functions utility
  # source("models/NivoPKTS_Final.R")  # PopPK model script 
  # source(paste0(script_path, "create_population.R"))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Define simulation inputs for mrgsolve
# Define sample times
  conc_times <- seq(from = 0, 364, by = 2)
  sample_times <- seq(from = 0, 364, by = 14)
  dose_times <- seq(from = 0, 350, by = 14)
  
# Residual unexplained variability
  input_dim <- nid*length(conc_times)
  EPS_df <- mrgsolve::smat(mod) %>%  # Omega block values from model
    as.matrix() %>%  # convert to matrix
    diag() %>%  # extract diagonal elements
    purrr::map_dfc(function(Z) rnorm(input_dim, mean = 0, sd = sqrt(Z)))  # allocate
  names(EPS_df) <- "EPS1"

# Create mrgsolve input data set
  input_mrgsim_df <- data.frame(
      ID = rep(ID, each = length(conc_times)),
      time = rep(conc_times, times = length(ID))
    ) %>% tibble::as_tibble() %>%
    dplyr::mutate(
      amt = if_else(time %in% dose_times, 240, 0),
      cmt = 1,
      evid = if_else(amt != 0, 1, 0),
      rate = if_else(amt != 0, -2, 0)
    )  %>%
    dplyr::bind_cols(EPS_df)
  
# Run mrgsolve
  output_mrgsim_df <- mod %>%
    mrgsolve::data_set(input_mrgsim_df) %>%  # set input data (observations/events)
    mrgsolve::idata_set(pop_df) %>%  # set individual data (sets tumour size)
    mrgsolve::carry_out(amt, evid, rate, cmt) %>%  # copy to simulated output
    mrgsolve::mrgsim() %>%  # simulate using mrgsolve
    tibble::as_tibble()
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Alter mrgsolve input to create NONMEM input data set
# Convert column names and add in population data
  temp_nonmem_df <- input_mrgsim_df %>%
    dplyr::group_by(ID) %>% tidyr::nest() %>%  # nest_by
    dplyr::bind_cols(pop_df) %>%
    tidyr::unnest() %>% 
    dplyr::select(ID, time, amt, cmt, evid, rate, dplyr::everything(), -ID1) %>%
    dplyr::rename(TIME = time, AMT = amt, CMT = cmt, EVID = evid, RATE = rate,
      PPVCL = ETA1, PPVV1 = ETA2, PPVV2 = ETA3, PPVEMAX = ETA4, PPVEC50 = ETA5, 
      PPVTG = ETA6, PPVR = ETA7, PPVHZ = ETA8, PPVTMAX = ETA9, ERRPRO = EPS1, 
      BASESLD = TUM_0, ADAPOS = ADApos, ADAUNK = ADAunk)
  
# Create dose data set
  dose_nonmem_df <- filter(temp_nonmem_df, TIME %in% dose_times)
  
# Create concentration data set and bind dose data set to it
  input_nonmem_df <- temp_nonmem_df %>%
    dplyr::mutate(
      AMT = 0,
      EVID = 0,
      RATE = 0
    ) %>%
    dplyr::bind_rows(dose_nonmem_df) %>%
    dplyr::arrange(ID, TIME, dplyr::desc(AMT)) %>%
    dplyr::rename(CID = ID) %>%
    dplyr::mutate(DV = ".") %>%
    tibble::add_column(DUM1 = 0, .after = "PPVCL") %>%
    tibble::add_column(DUM2 = 0, .after = "PPVEC50") %>%
    tibble::add_column(DUM3 = 0, .after = "PPVTMAX")
    
  readr::write_csv(input_nonmem_df, "output/input_nonmem_sim.csv")
  
# Tidy workspace
  rm(temp_nonmem_df)
  rm(dose_nonmem_df)
  