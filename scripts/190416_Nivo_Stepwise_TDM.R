# Nivolumab PK Model with Tumour Growth - TDM Step-wise Dosing
# -----------------------------------------------------------------------------
# Simulation of dosing 240 mg every 2 weeks initially, before using TDM with
#   proportional dosage changes.
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
  source("scripts/190320_Nivo_Induction.R")  # induction data
  # sourcing of induction data sources model and population
  # this will be swapped for one script sourcing all scripts for final run
  dose_interval <- 14
  dose_min <- 40
  dose_max <- 800
  dose_opts <- c(40, seq(80, 800, by = 20))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  TDMstep_fn <- function(induction_df) {
  # Set up a loop that will sample the individual's concentration, and determine
  #   the next dose.
  # Make all predicted concentrations and PK parameter values after 
  #   the first sample time equal to NA (aka NA_real_)
    tdmstep_df <- induction_df %>%
      dplyr::select(ID = ID2, dplyr::everything()) %>%
      dplyr::mutate(DV = dplyr::if_else(time > 14, NA_real_, DV))
  # Create tumour patient data for setting initial tumour size and assign 
  #   initial compartment value for model
    tdmstep_mod <- dplyr::summarise_at(tdmstep_df, "TUM", dplyr::first) %>%
      mrgsolve::init(.x = mod, .)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Loop until all doses have been optimised
    trough_target <- 2.5
    trough_upper <- 5
    repeat {
    # Determine latest sample time and dose
      last_sample <- tidyr::drop_na(tdmstep_df) %>%
        dplyr::summarise_at("time", max) %>%
        unlist()
      next_dose <- last_sample + dose_interval
      prev_dose <- last_sample - dose_interval
      last_dose <- dplyr::filter(tdmstep_df, time == prev_dose) %>%
        dplyr::pull(amt)
      this_dose <- dplyr::filter(tdmstep_df, time == last_sample) %>%
          dplyr::pull(amt)
    # Determine last sampled conc and last dose
      last_conc <- dplyr::filter(tdmstep_df, time == last_sample) %>%
        dplyr::pull(DV)
    # Determine next dose
      if (last_conc <= trough_upper & last_conc > trough_target) {
        dose_par <- this_dose
      } else if (last_conc > trough_upper & last_conc <= trough_upper + 5) {
        dose_par <- this_dose - 40
        if (dose_par < dose_min) dose_par <- dose_min
      } else if (last_conc <= trough_target & last_conc > trough_target - 1.25) {
        dose_par <- this_dose + 40
        if (dose_par > dose_max) dose_par <- dose_max
      } else if (last_conc > trough_upper + 5) {
        dose_par <- this_dose - 80
        if (dose_par < dose_min) dose_par <- dose_min
      } else if (last_conc <= trough_target - 1.25) {
        dose_par <- this_dose + 80
        if (dose_par > dose_max) dose_par <- dose_max
      }
    # Create simulation input
      input_tdmstep_df <- tdmstep_df %>%
        dplyr::mutate(amt = dplyr::if_else(
          time == next_dose, dose_par, amt
        ))
    # Simulate to represent time passing since last trough
      tdmstep_df <- tdmstep_mod %>%
        mrgsolve::data_set(data = input_tdmstep_df) %>%
        mrgsolve::idata_set(data = pop_df) %>%
        mrgsolve::carry_out(amt, evid, rate, cmt) %>% 
        mrgsolve::mrgsim() %>%
        tibble::as_tibble() %>% 
        dplyr::select(ID, time, amt, evid, rate, cmt, DV, AUC, TUM,
           AGE, ALB, BWT, GFR, SEX, ECOG, EPS1) %>%  # select important columns
        dplyr::mutate(Cavg = c(0, diff(AUC))) %>%  # calculate delta AUC (ddply .fun)
        dplyr::mutate(DV = dplyr::if_else(time > next_dose, NA_real_, DV))
    # End loop once a year of optimised dosing is complete
      if (last_sample == 350) break
    }  # brackets closing "repeat
    # browser()
    tdmstep_df
  }  # brackets closing "bayes_fn"
  
  tictoc::tic()
  output_tdmstep_df <- trough_flat_df %>%
    dplyr::filter(ID %in% 1:100) %>%
  { tibble::add_column(., ID2 = .$ID) } %>%  # so that ID is carried inside of the nest structure
    dplyr::group_by(ID) %>% tidyr::nest() %>%  # create list column for ID data
    dplyr::mutate(bayes = purrr::map(data, TDMstep_fn))  # create new list column using bayes_fn
  tictoc::toc() 
    
  readr::write_rds(output_tdmstep_df, path = "output/stepwise_tdm.rds")
  