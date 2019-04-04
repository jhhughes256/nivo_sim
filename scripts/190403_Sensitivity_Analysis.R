# Nivolumab PK Model with Tumour Growth - Flat Dosing Regimen/Induction
# -----------------------------------------------------------------------------
# Simulation of dosing 240 mg every 2 weeks
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Clear workspace
  # rm(list=ls(all=TRUE))
  # graphics.off()

# Set working directory
# if not working with RStudio project place working directory here
# setwd("C:/.../nivo_sim/")

# Load package libraries
  library(dplyr)  # Split and rearrange data - required for mrgsolve
  library(mrgsolve)  # Metrum differential equation solver for pharmacometrics
  # library(ggplot2)  # Graphical package

# # Source PopPK model script
#   source("scripts/functions_utility.R")  # functions utility
#   source("scripts/190320_Nivo_Population.R")  # population data (sources model too)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Create simulation input dataset
# Set all ETAs to zero
  pop_df <- filter(pop_df, ID == 1) %>%
    dplyr::mutate(
      ETA1 = 0, ETA2 = 0, ETA3 = 0, ETA4 = 0, ETA5 = 0, 
      ETA6 = 0, ETA7 = 0, ETA8 = 0, ETA9 = 0
    )
  
# Define time points
  conc_times <- seq(from = 0, 364, by = 14)  # 1 year of half daily data
  dose_times <- 0:25*14  # 26 doses separated by 14 days
  
# Residual unexplained variability
  input_dim <- length(conc_times)
  EPS_df <- mrgsolve::smat(mod) %>%  # Omega block values from model
    as.matrix() %>%  # convert to matrix
    diag() %>%  # extract diagonal elements
    purrr::map_dfc(function(Z) rnorm(input_dim, mean = 0, sd = sqrt(Z)))  # allocate
  names(EPS_df) <- "EPS1"
  
# Duplicate population data for number of times and alter dose times to include
#   dosage events for mrgsolve
  input_flat_df <- data.frame(
    ID = rep(1, each = length(conc_times)),
    time = conc_times) %>%
    dplyr::mutate(
      amt = if_else(time %in% dose_times, 240, 0),
      cmt = 1,
      evid = if_else(amt != 0, 1, 0),
      rate = if_else(amt != 0, -2, 0)
    )  %>%
    dplyr::bind_cols(EPS_df)
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Create sensitivity analysis options
# ETA ID 1:9; ETA value -1.2 - 1.2; times 14, 112, 364; 
  input_sa_lst <- expand.grid(1:9, seq(-1.2, 1.2, by = 0.4))
  
# Simulate and extract troughs
  output_sa_df <- input_sa_lst %>%
    purrr::pmap_dfr(function(Var1, Var2) {
      newpop_df <- pop_df
      newpop_df[paste0("ETA", Var1)] <- Var2
      output_flat_df <- mod %>%
        mrgsolve::data_set(input_flat_df) %>%  # set input data (observations/events)
        mrgsolve::idata_set(newpop_df) %>%  # set individual data (sets tumour size)
        mrgsolve::carry_out(amt, evid, rate, cmt) %>%  # copy to simulated output
        mrgsolve::mrgsim() %>%  # simulate using mrgsolve
        tibble::as_tibble()
      tibble::tibble(
        ETA = Var1,
        Val = Var2,
        t14 = filter(output_flat_df, time == 14)$DV,
        t112 = filter(output_flat_df, time == 112)$DV,
        t364 = filter(output_flat_df, time == 364)$DV
      )
    }) %>%  # for each ETA determine the sensitivity
    dplyr::group_by(ETA) %>% tidyr::nest() %>%  # nest_by
    dplyr::mutate(data = purrr::map(data, function(df) {  # tidyverse_ddply
    # Determine DV at time == 0
      DV_0 <- dplyr::filter(df, Val == .Machine$double.eps) %>%
        dplyr::select(-Val)
      
      df %>%
        dplyr::filter(Val != .Machine$double.eps) %>%
        dplyr::mutate( 
          d14 = t14 - DV_0$t14, 
          d112 = t112 - DV_0$t112, 
          d364 = t364 - DV_0$t364,
          s14 = d14/Val, 
          s112 = d112/Val, 
          s364 = d364/Val
        ) %>%
        dplyr::select(s14, s112, s364) %>%
        dplyr::summarise_all(mean) %>%
        dplyr::mutate_at(c("s14", "s112", "s364"), signif, digits = 3)
      
    })) %>%
    tidyr::unnest()
  