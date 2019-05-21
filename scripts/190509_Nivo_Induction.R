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
  # library(dplyr)  # Split and rearrange data - required for mrgsolve
  # library(mrgsolve)  # Metrum differential equation solver for pharmacometrics
  # library(ggplot2)  # Graphical package

# Source PopPK model script
  source("scripts/functions_utility.R")  # functions utility
  source("scripts/190509_Nivo_Population.R")  # population data (sources model too)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Create simulation input dataset
# Define time points
  conc_times <- seq(from = 0, 364, by = 0.5)  # 1 year of half daily data
  dose_times <- 0:25*14  # 26 doses separated by 14 days
  
# Residual unexplained variability
  input_dim <- nid*length(conc_times)
  EPS_df <- mrgsolve::smat(mod) %>%  # Omega block values from model
    as.matrix() %>%  # convert to matrix
    diag() %>%  # extract diagonal elements
    purrr::map_dfc(function(Z) rnorm(input_dim, mean = 0, sd = sqrt(Z)))  # allocate
  names(EPS_df) <- "EPS1"
  
# Duplicate population data for number of times and alter dose times to include
#   dosage events for mrgsolve
  input_flat_df <- data.frame(
    ID = rep(ID, each = length(conc_times)),
    time = conc_times) %>%
    dplyr::mutate(
      amt = if_else(time %in% dose_times, 240, 0),
      cmt = 1,
      evid = if_else(amt != 0, 1, 0),
      rate = if_else(amt != 0, -2, 0)
    )  %>%
    dplyr::bind_cols(EPS_df)
  
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
  # saveRDS(trough_flat_df, "rds/190326_NivoPKTS_Reg_240_PKPDpop.rds")
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot patient data
# Set ggplot2 theme
#   theme_bw2 <- theme_set(theme_bw(base_size = 14))
#   theme_update(plot.title = element_text(hjust = 0.5))
# 
# # Plot individual patient concentrations
#   p <- NULL
#   p <- ggplot(data = output_flat_df)
#   p <- p + geom_line(aes(x = time, y = IPRED), colour = "red", size = 1)
#   p <- p + geom_point(aes(x = time, y = DV), colour = "blue", shape = 1, alpha = 0.5)
#   p <- p + labs(x = "Time (days)", y = "Concentration (mg/mL)")
#   p <- p + coord_cartesian(xlim = c(84, 168), ylim = NULL)
#   p <- p + facet_wrap(~ID)
#   p
#     
# # Plot individual patient tumour size
#   p <- NULL
#   p <- ggplot(data = output_flat_df)
#   p <- p + geom_line(aes(x = time, y = TUM), colour = "red", size = 1)
#   p <- p + labs(x = "Time (days)", y = "Tumour Size (mm)")
#   p <- p + coord_cartesian(xlim = NULL, ylim = NULL)
#   p <- p + facet_wrap(~ID)
#   p