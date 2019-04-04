# Nivolumab Model and Bayes Validation - Process Bayes Data
# ------------------------------------------------------------------------------
# Process and compare output from R bayes function and NONMEM MAP estimates
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
# NONMEM .fit file reading
# Set up environment
  run_name <- "NivoPKTS_NM_Bayes"
  file_name <- paste0("output/", run_name, ".nm7/", run_name)
  
# Process data (only needs to be run once as it saves a .csv)
  nmbayes_fit <- read.table(file = paste0(file_name, ".fit"), 
    sep = "", skip = 1, header = T, 
    na.strings = c("NA","***********","1.#INFE+00")
  )
  write.csv(nmbayes_fit, file = paste0(file_name, ".fit.csv"))

# Read in .csv file
	output_nmbayes_df <- read.csv(paste0(file_name, ".fit.csv")) %>%
	  tibble::as_tibble() %>%
	  dplyr::filter(AMT == 0) %>%
	  dplyr::select(-X)
	
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Validation plots
# Create information dense validation dataset
	val_bayes_df <- output_mrgsim_df %>%
	  dplyr::filter(time %in% tail(sample_times, -1)) %>%
	  dplyr::select(ID, time, ETA1, ETA2, ETA3, ETA4, ETA9) %>%
	  dplyr::rename(  # rename columns
	    ETA_CL = ETA1, ETA_V1 = ETA2, ETA_V2 = ETA3, 
	    ETA_EMAX = ETA4, ETA_TMAX = ETA9
	  ) %>%
	  tibble::add_column(  # add output columns to dataset
	    ETA_CLr = output_rbayes_df$ETA1, ETA_CLn = output_nmbayes_df$ETA1,
	    ETA_V1r = output_rbayes_df$ETA2, ETA_V1n = output_nmbayes_df$ETA2,
	    ETA_V2r = output_rbayes_df$ETA3, ETA_V2n = output_nmbayes_df$ETA3,
	    ETA_EMAXr = output_rbayes_df$ETA4, ETA_EMAXn = output_nmbayes_df$ETA4,
	    ETA_TMAXr = output_rbayes_df$ETA9, ETA_TMAXn = output_nmbayes_df$ETA5
	  ) %>%
	  dplyr::mutate(  # determine absolute difference between method and reference
	    ETA_dCLr = abs(ETA_CLr - ETA_CL), ETA_dCLn = abs(ETA_CLn - ETA_CL), 
	    ETA_dV1r = abs(ETA_V1r - ETA_V1), ETA_dV1n = abs(ETA_V1n - ETA_V1), 
	    ETA_dV2r = abs(ETA_V2r - ETA_V2), ETA_dV2n = abs(ETA_V2n - ETA_V2), 
	    ETA_dEMAXr = abs(ETA_EMAXr - ETA_EMAX), 
	    ETA_dEMAXn = abs(ETA_EMAXn - ETA_EMAX), 
	    ETA_dTMAXr = abs(ETA_TMAXr - ETA_TMAX), 
	    ETA_dTMAXn = abs(ETA_TMAXn - ETA_TMAX)
	  ) %>%
	  dplyr::mutate(  # determine absolute difference in method error
	    ETA_CL_err = ETA_dCLr - ETA_dCLn, ETA_V1_err = ETA_dV1r - ETA_dV1n, 
	    ETA_V2_err = ETA_dV2r - ETA_dV2n, ETA_EMAX_err = ETA_dEMAXr - ETA_dEMAXn, 
	    ETA_TMAX_err = ETA_dTMAXr - ETA_dTMAXn 
	  )
	
# Simplify the validation dataset for plotting
	plot_bayes_df <- val_bayes_df %>% 
	  dplyr::select(
	    ID, ETA_CL_err, ETA_V1_err, ETA_V2_err, ETA_EMAX_err, ETA_TMAX_err
	  ) %>%
	  dplyr::group_by(ID) %>% tidyr::nest() %>%
	  dplyr::mutate(data = purrr::map(data, function(df) {
	    dplyr::summarise_all(df, mean)
	  })) %>%
	  tidyr::unnest() %>%
	  tidyr::gather(key = "Parameter", value = "dPAR",
	    ETA_CL_err, ETA_V1_err, ETA_V2_err, ETA_EMAX_err, ETA_TMAX_err
	  ) %>%
	  dplyr::group_by(Parameter) %>% tidyr::nest() %>%  # nest_by
	  dplyr::mutate(data = purrr::map(data, function(df) {
	    tibble::tibble(
	      medPAR = median(df$dPAR),
	      hiPAR = quantile(df$dPAR, prob = 0.75),
	      loPAR = quantile(df$dPAR, prob = 0.25)
	    )
	  })) %>%
	  tidyr::unnest()
	
# Box plot
	p <- NULL
	p <- ggplot(aes(x = Parameter, colour = Parameter), data = plot_bayes_df)
	p <- p + geom_point(aes(y = medPAR), 
	  size = 2)
	p <- p + geom_errorbar(aes(ymin = loPAR, ymax = hiPAR),
	  size = 1)
	p <- p + labs(x = NULL, y = "Absolute Error of R with Respect to NONMEM")
	p <- p + theme(legend.position = "none")
	# p <- p + coord_trans(y = "log10")
	p
	