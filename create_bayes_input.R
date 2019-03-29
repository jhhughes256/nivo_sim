# Nivolumab Model and Bayes Validation - Bayes Input Data
# ------------------------------------------------------------------------------
# Placeholder
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
# Define Bayes input data
# Define input data for Bayes function
  input_rbayes_df <- output_mrgsim_df %>%
    dplyr::filter(time %in% c(dose_times, 364)) %>%  # filter to trough times
    dplyr::select(ID, time, amt, evid, rate, cmt, DV, AGE, ALB, BWT, GFR, SEX, 
      ECOG, RCC, OTHERC, ADApos, ADAunk, SQNSQ) %>%  # select important columns
    tibble::add_column(TUM = rep(pop_df$TUM_0, each = length(sample_times)))

# Define input data for NONMEM Bayes estimation
  tmp_nmbayes_df <- input_rbayes_df %>%
    dplyr::rename(CID = ID, TIME = time, AMT = amt, EVID = evid, RATE = rate, 
      CMT = cmt, ADAPOS = ADApos, ADAUNK = ADAunk, BASESLD = TUM)
  
  input_nmbayes_df <- dplyr::filter(tmp_nmbayes_df, TIME != 0) %>%
    dplyr::mutate(AMT = 0, EVID = 0, RATE = 0) %>%
    dplyr::bind_rows(dplyr::filter(tmp_nmbayes_df, TIME != 364)) %>%
    tibble::add_column(MDV = 0, .after = "DV") %>%
    dplyr::mutate(
      DV = dplyr::if_else(AMT !=0, ".", as.character(DV)),
      MDV = dplyr::if_else(AMT !=0, 1, 0)
    ) %>%
    dplyr::arrange(CID, TIME, dplyr::desc(AMT))
  
  readr::write_csv(input_nmbayes_df, "output/input_nonmem_bayes.csv")
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Define and run R Bayes function
  bayes_fn <- function(input_df) {
    browser()
  # Set up mrgsolve model
    bayes_mod <- dplyr::summarise_at(input_df, "TUM", dplyr::first) %>%
      mrgsolve::init(.x = mod, .)
  # Loop until successful minimisation
    run_once <- FALSE
    repeat {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Bayes estimation
    # Initial estimates for Bayes parameters
      init_par <- exp(double(5))
      if (run_once) {
        init_par <- init_par*exp(runif(5, min = -0.01, max = 0.01))
      }
    # Previous dependent variable values
      prev_DV <- dplyr::filter(input_df, time != 0) %>% 
        dplyr::pull(DV)
    # Define bayesian estimation function
      bayes_estimate_fn <- function(par) {
      # Describe parameters to be optimised within mrgsolve data set
        ETA <- log(par)
      # Define mrgsolve dataset
        output_df <- cbind(input_df, data.frame(
            ID = 1,  # must provide ID
            ETA1 = ETA[1], ETA2 = ETA[2], ETA3 = ETA[3], ETA4 = ETA[4], ETA9 = ETA[5]
          )) %>%
        # Run data through mrgsolve, with idata using initial tumour size
          mrgsolve::data_set(x = bayes_mod, data = .) %>%
          mrgsolve::carry_out(amt, evid, rate, cmt) %>% 
          mrgsolve::mrgsim() %>%
          tibble::as_tibble() %>% dplyr::slice(-1) %>%
        # Ensure IPRED has finite values
          dplyr::mutate(IPRED = dplyr::if_else(
            !is.finite(IPRED) | IPRED < .Machine$double.eps,
            .Machine$double.eps,
            IPRED
          ))
      # Define yhat
        yhat <- dplyr::filter(output_df, time %in% dose_times[-1]) %>%
          dplyr::pull(IPRED)
      # Posterior log-likelihood
      # Error model: Y = IPRE*(1 + EPS*PERR), Y = IPRE + IPRE*EPS*PERR
        loglikpost_sd <- mod$PERR
        loglikpost <- dnorm(prev_DV, mean = yhat, sd = loglikpost_sd, log = T)
      # Prior log-likelihood
        loglikprior <- dnorm(ETA, mean = 0, sd = ETABSV[1:4], log = T)
      # Return objective function value to be minimised
        return(-1*sum(loglikpost, loglikprior))
      }  # end bayes_estimate_fn
    # Run bayes_estimate_fn using optim()
      browser()
      bayes_estimate <- optim(init_par, bayes_estimate_fn, method = "L-BFGS-B",
        lower = rep(0.001, times = 8), upper = rep(1000, times = 8),
        control = list(
          parscale = init_par, fnscale = bayes_estimate_fn(init_par),
          factr = 1e12, pgtol = 1e-8
        )
      )
      run.once <- TRUE
      minimised <- "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
  		if (bayes_estimate$message == minimised) break
    }  # end loop as successful minimisation has occurred
  # Record bayes estimate optim results
    prev_bayes_log <- TRUE
    bayes_estimate_lst <- list(bayes_estimate_lst, bayes_estimate)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate concentrations and population parameters according to Bayes 
  #   estimates for output
  # Convert new ETA values (estimated in the exp() domain)
    input_sim_df <- 
      cbind(input_df, data.frame(
        ID = 1,
        ETA1 = log(bayes_estimate$par[1]),
        ETA2 = log(bayes_estimate$par[2]),
        ETA3 = log(bayes_estimate$par[3]),
        ETA4 = log(bayes_estimate$par[4]),
        ETA9 = log(bayes_estimate$par[5])
      ))
    output_sim_df <- input_sim_df %>%
      mrgsolve::data_set(x = bayes_mod, data = .) %>%
      mrgsolve::carry_out(amt, evid, rate, cmt) %>% 
      mrgsolve::mrgsim() %>%
      tibble::as_tibble() %>% dplyr::slice(-1)
  }
  
  output_rbayes_df <- input_rbayes_df %>%
    dplyr::group_by(ID) %>% tidyr::nest() %>%  # create list column for ID data
    dplyr::mutate(bayes = purrr::map(data, bayes_fn))  # create new list column using bayes_fn