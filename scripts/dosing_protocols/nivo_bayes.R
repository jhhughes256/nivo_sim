# Nivolumab PK Model with Tumour Growth - Bayes Dosing
# -----------------------------------------------------------------------------
# Simulation of dosing 240 mg every 2 weeks initially, before using Bayes to
#   optimise dosing.
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
  dose_max <- 1200
  bayes_etas <- list(
    full = c(1:7, 9),
    base = c(1:4, 6, 9),
    notg = c(1:4, 9),
    noem = c(1:3, 6, 9),
    notm = c(1:4, 6),
    pktg = c(1:3, 6),
    pkem = 1:4,
    pk = 1:3
  )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  bayes_fn <- function(induction_df) {
  # Source model so that it is compiled on each node
    source("models/NivoPKTS_Final.R")
  # Set up a loop that will sample the individual's concentration, estimate 
  #   empirical Bayes parameters, optimise their dose and administer until 
  #   time == 350 days
  # Make all predicted concentrations and PK parameter values after 
  #   the first sample time equal to NA (aka NA_real_)
    bayes_df <- induction_df %>%
      dplyr::mutate(DV = dplyr::if_else(time > 14, NA_real_, DV))
  # Create tumour patient data for setting initial tumour size and assign 
  #   initial compartment value for model
    bayes_mod <- dplyr::summarise_at(bayes_df, "TUM", dplyr::first) %>%
      mrgsolve::init(.x = mod, .)
	# Define a variable telling the loop if previous bayes results are  
	#   present. If they are present, they will be used as initial estimates  
  #   for the next dosing interval.
    prev_bayes_log <- FALSE
    bayes_estimate_lst <- list("Initial" = NA)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Loop until all doses have been optimised
    # try(
    repeat {
    # Estimate individual parameter values using:
    #   1. Trough sample from the end of the previous intervals
  	#   2. Covariate values measured at beginning and end of previous 
    #      intervals
  	#   3. Known doses that were administered during the previous intervals
    # Determine latest sample time
      last_sample <- tidyr::drop_na(bayes_df) %>%
        dplyr::summarise_at("time", max) %>%
        unlist()
    # Filter bayes data frame for data up until the latest sample time,
    #   then add in mrgsim columns
      input_estim_df <- bayes_df %>%
        dplyr::select(-EPS1) %>%
        dplyr::filter(time <= last_sample)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Bayesian Estimation
    # Loop until successful minimisation
      if (last_sample >= 112) { n_mod <- 1 }
      else if (last_sample == 14) { n_mod <- 7 }
      else { n_mod <- 1 }
      run_num <- 0
      repeat {
      # Initial estimates for Bayes parameters
        mod_etas <- bayes_etas[[n_mod]]
        n_eta <- length(mod_etas)
        init_par <- exp(double(n_eta))
        if (run_num > 0) {
          init_par <- init_par*exp(runif(n_eta, min = -0.01, max = 0.01))
        }
      # Previous dependent variable values
        prev_DV <- dplyr::filter(input_estim_df, time != 0) %>% 
          dplyr::pull(DV)
      # Define bayesian estimation function
        bayes_estimate_fn <- function(par) {
        # Describe parameters to be optimised within mrgsolve data set
          ETA <- log(par)
        # Define mrgsolve dataset
          output_estim_df <- cbind(input_estim_df, data.frame(ID = 1)) %>%  # must provide ID
            dplyr::mutate_at(paste0("ETA", mod_etas), function(x) {
              eta <- as.numeric(substr(deparse(substitute(x)), 4, 4))
              which_eta <- which(mod_etas == eta)
              ETA[which_eta]
            }) %>%
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
          yhat <- dplyr::filter(output_estim_df, time %in% dose_times[-1]) %>%
            dplyr::pull(IPRED)
        # Posterior log-likelihood
        # Error model: Y = IPRE*(1 + EPS*PERR), Y = IPRE + IPRE*EPS*PERR
          loglikpost_sd <- mod$PERR
          loglikpost <- dnorm(prev_DV, mean = yhat, sd = loglikpost_sd, log = T)
        # Prior log-likelihood
          loglikprior <- dnorm(ETA, mean = 0, sd = ETABSV[mod_etas], log = T)
        # Return objective function value to be minimised
          return(-1*sum(loglikpost, loglikprior))
        }  # end bayes_estimate_fn
      # Run bayes_estimate_fn using optim()
        bayes_estimate <- optim(init_par, bayes_estimate_fn, method = "L-BFGS-B",
          lower = rep(0.001, times = n_eta), upper = rep(1000, times = n_eta),
          control = list(
            parscale = init_par, fnscale = bayes_estimate_fn(init_par),
            factr = 1e12, pgtol = 1e-8
          )
        )
        
        run_num <- run_num + 1
        minimised <- "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
        if (bayes_estimate$message == minimised) break
        if (run_num >= n_mod) n_mod <- n_mod + 1
        
      }  # end loop as successful minimisation has occurred
    # Record bayes estimate optim results
      prev_bayes_log <- TRUE
      bayes_estimate_lst[[which(dose_times == last_sample)]] <- bayes_estimate
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Calculate concentrations according to new Bayes estimates for time 
    #   until dose intervention can occur
    # Determine next dose and trough that can be affected by dose optimisation
      next_dose <- last_sample + dose_interval
      next_trough <- last_sample + dose_interval*2
    # Convert new ETA values (estimated in the exp() domain)
      input_sim_df <- dplyr::filter(bayes_df, time <= next_trough) %>%
        dplyr::select(-EPS1) %>% cbind(data.frame(ID = 1)) %>%
        dplyr::mutate_at(paste0("ETA", mod_etas), function(x) {
          eta <- as.numeric(substr(deparse(substitute(x)), 4, 4))
          which_eta <- which(mod_etas == eta)
          log(bayes_estimate$par[which_eta])
        })
      output_sim_df <- input_sim_df %>%
        mrgsolve::data_set(x = bayes_mod, data = .) %>%
        mrgsolve::carry_out(amt, evid, rate, cmt) %>% 
        mrgsolve::mrgsim() %>%
        tibble::as_tibble() %>% dplyr::slice(-1)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Optimise dose for the individual using individual Bayes predicted 
    #   concentrations (and compartment amounts) at time of last sampling.
    #   Only optimise if the last trough concentration from sim_bayes_df is
    #   out of target range.
    # Determine last sampled conc
      last_bayes_conc <- dplyr::filter(output_sim_df, time == next_trough) %>%
        dplyr::pull(IPRED)
    # Define trough target and upper bound, begin if statement
      trough_target <- 2.5
      trough_upper <- 5
      if (last_bayes_conc < trough_target | last_bayes_conc >= trough_upper) {
      # Modify model code ready for simulation
        optim_mod <- mrgsolve::init(bayes_mod, list(
          CMT1 = dplyr::filter(output_sim_df, time == last_sample) %>% 
            dplyr::pull(CMT1),
          CMT2 = dplyr::filter(output_sim_df, time == last_sample) %>% 
            dplyr::pull(CMT2),
          TUM = dplyr::filter(output_sim_df, time == last_sample) %>% 
            dplyr::pull(TUM),
          SRV = dplyr::filter(output_sim_df, time == last_sample) %>% 
            dplyr::pull(SRV),
          DRP = dplyr::filter(output_sim_df, time == last_sample) %>% 
            dplyr::pull(DRP),
          AUC = dplyr::filter(output_sim_df, time == last_sample) %>% 
            dplyr::pull(AUC)
        ))
      # Set initial dose and error estimates
        init_par <- c(240, 0.01)
      # Subset input dataset so only future concentrations are predicted
        input_optim_df <- dplyr::filter(input_sim_df, 
          time >= last_sample & time <= next_trough
        )
      # Find the doses that maximise the likelihood of trough concentrations
      #   being the target concentration
        optimise_dose_fn <- function(par) {
        # Add fitted parameters to the input data set, then simulate 
        #   concentration-time profiles with fitted doses
          output_optim_df <- input_optim_df %>%
            dplyr::mutate(amt = dplyr::if_else(
              time == next_dose, par[1], amt
            )) %>%
            mrgsolve::data_set(x = optim_mod, data = .) %>%
            mrgsolve::carry_out(amt, evid, rate, cmt) %>% 
            mrgsolve::mrgsim(start = last_sample, end = last_sample) %>%
            tibble::as_tibble() %>% dplyr::slice(-1) %>%
            dplyr::mutate(IPRED = dplyr::if_else(
              !is.finite(IPRED) | IPRED < .Machine$double.eps,
              .Machine$double.eps,
              IPRED
            ))
        # Define yhat and the residual
          yhat <- dplyr::filter(output_optim_df, time == next_trough) %>%
            dplyr::pull(IPRED)
          res <- dnorm(trough_target, yhat, yhat*par[2], log = T)
        # Objective function value to be minimised
          return(-1*sum(res))
        }
        optimise_dose <- optim(init_par, optimise_dose_fn, method = "L-BFGS-B",
          lower = c(dose_min, 0.0001),
          upper = c(dose_max, Inf),
          control = list(parscale = init_par, factr = 1e7)
        )
      # Administer the individual the optimised dose
        input_bayes_df <- bayes_df %>%
          dplyr::mutate(ID = 1, amt = dplyr::if_else(
          time == next_dose, optimise_dose$par[1], amt
        ))
      } else {
      # Give previous dose
        prev_dose <- dplyr::filter(input_sim_df, time == last_sample) %>%
          dplyr::pull(amt)
        input_bayes_df <- bayes_df %>%
          dplyr::mutate(ID = 1, amt = dplyr::if_else(
          time == next_dose, prev_dose, amt
        ))
      }  # end if at target trough
    # Simulate to represent time passing since last trough
      bayes_df <- bayes_mod %>%
        mrgsolve::data_set(data = input_bayes_df) %>%
        mrgsolve::idata_set(data = pop_df) %>%
        mrgsolve::carry_out(amt, evid, rate, cmt) %>% 
        mrgsolve::mrgsim() %>%
        tibble::as_tibble() %>% 
        dplyr::select(time, amt, evid, rate, cmt, DV, AUC, TUM,
           AGE, ALB, BWT, GFR, SEX, ECOG, EPS1) %>%  # select important columns
        dplyr::mutate(Cavg = c(0, diff(AUC))) %>%  # calculate delta AUC (ddply .fun)
        dplyr::mutate(DV = dplyr::if_else(time > next_dose, NA_real_, DV))
    # End loop once a year of optimised dosing is complete
      if (next_trough == 364) break
    }  # brackets closing "repeat"
    # )
    bayes_df
  }  # brackets closing "bayes_fn"
  
  future::plan(future::multiprocess, workers = 10)
 
  tictoc::tic()
  output_bayes_df <- trough_flat_df %>%
    dplyr::group_by(ID) %>% tidyr::nest() %>%  # create list column for ID data
    # dplyr::mutate(bayes = purrr::map(data, bayes_fn))  # create new list column using bayes_fn
    dplyr::mutate(bayes = furrr::future_map(data, bayes_fn, .progress = T))  # create new list column using bayes_fn
  tictoc::toc()  
  