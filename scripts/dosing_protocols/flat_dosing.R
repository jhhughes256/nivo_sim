# Nivolumab PK Model with Tumour Growth - Flat Dosing Regimen/Induction
# -----------------------------------------------------------------------------
# Simulation of dosing 240 mg every 2 weeks
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
      amt = dplyr::if_else(time %in% dose_times, 240, 0),
      cmt = 1,
      evid = dplyr::if_else(amt != 0, 1, 0),
      rate = dplyr::if_else(amt != 0, -2, 0)
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

# Check for subjects with NaN simulated values, these will be resampled
  repeat {
    if (any(!is.finite(output_flat_df$DV))) {
    # Identify bad subjects
      temp_pop <- pop_df
      ID <- output_flat_df %>%
        dplyr::slice(which(!is.finite(output_flat_df$DV))) %>%
        dplyr::pull(ID) %>%
        unique()
      nid <- length(ID)

    # Create subjects to replace bad subjects
      flag <- 1
      source("population.R")
      browser()

    # Simulate replacements
      input_dim <- nid*length(conc_times)
      EPS_df <- mrgsolve::smat(mod) %>%  # Omega block values from model
        as.matrix() %>%  # convert to matrix
        diag() %>%  # extract diagonal elements
        purrr::map_dfc(function(Z) rnorm(input_dim, mean = 0, sd = sqrt(Z)))  # allocate
      names(EPS_df) <- "EPS1"

    # Create input dataset
      input_flat_df <- data.frame(
        ID = rep(ID, each = length(conc_times)),
        time = conc_times) %>%
        dplyr::mutate(
          amt = dplyr::if_else(time %in% dose_times, 240, 0),
          cmt = 1,
          evid = dplyr::if_else(amt != 0, 1, 0),
          rate = dplyr::if_else(amt != 0, -2, 0)
        )  %>%
        dplyr::bind_cols(EPS_df)

    # Simulate concentrations
      temp_flat_df <- mod %>%
        mrgsolve::data_set(input_flat_df) %>%  # set input data (observations/events)
        mrgsolve::idata_set(pop_df) %>%  # set individual data (sets tumour size)
        mrgsolve::carry_out(amt, evid, rate, cmt) %>%  # copy to simulated output
        mrgsolve::mrgsim() %>%  # simulate using mrgsolve
        tibble::as_tibble()

    # Add replacements to final dataset
      id <- ID
      output_flat_df <- output_flat_df %>%
        dplyr::filter(!(ID %in% id)) %>%
        dplyr::bind_rows(temp_flat_df) %>%
        dplyr::arrange(ID, time, amt)

      pop_df <- temp_pop %>%
        dplyr::filter(!(ID %in% id)) %>%
        dplyr::bind_rows(pop_df) %>%
        dplyr::arrange(ID)
    } else {
      break
    }
  }

# Extract trough data from output
  trough_flat_df <- output_flat_df %>%
    dplyr::filter(time %in% c(dose_times, 364)) %>%  # filter to trough times
    dplyr::select(ID, time, amt, evid, rate, cmt, DV, AUC, TUM,
      AGE, ALB, BWT, GFR, SEX, ECOG, EPS1) %>%  # select important columns
    dplyr::group_by(ID) %>%  # group by ID (initiate ddply)
    dplyr::mutate(Cavg = c(0, diff(AUC))) %>%  # calculate delta AUC (ddply .fun)
    dplyr::ungroup()  # ungroup (collect ddply output)

# Save data as .RDS for trough_simulation.rmd
  readr::write_rds(pop_df, "pop_df.rds")
  readr::write_rds(trough_flat_df, "flat_dosing.rds")
