# Nivolumab PK Model with Tumour Growth - Validating Bayes Function
# -----------------------------------------------------------------------------
# Validation of different models as simulated in 190403_Bayes_Testing.R
# Files are:
#   Full - All ETAs (no IIVHAZ) [8 ETAs]
#   Base - All "influential" ETAs (no EC50, R or IIVHAZ) [6 ETAs]
#   NoTG - Base with no TG ETA [5 ETAs] - done
#   NoEMAX - Base with no EMAX ETA [5 ETAs]
#   NoTMAX - Base with no TMAX ETA [5 ETAs]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Clear workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set working directory
  master_dir <- "C:/Users/hugjh001/Documents/nivo_sim/"

# Load package libraries
  library(magrittr)
  # library(dplyr)  # Split and rearrange data - required for mrgsolve
  # library(mrgsolve)  # Metrum differential equation solver for pharmacometrics
  library(ggplot2)  # Graphical package
  library(cowplot)
  
# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))
  
# Read in simulation data
  output_mrgsim_lst <- readr::read_csv(
    paste0(master_dir, "output/output_mrgsim.csv"), 
    col_types = readr::cols()) %>% 
    list() %>% rep(times = 5)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load Bayes estimation output
# Define run suffixes, ETA lists and ETA names
  run_suff <- c("full", "base", "notg", "noemax", "notmax")
  run_etas <- list(
    full = paste0("ETA", c(1:7, 9)),
    base = paste0("ETA", c(1:4, 6, 9)),
    notg = paste0("ETA", c(1:4, 9)),
    noemax = paste0("ETA", c(1:3, 6, 9)),
    notmax = paste0("ETA", c(1:4, 6))
  )
  eta_names <- paste0("ETA_", 
    c("CL", "V1", "V2", "EMAX", "EC50", "TG", "R", "IIVHAZ", "TMAX")
  )
  eta_bsv <- c(
    ETA1 = 0.096721, ETA2 = 0.099225, ETA3 = 0.185761, 
    ETA4 = 1.000000, ETA5 = 0.010000, ETA6 = 0.100000, 
    ETA7 = 0.500000, ETA8 = 0.100000, ETA9 = 0.044521
  )
  mod_sig <- 0.194
  
# R Bayes output (.rds files)
# Define run names
  rbrun_name <- paste0("output_rbayes_", run_suff)
  
# Read in data and bind to NONMEM output
  output_rbayes_tbl <- rbrun_name %>%
    purrr::map(function(file_name) {
      readr::read_rds(paste0(master_dir, "output/", file_name, ".rds"))
    }) %>%
    tibble::tibble(rbdata = .) %>%
    tibble::add_column(refdata = output_mrgsim_lst)
  
# NONMEM output (.fit files)
# Define run names
  nmrun_name <- paste0("NivoPKTS_NM_Bayes_", run_suff)
  file_name <- paste0("output/", nmrun_name, ".nm7/", nmrun_name)
  
# Process data (only needs to be run once as it saves a .csv)
  # file_name %>%
  #   purrr::map(function(file_name) {
  #     readr::read_table(file = paste0(file_name, ".fit"), 
  #       skip = 1, col_names = T, col_type = readr::cols(),
  #       na = c("NA","***********","1.#INFE+00")
  #     ) %>%
  #     readr::write_csv(path = paste0(file_name, ".fit.csv"))
  #   })
  
# Read in data
  output_bayes_tbl <- file_name %>%
    purrr::map(function(file_name) {
      readr::read_csv(paste0(master_dir, file_name, ".fit.csv"), col_types = readr::cols())
    }) %>%
    tibble::tibble(nmdata = .) %>%
    tibble::add_column(run = run_suff, .before = "nmdata") %>%
    dplyr::bind_cols(output_rbayes_tbl)
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Determine R Bayes absolute error in relation to NONMEM
  val_bayes_lst <- output_bayes_tbl %>%
    purrr::pmap(function(run, nmdata, rbdata, refdata) {
    # Prepare environment
      sample_times <- unique(rbdata$time)
      etas <- run_etas[[run]]
    # Isolate ETAs and DVs from each dataset
    # dim(x_eta) == c(2600, length(etas))
      eta_tbl <- tibble::tibble(ID = unique(rbdata$ID),
        rb_eta = rbdata %>%
          dplyr::filter(time == 14) %>%
          dplyr::select(ID, etas) %>%
          dplyr::group_by(ID) %>% tidyr::nest() %>% dplyr::pull(data),
        nm_eta = nmdata %>%
          dplyr::filter(AMT == 0 & TIME == 14) %>%
          dplyr::select(ID, paste0("ETA", 1:length(etas))) %>%
          dplyr::rename_at(paste0("ETA", 1:length(etas)), ~ etas) %>%
          dplyr::group_by(ID) %>% tidyr::nest() %>% dplyr::pull(data),
        rf_eta = refdata %>%
          dplyr::filter(time == 14) %>%
          dplyr::select(ID, etas) %>%
          dplyr::group_by(ID) %>% tidyr::nest() %>% dplyr::pull(data),
        rb_ipred = rbdata %>%
          dplyr::filter(time %in% sample_times) %>%
          dplyr::select(ID, IPRED) %>%
          dplyr::group_by(ID) %>% tidyr::nest() %>% dplyr::pull(data),
        nm_ipred = nmdata %>%
          dplyr::filter(AMT == 0 & TIME %in% sample_times) %>%
          dplyr::select(ID, IPRED) %>%
          dplyr::group_by(ID) %>% tidyr::nest() %>% dplyr::pull(data),
        rf_ipred = refdata %>%
          dplyr::filter(time %in% sample_times) %>%
          dplyr::select(ID, IPRED) %>%
          dplyr::group_by(ID) %>% tidyr::nest() %>% dplyr::pull(data),
        rf_dv = refdata %>%
          dplyr::filter(time %in% sample_times) %>%
          dplyr::select(ID, DV) %>%
          dplyr::group_by(ID) %>% tidyr::nest() %>% dplyr::pull(data)
      )
    # Calculate objective function
      obj_tbl <- dplyr::mutate(eta_tbl,
        rb_obj = eta_tbl %>%
          dplyr::select(ID, rf_dv, rb_ipred, rb_eta) %>%
          dplyr::mutate_at(c("rf_dv", "rb_ipred", "rb_eta"), 
            ~ purrr::map(., unlist, use.names = F)
          ) %>%
          purrr::pmap(function(ID, rf_dv, rb_ipred, rb_eta) {
            loglikpost <- dnorm(rf_dv, mean = rb_ipred, sd = mod_sig, log = T)
            loglikprior <- dnorm(rb_eta, mean = 0, sd = eta_bsv[etas], log = T)
            -1*sum(loglikpost, loglikprior)
          }) %>% unlist(),
        nm_obj = eta_tbl %>%
          dplyr::select(ID, rf_dv, nm_ipred, nm_eta) %>%
          dplyr::mutate_at(c("rf_dv", "nm_ipred", "nm_eta"), 
            ~ purrr::map(., unlist, use.names = F)
          ) %>%
          purrr::pmap(function(ID, rf_dv, nm_ipred, nm_eta) {
            loglikpost <- dnorm(rf_dv, mean = nm_ipred, sd = mod_sig, log = T)
            loglikprior <- dnorm(nm_eta, mean = 0, sd = eta_bsv[etas], log = T)
            -1*sum(loglikpost, loglikprior)
          }) %>% unlist(),
        rf_obj = eta_tbl %>%
          dplyr::select(ID, rf_dv, rf_ipred, rf_eta) %>%
          dplyr::mutate_at(c("rf_dv", "rf_ipred", "rf_eta"), 
            ~ purrr::map(., unlist, use.names = F)
          ) %>%
          purrr::pmap(function(ID, rf_dv, rf_ipred, rf_eta) {
            loglikpost <- dnorm(rf_dv, mean = rf_ipred, sd = mod_sig, log = T)
            loglikprior <- dnorm(rf_eta, mean = 0, sd = eta_bsv[etas], log = T)
            -1*sum(loglikpost, loglikprior)
          }) %>% unlist()
      )
      # Determine error
      out_tbl <- dplyr::mutate(obj_tbl,
          rb_err = rb_eta %>%
            purrr::map2(rf_eta, `-`) %>% 
            purrr::map(abs),
          nm_err = nm_eta %>%
            purrr::map2(rf_eta, `-`) %>% 
            purrr::map(abs),
          rb_obj_err = abs(rb_obj - rf_obj)/rf_obj,
          nm_obj_err = abs(nm_obj - rf_obj)/rf_obj,
          OBJ = rb_obj_err - nm_obj_err,
          rel_err = purrr::map2(rb_err, nm_err, `-`)
        ) %>%
        dplyr::select(ID, rel_err, OBJ)
    })

# Create plot data and plots using ggplot2
  val_plot_lst <- purrr::map(val_bayes_lst, function(eta_tbl) {
    plot_df <- eta_tbl %>%
      tidyr::unnest() %>%
      dplyr::select(-ID) %>%
      dplyr::summarise_all(list(
        med = ~ median(.), 
        lo = ~ quantile(., prob = 0.25), 
        hi = ~ quantile(., prob = 0.75)
      )) %>%
      tidyr::gather() %>%
      tidyr::separate(col = key, into = c("ETA", "metric"), sep = "_") %>%
      tidyr::spread(key = metric, value = value) %>%
      dplyr::mutate(
        ETA = c(eta_names[as.numeric(substr(head(ETA, -1), 4, 4))], "OBJ")
      )
    
    p <- NULL
    p <- ggplot(aes(x = ETA, colour = ETA), data = plot_df)
    p <- p + geom_hline(yintercept = 0, linetype = "dashed")
    p <- p + geom_point(aes(y = med), size = 2)
    p <- p + geom_errorbar(aes(ymin = lo, ymax = hi), size = 1)
  	p <- p + labs(x = NULL, y = "Absolute Error of R with Respect to NONMEM")
  	p <- p + coord_cartesian(ylim = c(-0.25, 0.75))
  	p <- p + theme(legend.position = "none")
  	p
  })

  