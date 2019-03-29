# Nivolumab Model and Bayes Validation - Process Simulation Data
# ------------------------------------------------------------------------------
# Process and compare output from mrgsolve and NONMEM models
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
  # library(ggplot2)  # graphical package

# Source external scripts
  source("nonmem/processSIMdata.R")
  # script_path <- "scripts/model_bayes_validation/"
  # source("scripts/functions_utility.R")  # functions utility
  # source("models/NivoPKTS_Final.R")  # PopPK model script 
  # source(paste0(script_path, "create_population.R"))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# NONMEM Simulation data processing
# Set up environment for function
  proj.dir <- getwd()
  master.dir <- paste(getwd(), "output", sep = "/")
  run.name <- "NivoPKTS_NM_Sim"
  
# Process data (only needs to be run once as it saves a .csv)
  # processSIMdata(paste0(run.name, ".ctl"))
  # setwd(proj.dir)
  
# Read in .csv file
  output_nonmem_df <- read.csv(
      paste0("output/", run.name, ".nm7/", run.name, ".fit.csv"),
      stringsAsFactors = F, na.strings = "."
    ) %>% tibble::as_tibble() %>%
    filter(AMT == 0)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Validation plots
# Concentration-time profile
  p <- NULL
  p <- ggplot()
  p <- p + stat_summary(aes(x = TIME, y = DV), data = output_nonmem_df,
    geom = "ribbon", fun.ymin = CI95lo, fun.ymax = CI95hi, 
    fill = cbPalette$blue, alpha = 0.5, size = 1)
  p <- p + stat_summary(aes(x = time, y = DV), data = output_mrgsim_df,
    geom = "line", fun.y = CI95lo, colour = cbPalette$orange, linetype = "dashed",
    size = 1)
  p <- p + stat_summary(aes(x = time, y = DV), data = output_mrgsim_df,
    geom = "line", fun.y = CI95hi, colour = cbPalette$orange, linetype = "dashed",
    size = 1)
  p <- p + stat_summary(aes(x = TIME, y = DV), data = output_nonmem_df,
    geom = "line", fun.y = median, colour = cbPalette$blue, size = 1)
  p <- p + stat_summary(aes(x = time, y = DV), data = output_mrgsim_df,
    geom = "line", fun.y = median, colour = cbPalette$orange, linetype = "dashed",
    size = 1)
  p <- p + scale_x_continuous("Time (days)", breaks = 0:8*56)
  p <- p + ylab("Concentration (mg/mL)")
  # p <- p + coord_cartesian(xlim = c(0, 224))
  # p <- p + scale_y_log10()
  p_CvT <- p
  
# Tumour-time profile
  p <- NULL
  p <- ggplot()
  p <- p + stat_summary(aes(x = TIME, y = TUM), data = output_nonmem_df,
    geom = "ribbon", fun.ymin = CI95lo, fun.ymax = CI95hi, 
    fill = cbPalette$blue, alpha = 0.5, size = 1)
  p <- p + stat_summary(aes(x = time, y = TUM), data = output_mrgsim_df,
    geom = "line", fun.y = CI95lo, colour = cbPalette$orange, linetype = "dashed",
    size = 1)
  p <- p + stat_summary(aes(x = time, y = TUM), data = output_mrgsim_df,
    geom = "line", fun.y = CI95hi, colour = cbPalette$orange, linetype = "dashed",
    size = 1)
  p <- p + stat_summary(aes(x = TIME, y = TUM), data = output_nonmem_df,
    geom = "line", fun.y = median, colour = cbPalette$blue, size = 1)
  p <- p + stat_summary(aes(x = time, y = TUM), data = output_mrgsim_df,
    geom = "line", fun.y = median, colour = cbPalette$orange, linetype = "dashed",
    size = 1)
  p <- p + scale_x_continuous("Time (days)", breaks = 0:8*56)
  p <- p + ylab("Concentration (mg/mL)")
  # p <- p + coord_cartesian(xlim = c(0, 224))
  # p <- p + scale_y_log10()
  p_TvT <- p
