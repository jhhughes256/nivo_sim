# Process NONMEM simulation .fit file for validating NONMEM model
# -----------------------------------------------------------------------------
# Aim is to replicate some of the plots from the paper

# Designed to be used with three separate simulation files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Clear workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set working directory
# if not working with RStudio project place working directory here
# setwd("C:/.../nivo_sim/")

# Load package libraries
  library(plyr)	      # New plyr - required for mrgsolve
  library(ggplot2)    # Graphical package
  
# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))

# Source processSIMdata function
  source("nonmem/processSIMdata.R")
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Simulation data processing
# Set up environment for function
  proj.dir <- getwd()
  master.dir <- paste(getwd(), "output", sep = "/")
  run.name <- paste0("nm_", c("20mgq2w", "1280mgq2w", "placebo"), "_haz")
  
# Process data (only needs to be run once as it saves a .csv)
  # llply(run.name, function(x) {
  #   processSIMdata(paste0(x, ".ctl"))
  # })
  setwd(proj.dir)
  
# Read in simdata
  simdata.list <- llply(run.name, function(x) {
    read.csv(paste0("output/", x, ".nm7/", x, ".fit.csv"),
      stringsAsFactors = F, na.strings = ".")
  })
  
# Correct DOSEGRP column
  simdata.list[[1]]$DOSEGRP <- 2
  simdata.list[[2]]$DOSEGRP <- 3
  
# Bind together
  simdata <- do.call(rbind, simdata.list)
  
# ID represents Dose group
# SIM represents patient ID
# ID == 1 - placebo group
# ID == 2 - 20mg q2w
# ID == 3 - 1280mg q2w
# Want to match the patients from the dose groups with their respective
#   placebo individual
  match.data <- ddply(simdata[simdata$DOSEGRP != 1,], .(DOSEGRP), function(df) {
    cbind(df, 
      data.frame(
        HAZPLAC = simdata$HAZRATE[simdata$DOSEGRP == 1]
      )
    )
  })
  
# Calculate the hazard ratio
  match.data$HAZRAT <- with(match.data, HAZRATE/HAZPLAC)
  
# Calculate summary statistics
  hazdata <- ddply(match.data, .(DOSEGRP, TIME), function(df) {
    data.frame(
      MEAN = mean(df$HAZRAT, na.rm = T),
      MED = median(df$HAZRAT, na.rm = T),
      CI90LO = quantile(df$HAZRAT, prob = 0.05, na.rm = T),
      CI90HI = quantile(df$HAZRAT, prob = 0.95, na.rm = T)
    )
  })
  
# Find out what time the hazard ratio is at 0.4
  closest <- which(abs(hazdata$MEAN - 0.4) == min(abs(hazdata$MEAN - 0.4), na.rm = T))
  hazdata[hazdata$TIME == hazdata$TIME[closest],]
  
# And there it is!!! The time at which the hazard ratio was calculated was
#   16 weeks!
  
# Save data as .RDS for nm_model_validation.rmd
  saveRDS(simdata, "rds/190301_nonmem_process_sim.rds")
  