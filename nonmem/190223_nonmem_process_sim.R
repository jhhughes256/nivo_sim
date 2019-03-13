# Process NONMEM simulation .fit file from 190220_nonmem_pop.R dataset
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
  
  match.data <- ddply(simdata[simdata$DOSEGRP != 1,], .(DOSEGRP), function(df) {
    cbind(df, 
      data.frame(
        HAZPLAC = simdata$HAZRATE[simdata$DOSEGRP == 1]
      )
    )
  })
  
  match.data$HAZRAT <- with(match.data, HAZRATE/HAZPLAC)
  
  hazdata <- ddply(match.data, .(DOSEGRP, TIME), function(df) {
    data.frame(
      MED = mean(df$HAZRAT, na.rm = T),
      CI90LO = quantile(df$HAZRAT, prob = 0.05, na.rm = T),
      CI90HI = quantile(df$HAZRAT, prob = 0.95, na.rm = T)
    )
  })
  
  closest <- which(abs(hazdata$MED - 0.4) == min(abs(hazdata$MED - 0.4), na.rm = T))
  hazdata[hazdata$TIME == hazdata$TIME[closest],]
  
  tumdata <- ddply(simdata, .(DOSEGRP, TIME), function(df) {
    data.frame(
      MED = median(df$TUM_SLD, na.rm = T),
      CI90LO = quantile(df$TUM_SLD, prob = 0.05, na.rm = T),
      CI90HI = quantile(df$TUM_SLD, prob = 0.95, na.rm = T)
    )
  })
  
  condata <- ddply(simdata, .(DOSEGRP, TIME), function(df) {
    data.frame(
      MED = median(df$CONC, na.rm = T),
      CI90LO = quantile(df$CONC, prob = 0.05, na.rm = T),
      CI90HI = quantile(df$CONC, prob = 0.95, na.rm = T)
    )
  })
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  p <- NULL
  p <- ggplot(data = hazdata)
  p <- p + geom_line(aes(x = TIME, y = MED), size = 1)
  p <- p + geom_line(aes(x = TIME, y = CI90LO), 
    size = 1, linetype = "dashed")
  p <- p + geom_line(aes(x = TIME, y = CI90HI), 
    size = 1, linetype = "dashed")
  p <- p + facet_wrap(~DOSEGRP)
  p
  
# See if tumour growth goes down
  p <- NULL
  p <- ggplot(data = tumdata)
  p <- p + geom_line(aes(x = TIME, y = MED), size = 1)
  p <- p + geom_line(aes(x = TIME, y = CI90LO), 
    size = 1, linetype = "dashed")
  p <- p + geom_line(aes(x = TIME, y = CI90HI), 
    size = 1, linetype = "dashed")
  p <- p + facet_wrap(~DOSEGRP)
  p
  
  p <- NULL
  p <- ggplot(data = condata[condata$DOSEGRP == 3,])
  p <- p + geom_line(aes(x = TIME, y = MED), size = 1)
  p <- p + geom_line(aes(x = TIME, y = CI90LO), 
    size = 1, linetype = "dashed")
  p <- p + geom_line(aes(x = TIME, y = CI90HI), 
    size = 1, linetype = "dashed")
  p
