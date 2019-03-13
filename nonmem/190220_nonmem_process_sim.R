# Process NONMEM simulation .fit file from 190220_nonmem_pop.R dataset
# -----------------------------------------------------------------------------
# Aim is to replicate some of the plots from the paper

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
  run.name <- "nivo_base"
  
# Process data (only needs to be run once as it saves a .csv)
  # processSIMdata(paste0(run.name, ".ctl"))
  setwd(proj.dir)
  
# Read in simdata
  simdata <- read.csv(paste0("output/", run.name, ".nm7/", run.name, ".fit.csv"),
    stringsAsFactors = F, na.strings = ".")
  
# ID represents Dose group
# SIM represents patient ID
# ID == 1 - placebo group
# ID == 2 - 20mg q2w
# ID == 3 - 1280mg q2w
  
  match.data <- ddply(simdata[simdata$ID != 1,], .(ID), function(df) {
    cbind(df, 
      data.frame(
        HAZS = simdata$SURS[simdata$ID == 1]
      )
    )
  })
  
  tumdata <- ddply(simdata, .(ID, TIME), function(df) {
    data.frame(
      MED = median(df$TUM_SLD, na.rm = T),
      CI90LO = quantile(df$TUM_SLD, prob = 0.05),
      CI90HI = quantile(df$TUM_SLD, prob = 0.95)
    )
  })
  
  condata <- ddply(simdata, .(ID, TIME), function(df) {
    data.frame(
      MED = median(df$CONC, na.rm = T),
      CI90LO = quantile(df$CONC, prob = 0.05),
      CI90HI = quantile(df$CONC, prob = 0.95)
    )
  })
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# See if tumour growth goes down
  p <- NULL
  p <- ggplot(data = tumdata)
  p <- p + geom_line(aes(x = TIME, y = MED), size = 1)
  p <- p + geom_line(aes(x = TIME, y = CI90LO), 
    size = 1, linetype = "dashed")
  p <- p + geom_line(aes(x = TIME, y = CI90HI), 
    size = 1, linetype = "dashed")
  p <- p + facet_wrap(~ID)
  p
  
  p <- NULL
  p <- ggplot(data = condata[condata$ID == 3,])
  p <- p + geom_line(aes(x = TIME, y = MED), size = 1)
  p <- p + geom_line(aes(x = TIME, y = CI90LO), 
    size = 1, linetype = "dashed")
  p <- p + geom_line(aes(x = TIME, y = CI90HI), 
    size = 1, linetype = "dashed")
  p
