# Nivolumab PK Model with Tumour Growth - Dosing Regimen 
# -----------------------------------------------------------------------------
# Simulate representative patient population with test dosing regimen
# Standard dosing - Dose = Target AUC*(CRCL+25) [Target AUC = 6mg.min/mL]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Clear workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set working directory
# if not working with RStudio project place working directory here
# setwd("C:/.../nivo_sim/")

# Load package libraries
  library(plyr)	      # New plyr - required for mrgsolve
  library(dplyr)	    # Split and rearrange data - required for mrgsolve
  library(mrgsolve)	  # Metrum differential equation solver for pharmacometrics
  library(ggplot2)    # Graphical package

# Source PopPK model script
  source("scripts/190130_NivoPKTS_Pop.R")
  regimen.name <- "190130"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Replicate test population for concentration dataset
  TIME.conc <- seq(from = 0, to = 672, by = 1)
  test.conc <- lapply(pop.data, rep.int, times = length(TIME.conc))
  test.conc <- as.data.frame(test.conc)
  test.conc <- test.conc[with(test.conc, order(test.conc$ID)), ]
  test.conc$time <- TIME.conc
  test.conc$amt <- 0
  test.conc$cmt <- 1
  test.conc$evid <- 0
  test.conc$rate <- 0
  head(test.conc)

# Replicate test population for dose dataset
  TIME.dose <- c(0, 336)
  test.dose <- lapply(pop.data, rep.int, times = length(TIME.dose))
  test.dose <- as.data.frame(test.dose)
  test.dose <- test.dose[with(test.dose, order(test.dose$ID)), ]
  test.dose$time <- TIME.dose
  test.dose$amt <- 240
  test.dose$cmt <- 1
  test.dose$evid <- 1
  test.dose$rate <- 0
  head(test.dose)

# Combine into simulation dataset
  test.input <- rbind(test.conc, test.dose)
  test.input <- test.input[with(test.input, 
    order(test.input$ID, test.input$time, test.input$evid)
  ), ]
  head(test.input)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Simulate standard dosing regimen
# Pipe dataset to model
  test.data <- mod %>% 
    data_set(test.input) %>% 
    carry.out(amt, evid) %>% 
    mrgsim()
  test.data <- as.data.frame(test.data)
  
# Check simulated data
  head(test.data)

# Plot simulated data
  p1 <- ggplot(test.data, aes(x = time, y = IPRED))
  p1 <- p1 + geom_line(col = "blue", size = 0.8)
  p1 <- p1 + facet_wrap(~ID)
  p1
  
  p2 <- ggplot(test.data, aes(x = time, y = TUMSLD))
  p2 <- p2 + geom_line(col = "red", size = 0.8)
  p2 <- p2 + facet_wrap(~ID)
  p2

# Save simulated data
  output.name <- paste0("output/NivoPKTS_", regimen.name, "_Data.csv")
  write.csv(test.data, file = output.name)