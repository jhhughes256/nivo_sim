# Nivolumab PK Model with Tumour Growth - Dosing Regimen 
# -----------------------------------------------------------------------------
# Simulation regimen for comparison to NONMEM model
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
  source("scripts/190302_NivoPKPD_Pop.R")
  regimen.name <- "190302"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Replicate test population for concentration dataset
  TIME.conc <- c(seq(from = 0, to = 13, by = 1), seq(from = 14, to = 735, by = 7))
  test.conc <- lapply(pop.data, rep.int, times = length(TIME.conc))
  test.conc <- as.data.frame(test.conc)
  test.conc <- test.conc[with(test.conc, order(test.conc$ID)), ]
  test.conc$time <- TIME.conc
  test.conc$amt <- 0
  test.conc$cmt <- 1
  test.conc$evid <- 0
  test.conc$rate <- 0
  test.conc$ii <- 0
  test.conc$addl <- 0
  head(test.conc)

# Replicate test population for dose dataset
  TIME.dose <- 0
  test.dose <- lapply(pop.data, rep.int, times = length(TIME.dose))
  test.dose <- as.data.frame(test.dose)
  test.dose <- test.dose[with(test.dose, order(test.dose$ID)), ]
  test.dose$time <- TIME.dose
  test.dose$amt <- 0
  test.dose$cmt <- 1
  test.dose$evid <- 1
  test.dose$rate <- 0
  test.dose$ii <- 14  # dose every 14 days
  test.dose$addl <- 51  # give 51 additional doses
  head(test.dose)

# Create list of simulation datasets for three populations
  dose.pops <- c(0, 10, 20, 40, 80, 160, 320, 640, 1280, 2560)
  data.in.list <- llply(dose.pops, function(dose) {
    if (dose != 0) { 
      pop.dose <- test.dose
      pop.dose$amt <- dose
      test.input <- rbind(test.conc, pop.dose)
    } else {
      test.input <- test.conc
    }
    test.input <- test.input[with(test.input, 
      order(test.input$ID, test.input$time, test.input$evid)
    ), ]
  })

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Simulate standard dosing regimen
# Pipe dataset to model
  test.data.list <- llply(1:length(dose.pops), function(i) {
    out <- mod %>% 
    data_set(data.in.list[[i]]) %>% 
    carry.out(amt, evid) %>% 
    mrgsim() %>%
    as.data.frame()
    out$DOSEGRP <- i
    out
  })
  
# Bind data together
  mrgdata <- do.call(rbind, test.data.list)
  mrgdata <- mrgdata[mrgdata$evid == 0,]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# ID represents Dose group
# SIM represents patient ID
# ID == 1 - placebo group
# ID == 2 - 20mg q2w
# ID == 3 - 1280mg q2w
# Want to match the patients from the dose groups with their respective
#   placebo individual
  match.data <- ddply(mrgdata[mrgdata$DOSEGRP != 1,], .(DOSEGRP), function(df) {
    cbind(df, 
      data.frame(
        HAZPLAC = mrgdata$HAZRATE[mrgdata$DOSEGRP == 1]
      )
    )
  })

# Calculate the hazard ratio
  match.data$HAZRAT <- with(match.data, HAZRATE/HAZPLAC)
  
# Calculate summary statistics
  hazdata <- ddply(match.data, .(DOSEGRP, time), function(df) {
    data.frame(
      MEAN = mean(df$HAZRAT, na.rm = T),
      MED = median(df$HAZRAT, na.rm = T),
      CI90LO = quantile(df$HAZRAT, prob = 0.05, na.rm = T),
      CI90HI = quantile(df$HAZRAT, prob = 0.95, na.rm = T)
    )
  })
  
# Find out what time the hazard ratio is at 0.4
  deltaMEAN <- abs(hazdata$MEAN - 0.4)
  min.deltaMEAN <- min(deltaMEAN, na.rm = T)
  closest <- which(deltaMEAN == min.deltaMEAN)
  hazdata[hazdata$time == hazdata$time[closest],]
  
# Save data as .RDS for nm_model_validation.rmd
  saveRDS(mrgdata, "rds/190302_NivoPKPD_Reg.rds")
  