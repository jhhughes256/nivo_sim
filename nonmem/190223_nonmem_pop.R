# Recreating the NONMEM populations from Nivolumab PK/PD Model Paper
# -----------------------------------------------------------------------------
# Script that produces NONMEM-compatible populations that match those used for
#   the simulations of the following paper:
#     C Liu, J Yu, H Li et al. (2017) Association of Time‐Varying Clearance of 
#     Nivolumab With Disease Dynamics and Its Implications on Exposure Response 
#     Analysis. Clin. Pharmacol. Ther., 101: 657-666. doi:10.1002/cpt.656

# Saves the population as three separate .csv's so that each simulated patient
#   is the same between dosages

# Column names for simulation input are:
#   ID TIME DV BASESLD ECOG AMT II ADDL MDV EVID FLAG DOSEGRP
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Clear workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set working directory
# if not working with RStudio project place working directory here
# setwd("C:/.../nivo_sim/")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Create simulation dataset
# Set dose groups
  dose <- c(0, 20, 1280)
  
# Set ID numbers (these will correlate with dose groups)
  id.seq <- c(1)
  
# Set simulation times
  times <- c(0:13, 2:105*7)
  
# Create simulation dataset for each dose
  simprep.list <- llply(dose, function(dose) {
  # Create amt and dv datasets and merge together
    simprep_amt <- data.frame(
      ID = id.seq,
      BASESLD = 54.6,
      ECOG = 0,
      AMT = dose,
      FLAG = 0,
      DOSEGRP = 1
    )
    simprep_dv <- data.frame(
      ID = rep(id.seq, each = length(times)),
      TIME = times,
      DV = 1,
      II = 0,
      ADDL = 0,
      MDV = 0,
      EVID = 0
    )
    simprep <- merge(simprep_dv, simprep_amt)
    
  # Set dosing frequency and how many additional doses to be given
    simprep$II[simprep$TIME == 0] <- 14
    simprep$ADDL[simprep$TIME == 0] <- 51
  
  # Fix up columns for NONMEM
    simprep$EVID[simprep$TIME == 0] <- 1
    simprep$DV[simprep$TIME == 0] <- "."
    simprep$MDV[simprep$TIME == 0] <- 1
    simprep$AMT[simprep$TIME != 0] <- "."
    simprep$II[simprep$TIME != 0] <- "."
    simprep$ADDL[simprep$TIME != 0] <- "."
    if (dose == 0) { simprep[1, ] <- c(1, 0, 1, ".", ".", 0, 0, 54.6, 0, ".", 0, 1) }
    
  # Rearrange columns
    simprep <- simprep[, c("ID", "TIME", "DV", "BASESLD", "ECOG", "AMT", "II", 
      "ADDL", "MDV", "EVID", "FLAG", "DOSEGRP")]
    names(simprep)[1] <- "#ID"
    
  # Request simprep as output
    simprep
  })
  
# Save to file
  filename_out <- "output/nm_pop_placebo.csv"
  write.csv(simprep.list[[1]], file = filename_out, quote = F, row.names = F)
  
  filename_out <- "output/nm_pop_20mgq2w.csv"
  write.csv(simprep.list[[2]], file = filename_out, quote = F, row.names = F)
  
  filename_out <- "output/nm_pop_1280mgq2w.csv"
  write.csv(simprep.list[[3]], file = filename_out, quote = F, row.names = F)
