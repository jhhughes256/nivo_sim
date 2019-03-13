# Recreating the NONMEM populations from Nivolumab PK/PD Model Paper
# -----------------------------------------------------------------------------
# Script that produces NONMEM-compatible populations that match those used for
#   the simulations of the following paper:
#     C Liu, J Yu, H Li et al. (2017) Association of Time‐Varying Clearance of 
#     Nivolumab With Disease Dynamics and Its Implications on Exposure Response 
#     Analysis. Clin. Pharmacol. Ther., 101: 657-666. doi:10.1002/cpt.656

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
  dose <- c(0, 10, 20, 40, 80, 160, 320, 640, 1280, 2560)
  
# Set ID numbers (these will correlate with dose groups)
  id.seq <- c(1:10)
  
# Set simulation times
  times <- seq(0, 56, by = 1)
  
# Create amt and dv datasets and merge together
  simprep_amt <- data.frame(
    ID = id.seq,
    BASESLD = 60,
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
  simprep$II[simprep$TIME == 0 & simprep$ID != 1] <- 14
  simprep$ADDL[simprep$TIME == 0 & simprep$ID != 1] <- 3

# Fix up columns for NONMEM
  simprep$EVID[simprep$TIME == 0] <- 1
  simprep$DV[simprep$TIME == 0] <- "."
  simprep$MDV[simprep$TIME == 0] <- 1
  simprep$AMT[simprep$TIME != 0] <- "."
  simprep$II[simprep$TIME != 0] <- "."
  simprep$ADDL[simprep$TIME != 0] <- "."
  simprep <- simprep[simprep$AMT != 0, ]
  
# Rearrange columns
  simprep <- simprep[, c("ID", "TIME", "DV", "BASESLD", "ECOG", "AMT", "II", 
    "ADDL", "MDV", "EVID", "FLAG", "DOSEGRP")]
  
# Save to file
  names(simprep)[1] <- "#ID"
  filename_out <- "output/nonmem_pop.csv"
  write.csv(simprep, file = filename_out, quote = F, row.names = F)
  
  
  
