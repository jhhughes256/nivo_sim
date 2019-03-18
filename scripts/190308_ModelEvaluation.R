# Evaluation of Hazard (NivoPKPD) and Modified PK (NivoPKTS_wPD) Nivolumab Models
# -----------------------------------------------------------------------------
# To evaluate the hazard model, the manuscript took the point estimate of the 
#   hazard ratio in the quartiles of three different exposure metrics and compared 
#   them to the hazard rate of the data at those same exposure quartiles. 

# The three exposure metrics were:
#   CavgT - Average concentration from first dose to the time of event/dropout
#   Cavg1 - Average concentration at the 1st cycle (1st - 13th day)
#   Cmin1 - Trough concentration at the 1st cycle (13th day)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Ready work environment
# Clear workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set working directory
# if not working with RStudio project place working directory here
# setwd("C:/.../nivo_sim/")

# Load package libraries
  library(plyr)  # ddply
  library(reshape2)  # melt
  
# Load .rds files containing model simulations
  pkpd.sim <- readRDS("rds/190302_NivoPKPD_Reg.rds")
  pkts.sim <- readRDS("rds/190302_NivoPKTS_Reg.rds")
  
# Create matched data
  pkpd.df <- ddply(pkpd.sim[pkpd.sim$DOSEGRP != 1,], .(DOSEGRP), function(df) {
    cbind(df, 
      data.frame(
        HAZPLAC = pkpd.sim$HAZRATE[pkpd.sim$DOSEGRP == 1]
      )
    )
  })
  pkts.df <- ddply(pkts.sim[pkts.sim$DOSEGRP != 1,], .(DOSEGRP), function(df) {
    cbind(df, 
      data.frame(
        HAZPLAC = pkts.sim$HAZRATE[pkts.sim$DOSEGRP == 1]
      )
    )
  })

  pkpd.df$HAZRAT <- with(pkpd.df, HAZRATE/HAZPLAC)
  pkts.df$HAZRAT <- with(pkts.df, HAZRATE/HAZPLAC)
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Determine exposure metric: CavgT
# **** MAY REQUIRE HIGHER SIMULATION RESOLUTION ****
  CavgT.fn <- function(df) {
    nodrp.df <- df[df$HASEVT == 0 & df$HASDRP == 0, ]
    DRPt <- tail(nodrp.df$time, 1)
    CavgT <- mean(nodrp.df$IPRED)
    out <- cbind(df, data.frame(CavgT, DRPt))
  }
  pkpd.df <- ddply(pkpd.df, .(ID, DOSEGRP), CavgT.fn)
  pkts.df <- ddply(pkts.df, .(ID, DOSEGRP), CavgT.fn)
  
# Create quartile function
  quartile.fn <- function(df, metric.name) {
    newcol <- paste0("Q", metric.name)
    metcol <- paste0("C", metric.name)
    quart <- quantile(unlist(df[metcol]))
    df[newcol] <- 1
    df[newcol][df[metcol] > quart["25%"]] <- 2
    df[newcol][df[metcol] > quart["50%"]] <- 3
    df[newcol][df[metcol] > quart["75%"]] <- 4
    df
  }
  
# Split into quartiles
  pkpd.df <- ddply(pkpd.df, .(DOSEGRP), quartile.fn, metric.name = "avgT")
  pkts.df <- ddply(pkts.df, .(DOSEGRP), quartile.fn, metric.name = "avgT")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Determine exposure metric: Cavg1
  Cavg1.fn <- function(df) {
    cycle1.df <- df[df$time < 14, ]
    Cavg1 <- mean(cycle1.df$IPRED)
    out <- cbind(df, data.frame(Cavg1))
  }
  pkpd.df <- ddply(pkpd.df, .(ID, DOSEGRP), Cavg1.fn)
  pkts.df <- ddply(pkts.df, .(ID, DOSEGRP), Cavg1.fn)
  
# Split into quartiles
  pkpd.df <- ddply(pkpd.df, .(DOSEGRP), quartile.fn, metric.name = "avg1")
  pkts.df <- ddply(pkts.df, .(DOSEGRP), quartile.fn, metric.name = "avg1")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Determine exposure metric: Cmin1
  Cmin1.fn <- function(df) {
    Cmin1 <- df$IPRED[df$time == 13]
    out <- cbind(df, data.frame(Cmin1))
  }
  pkpd.df <- ddply(pkpd.df, .(ID, DOSEGRP), Cmin1.fn)
  pkts.df <- ddply(pkts.df, .(ID, DOSEGRP), Cmin1.fn)
  
# Split into quartiles
  pkpd.df <- ddply(pkpd.df, .(DOSEGRP), quartile.fn, metric.name = "min1")
  pkts.df <- ddply(pkts.df, .(DOSEGRP), quartile.fn, metric.name = "min1")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Determine hazard rates for each exposure metric
  Qhaz.fn <- function(df, metric) {
    sub.df <- df[df$time == 112, ]
    Qmetric <- sub.df[, metric]
    mean.hazrat <- c(
      median(sub.df$HAZRAT[Qmetric == 1], na.rm = T),
      median(sub.df$HAZRAT[Qmetric == 2], na.rm = T),
      median(sub.df$HAZRAT[Qmetric == 3], na.rm = T),
      median(sub.df$HAZRAT[Qmetric == 4], na.rm = T)
    )
    ci90lo.hazrat <- c(
      quantile(sub.df$HAZRAT[Qmetric == 1], prob = 0.05),
      quantile(sub.df$HAZRAT[Qmetric == 2], prob = 0.05),
      quantile(sub.df$HAZRAT[Qmetric == 3], prob = 0.05),
      quantile(sub.df$HAZRAT[Qmetric == 4], prob = 0.05)
    )
    ci90hi.hazrat <- c(
      quantile(sub.df$HAZRAT[Qmetric == 1], prob = 0.95),
      quantile(sub.df$HAZRAT[Qmetric == 2], prob = 0.95),
      quantile(sub.df$HAZRAT[Qmetric == 3], prob = 0.95),
      quantile(sub.df$HAZRAT[Qmetric == 4], prob = 0.95)
    )
    
    data.frame(
      METRIC = metric,
      Q = rep(1:4),
      MED.HAZ = mean.hazrat,
      CI90LO.HAZ = ci90lo.hazrat,
      CI90HI.HAZ = ci90hi.hazrat
    )
  }
  
  pkpd.haz <- rbind(
    ddply(pkpd.df, .(DOSEGRP), Qhaz.fn, metric = "QavgT"),
    ddply(pkpd.df, .(DOSEGRP), Qhaz.fn, metric = "Qavg1"),
    ddply(pkpd.df, .(DOSEGRP), Qhaz.fn, metric = "Qmin1")
  )
  pkts.haz <- rbind(
    ddply(pkts.df, .(DOSEGRP), Qhaz.fn, metric = "QavgT"),
    ddply(pkts.df, .(DOSEGRP), Qhaz.fn, metric = "Qavg1"),
    ddply(pkts.df, .(DOSEGRP), Qhaz.fn, metric = "Qmin1")
  )
  pkpd.haz$MODEL <- "Hazard"
  pkts.haz$MODEL <- "Modified PK"
  haz.data <- rbind(pkpd.haz, pkts.haz)
  
# Try to plot it... but it doesn't look, good, confidence intervals are wayyy
#   wider than those in figure 4 of the manuscript
  p <- NULL
  p <- ggplot(data = haz.data[haz.data$DOSEGRP == 3,])
  p <- p + geom_point(aes(x = Q, y = MED.HAZ, colour = MODEL))
  p <- p + geom_errorbar(aes(x = Q, ymin = CI90LO.HAZ, ymax = CI90HI.HAZ, 
    colour = MODEL), size = 1)
  p <- p + facet_wrap(~ METRIC)
  p
  
  