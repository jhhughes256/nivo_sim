# Exploring the copula package for simulation of multivariate distributions
# -----------------------------------------------------------------------------
# Worked out from this stack overflow question:
# https://stackoverflow.com/questions/47684693/copula-and-simulation-of-binary-and-continuous-variables
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare work environment
# Clear workspace
  rm(list=ls(all=TRUE))
  graphics.off()
  
# Load libraries
  library(copula)
  library(ggplot2)
  library(cowplot)

# Define colourblind palette and custom palette
  cbPalette <- data.frame(
    grey = "#999999",
    orange = "#E69F00",
    skyblue = "#56B4E9",
    green = "#009E73",
    yellow = "#F0E442",
    blue = "#0072B2",
    red = "#D55E00",
    pink = "#CC79A7",
    stringsAsFactors = F
  )

# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Simple example
  set.seed(123)
  myCop <- normalCopula(param = c(0.7), dim = 2, dispstr = "un")
  
  out <- rCopula(1e5, myCop)
  out[, 1] <- qnorm(out[, 1], mean = 79.09, sd = 19.28)
  out[, 2] <- qnorm(out[, 2], mean = 168.86, sd = 10)
  
  cor(out)
  colMeans(out)
  c(sd(out[, 1]), sd(out[, 2]))
  
  out.df <- as.data.frame(out)
  names(out.df) <- c("WT", "HT")
  
  p1 <- NULL
  p1 <- ggplot(data = out.df)
  p1 <- p1 + geom_point(aes(x = WT, y = HT), alpha = 0.25, shape = 1,
    colour = cbPalette$blue)
  p1
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Slightly more complex example
# mean(WT)  = 79.09, sd(WT)  = 19.28; normal distribution
# mean(AGE) = 61.12, sd(AGE) = 11.12; normal distribution
# mean(SEX) = 0.667; binary
# mean(GFR) = 78.49, sd(GFR) = 21.63; normal distribution
# cor
#        WT   AGE   SEX   GFR
#  WT  1.00  0.00  0.00 -0.46
# AGE  0.00  1.00  0.00 -0.28
# SEX  0.00  0.00  1.00  0.20
# GFR -0.46 -0.28  0.20  1.00
  
  set.seed(123)
  myCop <- normalCopula(param = c(0, 0, -0.46, 0, -0.28, 0.26), 
    dim = 4, dispstr = "un")
  
  out <- rCopula(1e5, myCop)
  out[, 1] <- qnorm(out[, 1], mean = 79.09, sd = 19.28)
  out[, 2] <- qnorm(out[, 2], mean = 61.12, sd = 11.12)
  out[, 3] <- qbinom(out[, 3], size = 1, prob = 0.667)
  out[, 4] <- qnorm(out[, 4], mean = 78.49, sd = 21.63)
  
  cor(out)
  colMeans(out)
  c(sd(out[, 1]), sd(out[, 2]), sd(out[, 4]))
  
  out.df <- as.data.frame(out)
  names(out.df) <- c("WT", "AGE", "SEX", "GFR")
  out.df$SEXf <- factor(out.df$SEX)
  
  p2 <- NULL
  p2 <- ggplot(data = out.df)
  p2 <- p2 + geom_point(aes(x = WT, y = GFR), alpha = 0.25, shape = 1,
    colour = cbPalette$blue)
  p2 <- p2 + facet_wrap(~SEX)
  p2  
  
  p3 <- NULL
  p3 <- ggplot(data = out.df)
  p3 <- p3 + geom_point(aes(x = AGE, y = GFR), alpha = 0.25, shape = 1,
    colour = cbPalette$blue)
  p3 <- p3 + facet_wrap(~SEX)
  p3
  
  p4 <- NULL
  p4 <- ggplot(data = out.df)
  p4 <- p4 + geom_boxplot(aes(x = SEXf, y = GFR), alpha = 0.25, shape = 1,
    colour = cbPalette$blue)
  p4