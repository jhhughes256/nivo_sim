# Nivolumab  -  Simulate representative patient population with test dosing regimen
# Standard dosing - Dose = Target AUC*(CRCL+25) [Target AUC = 6mg.min/mL]

# ------------------------------------------------------------------------------
# Clear workspace
rm(list=ls(all=TRUE))
graphics.off()

# Set working directory
setwd("C:/Users/Ugo/Desktop/TestPK/")

# Load pakage libraries
library(plyr)	      #New plyr - required for mrgsolve
library(dplyr)	    #Split and rearrange data - required for mrgsolve
library(mrgsolve)	  #Metrum differential equation solver for pharmacometrics
library(ggplot2)    #Graphical package

# Source PopPK model script
source("NivolumabPopPKPopulationNoTS.R")

regimen.name        <- "Regimen"

# ------------------------------------------------------------------------------
# Replicate test population for concentration dataset
TIME.conc           <-  seq(from = 0, to = 672, by = 1)
test.conc           <-  lapply(population.data,rep.int, times=length(TIME.conc))
test.conc           <-  as.data.frame(test.conc)
test.conc           <-  test.conc[with(test.conc, order(test.conc$ID)),]
test.conc$time      <-  TIME.conc
test.conc$amt       <-  0
test.conc$cmt       <-  1
test.conc$evid      <-  0
test.conc$rate      <-  0
head(test.conc)

# Replicate test population for dose dataset
TIME.dose           <-  c(0,336)
test.dose           <-  lapply(population.data,rep.int, times=length(TIME.dose))
test.dose           <-  as.data.frame(test.dose)
test.dose           <-  test.dose[with(test.dose, order(test.dose$ID)),]
test.dose$time      <-  TIME.dose
test.dose$amt       <-  240
test.dose$cmt       <-  1
test.dose$evid      <-  1
test.dose$rate      <-  0
head(test.dose)

# Combine into simulation dataset
test.input          <-  rbind(test.conc, test.dose)
test.input          <-  test.input[with(test.input, order(test.input$ID, test.input$time, test.input$evid)),]
head(test.input)

# ------------------------------------------------------------------------------
# Simulate standard dosing regimen
test.data           <-  mod %>% data_set(test.input) %>% carry.out(amt, evid) %>% mrgsim()
test.data           <-  as.data.frame(test.data)
head(test.data)

# Plot
p <- ggplot(test.data, aes(x = time, y = IPRED))
p <- p + geom_line(col = "blue", size = 0.8)
p <- p + facet_wrap(~ID)
p

#Save simulated dataset
write.csv(test.data,file=paste("NivolumabPopPK",regimen.name,"_Data.csv",sep=""))