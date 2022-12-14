

# Fit gjam for bird application
# C. Lane Scher
# clanescher@gmail.com


library(tidyverse)
library(gjam)

# load data
load("DATA/birds/covar.rdata")
load("DATA/birds/counts.rdata")

all <- counts %>%
  inner_join(covar, by = "idYear")

# separate x and y
y <- all[,2:18]
x <- all[,c(1, 19, 27:28, 20:25, 29, 46:53)]

# formula
formulaI <- as.formula(~ elevation + 
                         prcpMean + prcpAnom +
                         tminMean + tminAnom +
                         forest + developed + water + wetlands + planted + shrub + herbaceous)

# set up gjam
effort <- list(columns=1:ncol(y), values = 150)

ng <- 30000
burnin <- ng/2
ml   <- list(ng = ng, burnin = burnin, typeNames = "DA", 
             effort = effort)

# run gjam
print(paste0("Number of species: ", ncol(y)))
startTime <- Sys.time()
print(paste0("Started at ", startTime))

out <- gjam(formulaI, xdata = x, ydata = y, modelList = ml)

endTime <- Sys.time()
print(paste0("Ended at ", endTime))

save(out, file = paste0("OUT/birds/birds-gjamOutput.rdata"))



