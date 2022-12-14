

###

# Conditional Prediction 
# C. Lane Scher, Kevin Kraus, Sarah Roberts, James S. Clark

# Lane Scher: clanescher@gmail.com

###

setwd("X:/conditionalPrediction/")


# Part 1: Simulations ----

source("CODE/11.simulation.r")
# Runs simulations and makes summary figure



# Part 2: Birds ----

source("CODE/21.counts.r")
# Cleans raw BBS data for analysis

source("CODE/22.nlcd.r")
# Downloads nlcd at BBS locations
# calls "CODE/nlcdFunctions.r"

source("CODE/23.env.r")
# Downloads environmental covariates, combines them with nlcd

source("CODE/24.gjam.r")
# Fit GJAM

source("CODE/25.output.r")
# make some figures with gjam output



# Part 3: Fish ----

# Fish analysis is contained in:
# CODE/31.fish.rmd




# Extra ----

# This makes the figures used in my ESA presentation
source("CODE/ESA.r")


