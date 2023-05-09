# This script is used to simulate the data 
# calls for functions already written to simulate the data and loads the parameters
# simulated from the prior distibutions

#Parameters needed for this function
load("idmFunctions/parameters_50.RData")
source("idmFunctions/functionForSimulationWithMissing.R")

#Simulation of the data with replication at sites
#List of input parameters
input_list_na <- list(input20)

#Number of replicates
nreplicates <- 100

# Simulate data with shared formulation of joint likelihood in IDM
simulationsShared <- pblapply(input_list_na, function(x){
  pblapply(1:nreplicates, function(z){
    data <- sim(x, seed = z, model = "shared")
  }, cl=4)
}, cl=4)

# Simulate data with covariate formulation of joint likelihood in IDM
simulationsCovariate <- pblapply(input_list_na, function(x){
  pblapply(1:nreplicates, function(z){
    data <- sim(x, seed = z, model = "covariate")
  }, cl=4)
}, cl=4)


#Pu both datasets together
simulations_all_na <- flatten(c(simulationsShared,
                                simulationsCovariate))

#save the results
save(simulations_all_na, file="simInterractionsNA50.RData")
save(input_list_na, file="sim_input_na.RData")



