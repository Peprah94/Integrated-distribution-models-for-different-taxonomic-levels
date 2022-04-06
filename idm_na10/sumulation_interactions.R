#Parameters needed for this function
load("/idm_functions/parameters_10.RData")
source("/idm_functions/function_for_simulation_with_missing.R")

#List of parameters for simulation
input_list_na <- list(input10)

#Number of replicates
nreplicates <- 200

simulations_all <- pblapply(input_list_na, function(x){
  pblapply(1:nreplicates, function(z){
    data <- sim(x, seed = z)
  }, cl=4)
}, cl=4)

simulations_all_na <- flatten(simulations_all)

#simulations_all =sim(input)
save(simulations_all_na, file="sim_interractions_na.RData")
save(input_list_na, file="sim_input_na.RData")



