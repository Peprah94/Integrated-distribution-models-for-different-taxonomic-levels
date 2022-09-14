# Functions and data needed for this script
load("idm_functions/parameters_50.RData")
source("idm_functions/function_for_simulation.R")
#Simulation of the data with replication at sites

#List of input parameters
input_list_na <- list(input50)

#Number of replicates
nreplicates <- 100
simulations_all1 <- pblapply(input_list_na, function(x){
  pblapply(1:nreplicates, function(z){
    data <- sim(x, seed = z, shared = "all")
  }, cl=4)
}, cl=4)

simulations_all_na1 <- flatten(simulations_all1)
 
 simulations_all2 <- pblapply(input_list_na, function(x){
   pblapply(1:nreplicates, function(z){
     data <- sim(x, seed = z, shared = "covariate")
   }, cl=4)
 }, cl=4)
 simulations_all_na2 <- flatten(simulations_all2)
 
 simulations_all3 <- pblapply(input_list_na, function(x){
   pblapply(1:nreplicates, function(z){
     data <- sim(x, seed = z, shared = "interractions")
   }, cl=4)
 }, cl=4)
 simulations_all_na3 <- flatten(simulations_all3)
 
 simulations_all_na <- c(simulations_all_na1, simulations_all_na2, simulations_all_na3)
 
#save the results
save(simulations_all_na, file="sim_interractions_na.RData")
save(input_list_na, file="sim_input_na.RData")



