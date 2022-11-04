source("nimble.R")

estimates <-  run_nimble_model(simulations_all[[3]], 
                               "IG", 
                               crossvalidate = FALSE, 
                               model = "shared")


save(estimates, file="estimate_data_genus_sol1.RData")

