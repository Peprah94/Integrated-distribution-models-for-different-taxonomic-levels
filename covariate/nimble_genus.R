source("nimble.R")

estimates <-  run_nimble_model(simulations_all[[3]],
                               method =  "IG", 
                               crossvalidate = FALSE, 
                               model = "covariate")


save(estimates, file="estimate_data_genus_sol1.RData")

