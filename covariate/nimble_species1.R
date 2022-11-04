source("nimble.R")

estimates <-  run_nimble_model(simulations_all[[2]],
                               method =  "Spe", 
                               crossvalidate = FALSE, 
                               model = "covariate")


save(estimates, file="estimate_data_species_hov1.RData")

