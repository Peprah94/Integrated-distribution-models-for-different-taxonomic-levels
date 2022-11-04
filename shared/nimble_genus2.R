source("nimble.R")

   estimates <-  run_nimble_model(simulations_all[[1]], 
                                  method = "IG", 
                                  crossvalidate = FALSE, 
                                  model = "shared")

save(estimates, file="estimate_data_genus_bum1.RData")

