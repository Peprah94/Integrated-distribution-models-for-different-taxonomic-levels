source("nimble.R")

   estimates <-  run_nimble_model(simulations_all[[1]], 
                                  method = "IG", 
                                  crossvalidate = FALSE, 
                                  model = "correlation")

save(estimates, file="estimate_data_genus_bum1.RData")

