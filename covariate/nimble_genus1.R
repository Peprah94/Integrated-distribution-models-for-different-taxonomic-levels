source("nimble.R")

   estimates <-  run_nimble_model(simulations_all[[2]],
                                  method =  "IG", 
                                  crossvalidate = FALSE, 
                                  model = "covariate")


save(estimates, file="estimate_data_genus_hov1.RData")

