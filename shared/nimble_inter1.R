source("nimble.R")

  estimates <- run_nimble_model(simulations_all[[2]], 
                                method = "IDM", 
                                crossvalidate = FALSE, 
                                model = "shared")


save(estimates, file="estimate_data_inter_hov1.RData")

