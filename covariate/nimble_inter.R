source("nimble.R")
  estimates <- run_nimble_model(simulations_all[[3]],
                                method =  "IDM", 
                                crossvalidate = FALSE, 
                                model = "covariate")


save(estimates, file="estimate_data_inter_sol1.RData")

