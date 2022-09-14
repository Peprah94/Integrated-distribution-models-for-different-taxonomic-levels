source("nimble.R")
  estimates <- run_nimble_model(simulations_all[[1]], "IDM", crossvalidate = TRUE, shared = "all")


save(estimates, file="estimate_data_inter_bum.RData")

