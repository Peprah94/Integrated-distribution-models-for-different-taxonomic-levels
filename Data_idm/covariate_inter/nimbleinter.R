source("nimble.R")
  estimates <- run_nimble_model(simulations_all[[3]], "IDM", crossvalidate = FALSE, shared = "covariate_inter")


save(estimates, file="estimate_data_inter_sol.RData")

