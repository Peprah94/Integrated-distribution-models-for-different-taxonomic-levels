# Load data
load("data_idm.RData")
source("cross_validation.R")

#Packages needed to run this script
library(doParallel)
library(parallel)
library(foreach)
library(iterators)

#clusters for parallization
cl <- makeCluster(3)
setDefaultCluster(cl)

# function for the parallelisation
sim <- function(i, data, method){
  loss_method = "predictive"
  model = as.list(rep(c("shared", "covariate", "correlation"), 1))
  ret <- run_nimble_model(simulations_all = data, method = method,
                          loss_method =  loss_method, 
                          shared = shared[[i]],
                          covariance_prior = "LV")
  return_value = c(loss_method, model[[i]], method,ret)
  return(return_value)
}

#export functions to the cluster
clusterExport(cl, c("run_nimble_model","myrunCrossValidate",
                    "MSElossFunction", "mycalcCrossVal", "RMSElossFunction", "calcCrossValSD", 
                    "generateRandomFoldFunction", "mysum", "nimble_sum", "sim"))

#export packages to the cluster
clusterEvalQ(cl, {
  library(ggplot2)
  library(stringr)
  library(nimble)
  library(parallel)
  library(foreach)
  library(coda)
  library(MCMCglmm)
  library(pbapply)
  library(ggmcmc)
  library(mcmcplots)
})


#run and save results

bumblebees <- pblapply(as.list(1:3), function(x){
  load("data_idm.RData")
  sim(x, data = simulations_all[[1]], method = "IDM")
}, cl = cl)

save(bumblebees, file="bumblebeesIDMCV.RData")

bumblebees1 <- pblapply(as.list(1:3), function(x){
  load("data_idm.RData")
  sim(x, data = simulations_all[[1]], method = "IDM1")
}, cl = cl)

save(bumblebees1, file="bumblebeesIDM1CV.RData")

bumblebees2 <- pblapply(as.list(1:3), function(x){
  load("data_idm.RData")
  sim(x, data = simulations_all[[1]], method = "IDM2")
}, cl = cl)

save(bumblebees2, file="bumblebeesIDM2CV.RData")

bumblebees3 <- pblapply(as.list(1:3), function(x){
  load("data_idm.RData")
  sim(x, data = simulations_all[[1]], method = "Spe")
}, cl = cl)

save(bumblebees3, file="bumblebeesSpeCV.RData")

bumblebees4 <- pblapply(as.list(1:3), function(x){
  load("data_idm.RData")
  sim(x, data = simulations_all[[1]], method = "IG")
}, cl = cl)

save(bumblebees4, file="bumblebeesIGCV.RData")

stopCluster(cl)


