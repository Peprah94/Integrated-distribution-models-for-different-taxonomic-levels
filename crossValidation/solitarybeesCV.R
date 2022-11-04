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
solitarybees <- pblapply(as.list(1:3), function(x){
  load("data_idm.RData")
  sim(x, data = simulations_all[[3]], method = "IDM")
}, cl = cl)

save(solitarybees, file="solitarybeesIDMCV.RData")

solitarybees1 <- pblapply(as.list(1:3), function(x){
  load("data_idm.RData")
  sim(x, data = simulations_all[[3]], method = "IDM1")
}, cl = cl)

save(solitarybees1, file="solitarybeesIDM1CV.RData")

solitarybees2 <- pblapply(as.list(1:3), function(x){
  load("data_idm.RData")
  sim(x, data = simulations_all[[3]], method = "IDM2")
}, cl = cl)

save(solitarybees2, file="solitarybeesIDM2CV.RData")

solitarybees3 <- pblapply(as.list(1:3), function(x){
  load("data_idm.RData")
  sim(x, data = simulations_all[[3]], method = "Spe")
}, cl = cl)

save(solitarybees3, file="solitarybeesSpeCV.RData")

solitarybees4 <- pblapply(as.list(1:3), function(x){
  load("data_idm.RData")
  sim(x, data = simulations_all[[3]], method = "IG")
}, cl = cl)

save(solitarybees4, file="solitarybeesIGCV.RData")

stopCluster(cl)