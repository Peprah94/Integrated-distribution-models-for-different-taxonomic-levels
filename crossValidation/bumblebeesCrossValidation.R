#run cross-validation for bumblebees

source("crossValidation/cross_validation.R")
load("CaseStudy/data_idm.RData")
methodVars <- c("IG", "IG","Spe", "IDM", "IDM")
modelVars <- c("shared", "covariate", "shared","shared","covariate")
#estimatesCaseStudy <- list()
#for(iCount in 1:5){
library(doParallel)
library(foreach)
library(parallel)


cl <- parallel::makeCluster(5)
doParallel::registerDoParallel(cl)
#doParallel::registerDoSEQ(cl)
setDefaultCluster(cl)

clusterExport(cl, c("run_nimble_model",
                    "mysum",
                    # "formatMatrix",
                    "nimble_sum",
                    #"nimbleFormatMatrix",
                    "calcCrossValSD",
                    "mycalcCrossVal",
                    "myrunCrossValidate",
                    "MSElossFunction"))

estimatesCaseStudy <- foreach(iter = seq_along(modelVars), .packages = c("pbapply", "nimble", "MCMCglmm", "coda", "parallel", "foreach", "doParallel")) %dopar% {
  tryCatch({ run_nimble_model(simulations_all[[1]], model = methodVars[iter],covariance_prior = "LV", method = modelVars[iter]) }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

save(estimatesCaseStudy, file=paste0("crossValidation/bumblebees/estimateCrossValidate.RData"))
