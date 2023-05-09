

#load data and function to run posterior predictive checks
source("CaseStudy/postPredCheck.R")
load("CaseStudy/data_idm.RData")
load("CaseStudy/hoverflies/estimateCaseStudy.RData")
allResults <- estimatesCaseStudy



# Parameterisations for function to run
simulations_all <- simulations_all[[3]] # PoMS data set for solitarybees
methodVars <- c("IDM", "Spe", "IDM", "IG", "IG")
modelVars <- c("shared", "shared", "covariate", "shared", "covariate")

# run the posterior predictive checks
ret <- lapply(1:5,
              predChecks,
              allResults,
              methodVars,
              modelVars,
              simulations_all,
              parallel = FALSE,
              n.samples = 2000)

#save results
save(ret, file = "CaseStudy/solitarybees/posteriorPredCheck.RData")











