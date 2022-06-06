load("data_idm.RData")

source("fnx_for estimation.R")

require(coda)
require(nimble)
require(devtools)
require(mcmcplots)
require(MCMCglmm)
#require(furrr)
#require(purrr)
#require(future.apply)
library(parallel)
#library(foreach)
#library(progressr)
library(pbapply)

#handlers(global = TRUE)
op <- pboptions(type="timer")

run_nimble_model <- function(simulations_all){
  
  # Set up
  
  #### Load required packages ####
  require(coda)
  require(nimble)
  require(devtools)
  require(mcmcplots)
  require(MCMCglmm)
  #require(furrr)
  #require(purrr)
  #require(future.apply)
  library(parallel)
  library(ggmcmc)
  #library(foreach)
  start_time <- Sys.time()
  ############################################
  code <- nimbleCode({
    
    #############################################################
    #                 PRIOR DISTRIBUTIONS                       #
    #############################################################
    for(spe.tag in 1:n.species){
      beta0[spe.tag]~ dnorm(0, sd=5) 
      alpha0[spe.tag] ~  dnorm(0,sd=5)
      mu.eta[spe.tag] <- 0
    }
    sigma.obs ~ dgamma(1,0.1)
    sigma.det ~ dgamma(1,0.1)
    
    #Vague prior for variance covariance matrix
    omega[1:n.species,1:n.species] ~ dwish(R[1:n.species,1:n.species], df)
    Cov[1:n.species,1:n.species] <- inverse(omega[1:n.species,1:n.species])
    
    for(site.tag in 1:n.sites){
      site.obs.vars[site.tag] ~ dnorm(0, sd = sigma.obs)
      site.det.vars[site.tag] ~ dnorm(0, sd = sigma.det)
      eta.lam[site.tag, 1:n.species] ~ dmnorm(mean = mu.eta[1:n.species],prec= omega[1:n.species, 1:n.species])
    }  
    
    for (spe.tag in 1:n.species) {
      for (ss in 1:n.species) {
        CorrIn[spe.tag, ss] <- Cov[spe.tag, ss]/sqrt(Cov[spe.tag, spe.tag] * Cov[ss, ss])
      }
    }
    
    #precmatrix[1:n.species, 1:n.species] <- inverse(Cov[1:n.species, 1:n.species])
    #Link between the abundance and occupancy
    
    for(year.tag in 1:n.years){ #loop over years
      for(site.tag in 1:n.sites){ #loop over sites
        for(spe.tag in 1:n.species){#loop over species
          mu[site.tag,spe.tag, year.tag] <- beta0[spe.tag] + eta.lam[site.tag,spe.tag] + site.obs.vars[site.tag]#+ beta1[spe.tag]*ecological_cov[site.tag]
          log(lambda[site.tag, spe.tag, year.tag]) <-mu[site.tag,spe.tag, year.tag]
        }
      } 
    }
    
    #Detection probability
    for(site.tag in 1:n.sites){
      for(spe.tag in 1:n.species){
        logit(p.tag[site.tag, spe.tag]) <- alpha0[spe.tag]+ site.det.vars[site.tag] #+ alpha1[spe.tag]*detection_cov[site.tag]
      }
    }
    

    # Likelihood: key definitions in the likelihood
    #############################################################
    #                 Abundance Model                           #
    #############################################################
    
    for(year.tag in 1:n.years){ #loop over years
      for(site.tag in 1:n.sites){ # loop over sites
        lambda.g[site.tag, year.tag] <- sum(lambda[site.tag,1:n.species, year.tag])
        for(visit.tag in 1:n.visit){ #loop over visits
          Y[site.tag, visit.tag, year.tag] ~ dpois(lambda.g[site.tag, year.tag])
        }
      }
    }
    
    #############################################################
    #                pis                           #
    #############################################################
    
    for(year.tag in 1:n.years){
      for(site.tag in 1:n.sites){
        for(spe.tag in 1:n.species){
          pis[site.tag, spe.tag, year.tag] <- (lambda[site.tag,spe.tag, year.tag]/lambda.g[site.tag, year.tag])
        }
      }
    }
    #############################################################
    #               shannon Index                              #
    #############################################################
    for(year.tag in 1:n.years){
      for(site.tag in 1:n.sites){
        for(spe.tag in 1:n.species){
          shans[site.tag, spe.tag, year.tag] <- log(pis[site.tag, spe.tag, year.tag])*(pis[site.tag, spe.tag, year.tag])
        }
      }
    }
    
    for(year.tag in 1:n.years){
      for(site.tag in 1:n.sites){
        shan[site.tag, year.tag] <- - sum(shans[site.tag, 1:n.species, year.tag])
      }
    }
    
  })
  
  data <- simulations_all
  #dimensions of the data
  dim_data <- dim(data[[1]]) 
  data_dim <- dim_data[2] 
  
  
  
  # Constants for the model 
  const <- list(n.sites = dim_data[1], 
                n.species= (dim_data[2]), 
                n.years = dim_data[4],
                n.replicates = 5,
                n.visit=dim_data[3]
  )
  
  
  #R and df
  Rmat <- diag(const$n.species)
  df <- const$n.species + 1 
  #omega <- LaplacesDemon:: rwishart(df, Rmat)
  
  #Data for the model
  idm_data <- list(Y = simulations_all[[2]], 
                   X = simulations_all[[1]], 
                   R = Rmat,
                   df = df)
  
  #Retrieving the true Presence/Absence
  pa_data <- function(data){
    zst <- apply(data, c(1,2,4), max)
    zst[zst!= 0] <- 1
    zst[is.na(zst)] <- 0
    return(zst)
  }
  
  
  idm_inits <- function(){list(beta0= rnorm(const$n.species,0,1),
                               alpha0=rnorm(const$n.species,0,1), #detection probability
                               sigma.det = rgamma(1,1,0.1),
                               sigma.obs = rgamma(1,1,0.1),
                               z=pa_data(idm_data$X),
                               omega = LaplacesDemon:: rwishart(df, Rmat)
  )#True Presence/Absence
  }
  initsList <- idm_inits() 
  
  modelInfo <- list(
    code = code,
    constants = const,
    data = idm_data,
    inits = idm_inits
  )
  
  #############################################################
  #                 Compile Model                             #
  #############################################################
  
  mwtc <- nimbleModel(code,
                      name = "mwtc",
                      data = idm_data, 
                      constants = const, 
                      inits = initsList)
  
  Cmwtc <- compileNimble(mwtc)
  
  mcmcconf <- configureMCMC(Cmwtc, 
                            monitors = c("beta0", "alpha0", "sigma.det",
                                         "sigma.obs","shan", "Cov", "omega"))
  
  #assign block samplers to alpha0, beta0, site.obs.vars and site.det.vars
  mcmcconf$removeSamplers(c("alpha0", "beta0", "site.obs.vars", "site.det.vars"))
  mcmcconf$addSampler(c("alpha0"), "RW_block")
  mcmcconf$addSampler(c("beta0"), "RW_block")
  mcmcconf$addSampler(c("site.obs.vars"), "RW_block")
  mcmcconf$addSampler(c("site.det.vars"), "RW_block")
  
  Rmcmc <- buildMCMC(mcmcconf)
  cmcmc <- compileNimble(Rmcmc, project = Cmwtc,resetFunctions = TRUE) 
  
  estimate <- runMCMC(cmcmc,
                      niter = 300000,
                      nburnin = 100000,
                      inits = initsList,
                      nchains=4,
                      thin = 50,
                      summary=TRUE,
                      samples=TRUE,
                      setSeed = TRUE,
                      samplesAsCodaMCMC=TRUE)
  end.time <- Sys.time()
  time_taken <- as.numeric(round(end.time-start_time,2))
  
  mcmclist <- ggs(estimate$samples)
  #ggmcmc(mcmclist, file = "model_simple-diag3.pdf")
  
 # print(ggs_density(mcmclist, family = "beta0"))
 # print(ggs_density(mcmclist, family = "alpha0"))
  # ggs_traceplot(mcmclist, family = "beta0")
  # ggs_traceplot(mcmclist, family = "alpha0")
  
  
  subset_parameters <- c(paste0("beta0[",1:const$n.species,"]"),
                         paste0("alpha0[",1:const$n.species,"]")
  )
  
  aa <-mcmclist%>%
    filter(Parameter %in% subset_parameters)%>%
    ggs_Rhat()
  Rhat_data <- aa$data[,5]
  all_rhat <- all(Rhat_data < 1.05) #Rhat is less than 1.05
  N_over_rhat <-length(which(Rhat_data > 1.05))/length(subset_parameters) #Rhats over 1.04 
  
  estimation <- estimate$summary$all.chains
  #est_fixed <- estimation[((const$n.species*const$n.species)+1):((const$n.species*const$n.species)+(const$n.species*2)),1]
  #est_pis <- matrix(estimation[,1][((const$n.species*const$n.species)+((const$n.species)*4)+1) : ((const$n.species*const$n.species)+((const$n.species)*4)+(const$n.sites*const$n.species) )], nrow=const$n.sites, ncol=const$n.species, byrow=FALSE)
  return(list(const, estimation,N_over_rhat))
}

## Run the simulations
#genus_estimates_data <-run_nimble_model(simulations_all[[1]])
#cl <- makeCluster(2)
#setDefaultCluster(cl)
#setDefaultCluster(cl)
genus_estimates_data <- run_nimble_model(simulations_all[[3]])
#genus_estimates_data <- foreach(i=1:6, .export=c('run_nimble_model'))%dopar% {
#  run_nimble_model(simulations_all[[i]])
#}

save(genus_estimates_data, file="estimate_genus_data.RData")
