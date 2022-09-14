#Load the required data
load("data_idm.RData")
source("cross_valid_function.R")

require(coda)
require(nimble)
require(MCMCglmm)
library(parallel)
library(pbapply)

op <- pboptions(type="timer")

mysum <- function(x){

    ret <- sum(x[!is.infinite(x)], na.rm = TRUE)
  #}
  return(ret)
}

nimble_sum <- nimbleRcall(
  prototype = function(
    x=double(1)
    # beta is a vector
  ) {},
  returnType = double(0), # outcome is a vector
  Rfun = 'mysum'
)

run_nimble_model <- function(simulations_all, method, crossvalidate, shared){
  #### Load required packages ####
  require(coda)
  require(nimble)
  require(mcmcplots)
  require(MCMCglmm)
  require(ggmcmc)
  library(parallel)
  library(pbapply)

  start_time <- Sys.time()
  
  #write nimble code
  code <- nimbleCode({
    
    #############################################################
    #                 PRIOR DISTRIBUTIONS                       #
    #############################################################

    for(spe.tag in 1:n.species){
      beta0[spe.tag]~ dnorm(0, sd=5) 
      alpha0[spe.tag] ~  dnorm(0,sd=5)
    }
    
    intercept_lambda ~ dnorm(0, sd = 10) # for shared covariate and interractions
    #detection and observation variation
      sigma.obs ~ dgamma(1,0.1)
      sigma.det ~ dgamma(1,0.1)
      sigma.year ~ dgamma(1,0.1)
      sigma.year.det ~ dgamma(1, 0.1)
    
  #Vague prior for variance covariance matrix
      # We used the inverse Wishart for this
      omega[1:n.species,1:n.species] ~ dwish(R[1:n.species,1:n.species], df)
      Cov[1:n.species,1:n.species] <- inverse(omega[1:n.species,1:n.species])
      
    for(site.tag in 1:n.sites){
    site.obs.vars[site.tag] ~ dnorm(0, tau = sigma.obs) # random site effect of observation process
    site.det.vars[site.tag] ~ dnorm(0, tau = sigma.det) # random site effect of detection probability
    eta.lam[site.tag, 1:n.species] ~ dmnorm(mean = mu.eta[1:n.species],
                                            prec= omega[1:n.species, 1:n.species])
    }  
      
      #year random effect
      for(year.tag in 1:n.years){
      year.obs.vars[year.tag] ~ dnorm(0, tau = sigma.year) # random site effect of observation process
      year.det.vars[year.tag] ~ dnorm(0, tau = sigma.year.det) # random site effect of detection probability
      }
      
    for (spe.tag in 1:n.species) {
      for (ss in 1:n.species) {
        CorrIn[spe.tag, ss] <- Cov[spe.tag, ss]/sqrt(Cov[spe.tag, spe.tag] * Cov[ss, ss])
      }
    }
    
      #Link between the abundance and occupancy
      if(shared == "all"){
        for(year.tag in 1:n.years){ #loop over years
          for(site.tag in 1:n.sites){ #loop over sites
            for(spe.tag in 1:n.species){#loop over species
              mu[site.tag,spe.tag, year.tag] <- beta0[spe.tag] + eta.lam[site.tag,spe.tag] + site.obs.vars[site.tag] + year.obs.vars[year.tag]
              log(lambda[site.tag, spe.tag, year.tag]) <- mu[site.tag,spe.tag, year.tag]
              cloglog(psi[site.tag, spe.tag, year.tag]) <- mu[site.tag,spe.tag, year.tag]
            }
          } 
        }
      }
      
      if(shared == "covariate_inter"){
        for(year.tag in 1:n.years){ #loop over years
          for(site.tag in 1:n.sites){ #loop over sites
            for(spe.tag in 1:n.species){#loop over species
              mu[site.tag,spe.tag, year.tag] <- beta0[spe.tag] + eta.lam[site.tag,spe.tag] + site.obs.vars[site.tag] + year.obs.vars[year.tag]
            cloglog(psi[site.tag, spe.tag, year.tag]) <-mu[site.tag,spe.tag, year.tag]
            log(lambda[site.tag, spe.tag, year.tag]) <- intercept_lambda  + eta.lam[site.tag,spe.tag] + site.obs.vars[site.tag] + year.obs.vars[year.tag]
          }
        }
        } 
      }
      
      if(shared == "interractions"){
        sigma.obs.lambda ~ dgamma(1,0.1)
        sigma.year.lambda ~ dgamma(1,0.1) 
        #create extra site and detection_effect for lambda
        for(site.tag in 1:n.sites){
          site.obs.vars.lambda[site.tag] ~ dnorm(0, tau = sigma.obs.lambda) # random site effect of observation process
        } 
        
        for(year.tag in 1:n.years){
          year.obs.vars.lambda[year.tag] ~ dnorm(0, tau = sigma.year.lambda) # random site effect of observation process
        } 
        
        for(year.tag in 1:n.years){ #loop over years
          for(site.tag in 1:n.sites){ #loop over sites
            for(spe.tag in 1:n.species){#loop over species
              # mu[site.tag,spe.tag, year.tag] <- beta0[spe.tag] + eta.lam[site.tag,spe.tag] + site.obs.vars[site.tag] + year.obs.vars[year.tag]
              cloglog(psi[site.tag, spe.tag, year.tag]) <- beta0[spe.tag] + eta.lam[site.tag,spe.tag] + site.obs.vars[site.tag] + year.obs.vars[year.tag]
              log(lambda[site.tag, spe.tag, year.tag]) <- intercept_lambda  + eta.lam[site.tag,spe.tag] + site.obs.vars.lambda[site.tag] + year.obs.vars.lambda[year.tag]
            }
          }
        } 
      }
      
      
    #Detection probability
      for(year.tag in 1:n.years){ #loop over years
    for(site.tag in 1:n.sites){
      for(spe.tag in 1:n.species){
        logit(p.tag[site.tag, spe.tag, year.tag]) <- alpha0[spe.tag] + site.det.vars[site.tag] + year.det.vars[year.tag]
      }
    }
      }
 
    
    # Likelihood: key definitions in the likelihood
    #############################################################
    #                 Abundance Model                           #
    #############################################################
      if(method == "IDM"){
        #Model for Group count data 
        for(year.tag in 1:n.years){ #loop over years
          for(site.tag in 1:n.sites){ # loop over sites
            lambda.g[site.tag, year.tag] <- sum(lambda[site.tag,1:n.species, year.tag])
            for(visit.tag in 1:n.visit){ #loop over visits
              Y[site.tag, visit.tag, year.tag] ~ dpois(lambda.g[site.tag, year.tag])
            }
          }
        }
        
#   #Model for species occupancy data
        for(year.tag in 1:n.years){ #loop over years
          for(site.tag in 1:n.sites){ #loop over sites
            for(spe.tag in 1:n.species){ #loop over species
              z[site.tag,spe.tag,year.tag] ~ dbern(psi[site.tag, spe.tag, year.tag])
              for(visit.tag in 1:n.visit){ #loop over visit
                X[site.tag,spe.tag,visit.tag, year.tag] ~ dbin(z[site.tag,spe.tag, year.tag]*p.tag[site.tag, spe.tag, year.tag], n.replicates)
              }
            }
          }
        }
        
        # estimate the value to use as input of the hills indices
        for(year.tag in 1:n.years){
          for(site.tag in 1:n.sites){
            for(spe.tag in 1:n.species){ 
            #pps[site.tag, spe.tag,year.tag] <- psi[site.tag,spe.tag, year.tag]*(1-pow((1- p.tag[site.tag, spe.tag, year.tag]),n.visit)) #z[site.tag,spe.tag,year.tag]
              pps[site.tag, spe.tag,year.tag] <- -log(1- psi[site.tag, spe.tag,year.tag]) 
              }
          }
        }
      }
      
      if(method == "IG"){
        for(year.tag in 1:n.years){
          for(site.tag in 1:n.sites){
            for(spe.tag in 1:n.species){
              z[site.tag,spe.tag,year.tag] <- 1
            }
          }
        }
        #Model for Group counnt data 
        for(year.tag in 1:n.years){ #loop over years
          for(site.tag in 1:n.sites){ # loop over sites
            lambda.g[site.tag, year.tag] <- sum(lambda[site.tag,1:n.species, year.tag])
            for(visit.tag in 1:n.visit){ #loop over visits
              Y[site.tag, visit.tag, year.tag] ~ dpois(lambda.g[site.tag, year.tag])
            }
          }
        }
        
        for(year.tag in 1:n.years){
          for(site.tag in 1:n.sites){
            for(spe.tag in 1:n.species){ 
              pps[site.tag, spe.tag,year.tag] <- lambda[site.tag,spe.tag, year.tag] #z[site.tag,spe.tag,year.tag]
            }
          }
        }
      }
      
      if(method == "Spe"){
        #Model for Species occupancy data
        for(year.tag in 1:n.years){
          for(site.tag in 1:n.sites){
            for(spe.tag in 1:n.species){
              z[site.tag,spe.tag,year.tag] ~ dbern(psi[site.tag, spe.tag, year.tag])
              for(visit.tag in 1:n.visit){
                X[site.tag,spe.tag,visit.tag, year.tag] ~ dbin(z[site.tag,spe.tag, year.tag]*p.tag[site.tag, spe.tag, year.tag], n.replicates)
              }
            }
          }
        }
        
        for(year.tag in 1:n.years){
          for(site.tag in 1:n.sites){
            for(spe.tag in 1:n.species){ 
              pps[site.tag, spe.tag,year.tag] <- -log(1- psi[site.tag, spe.tag,year.tag]) #z[site.tag,spe.tag,year.tag]
            }
          }
        }
        
      }
    
    #############################################################
    #                pis                           #
    #############################################################

      
      for(year.tag in 1:n.years){
        for(site.tag in 1:n.sites){
          pis.g[site.tag, year.tag] <- sum(pps[site.tag, 1:n.species,year.tag])
        }
      }
      
      for(year.tag in 1:n.years){
        for(site.tag in 1:n.sites){
      for(spe.tag in 1:n.species){
        pis[site.tag, spe.tag, year.tag] <- pps[site.tag, spe.tag,year.tag]/pis.g[site.tag, year.tag]
      }
    }
    }
#       
#     #############################################################
#     #               shannon Index                              #
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
    shan[site.tag, year.tag] <- - nimble_sum(shans[site.tag, 1:n.species, year.tag])
    }
    }
      
      for(year.tag in 1:n.years){
        for(site.tag in 1:n.sites){
          richness[site.tag, year.tag] <- sum(z[site.tag, 1:n.species, year.tag])
        }
      }    
 
  })
  
  # Extracting dimensions of data parameters
  data <- simulations_all
  dim_data <- dim(data[[1]]) 
  data_dim <- dim_data[2] 
  
  # Constants for the model 
  const <- list(n.sites = dim_data[1], 
                n.species= (dim_data[2]), 
                n.years = dim_data[4],
                n.replicates = 5,
                n.visit=dim_data[3],
                mu.eta = rep(0, dim_data[2])
  )
  

  #R and df
  Rmat <- diag(const$n.species) # In the paper, this is initial covariance matrix
  df <- const$n.species + 1 
  
  #Data for the model
  idm_data <- list(Y = simulations_all[[2]], 
                   X = simulations_all[[1]], 
                   R = Rmat,
                   df = df)
  
  #Retrieving the true Presence/Absence
  pa_data <- function(data){
    zst <- apply(data, c(1,2,4), function(x){
      max(x, na.rm = TRUE)
    })
    zst[is.infinite(zst)] <- 0
    zst[zst!= 0 ] <- 1
    return(zst)
  }

#initial values 
  idm_inits <- function(){list(beta0= rnorm(const$n.species,0,0.01),
                               alpha0=rnorm(const$n.species,0,0.01), #detection probability
                               sigma.det =100 ,
                               sigma.obs = 100,
                               sigma.obs.lambda = 100,
                               sigma.year = 100,
                               sigma.year.lambda = 100,
                               sigma.year.det = 100,
                               intercept_lambda = rnorm(1,0,1),
                               site.obs.vars = rnorm(const$n.sites, 0, 0.01),
                               site.obs.vars.lambda = rnorm(const$n.sites, 0, 0.01),
                               site.det.vars = rnorm(const$n.sites, 0,0.01),
                               year.obs.vars = rnorm(const$n.years, 0,0.01),
                               year.obs.vars.lambda = rnorm(const$n.years, 0,0.01),
                               year.det.vars = rnorm(const$n.years, 0,0.01),
                               z=pa_data(idm_data$X),
                               eta.lam = array(0, dim = c(const$n.sites, const$n.species)),
                               omega = diag(1, const$n.species) #Note this is the precision matrix
  )
  }
  
  initsList <- idm_inits() 
  
  
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
                            print = TRUE,
                            monitors = c("beta0", "alpha0", "sigma.det",
                                         "sigma.obs", "Cov", "shan", "omega", 
                                         "psi", "lambda", "z", "p.tag","intercept_lambda",
                                         "richness"))
 
 if(crossvalidate == TRUE){
   set.seed(0)
   indx <- sample(1:74, 74, replace = FALSE)
   if( method == "Spe"){
     fold_function <- function(i){
       if(i ==1 ){
       nodes <- mwtc$expandNodeNames(paste0('X[', indx[i: 34],', ,  ,]') )
       }
         if(i == 2){
           nodes <- mwtc$expandNodeNames(paste0('X[', indx[(i+33): 74],', , , ]') ) 
         }
       
       foldNodes_i <-mwtc$isData(nodes)
       ret <- nodes[foldNodes_i]
       return(ret)
     }
   }

     if(method == "IG"){
       fold_function <- function(i){
         if(i ==1 ){
           nodes <- mwtc$expandNodeNames(paste0('Y[', indx[i: 34],',, ]') )
         }
           if(i == 2){
             nodes <- mwtc$expandNodeNames(paste0('Y[', indx[(i+33): 74],', , ]') ) 
           }
             #if(i == 3){
             #  nodes <- mwtc$expandNodeNames(paste0('Y[', indx[(i+28): 45],',  , ]') ) 
            # }
            # if(i == 4){
            #   nodes <- mwtc$expandNodeNames(paste0('Y[', indx[(i+42): 60],',  , ]') )  
            # }
       #if(i == 5){
         #    nodes <- mwtc$expandNodeNames(paste0('Y[',indx[(i+56): 74],',  , ]') )  
          # }
         
         foldNodes_i <- mwtc$isData(nodes)
         ret <- nodes[foldNodes_i]
         return(ret)
       }
     }
if(method == "IDM"){
     fold_function <- function(i){
       if(i ==1 ){
         nodes <- mwtc$expandNodeNames(paste0('X[', indx[i: 34],', ,  ,]'), paste0('Y[', indx[i: 34],',, ]')  )
       }
       if(i == 2){
         nodes <- mwtc$expandNodeNames(paste0('X[', indx[(i+33): 74],', , , ]'), paste0('Y[', indx[(i+33): 74],', , ]') ) 
       }
       
       foldNodes_i <-mwtc$isData(nodes)
       ret <- nodes[foldNodes_i]
       return(ret)
     }
   }
       

 ret <- myrunCrossValidate(MCMCconfiguration = mcmcconf, 
                         k = 2,
                         nCores = 1,
                         lossFunction = "predictive",
                         MCMCcontrol = list(niter = 10000, #10000
                                            nburnin = 9900
                                            ), #5000
                         foldFunction = "random")
 }else{
   ret <- NULL
 }


  
   Rmcmc <- buildMCMC(mcmcconf)
  cmcmc <- compileNimble(Rmcmc, 
                         project = Cmwtc,
                         resetFunctions = TRUE) 
  
  estimate <- runMCMC(cmcmc,
                      niter = 100000, #300000
                      nburnin = 50000,#100000
                      #inits = initsList,
                      nchains = 3,
                      thin = 10,
                      summary=TRUE,
                      samples=TRUE,
                      setSeed = TRUE,
                      samplesAsCodaMCMC=TRUE, 
                      WAIC = FALSE)
  

  end.time <- Sys.time()
  time_taken <- as.numeric(round(end.time-start_time,2))
  
  #save the result as ggmcmc list
  mcmclist <- ggs(estimate$samples)
  
  #subset the parameters and check convergence

  #Estimate the Gelman diagnostic value of the subset parameters
  aa <- mcmclist%>%
    filter(grepl('alpha|beta',Parameter))
  
  subset_parameters<- unique(aa$Parameter)
  
  # Rhat values
  subset_Rhat <- aa%>%
    ggs_Rhat()
  
  Rhat_data <- subset_Rhat$data[,5]
  all_rhat <- all(Rhat_data < 1.1) #Rhat is less than 1.1
  N_over_rhat <-length(which(Rhat_data > 1.1))/length(subset_parameters)
  
  estimation <- estimate$summary$all.chains
  
  return(list(const, estimation,N_over_rhat,ret))

}


