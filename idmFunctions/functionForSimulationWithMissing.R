#Function needed to run the simulation
source("idmFunctions/fnxForEstimation.R")
##############################################################################################
#                         NOTES                                                              #
##############################################################################################
### We assume n.genus  g = 1,2,...,n.genus
# For each g genera, we assume n.species in that geunus
# We assume the data was collected at n.sites
#We also assume the total number of visits to each site is n.visit.
# We also assume that  some of the data are collected with species ID
# Most of them are collected with genus ID
#n.id = number of sites that have ID at either species or group level
#p.tag is the detectiion probability for the occurence model (species data)
##############################################################################################
set.seed(2020)

# Packages needed to run the script
require(coda)
require(nimble)
require(devtools)
require(mcmcplots)
require(MCMCglmm)
require(furrr)
require(purrr)
require(dplyr)
require(stringr)
require(MASS)
require("LaplacesDemon")
library(pbapply)
library("pscl")

# Function that simulates the covariance matrix
sigma_sim <- function(n){
  cov <- rinvwishart(nu=(n+1),S= diag(1,n, n))
  return(cov)
}

sim <- function(input, seed, method){
  n.sites = input$constants$n.sites
  n.species = input$constants$n.species
  n.visit = input$constants$n.visit
  n.id = input$constants$n.id
  n.gen = input$constants$n.gen
  n.replicate= input$constants$n.replicate
  q = input$constants$q
  rho = input$constants$rho

  set.seed(seed)

  #### CHECKS
  if(is.null(n.sites)){stop("please specify number of sites (n.sites)")}
  if(is.null(n.species)){stop("please specify number of species (n.species)")}
  if(is.null(n.visit)){stop("please specify number of visits (n.visit)")}
  if(is.null(n.id)){stop("please specify number of sites that have ID at the species level (n.id)")}
  if(is.null(n.gen)){stop("please specify number of species that have ID at the species level(n.gen)")}
  #if(is.null(cov)){stop("please provide covariates (n.cov)")}
  if(is.null(n.replicate)){stop("please number of replicates for the observation process")}
  if(!method %in% c("covariate", "shared")) stop("Method can only be covariate or shared")

  # are all the elements of the dimensions correct?

  if((length(n.sites) == length(n.species)) == FALSE)
  {stop("inputs(n.sites, n.species, n.visit) are not of the same length")}

  if((length(n.id) == length(n.gen)) == FALSE)
  {stop("inputs( n.id, n.gen) are not of the same length")}

  if((length(input$parameters$ecological$betaSpecies) !=n.species))
  {stop("length of ecological parameters = number of species")}

  #if((length(input$parameters$detection$alpha0) !=n.species) & (length(input$parameters$detection$alpha1)!= n.species))
  #{stop("length of detection parameters = number of species")}

  if(dim(input$interaction)[1] != dim(input$interaction)[1])
  {stop("Dimension of interaction variance must be n.species * n.species")}

  # Dimension of matrix to store the various parameter values

  N <-log_lambda <- vector("numeric", n.sites)
  z <- mu <-  lambda.s <- epsilon <- mu.lambda <- p.tag <- psi.s <- hillsLambda <- speHillsLambda <- groHillsLambda<- mu.Psi <- mu.all <- mu.Lambda <- matrix(NA, nrow = n.sites, ncol = n.species)
  x <- eta <- p.tag <-  array(NA, dim = c(n.sites, n.species, n.visit))
  y <- matrix(NA, nrow=n.sites, ncol=n.visit)

  # Simulation of species interaction effect
  for(site.tag in 1:n.sites){
    epsilon[site.tag,] = mvrnorm(1,rep(0, n.species),input$interaction )
  }

  visit = input$visit

  #### Simulation of data
  #linear predictor for ecological process
  if(method == "covariate"){
    for(site.tag in 1:n.sites){
      for(spe.tag in 1:n.species){
        mu.Psi[site.tag, spe.tag] = (input$parameters$ecological$betaSpecies)[spe.tag] + (input$parameters$ecological$betaLatitude)[spe.tag]* (input$covariate$ecological)[site.tag] + epsilon[site.tag,spe.tag]
        mu.Lambda[site.tag, spe.tag] = (input$parameters$ecological$betaLambda) + (input$parameters$ecological$betaLatitude)[spe.tag]* (input$covariate$ecological)[site.tag] + epsilon[site.tag,spe.tag]
        mu.all[site.tag, spe.tag] = (input$parameters$ecological$betaLambda) + (input$parameters$ecological$betaSpecies)[spe.tag] + (input$parameters$ecological$betaLatitude)[spe.tag]* (input$covariate$ecological)[site.tag] + epsilon[site.tag,spe.tag] #+ input$parameters$ecological$betaVisits[visit[site.tag]]
      }
    }
    }else{
      for(site.tag in 1:n.sites){
        for(spe.tag in 1:n.species){
      mu.all[site.tag, spe.tag] =  mu.Lambda[site.tag, spe.tag] = mu.Psi[site.tag, spe.tag] = (input$parameters$ecological$betaSpecies)[spe.tag] + (input$parameters$ecological$betaLatitude)[spe.tag]* (input$covariate$ecological)[site.tag] + epsilon[site.tag,spe.tag]
        }
      }
      }

    # estimation of occupancy probability
    lambda.s <- exp(mu.Lambda)

    psi.s <- invcloglog(mu.Psi)

    # Intensities for Hills Indices
    hillsLambda <- exp(mu.all)
    speHillsLambda <- exp(mu.Psi)
    groHillsLambda <- exp(mu.Lambda)

  #mean abundance for genus
  lambda.g <- rowSums(lambda.s)

  #detection probability.
  for(site.tag in 1:n.sites){
    for(spe.tag in 1:n.species){
      for(visit.tag in 1:n.visit){
      eta[site.tag, spe.tag, visit.tag]= input$parameters$detection$alphaSites[site.tag] + input$parameters$detection$alphaSpecies[spe.tag]+ input$parameters$ecological$betaVisits[visit.tag] #+ input$parameters$detection$alphaVisits[visit.tag]
      }
         }
  }
  p.tag <- invlogit(eta)

  #Simulation of observation process
  for(site.tag in 1:n.sites){
    for(spe.tag in 1:n.species){
      #True presence or absence
      z[site.tag, spe.tag] <- rbinom(1,1, psi.s[site.tag, spe.tag])
      for(k in 1:n.visit){
        # Genus counts
        y[site.tag,k] <- rnegbin(1, mu = lambda.g[site.tag], theta = rho)
        # Species Presence and absence observations
        x[site.tag, spe.tag,k] <- rbinom(1,n.replicate, psi.s[site.tag, spe.tag]*p.tag[site.tag, spe.tag, k])
      }
    }
  }


  # Making provision for the missing data at sites and species
#We assume 20 random sites for each visit will have a missing Species
  #and group Count information
  for(k in 1:n.visit){
    index_missing_genus <- sample(1:n.sites, n.id, replace=FALSE)
    index_missing_site_species <- sample(1:n.sites, n.id, replace=FALSE)
    y[-(index_missing_genus),k] <- NA
    x[-(index_missing_site_species), ,k] <- NA
  }

  # speciesData <- rbind(x[ , , 1],
  #                      x[ , , 2],
  #                      x[ , , 3],
  #                      x[ , , 4],
  #                      x[ , , 5],
  #                      x[ , , 6],
  #                      x[ , , 7],
  #                      x[ , , 8])
  #
  # genusData <- c(y[ ,1],
  #                y[ ,2],
  #                y[ ,3],
  #                y[ ,4],
  #                y[ ,5],
  #                y[ ,6],
  #                y[ ,7],
  #                y[ ,8])




  #Estimate hills Indices
  hillsIDM0 <- hill_index(0, hillsLambda)
  hillsSpe0 <- hill_index(0, speHillsLambda)
  hillsGro0 <- hill_index(0, groHillsLambda)
  hillsIDM1 <- hill_index(1, hillsLambda)
  hillsSpe1 <- hill_index(1, speHillsLambda)
  hillsGro1 <- hill_index(1, groHillsLambda)
  hillsIDM2 <- hill_index(2, hillsLambda)
  hillsSpe2 <- hill_index(2, speHillsLambda)
  hillsGro2 <- hill_index(2, groHillsLambda)



  # Returning the results
  data <- list(mat.species = x,
               mat.genus = y,
               visit = visit,
               #pis = pis,
               ecological_cov = input$covariate$ecological ,
               detection_cov = input$covariate$detection,
               p.tag = p.tag,
               psi.s = psi.s,
               z = z,
               hillsIDM0 = hillsIDM0,
               hillsIDM1 = hillsIDM1,
               hillsIDM2 = hillsIDM2,
               hillsSpe0 = hillsSpe0,
               hillsSpe1 = hillsSpe1,
               hillsSpe2 = hillsSpe2,
               hillsGro0 = hillsGro0,
               hillsGro1 = hillsGro1,
               hillsGro2 = hillsGro2
  )
  return(data)
}


##############################
# Simulating data for simulation study
#############################
load("idmFunctions/parameters_50.RData")
#List of input parameters
input_list_na <- list(input50)

#Number of replicates
nreplicates <- 100
simulationsShared <- pblapply(input_list_na, function(x){
  pblapply(1:nreplicates, function(z){
    data <- sim(x, seed = z, method = "shared")
  }, cl=4)
}, cl=4)

simulationsCovariate <- pblapply(input_list_na, function(x){
  pblapply(1:nreplicates, function(z){
    data <- sim(x, seed = z, method = "covariate")
  }, cl=4)
}, cl=4)



simulations_all_na <- flatten(c(simulationsShared,
                                simulationsCovariate))

#save the results
save(simulations_all_na, file="idmSimulations/simInterractionsNA50.RData")
save(input_list_na, file="idmSimulations/simInputNa.RData")
