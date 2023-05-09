# Script that simulates the parameters needed for simulation
set.seed(1994)
require("LaplacesDemon")
#simulation of covariance matrix
sigma_sim <- function(n){
  cov <- rinvwishart(nu=(n+1),S= diag(1,n, n))
  return(cov)
}

#load data formatted for PoMS data analysis
# extract the altitude values

load("CaseStudy/data_idm.RData")
latitude <- simulations_all[[1]][[4]]
latitude[is.na(latitude)] <- 100000
latitude <- as.numeric(scale(latitude))

# 20 species
N.species = 20

# Writing the parameters needed for simulation as lists and use them as imputs to 
# functionForSimulationWithMissing.R script

input20 <- list(
  constants = list(
    n.sites = 74,
    n.species = N.species,
    n.visit = 8,
    n.id= 50,
    n.gen= 50,
    rho = 1.5,
    nspecies_no_missing = floor(N.species/5),
    n.replicate = 5,
    q = 2
  ) ,
  covariate =list(
    ecological = latitude,
    detection = rnorm(74,0,1)
  ),
  parameters = list(
    ecological=list(
      betaSpecies = rnorm(N.species,0,0.2),
     # betaVisits = rnorm(N.species,0,0.2),
      betaLatitude = rnorm(N.species,-0.5,1),
     betaLambda = rnorm(1,0, 0.2),
     betaVisits = rnorm(8, 0, 0.5)
    ),
    detection=list(
      alphaSites = rnorm(74, 0, 0.3),
      alphaSpecies = rnorm(N.species, 0,1)
    )
  ),
  interaction = sigma_sim(N.species),
  visit = rep(1:8, each = 74)
)

save(input20, file="idmFunctions/parameters_20.RData")
