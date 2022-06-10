# Script that simulates the parameters needed for simulation
set.seed(1994)
require("LaplacesDemon")
#simulation of covariance matrix
sigma_sim <- function(n){
  cov <- rinvwishart(nu=(n+1),S= diag(1,n, n))
  return(cov)
}

#Parameters for simulation for S= 10 species
N.species = 10
input10 <- list(
  constants = list(
    n.sites = 75,
    n.species = N.species,
    n.visit = 4,
    n.id= 40,
    n.gen= 75,
    nspecies_no_missing = floor(N.species/5),
    n.replicate = 5,
    q = 2
  ) ,
  covariate =list(
    ecological = runif(75, -1,1),
    detection = runif(75,-1,1)
  ),
  parameters = list(
    ecological=list(
      beta0 = rnorm(N.species,0,2),
      beta1=rnorm(N.species,0,2)
    ),
    detection=list(
      alpha0 = rnorm(N.species,0,1),
      alpha1 = rnorm(N.species,0,1)
    )
  ),
  interaction = sigma_sim(N.species)
)

save(input10, file="parameters_10.RData")

#Simulation of replicates for S= 20 species
N.species = 20
input20 <- list(
  constants = list(
    n.sites = 75,
    n.species = N.species,
    n.visit = 4,
    n.id= 40,
    n.gen= 75,
    nspecies_no_missing = floor(N.species/5),
    n.replicate = 5,
    q = 2
  ) ,
  covariate =list(
    ecological = runif(75, -1,1),
    detection = runif(75,-1,1)
  ),
  parameters = list(
    ecological=list(
      beta0 = rnorm(N.species,0,2),
      beta1=rnorm(N.species,0,2)
    ),
    detection=list(
      alpha0 = rnorm(N.species,0,1),
      alpha1 = rnorm(N.species,0,1)
    )
  ),
  interaction = sigma_sim(N.species)
)

save(input20, file="parameters_20.RData")

# 50 species
N.species = 50
input50 <- list(
  constants = list(
    n.sites = 75,
    n.species = N.species,
    n.visit = 4,
    n.id= 40,
    n.gen= 75,
    nspecies_no_missing = floor(N.species/5),
    n.replicate = 5,
    q = 2
  ) ,
  covariate =list(
    ecological = runif(75, -1,1),
    detection = runif(75,-1,1)
  ),
  parameters = list(
    ecological=list(
      beta0 = rnorm(N.species,0,2),
      beta1=rnorm(N.species,0,2)
    ),
    detection=list(
      alpha0 = rnorm(N.species,0,1),
      alpha1 = rnorm(N.species,0,1)
    )
  ),
  interaction = sigma_sim(N.species)
)

save(input50, file="parameters_50.RData")
