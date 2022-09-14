library(unix)
unix::rlimit_as(100*10^9)

#Mean bias
mbias <- function(true, estimate){
  ret <- mean(true-estimate) 
  return(ret)
}

#Mean square error
mse <- function(true, estimate){
  ret <- mean((true-estimate)^2)
  return(ret)
}

#Standardized Root mean square error
srmse <- function(true, estimate){
diff <- true - estimate
sd_diff <- sd(diff)
ret <- mean(((true-estimate)^2)/sd_diff)
  return(ret)
}

# Root mean square error
rmse <- function(true, estimate){
  ret <- mean((true-estimate)^2)
  return(ret)
}

#Mean absolute error
mae <- function(true, estimate){
  ret <- abs(true-estimate)
  ret <- mean(ret[!is.infinite(ret)], na.rm = TRUE)
  return(ret)
}

#incidence estimate
incidence <- function(lambdas, prob, K=4){ #K is the number of visits
  incidence <- lambdas*(1-(1-prob)^K)
  return(incidence)
}

# Hill index
hill_index <- function(q, lambda, prop = TRUE){
  if(prop == TRUE){
pis = lambda
  }else{
    pis <- proportions(lambda, margin = 1)
  }
  if(q != 1){
    hill <- (rowSums(pis^q, na.rm = TRUE))^(1/(1-q))
  }else{
      hill <- exp(-rowSums(log(pis)*(pis), na.rm = TRUE))  
      }
  return(hill)
}

#Richness
richness <- function(z){
  rowSums(z)
}

#Shannon index
shan_index <- function(pis){
 shan <-  -rowSums(log(pis)*(pis), na.rm = TRUE)
  return(shan)
}

#Estimation of eveness
eveness <- function(pis){
  shan <- shan_index(pis)
  richness <- log(hill_index(0,pis))
  even <- shan/richness
return(even)
  }


# Function to subset parameters from MCMC outputs
subsetting_parameters <- function(x,name, method){
  if(method %in% c("mean", "sd")){
  if(method == "mean"){
  ret <- x[grep(name, rownames(x)),1]
  }
  if(method == "sd"){
    ret <- x[grep(name, rownames(x)),3]
  }
  }else{
      stop("method should have name mean or sd")
    }
  return(ret)
}

#cloglog function
myclog <- function(psi){
  return(-log(1-psi))
}

#inverse cloglog function
myinvclog <- function(mu){
  return(1-exp(-exp(mu)))
}



hill_index_alt <- function(q, lambdas, zs){
  pis <- proportions(lambdas * zs, margin = 1)
  if(q != 0 & q != 1){
    hill <- rowSums((pis^q))^(1/(1-q)) 
  }else{
    if(q == 0){
      hill <- rowSums(zs)
    }else{if(q == 1){
      hill <- exp(-rowSums(log(pis)*(pis), na.rm = TRUE))  
    }
    }
    
  }
  return(hill)
}


