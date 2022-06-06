library(unix)
unix::rlimit_as(100*10^9)
mbias <- function(true, estimate){
  ret <- mean(true-estimate) 
  return(ret)
}


mse <- function(true, estimate){
  ret <- mean((true-estimate)^2)
  return(ret)
}


hill_index <- function(q, pis){
    hill <- rowSums((pis^q))^(1/(1-q)) 
  return(hill)
}

#hill_index(2,pis)
#estimation of shannon index
shan_index <- function(pis){
  #shan <- exp(-rowSums(pis * log(pis)))
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



mbias <- function(true, estimate){
  ret <- mean(true-estimate) 
  return(ret)
}
aed <- function(pis){
  h0 <- hill_index(0, pis) 
  h1 <- exp(shan_index(pis))
  h2 <- hill_index(2, pis) 
  return(h0+((h1^2)/(h2*2)))
}


mse <- function(true, estimate){
  ret <- sqrt(mean((true-estimate)^2))
  return(ret)
}

#Trials
#shan_index(pis)
#hill_index(2, pis)
#eveness(pis)