#Need to set the working directory to the folder with the simulations
setwd("/Volumes/kwakupa/IDM_simulations")

# Packages needed to run the code
library(ggplot2)
library(pbapply)
library(reshape2)
library(purrr)
library(dplyr)
library(gganimate)
library(ggpubr)
library(pbapply)
library(ggplot2)
source("fnx_for estimation.R")

nspecies = 3; nsite=75; niter= 100

###################
# NO MISSING DATA
###################

#Read the data for the N.species = 10 data set
load("idm_miss_na_10/estimate_inter_na.RData")
load("idm_miss_na_10/estimate_genus_na.RData")
load("idm_miss_na_10/estimate_species_na.RData")

## IDM approach used
# The pi's for the shannon index
inter_prop10 <- lapply(inter_estimates_na, function(x){
  if(class(x)!= "try-error"){x[[1]]}
})

#estimates of parameters
inter_est10<- lapply(inter_estimates_na, function(x){
  if(class(x)!= "try-error"){x[[2]]}
})

## Genus only approach used
# The pi's for the shannon index
genus_prop10 <- lapply(genus_estimates_na, function(x){
  if(class(x)!= "try-error"){x[[1]]}
})

#Estimates of parameters
genus_est10 <- lapply(genus_estimates_na, function(x){
  if(class(x)!= "try-error"){x[[2]]}
})

## Species only approach used
# The pi's for the shannon index
species_est10 <- lapply(species_estimates_na, function(x){
  if(class(x)!= "try-error"){x[[2]]}
})

#Estimate of parameters
species_prop10 <- lapply(species_estimates_na, function(x){
  if(class(x)!= "try-error"){x[[1]]}
})


#The exponent for the Hills index
q <- as.list(seq(0,1,0.2 ))

#Estimating the Hills index for the three approaches
hills_inter10<- pblapply(inter_prop10, function(x){
  pblapply(q, function(z){
    if(!is.null(x)){
      data <- hill_index(z, x)
    }
  }, cl=1)
}, cl=1)

hills_genus10<- pblapply(genus_prop10, function(x){
  pblapply(q, function(z){
    if(!is.null(x)){
      data <- hill_index(z, x)
    }
  }, cl=1)
}, cl=1)

hills_species10<- pblapply(species_prop10, function(x){
  pblapply(q, function(z){
    if(!is.null(x)){
      data <- hill_index(z, x)
    }
  }, cl=1)
}, cl=1)

#Loading data for N.species =20
load("idm_miss_na_20/estimate_inter_na.RData")
inter1 <- inter_estimates_na
load("idm_miss_na_20/estimate_inter_na1.RData")
inter2 <- inter_estimates_na
inter_estimates_na <- c(inter1, inter2)
load("idm_miss_na_20/estimate_genus_na.RData")
genus1 <- genus_estimates_na
load("idm_miss_na_20/estimate_genus_na1.RData")
genus2 <- genus_estimates_na
genus_estimates_na <- c(genus1, genus2)
load("idm_miss_na_20/estimate_species_na.RData")

##IDM approached used
# Pi's
inter_prop20 <- lapply(inter_estimates_na, function(x){
  if(class(x)!= "try-error"){x[[1]]}
})

#parameters
inter_est20<- lapply(inter_estimates_na, function(x){
  if(class(x)!= "try-error"){x[[2]]}
})

#Genus only approach
#pi's
genus_prop20 <- lapply(genus_estimates_na, function(x){
  if(class(x)!= "try-error"){x[[1]]}
})

#parameters
genus_est20 <- lapply(genus_estimates_na, function(x){
  if(class(x)!= "try-error"){x[[2]]}
})

#Species only 
#parameters
species_est20 <- lapply(species_estimates_na, function(x){
  if(class(x)!= "try-error"){x[[2]]}
})

#pi's
species_prop20 <- lapply(species_estimates_na, function(x){
  if(class(x)!= "try-error"){x[[1]]}
})


hills_inter20<- pblapply(inter_prop20, function(x){
  pblapply(q, function(z){
    if(!is.null(x)){
      data <- hill_index(z, x)
    }
  }, cl=1)
}, cl=1)

hills_genus20<- pblapply(genus_prop20, function(x){
  pblapply(q, function(z){
    if(!is.null(x)){
      data <- hill_index(z, x)
    }
  }, cl=1)
}, cl=1)

hills_species20<- pblapply(species_prop20, function(x){
  pblapply(q, function(z){
    if(!is.null(x)){
      data <- hill_index(z, x)
    }
  }, cl=1)
}, cl=1)

#50 species
load("idm_miss_na_50/estimate_inter_na.RData")
inter1 <- inter_estimates_na
load("idm_miss_na_50/estimate_inter_na1.RData")
inter2 <- inter_estimates_na
load("idm_miss_na_50/estimate_inter_na2.RData")
inter3 <- inter_estimates_na
load("idm_miss_na_50/estimate_inter_na3.RData")
inter4 <- inter_estimates_na
load("idm_miss_na_50/estimate_inter_na4.RData")
inter5 <- inter_estimates_na
load("idm_miss_na_50/estimate_inter_na5.RData")
inter6 <- inter_estimates_na
load("idm_miss_na_50/estimate_inter_na6.RData")
inter7 <- inter_estimates_na
load("idm_miss_na_50/estimate_inter_na7.RData")
inter8 <- inter_estimates_na
load("idm_miss_na_50/estimate_inter_na8.RData")
inter9 <- inter_estimates_na
load("idm_miss_na_50/estimate_inter_na9.RData")
inter10 <- inter_estimates_na
inter_estimates_na <- c(inter1, inter2, inter3, inter4, inter5, 
                        inter6, inter7, inter8, inter9, inter10)


load("idm_miss_na_50/estimate_genus_na.RData")
genus1 <- genus_estimates_na
load("idm_miss_na_50/estimate_genus_na1.RData")
genus2 <- genus_estimates_na
load("idm_miss_na_50/estimate_genus_na2.RData")
genus3 <- genus_estimates_na
load("idm_miss_na_50/estimate_genus_na3.RData")
genus4 <- genus_estimates_na
load("idm_miss_na_50/estimate_genus_na4.RData")
genus5 <- genus_estimates_na
load("idm_miss_na_50/estimate_genus_na5.RData")
genus6 <- genus_estimates_na
load("idm_miss_na_50/estimate_genus_na6.RData")
genus7 <- genus_estimates_na
load("idm_miss_na_50/estimate_genus_na7.RData")
genus8 <- genus_estimates_na
load("idm_miss_na_50/estimate_genus_na8.RData")
genus9 <- genus_estimates_na
load("idm_miss_na_50/estimate_genus_na9.RData")
genus10 <- genus_estimates_na
genus_estimates_na <- c(genus1, genus2, genus3, genus4, genus5, 
                        genus6, genus7, genus8, genus9, genus10)

load("idm_miss_na_50/estimate_species_na.RData")
species1 <- species_estimates_na
load("idm_miss_na_50/estimate_species_na1.RData")
species2 <- species_estimates_na
load("idm_miss_na_50/estimate_species_na2.RData")
species3 <- species_estimates_na
load("idm_miss_na_50/estimate_species_na3.RData")
species4 <- species_estimates_na
load("idm_miss_na_50/estimate_species_na4.RData")
species5 <- species_estimates_na
load("idm_miss_na_50/estimate_species_na5.RData")
species6 <- species_estimates_na
load("idm_miss_na_50/estimate_species_na6.RData")
species7 <- species_estimates_na
load("idm_miss_na_50/estimate_species_na7.RData")
species8 <- species_estimates_na
load("idm_miss_na_50/estimate_species_na8.RData")
species9 <- species_estimates_na
load("idm_miss_na_50/estimate_species_na9.RData")
species10 <- species_estimates_na
species_estimates_na <- c(species1, species2, species3, species4, species5, 
                        species6, species7, species8, species9, species10)


inter_prop50 <- lapply(inter_estimates_na, function(x){
  if(class(x)!= "try-error"){x[[1]]}
})


inter_est50 <- lapply(inter_estimates_na, function(x){
  if(class(x)!= "try-error"){x[[2]]}
})

genus_prop50 <- lapply(genus_estimates_na, function(x){
  if(class(x)!= "try-error"){x[[1]]}
})


genus_est50 <- lapply(genus_estimates_na, function(x){
  if(class(x)!= "try-error"){x[[2]]}
})


species_est50 <- lapply(species_estimates_na, function(x){
  if(class(x)!= "try-error"){x[[2]]}
})

species_prop50 <- lapply(species_estimates_na, function(x){
  if(class(x)!= "try-error"){x[[1]]}
})




hills_inter50<- pblapply(inter_prop50, function(x){
  pblapply(q, function(z){
    if(!is.null(x)){
      data <- hill_index(z, x)
    }
  }, cl=1)
}, cl=1)

hills_genus50 <- pblapply(genus_prop50, function(x){
  pblapply(q, function(z){
    if(!is.null(x)){
      data <- hill_index(z, x)
    }
  }, cl=1)
}, cl=1)

hills_species50 <- pblapply(species_prop50, function(x){
  pblapply(q, function(z){
    if(!is.null(x)){
      data <- hill_index(z, x)
    }
  }, cl=1)
}, cl=1)



####################
# TRUE VALUES
####################

# True estimates of hills numbers
load("idm_miss_na_10/sim_interractions_na.RData")
true_sim_10 <- simulations_all_na[1:100]
load("idm_miss_na_20/sim_interractions_na.RData")
true_sim_20 <- simulations_all_na[1:100]
load("idm_miss_na_50/sim_interractions_na.RData")
true_sim_50 <- simulations_all_na


true_sim_prop10 <- lapply(true_sim_10, function(x){
  if(class(x)!= "try-error"){x[[3]]}
})

true_sim_prop20 <- lapply(true_sim_20, function(x){
  if(class(x)!= "try-error"){x[[3]]}
})

true_sim_prop50 <- lapply(true_sim_50, function(x){
  if(class(x)!= "try-error"){x[[3]]}
})



#Estimating hills numbers
hills_true_sim10<- pblapply(true_sim_prop10, function(x){
  pblapply(q, function(z){
    if(!is.null(x)){
      data <- hill_index(z, x)
    }
  }, cl=1)
}, cl=1)

hills_true_sim20<- pblapply(true_sim_prop20, function(x){
  pblapply(q, function(z){
    if(!is.null(x)){
      data <- hill_index(z, x)
    }
  }, cl=1)
}, cl=1)

hills_true_sim50<- pblapply(true_sim_prop50, function(x){
  pblapply(q, function(z){
    if(!is.null(x)){
      data <- hill_index(z, x)
    }
  }, cl=1)
}, cl=1)



#Function to put the data togehter
data_melt1 <- function(hills_inter){
  data_new_inter <-  list()
  for(j in 1: length(hills_inter)){
    data_inter <- matrix(NA, nrow=75, ncol=length(q))
    for(i in 1:length(q)){
      if(is.null(hills_inter[[j]][[i]])){
        data_inter[,i] <- rep(NA, 75)
      }else{
        data_inter[,i] <- hills_inter[[j]][[i]]
      }
    }
    data_inter <- cbind(site=seq(1:75),as.data.frame(data_inter))
    data_new_inter[[j]] <- melt(data_inter, id.vars = c("site"))
  }
  return(data_new_inter)
}

#Function to put the data togehter


#Comnbining the all number of species data
hills_inter <- c(hills_inter10, hills_inter20,
                 hills_inter50)
hills_species <- c(hills_species10, hills_species20,
                   hills_species50)
hills_genus <- c(hills_genus10, hills_genus20,
                 hills_genus50)

hills_true <- c(hills_true_sim10, hills_true_sim20,
                hills_true_sim50)

data_melt <- function(hills_inter){
  data_new_inter <-list()
  for(j in 1: length(hills_inter)){
    data_inter <-  matrix(NA, nrow=75, ncol=length(q))
    for(i in 1:length(q)){
      if(is.null(hills_inter[[j]][[i]])){
        data_inter[,i] <- rep(NA, 75)
      }else{
        data_inter[,i] <- hills_inter[[j]][[i]]
      }
    }
    
    
    data_inter <- cbind(site=seq(1:75),as.data.frame(data_inter))
    
    
    data_new_inter[[j]] <- melt(data_inter, id.vars = c("site"))
    
  }
  return(data_new_inter)
}



data_inter <- data_melt(hills_inter)
data_species <- data_melt(hills_species)
data_genus <- data_melt(hills_genus)
data_true <- data_melt(hills_true)
data <- c(data_inter, data_species, data_genus, data_true)

data_flat <- do.call("rbind", data)
n <- (dim(data_flat)[1]/(nspecies+1)) #Each group has 240 simulated results with 75 sites and q hills exponents
n1 <- niter * nsite*length(q) #The numbers for each group
groups <- c(rep("IDM",n ), rep("Species Only",n ),rep("Insect Group only",n ), rep("Truth",n))
sim_group <- rep(c(rep("N.SP=10", n1), rep("N.SP=20", n1),
                   rep("N.SP= 50", n1)),(nspecies+1))

all_data <- cbind(data_flat, groups=groups, sim_group=sim_group)



all_data1 <- all_data%>%
  dplyr::filter(site==75)

ggplot()+
  geom_boxplot(all_data1,mapping = aes(x=variable, y=value, fill=groups))+
  theme_bw()+
  ylab("Hill's numbers")+
  xlab("q")+
  ylim(c(0,50))+
  theme(legend.position = c(0.95,1))+
  scale_fill_manual(values = c("#FC4E07","#00AFBB", "#E7B800", "black"),name="Method")+
  scale_x_discrete(labels=c("0","0.2","0.4",
                            "0.6","0.8","1"))+
  facet_wrap( ~site+sim_group, ncol=2)

#Putting all the pi's together
inter_prop <- c(inter_prop10, inter_prop20,inter_prop50) #NB: med if for the median and the prop is for the mean
genus_prop <- c(genus_prop10, genus_prop20, genus_prop50)
species_prop <- c(species_prop10, species_prop20, species_prop50)
true_sim <- c(true_sim_prop10 ,true_sim_prop20 ,
              true_sim_prop50)

#################
# AED
################
aed_inter<- pblapply(inter_prop, function(x){
  if(!is.null(x)){
    data <- aed( x)
  }
}, cl=1)

aed_genus<- pblapply(genus_prop, function(x){
  if(!is.null(x)){
    data <- aed(x)
  }
}, cl=1)

aed_species<- pblapply(species_prop, function(x){
  if(!is.null(x)){
    data <- aed(x)
  }
}, cl=1)

#AED of the true values
aed_true <-pblapply(true_sim, function(x){
  if(!is.null(x)){
    data <- aed(x)
  }
}, cl=1)

data_inter <- data_melt(aed_inter)
data_genus <- data_melt(aed_genus)
data_species <- data_melt(aed_species)
data_true <- data_melt(aed_true)
data <- c(data_inter, data_species, data_genus, data_true)
data_flat <- do.call("rbind", data)
n <-  (dim(data_flat)[1]/(nspecies+1))
n1 <- (niter)*nsite*length(q)
groups <- c(rep("IDM",n ), rep("Species Only",n ),rep("Insect Group only",n ), rep("Truth",n))
sim_group <- rep(c(rep("N.SP=10", n1), rep("N.SP=20", n1),
                   rep("N.SP= 50", n1)),(nspecies+1))
all_data1 <- cbind(data_flat, groups=groups, sim_group=sim_group)
all_data1%>%
  filter(site == 10)%>%
  ggplot()+
  geom_boxplot(mapping = aes(x=as.factor(sim_group), y=as.numeric(value), fill=as.factor(groups)))+
  theme_bw()+
  ylab("AED Index")+
  xlab("Number of species")+
  theme(legend.position = "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  #scale_fill_discrete(name="Method")+
  scale_fill_manual(values = c("#FC4E07","#00AFBB", "#E7B800","black"),name="Method")+
  facet_wrap(~as.factor(site))



#N.species =20
#all_data1%>%
#  filter(site <= 2)%>%

all_data_nspecies_20 <- all_data1%>%
  filter(site == 75)%>%
  filter(sim_group == "N.SP=20"| sim_group =="N.SP= 20 (MI)")

ggplot()+
  geom_boxplot(data=all_data_nspecies_20,mapping = aes(x=as.factor(sim_group), y=as.numeric(value), fill=as.factor(groups)))+
  theme_bw()+
  ylab("AED Index")+
  xlab("Number of species")+
  theme(legend.position = "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  #scale_fill_discrete(name="Method")+
  scale_fill_manual(values = c("#FC4E07","#00AFBB", "#E7B800","black"),name="Method")+
  facet_wrap(~as.factor(site))

#ggarrange(gg1,gg2, nrow=2,common.legend=TRUE, legend = "right" )

#################################
# Shannon index and eveness
##################################
###Shannon index and eveness
shannon_inter<- pblapply(inter_prop, function(x){
  if(!is.null(x)){
    data <- shan_index( x)
  }
}, cl=1)

shannon_genus<- pblapply(genus_prop, function(x){
  if(!is.null(x)){
    data <- shan_index(x)
  }
}, cl=1)

shannon_species<- pblapply(species_prop, function(x){
  if(!is.null(x)){
    data <- shan_index(x)
  }
}, cl=1)


data_melt <- function(hills_inter11, hills_species11, hills_genus11){
  data_new_inter <- data_new_species <- data_new_genus <- list()
  for(j in 1: length(hills_inter11)){
    data_inter <- data_species <-data_genus <- vector("numeric", length=75)
    
    if(is.null(hills_inter11[[j]])){
      data_inter <- rep(NA, 75)
    }else{
      data_inter<- hills_inter11[[j]]
    }
    
    #Species
    if(is.null(hills_species11[[j]])){
      data_species <- rep(NA, 75)
    }else{
      data_species <- hills_species11[[j]]
    }
    
    
    #genus
    if(is.null(hills_genus11[[j]])){
      data_genus <- rep(NA, 75)
    }else{
      data_genus<- hills_genus11[[j]]
    }
    
    
    data_inter1 <- data.frame(site=seq(1:75),value=data_inter)
    data_species1 <- data.frame(site=seq(1:75),value=data_species)
    data_genus1 <- data.frame(site=seq(1:75),value=data_genus)
    
    
    data_new_inter[[j]] <- data_inter1
    data_new_species[[j]] <- data_species1
    data_new_genus[[j]] <- data_genus1
    #gg[[j]] <- ggplot(data = data_new)+
    #geom_line(aes(x, value, col=variable))
  }
  return(list(data_new_inter, data_new_species, data_new_genus))
}

#Shannon index
data <- data_melt(shannon_inter, shannon_species, shannon_genus)
data_flat <- flatten(data)
n <- (length(data_flat)/(nspecies))*nsite
n1 <- (niter)*nsite
groups <- c(rep("IDM",n ), rep("Species Only",n ),rep("Insect Group only",n ))
sim_group <- rep(c(rep("N.SP=10", n1), rep("N.SP=20", n1),
                   rep("N.SP= 50", n1)),(nspecies))

all_data1 <- cbind(do.call("rbind", data_flat), groups=groups, sim_group=sim_group)

all_data1%>%
  filter(site <= 4)%>%
  ggplot()+
  geom_boxplot(mapping = aes(x=as.factor(sim_group), y=as.numeric(value), fill=as.factor(groups)))+
  theme_bw()+
  ylab("H'")+
  xlab("Number of species")+
  theme(legend.position = "right",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  #scale_fill_discrete(name="Method")+
  scale_fill_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),name="Method")+
  facet_wrap(~as.factor(site))

############################
#   MSE and MBIAS
#########################

data_mse <- function(hills_inter1, hills_species1, hills_genus1,true_est){
  data_new_inter <- data_new_species <- data_new_genus <- list()
  for(j in 1: length(hills_inter1)){
    data_inter <- data_species <-data_genus <- vector("numeric", length=75)
    
    if(is.null(hills_inter1[[j]])){
      data_inter <- rep(NA, 75)
    }else{
      data_inter<- hills_inter1[[j]]
    }
    
    #Species
    if(is.null(hills_species1[[j]])){
      data_species <- rep(NA, 75)
    }else{
      data_species <- hills_species1[[j]]
    }
    
    
    #genus
    if(is.null(hills_genus1[[j]])){
      data_genus <- rep(NA, 75)
    }else{
      data_genus<- hills_genus1[[j]]
    }
    
    
    #data_inter <- data.frame(site=seq(1:75),value=data_inter)
    #data_species <- data.frame(site=seq(1:75),value=data_species)
    #data_genus <- data.frame(site=seq(1:75),value=data_genus)
    
    
    data_new_inter[[j]] <- mse(true_est[[j]],data_inter)
    data_new_species[[j]] <- mse(true_est[[j]],data_species)
    data_new_genus[[j]] <- mse(true_est[[j]],data_genus)
    #gg[[j]] <- ggplot(data = data_new)+
    #geom_line(aes(x, value, col=variable))
  }
  return(list(data_new_inter, data_new_species, data_new_genus))
}


true_shan <-pblapply(true_sim, function(x){
  if(!is.null(x)){
    data <- shan_index( x)
  }
}, cl=1)


true_aed <- pblapply(true_sim, function(x){
  if(!is.null(x)){
    data <- aed( x)
  }
}, cl=1)

#Estimates of MSE and Mbias
shannon_MSE <- data_mse(shannon_inter, shannon_species, shannon_genus, true_shan)%>%
  flatten()%>%
  unlist()
aed_MSE <- data_mse(aed_inter, aed_species, aed_genus, true_aed)%>%
  flatten()%>%
  unlist()


n <- niter*3 #Number in each group
n1 <-niter
groups <- c(rep("IDM",n ), rep("Species only",n ),rep("Insect Group only",n ))
sim_group <- rep(c(rep("N.SP=10", n1), rep("N.SP=20", n1),
                   rep("N.SP= 50", n1)),nspecies)
#sim_group <- rep(c(rep("sim1", 20), rep("sim2", 20), rep("sim3", 20)),3)
est_data <- data.frame(groups, sim_group, shannon_MSE, aed_MSE)

est_data_melt <-est_data%>%
  select(groups,sim_group, shannon_MSE,aed_MSE)#%>%
#melt(id.vars=c("groups", "sim_group"))


gg_shan_mse <- ggplot(est_data_melt)+
  geom_boxplot(aes(x=sim_group, y=shannon_MSE, fill=groups))+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(legend.position = "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7))+
  scale_fill_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),name="Method")#+
#ylim(c(0,500))

gg_aed_mse <- ggplot(est_data_melt)+
  geom_boxplot(aes(x=sim_group, y=aed_MSE, fill=groups))+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(legend.position = "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7))+
  scale_fill_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),name="Method")#+
#ylim(c(0,500))

ggpubr::ggarrange( gg_shan_mse, gg_aed_mse,
           ncol=2, 
           nrow = 1, 
           common.legend = TRUE, 
           legend = "bottom",
           labels = c("(a) RMSE of H'", "(b) RMSE of AED"))


#####################################################
#   Estimates
####################################################
subsetting_parameters <- function(x,name){
  ret <- x[grep(name, names(x))]
  return(ret)
}

########
# ALPHA 0
########





#######################
#Estimates data
#####################

###################
# NO MISSING DATA
###################
load("idm_miss_na_10/sim_input_na.RData")
inputs <- list(c(input_list_na[[1]]$parameters$detection$alpha0))
#inputs <- list(inputs1)
true_est10 <- rep(inputs, each=niter)
inter_est_mse_10 <- Map(mse, true_est10, lapply(inter_est10, function(x){subsetting_parameters(x, "alpha0")}))
genus_est_mse_10 <- Map(mse, true_est10, lapply(genus_est10, function(x){subsetting_parameters(x, "alpha0")}))
species_est_mse_10 <- Map(mse, true_est10, lapply(species_est10, function(x){subsetting_parameters(x, "alpha0")}))
estimates_mse10 <- c(inter_est_mse_10,genus_est_mse_10,species_est_mse_10)

load("idm_miss_na_20/sim_input_na.RData")
inputs <- list(c(input_list_na[[1]]$parameters$detection$alpha0))
true_est20 <- rep(inputs, each=niter)
inter_est_mse_20 <- Map(mse, true_est20, lapply(inter_est20, function(x){subsetting_parameters(x, "alpha0")}))
genus_est_mse_20 <- Map(mse, true_est20, lapply(genus_est20, function(x){subsetting_parameters(x, "alpha0")}))
species_est_mse_20 <- Map(mse, true_est20, lapply(species_est20, function(x){subsetting_parameters(x, "alpha0")}))
estimates_mse20 <- c(inter_est_mse_20,genus_est_mse_20,species_est_mse_20)


load("idm_miss_na_50/sim_input_na.RData")
inputs <- list(c(input_list_na[[1]]$parameters$detection$alpha0))
true_est50 <- rep(inputs, each=niter)
inter_est_mse_50 <- Map(mse, true_est50, lapply(inter_est50, function(x){subsetting_parameters(x, "alpha0")}))
genus_est_mse_50 <- Map(mse, true_est50, lapply(genus_est50, function(x){subsetting_parameters(x, "alpha0")}))
species_est_mse_50 <- Map(mse, true_est50, lapply(species_est50, function(x){subsetting_parameters(x, "alpha0")}))
estimates_mse50 <- c(inter_est_mse_50,genus_est_mse_50,species_est_mse_50)

all_estimates_mse <- c(estimates_mse10,estimates_mse20,estimates_mse50)
n_est <- niter #Number of results for each inter etc. 
n1 <- niter*3
n <- niter
groups <- rep(c(rep("IDM",n), rep("Insect group only",n ),rep("Species only",n )),3)
sim_group <- rep(c(rep("N.SP=10", n1),  
                   rep("N.SP=20", n1),
                   rep("N.SP=50", n1) ),1)

all_data_est <- data.frame(cbind(do.call("rbind", all_estimates_mse),
                                 groups=groups, sim_group=sim_group))
colnames(all_data_est) <- c("mse", "groups", "sim_group")


gg2 <-all_data_est%>%
  filter(groups != "Insect group only")%>%
  ggplot()+
  geom_boxplot(mapping=aes(x=as.factor(sim_group), y= as.numeric(mse), fill=as.factor(groups)))+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(legend.position = "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7))+
  #scale_fill_discrete(name="Method")+
  scale_fill_manual(values = c("#00AFBB", "#E7B800"),name="Method")#+
#ylim(c(1.25,1.75))
#ggarrange(gg1,gg2, nrow=1, ncol=2, common.legend = TRUE, legend="bottom")

########
# ALPHA 1
########
load("idm_miss_na_10/sim_input_na.RData")
inputs <- list(c(input_list_na[[1]]$parameters$detection$alpha1))
#inputs <- list(inputs1)
true_est10 <- rep(inputs, each=niter)
inter_est_mse_10 <- Map(mse, true_est10, lapply(inter_est10, function(x){subsetting_parameters(x, "alpha1")}))
genus_est_mse_10 <- Map(mse, true_est10, lapply(genus_est10, function(x){subsetting_parameters(x, "alpha1")}))
species_est_mse_10 <- Map(mse, true_est10, lapply(species_est10, function(x){subsetting_parameters(x, "alpha1")}))
estimates_mse10 <- c(inter_est_mse_10,genus_est_mse_10,species_est_mse_10)

load("idm_miss_na_20/sim_input_na.RData")
inputs <- list(c(input_list_na[[1]]$parameters$detection$alpha1))
true_est20 <- rep(inputs, each=niter)
inter_est_mse_20 <- Map(mse, true_est20, lapply(inter_est20, function(x){subsetting_parameters(x, "alpha1")}))
genus_est_mse_20 <- Map(mse, true_est20, lapply(genus_est20, function(x){subsetting_parameters(x, "alpha1")}))
species_est_mse_20 <- Map(mse, true_est20, lapply(species_est20, function(x){subsetting_parameters(x, "alpha1")}))
estimates_mse20 <- c(inter_est_mse_20,genus_est_mse_20,species_est_mse_20)


load("idm_miss_na_50/sim_input_na.RData")
inputs <- list(c(input_list_na[[1]]$parameters$detection$alpha1))
true_est50 <- rep(inputs, each=niter)
inter_est_mse_50 <- Map(mse, true_est50, lapply(inter_est50, function(x){subsetting_parameters(x, "alpha1")}))
genus_est_mse_50 <- Map(mse, true_est50, lapply(genus_est50, function(x){subsetting_parameters(x, "alpha1")}))
species_est_mse_50 <- Map(mse, true_est50, lapply(species_est50, function(x){subsetting_parameters(x, "alpha1")}))
estimates_mse50 <- c(inter_est_mse_50,genus_est_mse_50,species_est_mse_50)

all_estimates_mse <- c(estimates_mse10,estimates_mse20,estimates_mse50)
n_est <- niter #Number of results for each inter etc. 
n1 <- niter*3
n <- niter
groups <- rep(c(rep("IDM",n), rep("Insect group only",n ),rep("Species only",n )),3)
sim_group <- rep(c(rep("N.SP=10", n1),  
                   rep("N.SP=20", n1),
                   rep("N.SP=50", n1) ),1)

all_data_est <- data.frame(cbind(do.call("rbind", all_estimates_mse),
                                 groups=groups, sim_group=sim_group))
colnames(all_data_est) <- c("mse", "groups", "sim_group")


gg4 <-all_data_est%>%
  filter(groups != "Insect group only")%>%
  ggplot()+
  geom_boxplot(mapping=aes(x=as.factor(sim_group), y= as.numeric(mse), fill=as.factor(groups)))+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(legend.position = "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7))+
  #scale_fill_discrete(name="Method")+
  scale_fill_manual(values = c("#00AFBB", "#E7B800"),name="Method")#+
#ylim(c(1.25,1.75))
#ggarrange(gg1,gg2, nrow=1, ncol=2, common.legend = TRUE, legend="bottom")


########
# BETA0
########
load("idm_miss_na_10/sim_input_na.RData")
inputs <- list(c(input_list_na[[1]]$parameters$ecological$beta0))
true_est10 <- rep(inputs, each=niter)

inter_est_mse_10 <- Map(mse, true_est10, lapply(inter_est10, function(x){subsetting_parameters(x, "beta0")}))
genus_est_mse_10 <- Map(mse, true_est10, lapply(genus_est10, function(x){subsetting_parameters(x, "beta0")}))
species_est_mse_10 <- Map(mse, true_est10, lapply(species_est10, function(x){subsetting_parameters(x, "beta0")}))
estimates_mse10 <- c(inter_est_mse_10,genus_est_mse_10,species_est_mse_10)

load("idm_miss_na_20/sim_input_na.RData")
inputs <- list(c(input_list_na[[1]]$parameters$ecological$beta0))
true_est20 <- rep(inputs, each=niter)
inter_est_mse_20 <- Map(mse, true_est20, lapply(inter_est20, function(x){subsetting_parameters(x, "beta0")}))
genus_est_mse_20 <- Map(mse, true_est20, lapply(genus_est20, function(x){subsetting_parameters(x, "beta0")}))
species_est_mse_20 <- Map(mse, true_est20, lapply(species_est20, function(x){subsetting_parameters(x, "beta0")}))
estimates_mse20 <- c(inter_est_mse_20,genus_est_mse_20,species_est_mse_20)


load("idm_miss_na_50/sim_input_na.RData")
inputs <- list(c(input_list_na[[1]]$parameters$ecological$beta0))
true_est50 <- rep(inputs, each=niter)
inter_est_mse_50 <- Map(mse, true_est50, lapply(inter_est50, function(x){subsetting_parameters(x, "beta0")}))
genus_est_mse_50 <- Map(mse, true_est50, lapply(genus_est50, function(x){subsetting_parameters(x, "beta0")}))
species_est_mse_50 <- Map(mse, true_est50, lapply(species_est50, function(x){subsetting_parameters(x, "beta0")}))
estimates_mse50 <- c(inter_est_mse_50,genus_est_mse_50,species_est_mse_50)

all_estimates_mse <- c(estimates_mse10,estimates_mse20,estimates_mse50)
n_est <- niter #Number of results for each inter etc. 
n1 <- niter*3
n <- niter
groups <- rep(c(rep("IDM",n), rep("Insect group only",n ),rep("Species only",n )),3)
sim_group <- rep(c(rep("N.SP=10", n1),  
                   rep("N.SP=20", n1),
                   rep("N.SP=50", n1) ),1)

all_data_est <- data.frame(cbind(do.call("rbind", all_estimates_mse),
                                 groups=groups, sim_group=sim_group))
colnames(all_data_est) <- c("mse", "groups", "sim_group")

gg6 <-ggplot(all_data_est)+
  geom_boxplot(mapping=aes(x=as.factor(sim_group), y= as.numeric(mse), fill=as.factor(groups)))+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(legend.position = "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7))+
  scale_fill_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),name="Method")#+
#  ylim(c(1.25,3))

########
# BETA1
########
load("idm_miss_na_10/sim_input_na.RData")
inputs <- list(c(input_list_na[[1]]$parameters$ecological$beta1))
true_est10 <- rep(inputs, each=niter)

inter_est_mse_10 <- Map(mse, true_est10, lapply(inter_est10, function(x){subsetting_parameters(x, "beta1")}))
genus_est_mse_10 <- Map(mse, true_est10, lapply(genus_est10, function(x){subsetting_parameters(x, "beta1")}))
species_est_mse_10 <- Map(mse, true_est10, lapply(species_est10, function(x){subsetting_parameters(x, "beta1")}))
estimates_mse10 <- c(inter_est_mse_10,genus_est_mse_10,species_est_mse_10)

load("idm_miss_na_20/sim_input_na.RData")
inputs <- list(c(input_list_na[[1]]$parameters$ecological$beta1))
true_est20 <- rep(inputs, each=niter)
inter_est_mse_20 <- Map(mse, true_est20, lapply(inter_est20, function(x){subsetting_parameters(x, "beta1")}))
genus_est_mse_20 <- Map(mse, true_est20, lapply(genus_est20, function(x){subsetting_parameters(x, "beta1")}))
species_est_mse_20 <- Map(mse, true_est20, lapply(species_est20, function(x){subsetting_parameters(x, "beta1")}))
estimates_mse20 <- c(inter_est_mse_20,genus_est_mse_20,species_est_mse_20)


load("idm_miss_na_50/sim_input_na.RData")
inputs <- list(c(input_list_na[[1]]$parameters$ecological$beta1))
true_est50 <- rep(inputs, each=niter)
inter_est_mse_50 <- Map(mse, true_est50, lapply(inter_est50, function(x){subsetting_parameters(x, "beta1")}))
genus_est_mse_50 <- Map(mse, true_est50, lapply(genus_est50, function(x){subsetting_parameters(x, "beta1")}))
species_est_mse_50 <- Map(mse, true_est50, lapply(species_est50, function(x){subsetting_parameters(x, "beta1")}))
estimates_mse50 <- c(inter_est_mse_50,genus_est_mse_50,species_est_mse_50)

all_estimates_mse <- c(estimates_mse10,estimates_mse20,estimates_mse50)
n_est <- niter #Number of results for each inter etc. 
n1 <- niter*3
n <- niter
groups <- rep(c(rep("IDM",n), rep("Insect group only",n ),rep("Species only",n )),3)
sim_group <- rep(c(rep("N.SP=10", n1),  
                   rep("N.SP=20", n1),
                   rep("N.SP=50", n1) ),1)

all_data_est <- data.frame(cbind(do.call("rbind", all_estimates_mse),
                                 groups=groups, sim_group=sim_group))
colnames(all_data_est) <- c("mse", "groups", "sim_group")
gg8 <-ggplot(all_data_est)+
  geom_boxplot(mapping=aes(x=as.factor(sim_group), y= as.numeric(mse), fill=as.factor(groups)))+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(legend.position = "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7))+
  scale_fill_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),name="Method")#+
#ylim(c(1.35,2))
ggpubr::ggarrange(gg2,gg4,gg6,gg8, nrow=2, ncol=2, common.legend = TRUE, legend="bottom")

###################
# All plots together
################
figure <- ggpubr::ggarrange(gg6,gg8,gg2, gg4, gg_shan_mse, gg_aed_mse, nrow=2, ncol=4, common.legend = TRUE, legend="top",
                    labels = c("(a) beta0", "(b) beta1",
                               "(c) ahpha0", "(d) alpha1",
                               "(e) H'", "(f) AED"),
                    font.label=list(color="black",size=9), 
                    hjust = -1, vjust = -0.3)
annotate_figure(figure,
                bottom = text_grob("Number of species", 
                                   hjust = 1),
                left = text_grob("Root Mean Square Error", rot = 90)
)
ggsave("metrics_hills.png")

##########
#Comnbining the two number of species


hills_mean <- function(hills_inter){
  data_new_inter <-list()
  for(j in 1: length(hills_inter)){
    data_inter <-  matrix(NA, nrow=75, ncol=length(q))
    for(i in 1:length(q)){
      if(is.null(hills_inter[[j]][[i]])){
        data_inter[,i] <- rep(NA, 75)
      }else{
        data_inter[,i] <- hills_inter[[j]][[i]]
      }
    }
    
    
    # data_inter <- cbind(site=seq(1:75),as.data.frame(data_inter))
    
    
    data_new_inter[[j]] <- colMeans(data_inter)
    
  }
  return(data_new_inter)
}



data_inter <- hills_mean(hills_inter)
data_species <- hills_mean(hills_species)
data_genus <- hills_mean(hills_genus)
data_true <- hills_mean(hills_true)
data <- c(data_inter, data_species, data_genus, data_true)

data_flat <- do.call("rbind", data)
n <- (dim(data_flat)[1]/(nspecies+1)) #Each group has 240 simulated results with 75 sites and q hills exponents
n1 <- niter #The numbers for each group
groups <- c(rep("IDM",n ), rep("Species Only",n ),rep("Insect Group only",n ), rep("Truth",n))
sim_group <- rep(c(rep("N.SP=10", n1), rep("N.SP=20", n1),
                   rep("N.SP=50", n1)),(nspecies+1))
colnames(data_flat) <- c("Q0", "Q0.2","Q0.4","Q0.6","Q0.8","Q1")
iterations <- rep(seq(1, niter, 1),16)

#all_data <- cbind(data_flat, groups=groups, sim_group=sim_group, iterations=iterations)
#colnames(all_data) <- c("Q0", "Q0.2","Q0.4","Q0.6","Q0.8","Q1", groups, sim_group)

all_data_melt <- reshape2::melt(data_flat, id.vars=c("iterations"))

all_data_plot <- cbind(all_data_melt,
                       groups=rep(groups, 6), 
                       sim_group=rep(sim_group,6))

ggplot(data=all_data_plot)+
  geom_boxplot(mapping=aes(x=as.factor(Var2), y= value, fill=as.factor(groups)))+  
  ylab("Hill's numbers")+
  xlab("q")+
  theme_bw()+
  #theme(legend.position = c(0.95,1))+
  theme(legend.position = "right")+
  scale_fill_manual(values = c("#FC4E07","#00AFBB", "#E7B800","black"),name="Method")+
  # scale_color_manual(values = c("#FC4E07","#00AFBB", "#E7B800","black"),name="Method")+
  scale_x_discrete(labels=c("0","0.2","0.4",
                            "0.6","0.8", "1"))+
  facet_wrap(~factor(sim_group, levels=c("N.SP=10","N.SP=20","N.SP=50")),
             nrow=2, scales = "free_y")

ggsave("average_hills.png")


