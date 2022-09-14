#Need to set the working directory to the folder with the simulations
setwd("/Volumes/kwakupa/IDM_new/IDM_simulations_new")

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
library(scales)
require("LaplacesDemon")
source("/Volumes/kwakupa/IDM_new/idm_functions/fnx_for_estimation.R")

nspecies = 3; nsite=75; niter= 300

#load data
load("idm_miss_na_50/estimate_inter_na.RData")
inter1 <- inter_estimates_na
load("idm_miss_na_50/estimate_inter_na1.RData")
inter2 <- inter_estimates_na
load("idm_miss_na_50/estimate_inter_na2.RData")
 inter3 <- inter_estimates_na
 load("idm_miss_na_50/estimate_inter_na3.RData")
 inter4 <- inter_estimates_na
# load("idm_miss_na_50/estimate_inter_na4.RData")
# inter5 <- inter_estimates_na
inter_estimates_na <- c(inter1, inter2, inter3, inter4)

load("idm_miss_na_50/estimate_genus_na.RData")
genus1 <- genus_estimates_na
load("idm_miss_na_50/estimate_genus_na1.RData")
genus2 <- genus_estimates_na
 load("idm_miss_na_50/estimate_genus_na2.RData")
 genus3 <- genus_estimates_na
 load("idm_miss_na_50/estimate_genus_na3.RData")
 genus4 <- genus_estimates_na
# load("idm_miss_na_50/estimate_genus_na4.RData")
# genus5 <- genus_estimates_na
genus_estimates_na_nomiss <- c(genus1, genus2, genus3, genus4)#, genus5)

load("idm_miss_na_50/estimate_species_na.RData")
species1 <- species_estimates_na
load("idm_miss_na_50/estimate_species_na1.RData")
species2 <- species_estimates_na
load("idm_miss_na_50/estimate_species_na2.RData")
species3 <- species_estimates_na
 load("idm_miss_na_50/estimate_species_na3.RData")
 species4 <- species_estimates_na
# load("idm_miss_na_50/estimate_species_na4.RData")
# species5 <- species_estimates_na
species_estimates_na <- c(species1, species2, species3, species4)#, species5)

#all data together
all_data_nomiss <- c(inter_estimates_na, species_estimates_na)



# #load data
# load("idm_miss_50/estimate_inter_na.RData")
# inter1 <- inter_estimates_na
# load("idm_miss_50/estimate_inter_na1.RData")
# inter2 <- inter_estimates_na
# #load("idm_miss_50/estimate_inter_na2.RData")
# # inter3 <- inter_estimates_na
# # load("idm_miss_50/estimate_inter_na3.RData")
# # inter4 <- inter_estimates_na
# # load("idm_miss_50/estimate_inter_na4.RData")
# # inter5 <- inter_estimates_na
# inter_estimates_miss_na <- c(inter1, inter2)
# 
# load("idm_miss_50/estimate_genus_na.RData")
# genus1 <- genus_estimates_na
# load("idm_miss_50/estimate_genus_na1.RData")
# genus2 <- genus_estimates_na
# # load("idm_miss_50/estimate_genus_na2.RData")
# # genus3 <- genus_estimates_na
# # load("idm_miss_50/estimate_genus_na3.RData")
# # genus4 <- genus_estimates_na
# # load("idm_miss_50/estimate_genus_na4.RData")
# # genus5 <- genus_estimates_na
# genus_estimates_miss_na <- c(genus1, genus2)#, genus3, genus4, genus5)
# 
# load("idm_miss_50/estimate_species_na.RData")
# species1 <- species_estimates_na
# load("idm_miss_50/estimate_species_na1.RData")
# species2 <- species_estimates_na
# # load("idm_miss_50/estimate_species_na2.RData")
# # species3 <- species_estimates_na
# # load("idm_miss_50/estimate_species_na3.RData")
# # species4 <- species_estimates_na
# # load("idm_miss_50/estimate_species_na4.RData")
# # species5 <- species_estimates_na
# species_estimates_miss_na <- c(species1, species2)#, species3, species4, species5)
# 
# all_data_miss <- c(inter_estimates_miss_na, species_estimates_miss_na)
# 
# all_data <- c(all_data_nomiss, all_data_miss)
all_data <- c(all_data_nomiss, genus_estimates_na_nomiss)

#retrieve parameters
mean_richness <- lapply(all_data, function(x){
  if(class(x[[1]])[1]!= "try-error"){
    ret <- subsetting_parameters(x[[1]], "richness", "mean")
    #psi <- 1 - exp(-ret)
    #matrix(ret, nrow= nsite, ncol = 50, byrow =  FALSE)
    }
})

sd_richness <- lapply(all_data, function(x){
  if(class(x[[1]])[1]!= "try-error"){
    ret <- subsetting_parameters(x[[1]], "richness", "sd")
  }
})

mean_hills_index1 <- lapply(all_data, function(x){
  if(class(x[[1]])[1]!= "try-error"){
    ret <- subsetting_parameters(x[[1]], "hills_index1", "mean")
    #psi <- 1 - exp(-ret)
    #matrix(ret, nrow= nsite, ncol = 50, byrow =  FALSE)
  }
})

sd_hills_index1 <- lapply(all_data, function(x){
  if(class(x[[1]])[1]!= "try-error"){
    ret <- subsetting_parameters(x[[1]], "hills_index1", "sd")
  }
})

mean_hills_index2 <- lapply(all_data, function(x){
  if(class(x[[1]])[1]!= "try-error"){
    ret <- subsetting_parameters(x[[1]], "hills_index2", "mean")
    #psi <- 1 - exp(-ret)
    #matrix(ret, nrow= nsite, ncol = 50, byrow =  FALSE)
  }
})

sd_hills_index2 <- lapply(all_data, function(x){
  if(class(x[[1]])[1]!= "try-error"){
    ret <- subsetting_parameters(x[[1]], "hills_index2", "sd")
  }
})

estimated_hills_indices <- c(mean_richness, mean_hills_index1, mean_hills_index2)
sd_hills_indices <- c(sd_richness, sd_hills_index1, sd_hills_index2)


#Extract true values
load("/Volumes/kwakupa/IDM_new/IDM_simulations_new/idm_miss_na_50/sim_interractions_na.RData")
true_values_nomiss <- simulations_all_na[1:niter]  
true_values <- c(true_values_nomiss)

true_richness <- lapply(true_values, function(x){
  x$richness
})

true_incidence_hills1 <- lapply(true_values, function(x){
  x$incidence_hills1
})

true_incidence_hills2 <- lapply(true_values, function(x){
  x$incidence_hills2
})

true_abundance_hills0 <- lapply(true_values, function(x){
  x$abundance_hills0 
})

true_abundance_hills1 <- lapply(true_values, function(x){
  x$abundance_hills1
})

true_abundance_hills2 <- lapply(true_values, function(x){
  x$abundance_hills2
})

true_hills_indices <- c(true_richness, true_richness, true_abundance_hills0, 
                        true_incidence_hills1, true_incidence_hills1, true_abundance_hills1,
                        true_incidence_hills2 ,true_incidence_hills2 ,true_abundance_hills2)



#estimating mae

hills_mae <- Map(mae, true_hills_indices, estimated_hills_indices)%>%
  unlist()%>%
  data.frame%>%
  mutate(method = rep(rep(c("IDM", "Species Only", "Group Only"), each = niter),3), 
         index = rep(rep(0:2, each = niter*3)))%>%
  ggplot()+
  geom_boxplot(mapping=aes(x=as.factor(index), y= as.numeric(.), fill=as.factor(method)))+
  theme_bw()+
  ylab("MAE")+
  xlab("q")+
  theme(legend.position = "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7),plot.title = element_text(size = 9),
        axis.title.x = element_text(size = 8), axis.title.y =element_text(size = 8) )+
  ylim(c(0,10))+
  #scale_y_continuous(trans="log",breaks=seq(-0.1,10,0.5))+
  ggtitle( "a) MAE of Hills indices")+
  scale_fill_manual(values = c("#FC4E07","#00AFBB", "#E7B800" ),name="Method")



#Average hills indices across all sites
average_estimated_hills_index <- lapply(estimated_hills_indices, function(x){
  mean(x)
})

average_true_hills_index <- lapply(true_hills_indices, function(x){
  mean(x)
})

average_hills_indices <- c(average_estimated_hills_index, average_true_hills_index)%>%
  unlist()%>%
  data.frame()%>%
  mutate(method = c(rep(rep(c("IDM", "Species Only", "Group Only"), each = niter),3),
                    rep(rep(c("Truth_incidence", "Truth_incidence", "Truth_abundance"), each = niter),3)), 
         index = c(rep(rep(0:2, each = niter*3)), 
         rep(rep(0:2, each = niter*3))))%>%
  ggplot()+
  geom_boxplot(mapping=aes(x=as.factor(index), y= as.numeric(.), fill=as.factor(method)))+
  theme_bw()+
  ylab("Hills index")+
  xlab("q")+
  theme(legend.position = "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7),plot.title = element_text(size = 9),
        axis.title.x = element_text(size = 8), axis.title.y =element_text(size = 8) )+
  ylim(c(10,50))+
  ggtitle( "b) Average Hills indices")+
  scale_fill_manual(values = c("#FC4E07","#00AFBB", "#E7B800", "black", "blue" ),name="Method")

#mean_precision of estimates
mean_precision_hill_indices <- lapply(sd_hills_indices, function(x){
  mean(x)
})%>%
  unlist()%>%
  data.frame(.)%>%
  mutate(precision = 1/.,
         method = rep(rep(c("IDM", "Species Only", "Group Only"), each = niter),3), 
         index = rep(rep(0:2, each = niter*3)))%>%
  ggplot()+
  geom_boxplot(mapping=aes(x=as.factor(index), y= as.numeric(.), fill=as.factor(method)))+
  theme_bw()+
  ylab("Mean precision")+
  xlab("q")+
  theme(legend.position = "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7),plot.title = element_text(size = 9),
        axis.title.x = element_text(size = 8), axis.title.y =element_text(size = 8) )+
  ylim(c(0,4))+
  ggtitle( "c) Mean precion of \n Hill indices")+
  scale_fill_manual(values = c("#FC4E07","#00AFBB", "#E7B800" ),name="Method")
  

## beta0, beta1

load("idm_miss_na_50/sim_input_na.RData")
beta0 <- list(c(input_list_na[[1]]$parameters$ecological$beta0))
beta1 <- list(c(input_list_na[[1]]$parameters$ecological$beta1))
alpha0 <- list(c(input_list_na[[1]]$parameters$detection$alpha0))
alpha1 <- list(c(input_list_na[[1]]$parameters$detection$alpha1))
true_beta0 <- rep(beta0, each=niter*3)
true_beta1 <- rep(beta1, each=niter*3)
true_alpha0 = rep(alpha0, each=niter*3)
true_alpha1 = rep(alpha1, each=niter*3)


mae_beta0 <- Map(mae, true_beta0, lapply(all_data, function(x){subsetting_parameters(x[[1]], "beta0", "mean")}))%>%unlist()
mae_beta1 <- Map(mae, true_beta1, lapply(all_data, function(x){subsetting_parameters(x[[1]], "beta1", "mean")}))%>%unlist()
mae_intercept_lambda <- Map(mae, -4, lapply(all_data, function(x){subsetting_parameters(x[[1]], "intercept_lambda", "mean")}))%>%unlist()
mae_intercept_lambda[(niter+1): (niter*2)] <-   NA
mae_alpha0 <- Map(mae, true_alpha0, lapply(all_data, function(x){subsetting_parameters(x[[1]], "alpha0", "mean")}))%>%
  unlist()
  mae_alpha0[(niter*2+1): (niter*3)] <-   NA #take out the values for Group only model. The values returned for the Group Only are the prior values
mae_alpha1 <- Map(mae, true_alpha1, lapply(all_data, function(x){subsetting_parameters(x[[1]], "alpha1", "mean")}))%>%unlist()
mae_alpha1[(niter*2+1): (niter*3)] <-   NA #take out the values for Group only model. The values returned for the Group Only are the prior values


mae_estimates_parameters <- data.frame(estimates =  c(mae_beta0, mae_beta1, mae_intercept_lambda,
                                                      mae_alpha0, mae_alpha1),
                                       variables = c(rep("beta0", niter*3),
                                                     rep("beta1", niter*3),
                                                     rep("gamma", niter*3),
                                                     rep("alpha0", niter*3),
                                                     rep("alpha1", niter*3)),
                                       method = rep(rep(c("IDM", "Species Only", "Group Only"), each = niter),5))%>%
  ggplot()+
  geom_boxplot(mapping=aes(x=as.factor(variables), y= as.numeric(estimates), fill=as.factor(method)))+
  theme_bw()+
  ylab("MAE")+
  xlab("Parameters")+
  ggtitle( "d) MAE of parameters")+
  theme(legend.position = "bottom",axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5, size = 5),plot.title = element_text(size = 9),
        axis.title.x = element_text(size = 8), axis.title.y =element_text(size = 8) )+
  #scale_y_continuous(trans="log",breaks=seq(-0.1,10,0.5))+
  #ylim(c(0,2))+
  scale_fill_manual(values = c("#FC4E07","#00AFBB", "#E7B800" ),name="Method")
  

#Put all the plots together
sim_plot <- ggpubr::ggarrange(hills_mae,
                  average_hills_indices, 
                  mean_precision_hill_indices,
                  mae_estimates_parameters, nrow = 1, ncol = 4,
                  common.legend = TRUE)

ggsave("simulations.png", plot = sim_plot, dpi = 320, 
       height = 10, width = 18, units = "cm")


