setwd("/Volumes/kwakupa/Data_idm")

require(dplyr)
require(ggcorrplot)
require(purrr)
require(ggpubr)
require(heatmaply)
require(reshape2)
library(ggplot2)
require(egg)
require(readr)
if(!require(devtools)) install.packages("devtools")
if(!require(ggcorrplot)) devtools::install_github("kassambara/ggcorrplot")
if(!require(ggcorrplot)) install.packages("heatmaply")




#load data
load("bumblebees/data_idm.RData")

#load bumblebees results
load("bumblebees/estimate_inter_data.RData")
load("bumblebees/estimate_species_data.RData")
load("bumblebees/estimate_genus_data.RData")
bb_inter <- inter_estimates_data
bb_genus <- genus_estimates_data
bb_species <- species_estimates_data

#load solitary bees results
load("solitary bees/estimate_inter_data.RData")
load("solitary bees/estimate_species_data.RData")
load("solitary bees/estimate_genus_data.RData")
sb_inter <- inter_estimates_data
sb_genus <- genus_estimates_data
sb_species <- species_estimates_data

#load hoverflies results
load("hoverflies/estimate_inter_data.RData")
load("hoverflies/estimate_species_data.RData")
load("hoverflies/estimate_genus_data.RData")
hv_inter <- inter_estimates_data
hv_genus <- genus_estimates_data
hv_species <- species_estimates_data

all_data <- list(bb_inter, bb_species, bb_genus,
                 sb_inter, sb_species, sb_genus,
                 hv_inter, hv_species, hv_genus)

extract_covariance <- function(x,y){
  constants <- x[[1]]
  correlation_matrix <- matrix(x[[2]][grep("Cov", rownames(x[[2]])), 1], 
                              nrow=constants$n.species,
                              ncol=constants$n.species,
                              byrow = F)%>%
    cov2cor()
  
  rownames(correlation_matrix) <- y
  colnames(correlation_matrix) <- y
  #names_shan_indx <- paste("shan")
  #mtcars[grep("Merc", rownames(mtcars)), ]
  return(correlation_matrix)

}

extract_shan <- function(x){
  constants <- x[[1]]
  shan_index <- matrix(x[[2]][ grep("shan", rownames(x[[2]])), 1],
                       nrow=constants$n.sites,
                       ncol=2,
                       byrow = T)
  return(shan_index)
}

extract_sigma <- function(x){
  #constants <- x[[1]]
  sigma <- matrix(x[[2]][ grep("sigm", rownames(x[[2]])), 1],
                       nrow=1,
                       ncol=2,
                       byrow = T)
  return(sigma)
}

# Getting Shannon Index from data
shans <- lapply(all_data, extract_shan)%>%
  do.call("cbind", .)%>%
  as.data.frame()%>%
  mutate(site = simulations_all[[1]]$sites$X1km_square)

 colnames(shans) <- c("bbIDM17", "bbIDM18", 
           "bbSPE17", "bbSPE18",
           "bbGRO17", "bbGRO18",
           "sbIDM17", "sbIDM18", 
           "sbSPE17", "sbSPE18",
           "sbGRO17", "sbGRO18",
           "hvIDM17", "hvIDM18", 
           "hvSPE17", "hvSPE18",
           "hvGRO17", "hvGRO18",
           "sites")
 
 #write.csv(shans,"Results/shannon_index.csv", row.names = FALSE)

 #
 names_list <- c(rep(list(simulations_all[[1]]$species_names),3),
                 rep(list(simulations_all[[3]]$species_names),1), #change this to the node with 1 changed to 2
                   rep(list(simulations_all[[3]]$species_names),2),
                 rep(list(simulations_all[[2]]$species_names),2), #change this to the node with 1 changed to 2
                 rep(list(simulations_all[[2]]$species_names),1)
                   #rep(list(simulations_all[[3]]$species_names),3)
                 )

 #Abbreviate the names of the species
abb_names_list <- lapply(names_list, function(x){
  abbreviate(x, minlength = 8)
})
                   
 correlation_matrices <- pmap(list(all_data, abb_names_list), extract_covariance)
   
#plot correlation plots
 cor_plots <- list()
 for(i in 1:9){
   cor_plots[[i]] <- ggcorrplot(correlation_matrices[[i]], hc.order = TRUE)
 }

 ggarrange(cor_plots[[1]],cor_plots[[2]],cor_plots[[3]],
           cor_plots[[4]],cor_plots[[5]],cor_plots[[6]],
           cor_plots[[7]],cor_plots[[8]],cor_plots[[9]],
           ncol=3, nrow=3)
 #cor_plots[[9]]

 

 
 # plotting corr heatmap
 par(mfrow=c(3,3))
 cor_plots1 <- list()

 for(i in 1:9){
   cor_plots1[[i]] <-  heatmaply_cor(x = correlation_matrices[[i]], xlab = "Species list",
               ylab = "Species list", k_col = 2, k_row = 2,
               row_text_angle = 45,
               fontsize_row = 4,
               fontsize_col = 4,
               plot_method = "ggplot") 
   print(cor_plots1[[i]])
 }
 
 
 # Standard deviation
#method_names <- c("bbIDM", "bbSPE", "bbGO",
  #                "sbIDM", "sbSPE", "sbGO",
  #                "hvIDM", "hvSPE", "hvGO")
 
 sigmas <- lapply(all_data, extract_sigma)%>%
   do.call("rbind", .)%>%
   as.data.frame()%>%
   mutate(method = rep(c("IDM" , "SPE", "GO"), 3),
          insects = rep(c("Bumblebees", "Solitary bees", "Hoverflies"),each= 3))

colnames(sigmas) <- c("Detection", "Observed", "Methods", "Insects")
write.csv(sigmas,"Results/sigmas.csv", row.names = FALSE)

# Plot variances
plot_sigmas <- melt(sigmas, id.vars = c("Methods", "Insects"))%>%
  ggplot(aes(x= Methods, y = value, group = Insects, col= Insects))+
  geom_point()+
  geom_line(linetype="dashed")+
  theme_article()+
  facet_wrap(~variable, nrow=1)
plot_sigmas 

##################################
# Alpha diversity at locations
###################################
gr_ref <- read_csv("Results/gr_ref.csv")
names(gr_ref)[1] <- "sites"
sites <- data.frame(sites = simulations_all[[1]]$sites$X1km_square)
merged_data_sites <- merge(x=sites, y=gr_ref, by="sites", all.x = TRUE, all.y = FALSE)%>%
  select(sites, region)%>%
  arrange(sites)

#Fill in the missing region for a location
merged_data_sites[62,2] <- "ENG"

all_data_alpha <- merge(shans, merged_data_sites, by = "sites", all.x = TRUE, all.y = FALSE)
  
write.csv(all_data_alpha,"Results/shannon_index.csv", row.names = FALSE)