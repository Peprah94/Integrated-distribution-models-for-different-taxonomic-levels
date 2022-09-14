setwd("/Volumes/kwakupa/IDM_new/Data")
load("/Volumes/kwakupa/IDM_new/Data/data_idm.RData")
source("/Volumes/kwakupa/IDM_new/idm_functions/fnx_for_estimation.R")

require(dplyr)
require(ggcorrplot)
require(purrr)
require(ggpubr)
require(heatmaply)
require(reshape2)
library(ggplot2)
require(egg)
require(readr)
library(stringr)
if(!require(devtools)) install.packages("devtools")
theme_set(theme_bw())


#load data


#load estimates of parameters
load("estimate_data_inter_bum1.RData")
bb_inter <- estimates
load("estimate_data_species_bum1.RData")
bb_spe <- estimates
load("estimate_data_genus_bum1.RData")
bb_gen <- estimates
load("estimate_data_inter_hov1.RData")
hov_inter <- estimates
load("estimate_data_species_hov1.RData")
hov_spe <- estimates
load("estimate_data_genus_hov1.RData")
hov_gen <- estimates
load("estimate_data_inter_sol1.RData")
sol_inter <- estimates
load("estimate_data_species_sol1.RData")
sol_spe <- estimates
load("estimate_data_genus_sol1.RData")
sol_gen <- estimates

all_data <- list(bb_inter , bb_spe , bb_gen ,
              hov_inter , hov_spe , hov_gen,
              sol_inter ,  sol_spe , sol_gen )


############
# Plot of estimates
########
#Extract the covariance estimates
extract_covariance <- function(x,y){
  constants <- x[[1]]
  correlation_matrix <- matrix(x[[2]][grep("Cov", rownames(x[[2]])), 1], 
                              nrow=constants$n.species,
                              ncol=constants$n.species,
                              byrow = F)%>%
    cov2cor()
  
  rownames(correlation_matrix) <- y
  colnames(correlation_matrix) <- y
  return(correlation_matrix)

}

# function to extract shannon index for all estimates
extract_shan <- function(x){
  constants <- x[[1]]
  shan_index <- matrix(x[[2]][ grep("shan", rownames(x[[2]])), 1],
                       nrow=constants$n.sites,
                       ncol=2,
                       byrow = F)
  return(shan_index)
}

extract_richness <- function(x){
  constants <- x[[1]]
  richness <- matrix(x[[2]][ grep("richness", rownames(x[[2]])), 2],
                       nrow=constants$n.sites,
                       ncol=2,
                       byrow = F)
  return(richness)
}

#function to extract the sd of the site effects
extract_sigma <- function(x){
  sigma <- matrix(x[[2]][ grep("sigm", rownames(x[[2]])), 1],
                       nrow=1,
                       ncol=4,
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
           "hvIDM17", "hvIDM18", 
           "hvSPE17", "hvSPE18",
           "hvGRO17", "hvGRO18",
           "sbIDM17", "sbIDM18", 
           "sbSPE17", "sbSPE18",
           "sbGRO17", "sbGRO18",
  
           "sites")

 names_list <- c(rep(list(simulations_all[[1]]$species_names),3),
                   rep(list(simulations_all[[2]]$species_names),3),
                 rep(list(simulations_all[[3]]$species_names),3)
                 )
 
 
 # Getting species richness from the data
richness <- lapply(all_data, extract_richness)%>%
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
 
 names_list <- c(rep(list(simulations_all[[1]]$species_names),3),
                 rep(list(simulations_all[[2]]$species_names),3),
                 rep(list(simulations_all[[3]]$species_names),3)
 )
 

# Species interractions
 #Abbreviate the names of the species
abb_names_list <- lapply(names_list, function(x){
  ret <- abbreviate(x, minlength = 8)
  #word(x, 1,2, sep=" ")
  sub("(\\w+\\s+\\w+).*", "\\1", ret)
})
                   
 correlation_matrices <- pmap(list(all_data, abb_names_list), extract_covariance)
   title_text <- c("a) bumblebees \n IDM", "b) bumblebees \n Speices only", "c) bumblebees \n Group only",
                   "d) hoverflies \n IDM", "e) hoverflies \n Speices only", "f) hoverflies \n Group only",
                   "g) solitarybees \n IDM", "h) solitarybees \n Speices only", "i) solitarybees \n Group only")
 
  
  gg_corrplots <-  lapply(seq_along(correlation_matrices), function(x){ 
     ggcorrplot(correlation_matrices[[x]],
              hc.order = TRUE,
              type = "upper",
              #outline.color = "white",
              #ggtheme = ggplot2::theme_minimal(plot.title = element_text(size = 9)),
              #lab_size = 2,
              tl.cex = 2,
              tl.srt = 90,
              lab_size = 2,
              title = title_text[[x]])+
      theme(plot.title = element_text(size = 9) )
              }) 

   
   #gg1 <- ggpubr::ggarrange(gg_corrplots[[1]],gg_corrplots[[2]], gg_corrplots[[3]],
                     #gg_corrplots[[4]],gg_corrplots[[5]], gg_corrplots[[6]],
                     #gg_corrplots[[7]],gg_corrplots[[8]], gg_corrplots[[9]],
        #             ncol = 3, nrow=1, common.legend = TRUE, legend = "right") 
   
   gg2 <- ggpubr::ggarrange(gg_corrplots[[1]],gg_corrplots[[2]], gg_corrplots[[3]],
                            gg_corrplots[[4]],gg_corrplots[[5]], gg_corrplots[[6]],
                            gg_corrplots[[7]],gg_corrplots[[8]], gg_corrplots[[9]],
                            ncol = 3, nrow = 3, common.legend = TRUE, legend = "right",
                            heights = c(1.5,2,2)) 
   
   figure <- ggpubr::annotate_figure(gg2, left = "Species list", bottom = "Species list")
   
   ggsave("Figures/species_association.png", plot = figure, dpi = 320, 
          height = 18, width = 18, units = "cm")
   
 
 #extracting the sd of the site effects
 sigmas <- lapply(all_data, extract_sigma)%>%
   do.call("rbind", .)%>%
   as.data.frame()%>%
   mutate(method = rep(c("IDM" , "Species only", "Group Only"), 3),
          insects = rep(c("Bumblebees", "Solitary bees", "Hoverflies"),each= 3))

colnames(sigmas) <- c("Site detection", "Site Observed", "Year detection", "Year Observed","Methods", "Insect Group")
write.csv(sigmas,"Figures/sigmas.csv", row.names = FALSE)

# Plot variances
plot_sigmas <- melt(sigmas, id.vars = c("Methods", "Insect Group"))%>%
  ggplot(aes(x= Methods, y = value, group = `Insect Group`, col= `Insect Group`))+
  geom_point()+
  geom_line(linetype="dashed")+
  theme_bw()+
  ylab("Variance")+
  facet_wrap(~variable, nrow=2)
ggsave(filename = "Figures/site_variance.png",plot_sigmas)


#saving the shannon diversity with the region and site ID
gr_ref <- read_csv("gr_ref.csv")
names(gr_ref)[1] <- "sites"
sites <- data.frame(sites = simulations_all[[1]]$sites$X1km_square)
merged_data_sites <- merge(x=sites, 
                           y=gr_ref, 
                           by="sites", 
                           all.x = TRUE, 
                           all.y = FALSE)%>%
  select(sites, region)%>%
  arrange(sites)

#Fill in the missing region for a location
merged_data_sites[62,2] <- "ENG"
all_data_estimates <- merge(shans, 
                            merged_data_sites, 
                            by = "sites", 
                            all.x = TRUE, 
                            all.y = FALSE)
  
write.csv(all_data_estimates,"Figures/shannon_index.csv", row.names = FALSE)




#################################################################


###############
# crossvalidation_results
##############

cross_validation <- lapply(all_data, function(x){
  mean <- x[[4]]$CVvalue
  sd <- x[[4]]$CVstandardError
  ret <- c(mean, sd)
})%>%
 do.call("rbind", .)%>%
  data.frame()%>%
  mutate(group = rep(c("bumblebees", "hoverflies", "solitarybees"), each = 3),
         method = rep(c("IDM", "Species only", "Group only"), 3))
colnames(cross_validation) <- c("CV_value", "CV_standError", "group", "method")

cross_validation_results <- cross_validation%>%
  tidyr::pivot_longer(., c(CV_value,CV_standError))%>%
  tidyr::pivot_wider(., names_from = method)

write.csv(cross_validation_results,"Figures/cross_validation.csv", row.names = FALSE)

#########################
# Parameter estimates
########################

extract_beta0 <- function(x, y){
ret<- x[[2]][ grep("beta0", rownames(x[[2]])), c(1,4,5)]%>%
    data.frame()%>%
mutate(species_name = y)
colnames(ret) <- c("mean", "low", "upper", "species")
    return(ret)

}

beta0_matrices <- pmap(list(all_data, abb_names_list), extract_beta0)
beta0_matrices <- list(beta0_matrices[[1]], beta0_matrices[[2]],
                    beta0_matrices[[4]], beta0_matrices[[5]],
                    beta0_matrices[[7]], beta0_matrices[[8]]
                    )

title_text_alt <- c("a) bumblebees \n IDM", "b) bumblebees \n Speices only", 
                "c) hoverflies \n IDM", "d) hoverflies \n Speices only",
                "e) solitarybees \n IDM", "f) solitarybees \n Speices only")

gg_beta0 <-  lapply(seq_along(beta0_matrices), function(x){ 
beta0_matrices[[x]]%>%
    arrange(mean)%>%
    mutate(id = row_number(),
           species = factor(species, levels = unique(species)))%>%
  ggplot()+
    geom_point(aes(y = species, x = mean))+
    #scale_y_discrete(breaks = x$id, labels = x$species)+
    geom_segment(aes(y = species,yend = species,
                     x = low, xend = upper))+
    labs(title = title_text_alt[[x]])+
    xlab("")+
    ylab("")+
    theme_bw()+
    theme(plot.title = element_text(size = 7),
          axis.text.x = element_text(size = 8), 
          axis.text.y =element_text(size = 3) )
}) 


gg1 <- ggpubr::ggarrange(gg_beta0[[1]],gg_beta0[[2]], gg_beta0[[3]],
                         gg_beta0[[4]],gg_beta0[[5]], gg_beta0[[6]],
                         ncol = 3, nrow = 2, common.legend = TRUE, legend = "right",
                         heights = c(2,2)) 

figure_beta0 <- ggpubr::annotate_figure(gg1, left = "Species", bottom = "Beta0 estimate")

ggsave("Figures/data_beta0.png", plot = figure_beta0, dpi = 320, 
       height = 18, width = 18, units = "cm")


#########
# alpha0
#######
extract_alpha0 <- function(x, y){
  ret<- x[[2]][ grep("alpha0", rownames(x[[2]])), c(1,4,5)]%>%
    data.frame()%>%
    mutate(species_name = y)
  colnames(ret) <- c("mean", "low", "upper", "species")
  return(ret)
  
}

alpha0_matrices <- pmap(list(all_data, abb_names_list), extract_alpha0)
alpha0_matrices <- list(alpha0_matrices[[1]], alpha0_matrices[[2]],
                        alpha0_matrices[[4]], alpha0_matrices[[5]],
                        alpha0_matrices[[7]], alpha0_matrices[[8]]
)
gg_alpha0 <-  lapply(seq_along(alpha0_matrices), function(x){ 
  alpha0_matrices[[x]]%>%
    arrange(mean)%>%
    mutate(id = row_number(),
           species = factor(species, levels = unique(species)))%>%
    ggplot()+
    geom_point(aes(y = species, x = mean))+
    #scale_y_discrete(breaks = x$id, labels = x$species)+
    geom_segment(aes(y = species,yend = species,
                     x = low, xend = upper))+
    labs(title = title_text_alt[[x]])+
    xlab("")+
    ylab("")+
    theme_bw()+
    theme(plot.title = element_text(size = 7),
          axis.text.x = element_text(size = 8), 
          axis.text.y =element_text(size = 3) )
}) 


gg3 <- ggpubr::ggarrange(gg_alpha0[[1]],gg_alpha0[[2]], gg_alpha0[[3]],
                         gg_alpha0[[4]],gg_alpha0[[5]], gg_alpha0[[6]],
                         ncol = 3, nrow = 2, common.legend = TRUE, legend = "right",
                         heights = c(2,2)) 

figure_alpha0 <- ggpubr::annotate_figure(gg3, left = "Species", bottom = "Intercept estimate")

ggsave("Figures/data_alpha0.png", plot = figure_alpha0, dpi = 320, 
       height = 18, width = 18, units = "cm")


#############
# beta1
#############3
extract_intercept_lambda <- function(x){
  ret<- x[[2]][ grep("intercept_lambda", rownames(x[[2]])), c(1,4,5)]

  #colnames(ret) <- c("mean", "low", "upper")
  return(ret)
  
}

intercept_lambda <- lapply(all_data, function(x){
  extract_intercept_lambda(x)})%>%
  do.call("rbind", .)

#take out data that belongs to the species 
idx <- c(2,5,8)
intercept_lambda <- intercept_lambda[-idx,]%>%
  data.frame()%>%
  mutate(method = rep(c("IDM", "Group only"), 3),
         group = rep(c("bumblebees", "hoverflies", "solitarybees"), each = 2))

colnames(intercept_lambda) <- c("mean", "low", "high", "method", "group")

gg_intercept_lambda <-  intercept_lambda%>%
  mutate(group_method = c("bumblebees-IDM","bumblebees-Group Only" ,
                          "hoverflies-IDM","hoverflies-Group Only",
                          "solitarybees-IDM","solitarybees-Group Only"))%>%
# tidyr::pivot_longer(cols = c("group", "method"))
  arrange(mean)%>%
  mutate(id = row_number(),
         species = factor(group_method, levels = unique(group_method)))%>%
  ggplot()+
  geom_point(aes(y = species, x = mean))+
  #scale_y_discrete(breaks = x$id, labels = x$species)+
  geom_segment(aes(y = species,yend = species,
                   x = low, xend = high))+
  #labs(title = title_text_alt[[x]])+
  xlab("gamma estimate")+
  ylab("Method & Insect Group")+
  theme_bw()+
  theme(plot.title = element_text(size = 7),
        axis.text.x = element_text(size = 8), 
        axis.text.y =element_text(size = 5) )



ggsave("Figures/data_intercept_lambda.png", plot = gg_intercept_lambda, dpi = 320, 
       height = 8, width = 10, units = "cm")


library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

shannon <- read.csv("Figures/shannon_index.csv")

names(shannon) <- c("sites","bb_IDM_17", "bb_IDM_18", "bb_SPE_17", "bb_SPE_18",
                    "bb_GRO_17", "bb_GRO_18", "sb_IDM_17", "sb_IDM_18",
                    "sb_SPE_17", "sb_SPE_18", "sb_GRO_17", "sb_GRO_18",
                    "hv_IDM_17", "hv_IDM_18", "hv_SPE_17", "hv_SPE_18",
                    "hv_GRO_17", "hv_GRO_18", "region")
#sites
site_id <- shannon$sites

shannon_long <- pivot_longer(shannon, cols = c(-sites, -region),
                             names_to = c("group", "model", "year"),
                             names_sep = "_")


# new figure

Shannon_boxplot <- ggplot(data = shannon_long) +
  geom_boxplot(aes(x = group, y = value,
                   fill = model,
                   color = model),
               alpha = 0.2) +
  facet_wrap(~year) +
  ylab("Mean Shannon index across sites") +
  xlab("Insects") +
  theme_bw()+
  scale_fill_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),
                    name="Method", labels = c("Insect group", "IDM", "Species only")) +
  scale_color_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),
                     name="Method", labels = c("Insect group", "IDM", "Species only")) +
  scale_x_discrete(labels= c("Bumblebees", "Hoverflies", "Solitary bees"))



png("Figures/shannon_boxplot.png", 
    height = 8, width = 16, units = "cm", res = 320)

Shannon_boxplot

dev.off()


# map of sites (represented at 10Km because too small to see at 1Km)
# devtools::install_github("https://github.com/colinharrower/BRCmap")
library(BRCmap)

# par(mar = c(0.1,0.1,1,0.1))
plot_GIS(UK[UK$REGION == "Great Britain",], xlab="", ylab="", show.axis = FALSE, show.grid = FALSE, 
         fill.col = "lightgrey", no.margin = TRUE)

plotUK_gr(shannon_long$sites, gr_prec = 10000, border = "black")


# map of shannon

coords <- gr_let2num(shannon_long$sites)

all_data <- cbind(shannon_long, coords)

UK_countries_lowres@data$ID = rownames(UK_countries_lowres@data)
UK_countries_lowres_df <- ggplot2::fortify(UK_countries_lowres)

all_data <- all_data %>%
  mutate(model = case_when(model == "GRO" ~ "Group only",
                           model == "SPE" ~ "Species only",
                           model == "IDM" ~ "IDM"),
         year = case_when(year == 17 ~ 2017,
                          year == 18 ~ 2018))

bb_shannon_map <- ggplot(data = all_data %>%
                           filter(group == "bb" & year == "2017")) +
  geom_polygon(data = UK_countries_lowres_df, 
               aes(x = long, y = lat, group = group),
               fill = "gray93", colour = "black") +
  geom_point(aes(x = EASTING, y = NORTHING,
                 fill = value), shape = 22,
             color = "gray") +
  facet_wrap(~year+model, nrow = 1, ncol = 3) +
  theme_map() +
  coord_equal() +
  scale_fill_viridis_c(name = "Shannon", limits = c(0, 4))

#png("Figures/bb_shannon_map.png", 
 #   height = 16, width = 14, units = "cm", res = 320)

#bb_shannon_map

#dev.off()



sb_shannon_map <- ggplot(data = all_data %>%
                           filter(group == "sb" & year == "2017")) +
  geom_polygon(data = UK_countries_lowres_df, 
               aes(x = long, y = lat, group = group),
               fill = "gray93", colour = "black") +
  geom_point(aes(x = EASTING, y = NORTHING,
                 fill = value), shape = 22,
             color = "gray") +
  facet_wrap(~year+model, nrow = 1, ncol = 3) +
  theme_map() +
  coord_equal() +
  scale_fill_viridis_c(name = "Shannon", limits = c(0, 4))

#png("Figures/sb_shannon_map.png", 
#    height = 16, width = 14, units = "cm", res = 320)

#sb_shannon_map

#dev.off()


hv_shannon_map <- ggplot(data = all_data %>%
                           filter(group == "hv" & year == "2017")) +
  geom_polygon(data = UK_countries_lowres_df, 
               aes(x = long, y = lat, group = group),
               fill = "gray93", colour = "black") +
  geom_point(aes(x = EASTING, y = NORTHING,
                 fill = value), shape = 22,
             color = "gray") +
  facet_wrap(~year+model, nrow = 1, ncol = 3) +
  theme_map() +
  coord_equal() +
  scale_fill_viridis_c(name = "Shannon", limits = c(0, 4))

#png("Figures/hv_shannon_map.png", 
#    height = 16, width = 14, units = "cm", res = 320)

#hv_shannon_map

#dev.off()
estimates_18 <- ggpubr::ggarrange(bb_shannon_map, 
                  sb_shannon_map, 
                  hv_shannon_map, 
                  nrow = 3, 
                  ncol = 1, 
                  common.legend = TRUE,
                  legend = "right")

ggsave("Figures/estimate_shannon.png", 
       plot = estimates_18, 
       dpi = 320, 
       height = 18, 
       width = 18, 
       units = "cm")


bb_shannon_map <- ggplot(data = all_data %>%
                           filter(group == "bb")) +
  geom_polygon(data = UK_countries_lowres_df, 
               aes(x = long, y = lat, group = group),
               fill = "gray93", colour = "black") +
  geom_point(aes(x = EASTING, y = NORTHING,
                 fill = value), shape = 22,
             color = "gray") +
  facet_wrap(~year+model, nrow = 2, ncol = 3) +
  theme_map() +
  coord_equal() +
  scale_fill_viridis_c(name = "Shannon")

ggsave("Figures/bb_shannon_map.png", 
       plot = bb_shannon_map, 
       dpi = 320, 
       height = 18, 
       width = 18, 
       units = "cm")

sb_shannon_map <- ggplot(data = all_data %>%
                           filter(group == "sb")) +
  geom_polygon(data = UK_countries_lowres_df, 
               aes(x = long, y = lat, group = group),
               fill = "gray93", colour = "black") +
  geom_point(aes(x = EASTING, y = NORTHING,
                 fill = value), shape = 22,
             color = "gray") +
  facet_wrap(~year+model, nrow = 2, ncol = 3) +
  theme_map() +
  coord_equal() +
  scale_fill_viridis_c(name = "Shannon")
ggsave("Figures/sb_shannon_map.png", 
       plot = sb_shannon_map, 
       dpi = 320, 
       height = 18, 
       width = 18, 
       units = "cm")

hv_shannon_map <- ggplot(data = all_data %>%
                           filter(group == "hv")) +
  geom_polygon(data = UK_countries_lowres_df, 
               aes(x = long, y = lat, group = group),
               fill = "gray93", colour = "black") +
  geom_point(aes(x = EASTING, y = NORTHING,
                 fill = value), shape = 22,
             color = "gray") +
  facet_wrap(~year+model, nrow = 2, ncol = 3) +
  theme_map() +
  coord_equal() +
  scale_fill_viridis_c(name = "Shannon")

ggsave("Figures/hv_shannon_map.png", 
       plot = hv_shannon_map, 
       dpi = 320, 
       height = 18, 
       width = 18, 
       units = "cm")
