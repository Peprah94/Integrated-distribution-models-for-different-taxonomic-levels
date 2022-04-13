#load data
#setwd("~/Documents/GitHub/Thesis/idm project/Data")
library(unix)
unix::rlimit_as(100*10^9)
library(readr)
library(GGally) 
library(tidyverse)
library(lubridate)
library(reshape2)
library(parallel)
library(doParallel)
library(foreach)
library(MCMCglmm)
library(pbapply)

#source('sim_rep.R')
#source('~/Documents/GitHub/integrated_dist_model/nimble.R')
#source('~/Documents/GitHub/integrated_dist_model/estimation.R')
#source('shanestimate.R')

#source('~/Documents/GitHub/Thesis/idm project/new/sim_rep.R')
#source('~/Documents/GitHub/Thesis/idm project/new/nimble.R')
#source('~/Documents/GitHub/Thesis/idm project/new/estimation.R')
#source('~/Documents/GitHub/Thesis/idm project/new/all.R')

#Reading the data
#genus_counts <- read_csv("~/Documents/GitHub/Thesis/idm project/Data/FIT_counts.csv")
#pan_traps <- read_csv("~/Documents/GitHub/Thesis/idm project/Data/pan_traps.csv")
#species_lookup <- read_csv("~/Documents/GitHub/Thesis/idm project/Data/species_lookup.csv")
#source('nimble.R')
#source('estimation.R')

genus_counts <- read_csv("FIT_counts.csv")
pan_traps <- read_csv("pan_traps.csv")
species_lookup <- read_csv("species_lookup.csv")

#Formating the date into year, month and day
pan_traps$date <- dmy_hms(pan_traps$date)
genus_counts$date <- dmy_hms(genus_counts$date)

###joining the insect group data to the species data
#Rename the column names as taxon aggregated to species
names(pan_traps)[names(pan_traps) == "taxon_aggregated"] <- "species" 

#Merge the species pantrap data to the species_lookup and match them to their taxon_group
update_pantraps <- merge(pan_traps, species_lookup, by= "species")%>%
  dplyr::mutate(year = lubridate::year(date), 
                month = lubridate::month(date), 
                day = lubridate::day(date))

# Format the group taxon data (FSV counts)
update_genus <-genus_counts %>%
  dplyr::mutate(year = lubridate::year(date), 
                month = lubridate::month(date), 
                day = lubridate::day(date))%>%
  filter(insect_group != "all_insects_total")

# grouping variables
groups <- update_pantraps%>%
  group_by(insect_group)%>%
  group_rows()
#unique id for the species
unique_id <- update_pantraps%>%
  group_indices(update_pantraps$species)
#unique id for the groups
group_id <- update_pantraps%>%
  group_indices(update_pantraps$insect_group)
#unique id for the sites
site_id <- update_pantraps%>%
  group_indices(update_pantraps$X1km_square)

x <- sort(unique(update_pantraps$X1km_square))
id <- data.frame(X1km_square = x)
#Data with group id
update_pantraps <- cbind(update_pantraps, unique_id, group_id, site_id)


includecov = FALSE

unique_visit <- function(x){ 
  p <- which((!is.na(x)))
  q <- x[p]
  q <- as.numeric(q[order(p)])
  X<-rep(NA,5)
  X[1:length(q)] <- q
  return(X)
}

spe_visit <- function(x){ p <- which((!is.na(x)))
q <- x[p]
q <- as.numeric(q[order(p)])
X<-rep(NA,4)
X[1:length(q)] <- q
return(X)
}



#Function that picks the year and the taxon group and outputs the formatted data to be sent to jags
data_est <- function(year_id, group_name){
  all_year_spe <- update_pantraps%>%
    filter(insect_group== group_name)%>%
    select(c(11))
  all_year_indices <- data.frame(unique_id = sort(unique(all_year_spe$unique_id)))
  
  pat_year_species <- update_pantraps%>%
    filter(insect_group== group_name & year == year_id)%>%
    select(c(11))
  
  part_year_indices <- data.frame(unique_id = sort(unique(pat_year_species$unique_id)))
  
  species<- update_pantraps%>%
    filter(insect_group== group_name & year == year_id)
  
  unique_visits <- update_pantraps%>%
     filter(insect_group== group_name & year == year_id)%>%
    group_by("site_id", "date")%>%
    arrange(site_id)
    
  
  genus <- update_genus %>%
    filter(insect_group== group_name & year == year_id)
  
  
  all_data <-full_join(species,genus, by = c("date", "X1km_square"),all = TRUE )%>%
    ungroup()%>%
    dplyr::arrange(c("date"))%>% 
    dplyr::group_by(X1km_square)%>%
    mutate(rank=as.numeric(as.factor((date))))
    #group_indices()
    #matsindf::index_column(var_to_index = "X1km_square", time_var = "date")
  
  n.visits <- sort(unique(all_data$rank))
  n.species <- sum(!is.na(unique(all_data$unique_id)))
  spe_unique_id <- data.frame(unique_id=unique(all_data$unique_id)[!is.na(unique(all_data$unique_id))])
 species_data <- array(NA, dim=c(74, n.species,length(n.visits) ))
  for(i in 1:length(n.visits)){
      subset_data <- all_data%>%
      filter(rank==i)%>%
      select(-c(1,3:4,6:10,12:19))%>%
      melt(id=c("X1km_square", "unique_id"))%>%
      filter(!is.na(unique_id))%>%
      distinct(X1km_square,unique_id, .keep_all = TRUE)
      
      visit_species <- data.frame(unique_id = rep(unique(all_data$unique_id)[!is.na(unique(all_data$unique_id))], nrow(id)), X1km_square=rep(id$X1km_square,n.species))
      
      species_data[,,i]<-subset_data%>%
        full_join(visit_species, by = c("unique_id","X1km_square"), all = TRUE)%>%
        dcast( X1km_square~unique_id ,value.var ="value",fun.aggregate = mean, na.rm = FALSE, drop = FALSE)%>%
      mutate_all(~replace(.,is.na(.),0))%>%
      arrange(X1km_square)%>%
        select(-c(1))%>%
        as.matrix()
       
       # <-as.matrix(xx)
  }
  

  #################################################################
  #             GENUS
  #################################################################
 
 genus_data <- all_data%>%
   select(-c(1,3:6,6:14, 16:18))%>%
   dcast(X1km_square~ rank ,value.var ="count",fun.aggregate = mean, na.rm = TRUE, drop = FALSE)%>%
   full_join(id, by="X1km_square", all = TRUE)%>%
   mutate_all(~replace(.,is.na(.),0))%>%
   select(-c(1))%>%
   as.matrix()
 
  
  ##################################################################
  #       RETURNED DATA
  ##################################################################
  #dataspp <- array(unlist_kk,dim=c(74, (n.species), length(n.visits)))
  #dataspp <- round(apply(kk5_species[,-1,-5],c(1,2,3), as.numeric),0)
  #datagen <-  round(apply(visit_data[,-1], c(1,2), as.numeric),0)
  ecological_cov =rep(0, nrow(genus_data))
  detection_cov=rep(0, nrow(genus_data))
  shan.index = NA
  n.replicates = rep(5, 74)
  data <- list(mat.species=species_data,
               mat.genus =  genus_data,
               ecological_cov =ecological_cov,
               detection_cov = detection_cov)
  return(data)
}

years = list(2017, 2018)
group= list("bumblebees",  "hoverflies", "solitary_bees")
#nreplicates <- 10
simulations_all <- pblapply(group, function(x){
  pblapply(years, function(z){
    data <- data_est(year_id=z, group_name=x)
  }, cl=4)
}, cl=4)

simulations_all <- flatten(simulations_all)

save(simulations_all, file="data_idm.RData")
