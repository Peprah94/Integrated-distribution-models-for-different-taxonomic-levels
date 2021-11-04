#load data

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


#Function that picks the year and the taxon group and outputs the formatted data to be sent to NIMBLE
data_est <- function(year_id, group_name){
 
  #select the particular group data for all the years
all_year_spe <- update_pantraps%>%
  filter(insect_group== group_name)%>%
  select(c(11))
  
  #get the indices for the years
all_year_indices <- data.frame(unique_id = sort(unique(all_year_spe$unique_id)))

  #get the data for each year
pat_year_species <- update_pantraps%>%
  filter(insect_group== group_name & year == year_id)%>%
  select(c(11))
  
#get the indices for each year
part_year_indices <- data.frame(unique_id = sort(unique(pat_year_species$unique_id)))
  
#select the species data only
species<- update_pantraps%>%
  filter(insect_group== group_name & year == year_id)

 # select the genus data only
genus <- update_genus %>%
  filter(insect_group== group_name & year == year_id)

  #number of visit for each year
n.visits <- sort(unique(update_pantraps$month))
  
  #number of species in each year and each group
n.species <- length(unique(species$unique_id))
  
  #put the data into each group
kk_species <-  kk1_species <- list(c())
for(i in 1:length(n.visits)){
  kk1_species[[i]]<- pivot_wider(data = update_pantraps, values_from="n_traps_present", 
                                          names_from="month")%>%
   filter(insect_group== group_name & year == year_id)%>%
   select(-c(1,3:8, 10,11))%>%
   melt(id=c("X1km_square", "unique_id"))%>%
   filter(variable == n.visits[i])%>%
  dcast( X1km_square~unique_id ,value.var ="value",fun.aggregate = mean, na.rm = TRUE, drop = FALSE)%>%
    mutate_all(~replace(.,is.nan(.),NA))%>%
    select(-1)
 
species_sites <- pivot_wider(data = update_pantraps, values_from="n_traps_present", 
                            names_from="month")%>%
  filter(insect_group== group_name & year == year_id)%>%
  select(-c(1,3:8, 10,11))%>%
  melt(id=c("X1km_square", "unique_id"))%>%
  filter(variable == n.visits[i])%>%
  dcast( X1km_square~unique_id ,value.var ="value",fun.aggregate = mean, na.rm = TRUE, drop = FALSE)%>%
  mutate_all(~replace(.,is.nan(.),NA))%>%
  select(1)
}

species_dim <- dim(kk1_species[[1]])
species_rows_no <- length(all_year_indices$unique_id)
kk2_species <- array(unlist(kk1_species), dim = c(species_dim[1],species_dim[2],5))
visit_species_data <- array(NA, dim=c(species_dim[1],species_dim[2],5))
for(i in 1:species_dim[1]){
  for(j in 1:species_dim[2]){
    if(sum(is.na(kk2_species[i,j,])) < 5){
    visit_species_data[i,j,] <- unique_visit(kk2_species[i,j,])
    }else{
      visit_species_data[i,j,] <- rep(NA,5)   
     }
}
}

##
#kk4_species <- array(NA, dim=c(74,(species_dim[2]+1),5))
  data_sp <- array(NA, dim=c(74,(species_rows_no+1),5))
kk5_species <-array(NA, dim=c(74,species_rows_no,5))
for(nvisit in 1:5){
  #as.array(cbind(species_list,visit_species_data[ , ,nvisit]))
  kkk <- cbind(species_sites,visit_species_data[ , ,nvisit])
  kk4_species <-kkk%>%
    full_join(id, by = "X1km_square",all = TRUE)

  #t((kk4_species[,-1,nvisit]))
 kk5_species[,,nvisit] <- t(kk4_species[,-1])%>%
    cbind(part_year_indices)%>%    
    full_join(all_year_indices, by = "unique_id", all = TRUE)%>%  
    select(-75)%>%  
    t()
 
 data_sp[,,nvisit] <- as.matrix(cbind(kk4_species[,1,1], kk5_species[,,nvisit]))
}


#################################################################
#             GENUS
#################################################################
gen_visit <- sort(unique(genus$month))
update_genus <- update_genus%>%arrange(month)
  kk_genus<- pivot_wider(data = update_genus, values_from="count", 
                                names_from="month")%>%
    filter(insect_group== group_name & year == year_id)%>%
    select(-c(2,3,4,5,))%>%
    melt(id=c("X1km_square"))%>%
    #distinct(X1km_square,variable, .keep_all = TRUE)%>%
    dcast(X1km_square~ variable ,value.var ="value",fun.aggregate = mean, na.rm = TRUE, drop = FALSE)%>%
    mutate_all(~replace(.,is.nan(.),NA))


 #Sort out the data into visits 
  visit_data = matrix(NA, nrow=nrow(kk_genus), ncol=4)
  for(i in 1: nrow(kk_genus)){
    visit_data[i,] <- apply(kk_genus[i,-1],1, spe_visit)
  }
  
  visit_data <- as.data.frame(cbind(kk_genus[,1], visit_data))
  colnames(visit_data) <- c("X1km_square", "visit1", "visit2", "visit3", "visit4")
  
 visit_data <- visit_data %>%
    full_join(id, by = "X1km_square",all = TRUE )


##################################################################
#       RETURNED DATA
##################################################################
#dataspp <- array(unlist_kk,dim=c(74, (n.species), length(n.visits)))
 dataspp <- round(apply(kk5_species[,,-5],c(1,2,3), as.numeric),0)
datagen <-  round(apply(visit_data[,-1], c(1,2), as.numeric),0)
ecological_cov =rep(0, nrow(datagen))
detection_cov=rep(0, nrow(datagen))
shan.index = NA
n.replicates = rep(5, 74)
data <- list(mat.species=dataspp,
             mat.genus = datagen,
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
