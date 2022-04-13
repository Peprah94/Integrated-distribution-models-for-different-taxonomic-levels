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
load("/Volumes/kwakupa/idm_data/data_idm.RData")
source("fnx_for estimation.R")

#sites
site_id <- simulations_all[[1]]$sites$X1km_square
#BUMBLEBEES
load("/Volumes/kwakupa/idm_data/estimate_inter_data.RData")
load("/Volumes/kwakupa/idm_data/estimate_species_data.RData")
load("/Volumes/kwakupa/idm_data/estimate_genus_data.RData")
bb_inter_17 <- inter_estimates_data
bb_species_17 <- species_estimates_data
bb_genus_17 <- genus_estimates_data

load("/Volumes/kwakupa/idm_data1/estimate_inter_data.RData")
load("/Volumes/kwakupa/idm_data1/estimate_species_data.RData")
load("/Volumes/kwakupa/idm_data1/estimate_genus_data.RData")
bb_inter_18 <- inter_estimates_data
bb_species_18 <- species_estimates_data
bb_genus_18 <- genus_estimates_data

# HOVERFLIES
load("/Volumes/kwakupa/idm_data2/estimate_inter_data.RData")
load("/Volumes/kwakupa/idm_data2/estimate_species_data.RData")
load("/Volumes/kwakupa/idm_data2/estimate_genus_data.RData")
hv_inter_17 <- inter_estimates_data
hv_species_17 <- species_estimates_data
hv_genus_17 <- genus_estimates_data

load("/Volumes/kwakupa/idm_data3/estimate_inter_data.RData")
load("/Volumes/kwakupa/idm_data3/estimate_species_data.RData")
load("/Volumes/kwakupa/idm_data3/estimate_genus_data.RData")
hv_inter_18 <- inter_estimates_data
hv_species_18 <- species_estimates_data
hv_genus_18 <- genus_estimates_data

# SOLITARY BEES
load("/Volumes/kwakupa/idm_data4/estimate_inter_data.RData")
load("/Volumes/kwakupa/idm_data4/estimate_species_data.RData")
load("/Volumes/kwakupa/idm_data4/estimate_genus_data.RData")
sb_inter_17 <- inter_estimates_data
sb_species_17 <- species_estimates_data
sb_genus_17 <- genus_estimates_data

load("/Volumes/kwakupa/idm_data5/estimate_inter_data.RData")
load("/Volumes/kwakupa/idm_data5/estimate_species_data.RData")
load("/Volumes/kwakupa/idm_data5/estimate_genus_data.RData")
sb_inter_18 <- inter_estimates_data
sb_species_18 <- species_estimates_data
sb_genus_18 <- genus_estimates_data


# estimating_shan <- function(x){
#   p <- matrix(x[[2]][(2*x[[1]]$n.species +1): (2*x[[1]]$n.species + (x[[1]]$n.sites*x[[1]]$n.species)),1],
#               nrow = x[[1]]$n.sites,
#               ncol=x[[1]]$n.species,
#               byrow=FALSE)
#   shan <- shan_index(p)
#   return(shan)
# }

estimating_shan <- function(x){
  p <- tail(x[[2]], x[[1]]$n.sites)
  #shan <- shan_index(p)
  return(p)
}


bb_shan_inter_17 <- estimating_shan(bb_inter_17)
bb_shan_species_17 <- estimating_shan(bb_species_17)
bb_shan_genus_17 <- estimating_shan(bb_genus_17)
bb_shan_inter_18 <- estimating_shan(bb_inter_18)
bb_shan_species_18 <- estimating_shan(bb_species_18)
bb_shan_genus_18 <- estimating_shan(bb_genus_18)

hv_shan_inter_17 <- estimating_shan(hv_inter_17)
hv_shan_species_17 <- estimating_shan(hv_species_17)
hv_shan_genus_17 <- estimating_shan(hv_genus_17)
hv_shan_inter_18 <- estimating_shan(hv_inter_18)
hv_shan_species_18 <- estimating_shan(hv_species_18)
hv_shan_genus_18 <- estimating_shan(hv_genus_18)

sb_shan_inter_17 <- estimating_shan(sb_inter_17)
sb_shan_species_17 <- estimating_shan(sb_species_17)
sb_shan_genus_17 <- estimating_shan(sb_genus_17)
sb_shan_inter_18 <- estimating_shan(sb_inter_18)
sb_shan_species_18 <- estimating_shan(sb_species_18)
sb_shan_genus_18 <- estimating_shan(sb_genus_18)

#Putting all data together

# all_data_17 <- data.frame(bb_shan_inter_17=bb_shan_inter_17,
#                        bb_shan_species_17=bb_shan_species_17,
#                        bb_shan_genus_17=bb_shan_genus_17,
#                        sites = site_id)

all_data_17 <- data.frame(rbind(bb_shan_inter_17,
                                bb_shan_species_17,
                                bb_shan_genus_17),
                          sites = rep(site_id, 3),
                          group= rep(c("IDM","Species Only", "Genus Only"), 
                                     each=nrow(bb_shan_inter_17)))
all_data_18 <- data.frame(rbind(bb_shan_inter_18,
                                bb_shan_species_18,
                                bb_shan_genus_18),
                          sites = rep(site_id, 3),
                          group= rep(c("IDM","Species Only", "Genus Only"), 
                                     each=nrow(bb_shan_inter_18)))



gg1 <- ggplot(all_data_17, aes(sites))+
  geom_line(mapping=aes(x=sites, y= Mean, group=as.factor(group),color=as.factor(group)))+
  geom_ribbon(aes(ymin=X95.CI_low, ymax=X95.CI_upp,group=as.factor(group), fill=as.factor(group)), alpha=0.2)+
  theme_bw()+
  theme(legend.position = "bottom",axis.text.x = element_blank(), axis.ticks = element_blank())+
  ylab("")+
  xlab("")+
  scale_color_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),name="Method")+
  scale_fill_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),name="Method")+
  ylim(c(0,3))

gg2 <- ggplot(all_data_18, aes(sites))+
  geom_line(mapping=aes(x=sites, y= Mean, group=as.factor(group),color=as.factor(group)))+
  geom_ribbon(aes(ymin=X95.CI_low, ymax=X95.CI_upp,group=as.factor(group), fill=as.factor(group)), alpha=0.2)+
  theme_bw()+
  theme(legend.position = "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7))+
  ylab("")+
  xlab("")+
  scale_color_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),name="Method")+
  scale_fill_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),name="Method")+
  ylim(c(0,3))


#Hoverflies
all_data_17 <- data.frame(rbind(hv_shan_inter_17,
                                hv_shan_species_17,
                                hv_shan_genus_17),
                          sites = rep(site_id, 3),
                          group= rep(c("IDM","Species Only", "Genus Only"), 
                                     each=nrow(hv_shan_inter_17)))
all_data_18 <- data.frame(rbind(hv_shan_inter_18,
                                hv_shan_species_18,
                                hv_shan_genus_18),
                          sites = rep(site_id, 3),
                          group= rep(c("IDM","Species Only", "Genus Only"), 
                                     each=nrow(hv_shan_inter_18)))



gg3 <- ggplot(all_data_17, aes(sites))+
  geom_line(mapping=aes(x=sites, y= Mean, group=as.factor(group),color=as.factor(group)))+
  geom_ribbon(aes(ymin=X95.CI_low, ymax=X95.CI_upp,group=as.factor(group), fill=as.factor(group)), alpha=0.2)+
  theme_bw()+
  theme(legend.position = "bottom",axis.text.x = element_blank(), axis.ticks = element_blank())+
  ylab("")+
  xlab("")+
  scale_color_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),name="Method")+
  scale_fill_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),name="Method")+
  ylim(c(0,4))

gg4 <- ggplot(all_data_18, aes(sites))+
  geom_line(mapping=aes(x=sites, y= Mean, group=as.factor(group),color=as.factor(group)))+
  geom_ribbon(aes(ymin=X95.CI_low, ymax=X95.CI_upp,group=as.factor(group), fill=as.factor(group)), alpha=0.2)+
  theme_bw()+
  theme(legend.position = "bottom",axis.text.x = element_blank(), axis.ticks = element_blank())+
  ylab("")+
  xlab("")+
  scale_color_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),name="Method")+
  scale_fill_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),name="Method")+
  ylim(c(0,4))


#Solitary bees
all_data_17 <- data.frame(rbind(sb_shan_inter_17,
                                sb_shan_species_17,
                                sb_shan_genus_17),
                          sites = rep(site_id, 3),
                          group= rep(c("IDM","Species Only", "Genus Only"), 
                                     each=nrow(sb_shan_inter_17)))
all_data_18 <- data.frame(rbind(sb_shan_inter_18,
                                sb_shan_species_18,
                                sb_shan_genus_18),
                          sites = rep(site_id, 3),
                          group= rep(c("IDM","Species Only", "Genus Only"), 
                                     each=nrow(sb_shan_inter_18)))



gg5 <- ggplot(all_data_17, aes(sites))+
  geom_line(mapping=aes(x=sites, y= Mean, group=as.factor(group),color=as.factor(group)))+
  geom_ribbon(aes(ymin=X95.CI_low, ymax=X95.CI_upp,group=as.factor(group), fill=as.factor(group)), alpha=0.2)+
  theme_bw()+
  theme(legend.position = "bottom",axis.text.x = element_blank(), axis.ticks = element_blank())+
  ylab("")+
  xlab("")+
  scale_color_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),name="Method")+
  scale_fill_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),name="Method")+
  ylim(c(0,4))

gg6 <- ggplot(all_data_18, aes(sites))+
  geom_line(mapping=aes(x=sites, y= Mean, group=as.factor(group),color=as.factor(group)))+
  geom_ribbon(aes(ymin=X95.CI_low, ymax=X95.CI_upp,group=as.factor(group), fill=as.factor(group)), alpha=0.2)+
  theme_bw()+
  theme(legend.position = "bottom",axis.text.x = element_blank(), axis.ticks = element_blank())+
  ylab("")+
  xlab("")+
  scale_color_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),name="Method")+
  scale_fill_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),name="Method")+
  ylim(c(0,4))


figure <- ggarrange(gg5, gg6, gg3, gg4, gg1, gg2, ncol=1,nrow=6, common.legend = TRUE, 
          labels = c("(a) Hoverflies - 2017", "(b) Hoverflies - 2018",
                     "(c) Solitary bees - 2017", "(d) Solitary bees - 2018",
                     "(e) Bublebees - 2017", "(f) Bumblebees - 2018"),
          font.label=list(color="black",size=9), 
          hjust = -1, vjust = -0.3)
annotate_figure(figure,
                bottom = text_grob("Site ID", 
                                   hjust = 1),
                left = text_grob("Shannon Index", rot = 90)
)

##################################
# Alpha diversity at locations
###################################
gr_ref <- read_csv("gr_ref.csv")
names(gr_ref)[1] <- "site"
sites <- data.frame(site= site_id)
merged_data_sites <- merge(x=sites, y=gr_ref, by="site", all.x = TRUE, all.y = FALSE)%>%
  select(site, region)%>%
  arrange(site)

#Fill in the missing region for a location
merged_data_sites[62,2] <- "ENG"


# bumblebees
all_data_17 <- data.frame(rbind(bb_shan_inter_17,
                                bb_shan_species_17,
                                bb_shan_genus_17),
                          sites = rep(site_id, 3),
                          group= rep(c("IDM","Species Only", "Genus Only"), 
                                     each=nrow(bb_shan_inter_17)),
                          region = rep(merged_data_sites$region,3))%>%
  group_by(group, region)%>%
  summarize(Mean = mean(Mean, na.rm=TRUE),
            sd = sd(St.Dev.,na.rm = TRUE),
            lower = mean(X95.CI_low,na.rm = TRUE),
            upper = mean(X95.CI_upp,na.rm = TRUE))

all_data_18 <- data.frame(rbind(bb_shan_inter_18,
                                bb_shan_species_18,
                                bb_shan_genus_18),
                          sites = rep(site_id, 3),
                          group= rep(c("IDM","Species Only", "Genus Only"), 
                                     each=nrow(bb_shan_inter_18)),
                          region = rep(merged_data_sites$region,3))%>%
  group_by(group, region)%>%
  summarize(Mean = mean(Mean, na.rm=TRUE),
            sd = sd(St.Dev.,na.rm = TRUE),
            lower = mean(X95.CI_low,na.rm = TRUE),
            upper = mean(X95.CI_upp,na.rm = TRUE))

ggg1 <- ggplot(all_data_17, mapping=aes(x=region, y=Mean, group=group,col=group))+
  geom_point(position=position_dodge(width=0.2))+
geom_errorbar(aes(ymin = (Mean-1.96*sd), ymax = (Mean+1.96*sd)),position = "dodge",width = 0.2)+
  theme_bw()+
  xlab("")+
  ylab("")+
  ylim(c(2.1,3.7))+
  scale_color_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),name="Method")+
  scale_x_discrete(labels=c("ENG" = "England", "SCO" = "Scotland","WAL" = "Wales"))


ggg2 <- ggplot(all_data_18, mapping=aes(x=region, y=Mean, group=group,col=group))+
  geom_point(position=position_dodge(width=0.2))+
  geom_errorbar(aes(ymin = (Mean-1.96*sd), ymax = (Mean+1.96*sd)),position = "dodge",width = 0.2)+
  theme_bw()+
  xlab("")+
  ylab("")+
  ylim(c(2.1,3.7))+
  scale_color_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),name="Method")+
  scale_x_discrete(labels=c("ENG" = "England", "SCO" = "Scotland","WAL" = "Wales"))


# Hoverflies
all_data_17 <- data.frame(rbind(hv_shan_inter_17,
                                hv_shan_species_17,
                                hv_shan_genus_17),
                          sites = rep(site_id, 3),
                          group= rep(c("IDM","Species Only", "Genus Only"), 
                                     each=nrow(bb_shan_inter_17)),
                          region = rep(merged_data_sites$region,3))%>%
  group_by(group, region)%>%
  summarize(Mean = mean(Mean, na.rm=TRUE),
            sd = sd(St.Dev.,na.rm = TRUE),
            lower = mean(X95.CI_low,na.rm = TRUE),
            upper = mean(X95.CI_upp,na.rm = TRUE))

all_data_18 <- data.frame(rbind(hv_shan_inter_18,
                                hv_shan_species_18,
                                hv_shan_genus_18),
                          sites = rep(site_id, 3),
                          group= rep(c("IDM","Species Only", "Genus Only"), 
                                     each=nrow(bb_shan_inter_18)),
                          region = rep(merged_data_sites$region,3))%>%
  group_by(group, region)%>%
  summarize(Mean = mean(Mean, na.rm=TRUE),
            sd = sd(St.Dev.,na.rm = TRUE),
            lower = mean(X95.CI_low,na.rm = TRUE),
            upper = mean(X95.CI_upp,na.rm = TRUE))

ggg3 <- ggplot(all_data_17, mapping=aes(x=region, y=Mean, group=group,col=group))+
  geom_point(position=position_dodge(width=0.2))+
  geom_errorbar(aes(ymin = (Mean-1.96*sd), ymax = (Mean+1.96*sd)),position = "dodge",width = 0.2)+
  theme_bw()+
  xlab("")+
  ylab("")+
  ylim(c(1,3.7))+
  scale_color_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),name="Method")+
  scale_x_discrete(labels=c("ENG" = "England", "SCO" = "Scotland","WAL" = "Wales"))


ggg4 <- ggplot(all_data_18, mapping=aes(x=region, y=Mean, group=group,col=group))+
  geom_point(position=position_dodge(width=0.2))+
  geom_errorbar(aes(ymin = (Mean-1.96*sd), ymax = (Mean+1.96*sd)),position = "dodge",width = 0.2)+
  theme_bw()+
  xlab("")+
  ylab("")+
  ylim(c(1,3.7))+
  scale_color_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),name="Method")+
  scale_x_discrete(labels=c("ENG" = "England", "SCO" = "Scotland","WAL" = "Wales"))


# Solitary bees
all_data_17 <- data.frame(rbind(sb_shan_inter_17,
                                sb_shan_species_17,
                                sb_shan_genus_17),
                          sites = rep(site_id, 3),
                          group= rep(c("IDM","Species Only", "Genus Only"), 
                                     each=nrow(bb_shan_inter_17)),
                          region = rep(merged_data_sites$region,3))%>%
  group_by(group, region)%>%
  summarize(Mean = mean(Mean, na.rm=TRUE),
            sd = sd(St.Dev.,na.rm = TRUE),
            lower = mean(X95.CI_low,na.rm = TRUE),
            upper = mean(X95.CI_upp,na.rm = TRUE))

all_data_18 <- data.frame(rbind(sb_shan_inter_18,
                                sb_shan_species_18,
                                sb_shan_genus_18),
                          sites = rep(site_id, 3),
                          group= rep(c("IDM","Species Only", "Genus Only"), 
                                     each=nrow(bb_shan_inter_18)),
                          region = rep(merged_data_sites$region,3))%>%
  group_by(group, region)%>%
  summarize(Mean = mean(Mean, na.rm=TRUE),
            sd = sd(St.Dev.,na.rm = TRUE),
            lower = mean(X95.CI_low,na.rm = TRUE),
            upper = mean(X95.CI_upp,na.rm = TRUE))

ggg5 <- ggplot(all_data_17, mapping=aes(x=region, y=Mean, group=group,col=group))+
  geom_point(position=position_dodge(width=0.2))+
  geom_errorbar(aes(ymin = (Mean-1.96*sd), ymax = (Mean+1.96*sd)),position = "dodge",width = 0.2)+
  theme_bw()+
  xlab("")+
  ylab("")+
  ylim(c(2.1,3.65))+
  scale_color_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),name="Method")+
  scale_x_discrete(labels=c("ENG" = "England", "SCO" = "Scotland","WAL" = "Wales"))


ggg6 <- ggplot(all_data_18, mapping=aes(x=region, y=Mean, group=group,col=group))+
  geom_point(position=position_dodge(width=0.2))+
  geom_errorbar(aes(ymin = (Mean-1.96*sd), ymax = (Mean+1.96*sd)),position = "dodge",width = 0.2)+
  theme_bw()+
  xlab("")+
  ylab("")+
  ylim(c(2.1,3.65))+
  scale_color_manual(values = c("#FC4E07","#00AFBB", "#E7B800"),name="Method")+
  scale_x_discrete(labels=c("ENG" = "England", "SCO" = "Scotland","WAL" = "Wales"))

figure <- ggarrange(ggg5, ggg6, ggg3, ggg4, ggg1, ggg2, ncol=2,nrow=3, common.legend = TRUE, 
                    labels = c("(a) Hoverflies - 2017", "(b) Hoverflies - 2018",
                               "(c) Solitary bees - 2017", "(d) Solitary bees - 2018",
                               "(e) Bublebees - 2017", "(f) Bumblebees - 2018"),
                    font.label=list(color="black",size=9), 
                    hjust = -1, vjust = -0.3)
annotate_figure(figure,
                bottom = text_grob("Region", 
                                   hjust = 1),
                left = text_grob("Alpha Shannon Diversity", rot = 90)
)



