library(dplyr)

load("/Volumes/kwakupa-1/IDM_new/Data/cross_validation/bumblebeesIDMCV.RData")
load("/Volumes/kwakupa-1/IDM_new/Data/cross_validation/bumblebeesIDM1CV.RData")
load("/Volumes/kwakupa-1/IDM_new/Data/cross_validation/bumblebeesIDM2CV.RData")
load("/Volumes/kwakupa-1/IDM_new/Data/cross_validation/bumblebeesIGCV.RData")
load("/Volumes/kwakupa-1/IDM_new/Data/cross_validation/bumblebeesSpeCV.RData")

bb <- do.call("rbind", bumblebees)[,1:4]%>%
  as.data.frame()%>%
  select(CVvalue)
bb <- rep(unlist(bb$CVvalue), 5)
all_bumblebees <- do.call("rbind", c(bumblebees,
                                     bumblebees1,
                                     bumblebees2,
                                     bumblebees3,
                                     bumblebees4))[,1:4]%>%
  as.data.frame()%>%
  mutate(CVvalue = unlist(CVvalue),
    change = bb - CVvalue )


# hoverflies
load("/Volumes/kwakupa-1/IDM_new/Data/cross_validation/hoverfliesIDMCV.RData")
load("/Volumes/kwakupa-1/IDM_new/Data/cross_validation/hoverfliesIDM1CV.RData")
load("/Volumes/kwakupa-1/IDM_new/Data/cross_validation/hoverfliesIDM2CV.RData")
load("/Volumes/kwakupa-1/IDM_new/Data/cross_validation/hoverfliesIGCV.RData")
load("/Volumes/kwakupa-1/IDM_new/Data/cross_validation/hoverfliesSpeCV.RData")
hv <- do.call("rbind", hoverflies)[,1:4]%>%
  as.data.frame()%>%
  select(CVvalue)
hv <- rep(unlist(hv$CVvalue), 5)

all_hoverflies <- do.call("rbind", c(hoverflies,
                                     hoverflies1,
                                     hoverflies2,
                                     hoverflies3,
                                     hoverflies4))[,1:4]%>%
  as.data.frame()%>%
  mutate(CVvalue = unlist(CVvalue),
         change = CVvalue - hv)

# solitarybees
load("/Volumes/kwakupa-1/IDM_new/Data/cross_validation/solitarybeesIDMCV.RData")
load("/Volumes/kwakupa-1/IDM_new/Data/cross_validation/solitarybeesIDM1CV.RData")
load("/Volumes/kwakupa-1/IDM_new/Data/cross_validation/solitarybeesIDM2CV.RData")
load("/Volumes/kwakupa-1/IDM_new/Data/cross_validation/solitarybeesIGCV.RData")
load("/Volumes/kwakupa-1/IDM_new/Data/cross_validation/solitarybeesSpeCV.RData")

sb <- do.call("rbind", solitarybees)[,1:4]%>%
  as.data.frame()%>%
  select(CVvalue)
sb <- rep(unlist(sb$CVvalue), 5)

all_solitarybees <- do.call("rbind", c(solitarybees,
                                     solitarybees1,
                                     solitarybees2,
                                     solitarybees3,
                                     solitarybees4))[,1:4]%>%
  as.data.frame()%>%
  mutate(CVvalue = unlist(CVvalue),
         change = CVvalue - sb)
