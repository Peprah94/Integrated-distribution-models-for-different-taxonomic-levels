library(readr)

shanEstimates0 <- read_csv("simEstimates/shanSimulationsEstimates0.csv")
shanEstimates1 <- read_csv("simEstimates/shanSimulationsEstimates1.csv")
shanEstimates2 <- read_csv("simEstimates/shanSimulationsEstimates2.csv")

shanEstimates <- rbind(shanEstimates0,
                       shanEstimates1,
                       shanEstimates2)

shanEst <- shanEstimates%>%
  ggplot()+
  geom_boxplot(aes(x = model, y = V1, fill = truth))+
 # facet_grid(~truth)+
  ylab("MSE(H')")+
  theme_classic()+
  theme(legend.title = element_blank(),
        legend.position = "bottom")

#text = element_text(size = 20),

# save shan est
ggsave(filename = "Figures/simHills.svg",
       shanEst,
       width = 10,
       height = 8,
       units = "in")

parEstimates <- read_csv("Figures/parEstimates.csv")

colnames(parEstimates)

parsEst <- parEstimates%>%
  mutate(muLatitude = muLatitude + 0.5,
         sigmaLatitude = sigmaLatitude -1,
         sigmaLambda = sigmaLambda - 0.2,
         sigmaAlphaSites = sigmaAlphaSites - 0.3,
         sigmaAlphaSpecies = sigmaAlphaSpecies - 1)%>%
  melt(id.vars = c("model", "truth"))%>%
  ggplot()+
  geom_boxplot(mapping = aes(x = variable, y = value, fill = model))+
  facet_grid(~truth) + ylab("Bias of Parameters")+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  theme_classic()+
  theme(text = element_text(size = 20),
        legend.title = element_blank(),
        legend.position = "bottom")

ggsave(filename = "Figures/simBiasParsEsts.png",
       shanEst,
       width = 10,
       height = 8,
       units = "in")
