#Figure S2. The number of foraging patches used is predicted by A) the distance that individuals moved each night and B) the total time that individuals were tracked each night. The number of feeding clusters identified was not predicted by either C) the distance moved per night or D) the number of foraging patches used. 

library(tidyverse)
library(lme4)

load("./data/13_NightSumValues.Rdata")

ptchTdis.m <- lmer(nPatches~totalDistance + (1|batID), data=nightNets)
plot.ptchTdis.m <- data.frame(ptchTdis.m@frame, fitted.re = fitted(ptchTdis.m))

patchDist <- ggplot(nightNets, aes(x = totalDistance, y = nPatches, color = groupID))+
  geom_point(size=3)+
  scale_color_manual(values = mycols, 
                     name = "", 
                     breaks = c("blue", "brown", "yellow"),
                     labels = c("Group 1", "Group 2", "Group 3"))+
  theme(legend.position = c(0.02, 0.94), 
        legend.background = element_blank(),
        legend.title = element_blank())+
  geom_smooth(data = plot.ptchTdis.m,
              aes(x =totalDistance, y = fitted.re), method = "lm", se = F, color = "black")+
  labs(x = "Distance moved per night (km)", 
       y = "Number of patches")

cp.m <- lmer(nPatches~timeTrack.min+(1|batID), data=nightNets)
plot.cp.m <- data.frame(cp.m@frame, fitted.re = fitted(cp.m))

patchTime <- ggplot(nightNets, aes(x = timeTrack.min, y = nPatches, color = groupID))+
  geom_point(size=3)+
  scale_color_manual(values = mycols, name = "Group")+
  geom_smooth(data = plot.cp.m,
              aes(x =timeTrack.min, y = fitted.re), method = "lm", se = F, color = "black")+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(x = "Time tracked per night (min)", 
       y = "Number of foraging patches")

feedDist <- ggplot(nightNets, aes(x = totalDistance, y = nFeedingClusters, color = groupID))+
  geom_point(size=3)+
  scale_color_manual(values = mycols, name = "Group")+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(x = "Distance moved per night (km)", 
       y = "Number of feeding clusters")
cl.m <- lmer(nFeedingClusters~totalDistance+(1|batID), data=nightNets)

clusPatch <- ggplot(nightNets, aes(x = nPatches, y = nFeedingClusters, color = groupID))+
  geom_point(size=3)+
  scale_color_manual(values = mycols, name = "Group")+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(x = "Number of foraging patches", 
       y = "Number of feeding clusters")
clusPat.m <- lmer(nFeedingClusters~nPatches+(1|batID), data=nightNets)

pdf("./output/Fig S2 - PatchClusterTimeDist.pdf", width = 8, height = 8)
plot_grid(patchDist, patchTime, feedDist, clusPatch,labels = c('A', 'B', 'C', 'D'), label_size = 12,  ncol=2)
dev.off()