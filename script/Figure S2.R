#Figure S2. The number of foraging patches used is predicted by A) the distance that individuals moved each night and B) the total time that individuals were tracked each night. The number of feeding clusters identified was not predicted by either C) the distance moved per night or D) the number of foraging patches used. 

library(tidyverse)
library(lme4)
library(MuMIn)

load("./data/13_NightSumValues.Rdata")
mycols <- viridisLite::plasma(9)[c(1, 5, 8)]


ptchTdis.m <- lmer(nPatches~totalDistance + (1|batID:groupID), data=nightNets)
ptchTdis.m2 <- lmer(nPatches~poly(totalDistance,2) + (1|batID:groupID), data=nightNets)
AICc(ptchTdis.m, ptchTdis.m2) #poly wins
plot.ptchTdis.m2 <- data.frame(ptchTdis.m2@frame, fitted.re = fitted(ptchTdis.m2))
Anova(ptchTdis.m)
anova(ptchTdis.m, type = 2, ddf = "Satterthwaite")
summary(ptchTdis.m)
r.squaredGLMM(ptchTdis.m)

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
              aes(x =totalDistance, y = fitted.re), method = "lm", formula = y~poly(x,2), se = T, color = "black")+
  labs(x = "Distance moved per night (km)", 
       y = "Number of patches")

cp.m <- lmer(nPatches~timeTrack.min+(1|batID), data=nightNets)
plot.cp.m <- data.frame(cp.m@frame, fitted.re = fitted(cp.m))

patchTime <- ggplot(nightNets, aes(x = timeTrack.min, y = nPatches, color = groupID))+
  geom_point(size=3)+
  scale_color_manual(values = mycols, name = "Group")+
  geom_smooth(data = plot.cp.m,
              aes(x =timeTrack.min, y = fitted.re), method = "lm", se = T, color = "black")+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(x = "Time tracked per night (min)", 
       y = "Number of foraging patches")

flwrDis.m2 <- lmer(nFeedingClusters~poly(totalDistance,2) + (1|batID:groupID), data=nightNets)
AICc(ptchTdis.m, ptchTdis.m2) #poly wins
Anova(flwrDis.m2)
summary(flwrDis.m2)
r.squaredGLMM(flwrDis.m2)

feedDist <- ggplot(nightNets, aes(x = totalDistance, y = nFeedingClusters, color = groupID))+
  geom_point(size=3)+
  scale_color_manual(values = mycols, name = "Group")+
  theme(legend.position = "none",
        legend.background = element_blank())+
  geom_smooth(method = "lm", formula = y~poly(x,2), se = T, color = "black", lty=2)+
  labs(x = "Distance moved per night (km)", 
       y = "Number of feeding clusters")

cl.m <- lmer(nFeedingClusters~totalDistance+(1|batID), data=nightNets)

clusPatch <- ggplot(nightNets, aes(x = nPatches, y = nFeedingClusters, color = groupID))+
  geom_point(size=3)+
  scale_color_manual(values = mycols, name = "Group")+
  theme(legend.position = "none",
        legend.background = element_blank())+
  geom_smooth(method = "lm", se = T, color = "black", lty=2)+
  labs(x = "Number of foraging patches", 
       y = "Number of feeding clusters")
clusPat.m <- lmer(nFeedingClusters~nPatches+(1|batID), data=nightNets)

pdf("./output/Fig S2 - PatchClusterTimeDist.pdf", width = 8, height = 8)
plot_grid(patchDist, patchTime, feedDist, clusPatch,labels = c('A', 'B', 'C', 'D'), label_size = 12,  ncol=2)
dev.off()