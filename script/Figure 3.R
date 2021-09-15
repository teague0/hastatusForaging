#Figure 3
library(tidyverse)
library(lme4)
library(car)
library(cowplot)
theme_set(theme_cowplot())

load("./processed data/NightlySummaries 2020-08-04.Rdata")

actNet.r <- lmer(state1~strength + (1|batID), data = nightNets)
Anova(actNet.r)
# Response: state1
# Chisq Df Pr(>Chisq)  
# strength 4.4804  1    0.03429 *
plot.actNet.r <- data.frame(actNet.r@frame, fitted.re = fitted(actNet.r))
actNet.m <- lmer(state3~strength + (1|batID), data = nightNets)
Anova(actNet.m)
# Response: state3
# Chisq Df Pr(>Chisq)  
# strength 4.5993  1    0.03198 *
plot.actNet.m <- data.frame(actNet.m@frame, fitted.re = fitted(actNet.m))

rest <- ggplot()+
  geom_point(data = nightNets, 
             aes(x = strength, y = state1, color=groupID), size = 3)+
  geom_smooth(data = plot.actNet.r,
              aes(x = strength, y = fitted.re), method = "lm", se = F, color = "black")+
  theme(legend.position = "none",
        legend.background = element_blank())+
  scale_color_manual(values = mycols, 
                     name = "", 
                     breaks = c("blue", "brown", "yellow"),
                     labels = c("Group 1", "Group 2", "Group 3"))+
  theme(legend.position = c(0.03, 0.9), 
        legend.background = element_blank(),
        legend.title = element_blank())+ 
  labs(y = "Proportion in rest", 
       x = "Network strength")
forage <- ggplot()+
  geom_point(data = nightNets, 
             aes(x = strength, y = state2, color=groupID), size = 3)+
  scale_color_manual(values = mycols, name="Group")+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(y = "Proportion in slow flight", 
       x = "Network strength")
moveF <- ggplot()+
  geom_point(data = nightNets, 
             aes(x = strength, y = state3, color=groupID), size = 3)+
  scale_color_manual(values = mycols, name="Group")+
  geom_smooth(data = plot.actNet.m,
              aes(x = strength, y = fitted.re), method = "lm", se = F, color = "black")+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(y = "Proportion in move", 
       x = "Network strength")
commute <- ggplot()+
  geom_point(data = nightNets, 
             aes(x = strength, y = state4, color=groupID), size = 3)+
  scale_color_manual(values = mycols, name="Group")+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(y = "proportion in commute", 
       x = "network strength")