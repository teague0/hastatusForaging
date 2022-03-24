#Figure 3. Relationships between proximity network strength centrals of individual bats and the proportion of time that they spend per night A) Resting, B) Moving between patch, C) slow foraging and feeding flight, and D) Commuting.

library(tidyverse)
library(lme4)
library(car)
library(cowplot)
library(effects)
library(MuMIn)
theme_set(theme_cowplot())

load("./data/13_NightSumValues.Rdata") #nightNets
#mycols <- viridisLite::viridis(6)[c(3, 4, 6)]
mycols <- viridisLite::plasma(9)[c(1, 5, 8)]

#state 1:rest ; state 2:slow flight, state 3:move;  state 4: commute
actNet.r <- glmer(state1~strength + (1|batID:groupID), data = nightNets, family = binomial)

actNet.r.d <- glmer(state1~degree + (1|batID:groupID), data = nightNets, family = binomial)
summary(actNet.r)
Anova(actNet.r) #There is no effect here
anova(actNet.r, type = 2, ddf = "Satterthwaite")
plot.actNet.r <- data.frame(actNet.r@frame, fitted.re = fitted(actNet.r))
r.squaredGLMM(actNet.r) #0.28528124 0.28528124

actNet.s <- glmer(state2~strength + (1|batID), data = nightNets, family = binomial)
plot.actNet.s <- data.frame(actNet.s@frame, fitted.re = fitted(actNet.s))
summary(actNet.s)
Anova(actNet.s) #There is no effect here
r.squaredGLMM(actNet.s) #0.05285176 0.05285177

actNet.m <- glmer(state3~strength + (1|batID), data = nightNets, family = binomial)
plot.actNet.m <- data.frame(actNet.m@frame, fitted.re = fitted(actNet.m))
summary(actNet.m)
Anova(actNet.m) #There is no effect here
r.squaredGLMM(actNet.m) #0.045154005 0.045154005

actNet.c <- glmer(state4~strength + (1|batID), data = nightNets, family = binomial)
plot.actNet.c <- data.frame(actNet.c@frame, fitted.re = fitted(actNet.c))
summary(actNet.c)
Anova(actNet.c) #There is no effect here
r.squaredGLMM(actNet.c) #0.9999959 0.9999959

rest <- ggplot()+
  geom_point(data = nightNets, 
             aes(x = strength, y = state1, color=groupID), size = 3)+
  geom_smooth(data = plot.actNet.r,
              aes(x = strength, y = fitted.re), 
              method = "lm", 
              se = TRUE, level = 0.95, 
              color = "black", lty = 2)+
    theme(legend.position = "none",
        legend.background = element_blank())+
  scale_color_manual(values = mycols, 
                     name = "", 
                     breaks = c("blue", "brown", "yellow"),
                     labels = c("Group 1", "Group 2", "Group 3"))+
  ylim(0,0.6)+
  theme(legend.position = c(0.6, 0.9), 
        legend.background = element_blank(),
        legend.title = element_blank())+ 
  labs(y = "proportion in rest", 
       x = "strength centrality")

forage <- ggplot()+
  geom_point(data = nightNets, 
             aes(x = strength, y = state2, color=groupID), size = 3)+
  geom_smooth(data = plot.actNet.s,
              aes(x = strength, y = fitted.re), 
              method = "lm", 
              se = TRUE, level = 0.95, 
              color = "black", lty = 2)+
  scale_color_manual(values = mycols, name="Group")+
  ylim(0,0.6)+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(y = "proportion in slow flight", 
       x = "strength centrality")

moveF <- ggplot()+
  geom_point(data = nightNets, 
             aes(x = strength, y = state3, color=groupID), size = 3)+
  scale_color_manual(values = mycols, name="Group")+
  geom_smooth(data = plot.actNet.m,
              aes(x = strength, y = fitted.re), 
              method = "lm", se = F, 
              color = "black", lty = 2)+
  ylim(0,0.6)+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(y = "proportion in move", 
       x = "strength centrality")

commute <- ggplot()+
  geom_point(data = nightNets, 
             aes(x = strength, y = state4, color=groupID), size = 3)+
  geom_smooth(data = plot.actNet.c,
              aes(x = strength, y = fitted.re), 
              method = "lm", 
              se = TRUE, level = 0.95, 
              color = "black", lty = 2)+
  scale_color_manual(values = mycols, name="Group")+
  ylim(0,0.6)+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(y = "proportion in commute", 
       x = "strength centrality")

pdf("./output/Fig S5 - ActBudgNetStr_wide_NoEffects.pdf",  width = 12, height = 4)
plot_grid(rest, moveF, forage, commute, labels = c('A', 'B', 'C', 'D'), label_size = 12, rows = 1)
dev.off()


