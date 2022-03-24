#Figure 3. Relationships between proximity network strength centrality of individual bats and A) time tracked per night, B) foraging efficiency, C) DEE.

library(tidyverse)
library(lme4)
library(car)
library(cowplot)
library(effects)
theme_set(theme_cowplot())

load("./data/13_NightSumValues.Rdata")
#mycols <- viridisLite::viridis(6)[c(3, 4, 6)]
mycols <- viridisLite::plasma(9)[c(1, 5, 8)]

#state 1:rest ; state 2:slow flight, state 3:move;  state 4: commute
actNet.r <- glmer(state1~strength + (1|batID:groupID), data = nightNets, family = binomial)
summary(actNet.r)
Anova(actNet.r) #There is no effect here
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
  # geom_ribbon(data = df.ef,
  #             aes(x = strength, 
  #                  ymax = upper,
  #                  ymin = lower),
  #             alpha = 0.4, fill = 'blue')+
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
  labs(y = "Proportion in rest", 
       x = "Network strength")

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
  labs(y = "Proportion in slow flight", 
       x = "Network strength")

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
  labs(y = "Proportion in move", 
       x = "Network strength")

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
       x = "network strength")

pdf("./output/Fig 3 - ActBudgNetStr_wide_NoEffects.pdf",  width = 12, height = 4)
plot_grid(rest, moveF, forage, commute, labels = c('A', 'B', 'C', 'D'), label_size = 12, rows = 1)
dev.off()


