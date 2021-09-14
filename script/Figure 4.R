#Figure 4

load("./processed data/NightlySummaries 2020-08-04.Rdata")

kunz <- read.csv("./raw data/Kunz1998_DEE.csv")

bocasPlot <- nightNets %>% dplyr::select(groupID, batID, timeTrack.min, dee.kJ)
kunzPlot <- kunz %>% dplyr::select("groupID" = Harem, "batID"= bat1D, "dee.kJ"= dee.kJday, "timeTrack.min" = mins.out)
kunzPlot$groupID <- rep("Kunz")
bocasKunzPlot <- bocasPlot %>% bind_rows(kunzPlot)

t.m <- lmer(dee.kJ~timeTrack.min+(1|batID), data=nightNets)
plot.t.m <- data.frame(t.m@frame, fitted.re = fitted(t.m))

kunzPlotCols <- c(mycols, "#FF0000")
timeDEE <- ggplot(bocasKunzPlot, aes(x = timeTrack.min, y = dee.kJ, color = groupID))+
  geom_point(size=3)+
  scale_color_manual(values = kunzPlotCols,
                     breaks = c("blue", "brown", "yellow", "Kunz"),
                     labels = c("Group 1", "Group 2", "Group 3", "Kunz et al\n1998" ))+
  geom_smooth(data = plot.t.m, 
              aes(x = timeTrack.min, y = fitted.re), method = "lm", se = F, color = "black")+
  theme(legend.position = c(0.07, 0.9),
        legend.background = element_blank(), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12))+
  labs(x = "Time tracked per night (min)", 
       y = "Daily energy\nexpenditure (kJ)")


flwr.m <- lmer(nFeedingClusters~NFlowers+(1|batID), data=nightNets)
summary(cl.m)
Anova(flwr.m)
# Response: nFeedingClusters
# Chisq Df Pr(>Chisq)    
# NFlowers 22.22  1  2.432e-06 ***
r.squaredGLMM(flwr.m)

plot.flwr.m <- data.frame(flwr.m@frame, fitted.re = fitted(flwr.m))

metFlwrNClust <- ggplot()+
  geom_point(data = nightNets, 
             aes(x = NFlowers, y = nFeedingClusters, color = groupID),size=3)+
  scale_color_manual(values = mycols, name = "Group")+
  geom_smooth(data = plot.flwr.m,
              aes(x = NFlowers, y = fitted.re), method = "lm", se = F, color = "black")+
  #geom_abline(slope = 1, intercept = 0)+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(x = "Number of full flowers needed", 
       y = "Number of feeding clusters used")


strDEE <- ggplot()+
  geom_point(data= nightNets, aes(x = strength, y = dee.kJ, color = groupID), size=3)+
  scale_color_manual(values = mycols, name = "Group")+
  theme(legend.position = "none",
        legend.background = element_blank(),
        legend.text = element_text(size = 8))+
  labs(x = "Network strength", 
       y = "Daily energy\n expenditure (kJ)")


timeI <- lmer(timeTrack.min~strength+(1|batID), data=nightNets)
Anova(timeI)
# Response: timeTrack.min
# Chisq Df Pr(>Chisq)  
# strength 6.0706  1    0.01374 *

plot.timeI <- data.frame(timeI@frame, fitted.re = fitted(timeI))
head(plot.timeI)

strTime <- ggplot()+
  geom_point(data=nightNets, aes(y = timeTrack.min, x = strength, color = groupID), size=3)+
  scale_color_manual(values = mycols)+
  geom_smooth(aes(y = fitted.re, x=strength), 
              method = "lm", data = plot.timeI, se=FALSE, color = "black")+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(y = "Time tracked per night (min)", 
       x = "Network strength")

pdf("./output/Fig 4 - DEE_NetStr_Flowers.pdf", width = 8, height = 8)
plot_grid(timeDEE, metFlwrNClust, strDEE, strTime, labels = c('A', 'B', 'C', 'D'), label_size = 12, ncol=2)
dev.off()
