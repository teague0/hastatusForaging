# Figure 4. Estimated daily energy expenditure and the foraging returns to this DEE. A) Daily energy expenditure relative to the amount of time the bat was out of the cave, with comparative data from [27]. The line shows the relationship derived from the bats tracked in this study. B) The number of flowers needed based on average energy content to support the energy requirements of bats predict the number of feeding clusters used. C) Social network strength does not predict daily energy expenditure, but D) increases with the amount of time tracked per night. 

library(tidyverse)
library(lme4)
library(cowplot)
theme_set(theme_cowplot())

load("./data/13_NightSumValues.Rdata")

#Flower power: calculate energy potential from balsa trees. This is the rationale for information calculated for Figure 4. 

#Derive the energy content per ml of balsa nectar from Kays et al 2012
# - Total nectar produced by a flower is *estimated* at 25.5 ml
# - Each flower produced 11.6 kcal per night
# - Large peak is 60 flowers per patch, but there is a normal distribution of the number of flowers open per patch over the season.
# - Flowers open with 4.9 ± 1.3 ml of nectar. There is a sharp decline in nectar production over the night. It starts at 2.6 ml / hour and then declines to 0.4 ml per hour.
# - Nectar sugar concentration decreases over the night from 13.3% at 18h to 7.9 % at 06 h. This an average of 12.4% sugar. This is 0.124 g sucrose per ml nectar
# - Roland estimates 11.75 kcal over 1 night (but we don't care about the full night)
# - 3.87 kcal/g sugar
# - 0.47988 kcal per ml (3.87*0.124)
# - 1 kcal = 4.184 kJ
# - 2.007818 kJ / ml nectar
# - peak flower availability is about 60 per patch or 5 per m2. Mean is more like 20 & 2.5 flowers per m2
# - 602.3454 kJ per crown (2.007818 kJ   / ml * 5 ml / flower * 60 flowers per crown) 

#If we assume that bats drain 5 ml of the most profitable nectar

#How many ml of nectar & # of flowers would a bat need to power RMR + flight power for the day
#This can be added to the nightNets summary table

#Totals needed per night.
nightNets$mlRequired <- nightNets$dee.kJ / 2.007818
nightNets$NFlowers <-  nightNets$mlRequired / 5


#Calculate Cumulative energy expenditure per bat day.
nas <- which(is.na(hastMorph$PmetTotal))
hastMorph$PmetTotal[nas] <- 0
tlagNas <- which(is.na(hastMorph$tlag))
hastMorph$tlag[tlagNas] <- 0

batDay <- split(hastMorph, f=hastMorph$batIDday)

cumEE <- lapply(batDay, function(x){
  x$cumEE.kJ <- cumsum(x$PmetTotal)/1000
  x$timeElapsed.s <-  cumsum(x$tlag)
  return(x)
})

hastMorph <- do.call("rbind", cumEE)
mlNeeded <- hastMorph$cumEE.kJ / 2.007818
flwrNeeded <- mlNeeded / 5

kunz <- read.csv("./data/Kunz1998_DEE.csv")
mycols <- viridisLite::viridis(6)[c(3, 4, 6)]

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
plot.timeI <- data.frame(timeI@frame, fitted.re = fitted(timeI))


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
