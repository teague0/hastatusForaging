#This version of night summaries builds on #8 (recalculation of time lags, speeds, etc) and 10 to estimate energy consumption. 10 also reclassified some points into state1 based on the long time lags.

library(tidyverse)
library(sp)
library(sf)
library(lubridate)
library(geosphere)
library(raster)
library(maptools)


load("./processed data/HastatusSegStateBiodatPwr 2020-08-09.Rdata") #hastMorph
load("./processed data/LaGrutaCoords.Rdata")


#### Data Product from this script:
load("./processed data/NightlySummaries 2020-08-04.Rdata") #nightNets
load("./processed data/NightSumValues 2021-05-27.Rdata") #nightNets


#Recalculate time budgets ####
commuteForageTime <- hast_df %>% group_by(groupID, batID, date, behaviour) %>% 
  summarise(nlocsBehav = n(),
            sumTimeBehav = sum(tlag, na.rm = TRUE)) %>% 
  mutate(freqComFor = nlocsBehav / sum(nlocsBehav))
ggplot(actBudg, aes(x = behaviour, y = freqBehav))+geom_point()

actBudg <- hast_df %>% group_by(groupID, batID, batID_day, date, newState) %>% 
  summarise(nlocsBehav = n(),
            sumTimeBehav = sum(tlag, na.rm = TRUE)) %>% 
  mutate(freqAct = nlocsBehav / sum(nlocsBehav))

actBudg_p <- actBudg %>% filter(!is.na(newState)) %>% 
  ggplot( )+
  geom_boxplot(aes(x = newState, y = freqAct, fill = groupID))+
  scale_x_discrete(breaks = c("state1", "state2", "state3", "state4"),
                   labels = c("rest", "slow flight" ,"move", "commute"))+
  scale_fill_manual(values = rep(mycols, 3), 
                     name = "Group")+
  theme_bw()+
  theme(legend.position = c(0.15, 0.8),
        legend.background = element_blank())+
  labs(x = "", 
       y = "frequency", 
       subtitle =  "Nightly activity budgets from HMM")

behaviorStates <- c("rest", "slow flight" ,"move", "commute")
names(behaviorStates) <- c("state1", "state2", "state3", "state4")

library(ggridges)
speedBehav <- hast_df %>% filter(!is.na(newState)) %>% 
  ggplot()+
  geom_density_ridges(aes(x = ground_speedRecalc, y = newState, fill = newState))+
  scale_y_discrete(breaks = c("state1", "state2", "state3", "state4"),
                   labels = c("rest", "slow flight" ,"move", "commute"))+
  xlim(0,20)+
  theme_ridges()+
  scale_fill_viridis_d(option = "C", alpha = 0.7)+
  theme(legend.position = "none")+
  labs(x = expression(paste("ground speed m ", s^-1, sep="")),
       y = "behavioral state", 
       subtitle = "ground speeds per activity")+
  facet_wrap(~groupID)

airspeedBehav <- hast_df %>% filter(!is.na(newState)) %>% 
  ggplot()+
  geom_density_ridges(aes(x = airspeedRecalc, y = newState, fill = newState))+
  scale_y_discrete(breaks = c("state1", "state2", "state3", "state4"),
                   labels = c("rest", "slow flight" ,"move", "commute"))+
  xlim(0,20)+
  theme_ridges()+
  scale_fill_viridis_d(option = "C", alpha = 0.7)+
  theme(legend.position = "none")+
  labs(x = expression(paste("airspeed m ", s^-1, sep="")),
       y = "behavioral state", 
       subtitle = "airspeeds per activity")+
  facet_wrap(~groupID)


pdf("./output/ActBudgSpeed.pdf",  width = 8, height = 4)
plot_grid(actBudg_p, speedBehav, rel_widths = c(1,1.5))
dev.off()

pdf("./output/Speeds.pdf",  width = 8, height = 4)
plot_grid(speedBehav, airspeedBehav)
dev.off()



##Night sums ####
nightSums <- hast_df %>% group_by(groupID, batID, date, batID_day) %>% 
  dplyr::summarise(minTimeUTC = min(timestamps),
            maxTimeUTC = max(timestamps),
            minTimeLocal = min(with_tz(timestamps, tz="America/Panama")),
            maxTimeLocal = max(with_tz(timestamps, tz="America/Panama")),
            timeTrack.min = round(as.numeric(as.duration(maxTimeLocal-minTimeLocal))/60, 2),
            timeLeftDay.min = 1440 - timeTrack.min,
            nlocs = n(),
            meanTimeLag.s = mean(tlag, na.rm = TRUE),
            maxTimeLag.s = max(tlag, na.rm = TRUE),
            totalDistance = sum(stepLength, na.rm = TRUE)/1000,
            nPatches = length(unique(patchID)), 
            nClus1 = length(unique(dbClus, na.rm=TRUE)),
            nClus2 = length(unique(dbClus_feed2, na.rm=TRUE)),
            nClus3 = length(unique(dbClus_feed3, na.rm=TRUE)),
            nClus4 = length(unique(dbClus_feed4, na.rm=TRUE)),
            Pmech.kJ = sum(pwrJ, na.rm = TRUE)/1000, 
            Pmet.kJ = sum(PmetTotal, na.rm = TRUE)/1000,
            tlagSum.min = sum(tlag, na.rm = TRUE)/60,
            bodymass = unique(mass.deploy)) %>% 
  mutate(Pmet.unTracked.kJ = (timeLeftDay.min * 23.8/60 * bodymass)/1000, 
         dee.kJ = Pmet.unTracked.kJ + Pmet.kJ, 
         nFeedingClusters = nClus2 + nClus3 + nClus4)


#can I also add in time to first patch & first flower?

firstPtchPt <- hast_df %>% 
  filter(!is.na(patchID)) %>% 
  group_by(groupID, batID, date, batID_day) %>%
  slice(which.min(timestamp)) %>% 
  dplyr::select(batID_day, firstPatch = patchID, firstPatchTimeUTC = timestamp, firstPatchDist = PtdistToCave)

firstFlwPt <- hast_df %>% 
  group_by(groupID, batID, date, batID_day, patchID) %>%
  gather("clustNum", "value", dbClus_feed2:dbClus_feed4 ) %>% 
  filter(value == 1) %>% 
  slice(which.min(timestamp)) %>% 
  dplyr::select(batID_day, firstFlwPatch = patchID, firstFlwrID = value, firstFlwrTimeUTC = timestamp, firstFlwrDist = PtdistToCave)

firstPatchFlw <- left_join(firstPtchPt, firstFlwPt)
nightSums <- left_join(nightSums, firstPatchFlw)

#Add in the social network metrics to the daily summaries
load("./processed data/ObservedSocNetworkMetrics.Rdata") #observed
nightNets <- nightSums %>% left_join(observed, by = c("batID" = "ID"))

startDist <- hast_df %>% group_by(batID_day) %>%
  filter(timestamps == min(with_tz(timestamps, tz="America/Panama"))) %>%
  dplyr::select(batID_day, firstFixDist = PtdistToCave)
lastDist <- hast_df %>% group_by(batID_day) %>%
  filter(timestamps == max(with_tz(timestamps, tz="America/Panama"))) %>%
  dplyr::select(batID_day, lastFixDist = PtdistToCave)
firstLast <- startDist %>% left_join(lastDist) 

nightNets <- nightNets %>% left_join(firstLast)
nightNets <- nightNets %>% mutate(firstFixDist = firstFixDist/1000,
                                  lastFixDist = lastFixDist/1000,
                                  elapsedTimePatch = firstPatchTimeUTC - minTimeUTC,
                                  elapsedTimeFlw = firstFlwrTimeUTC - minTimeUTC,
                                  trackDistPatch = firstFixDist - firstPatchDist/1000,
                                  trackDistFlw = firstFixDist - firstFlwrDist/1000)
save(nightNets, file = "./processed data/NightSumValues 2021-05-27.Rdata")


#How do the energy estimates compare to the Kunz 1998 data?
kunz <- read.csv("./raw data/Kunz1998_DEE.csv")
mycols <- viridisLite::viridis(6)[c(3, 4, 6)]
shape <- c("square", "circle", "triangle")

library(lme4)
library(car)
library(cowplot)
theme_set(theme_cowplot())

###Explore effects of time tracking ####
tDist.m <- lmer(totalDistance~timeTrack.min+(1|batID), data=nightNets) #a curve needs to be fit 
Anova(tDist.m)
# Response: totalDistance
# Chisq Df Pr(>Chisq)    
# timeTrack.min 49.601  1  1.884e-12 ***
summary(tDist.m)

plot.tDist.m <- data.frame(tDist.m@frame, fitted.re = fitted(tDist.m))
head(plot.tDist.m)

timeDist <- ggplot()+
  geom_point(data=nightNets, 
             aes(x = timeTrack.min, y = totalDistance, color = groupID), size=3)+
  scale_color_manual(values = mycols, name = "Group")+
  #geom_point(data=kunz, aes(x = mins.out, y = dee.kJday), color = "red")+
  geom_smooth(data = plot.tDist.m,
              aes(x = timeTrack.min, y = fitted.re),method = "lm",se = FALSE, color = "black")+
  theme_bw()+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(x = "time tracked per night (min)", 
       y = "total distance moved\nper night (km)", 
       subtitle =  "Total distance increases\nwith tracking time")


firstTime <- ggplot(nightNets, aes(y = timeTrack.min, x = firstFixDist, color = groupID))+
  geom_point(size=3)+
  scale_color_manual(values = mycols, name = "Group")+
  theme_bw()+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(y = "time tracked per night (min)", 
       x = "first fix distance\nfrom roost (km)", 
       subtitle =  "Distance of first fix of the night\nis independent of time tracked")
fd.m <-  lmer(timeTrack.min~firstFixDist+(1|batID), data=nightNets)
Anova(fd.m) #p=0.07


ld.m <-  lmer(timeTrack.min~lastFixDist+(1|batID), data=nightNets)
Anova(ld.m) #p<0.0017
plot.ld.m <- data.frame(ld.m@frame, fitted.re = fitted(ld.m))
head(plot.ld.m)

lastTime <- ggplot()+
  geom_point(data = nightNets, 
             aes(y = timeTrack.min, x = lastFixDist, color = groupID), size=3)+
  scale_color_manual(values = mycols, name = "Group")+
  geom_smooth(data = plot.ld.m, 
              aes(x = lastFixDist, y = fitted.re), method = "lm", se = F, color = "black")+
  theme_bw()+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(y = "time tracked per night (min)", 
       x = "last fix distance\nfrom roost (km)", 
       subtitle =  "Animals with far last fix\nare tracked for less time")

t.m <- lmer(dee.kJ~timeTrack.min+(1|batID), data=nightSums)
t.m2 <- lmer(dee.kJ~timeTrack.min^2+(1|batID), data=nightSums)
summary(t.m)
Anova(t.m2)
# Response: dee.kJ
# Chisq Df Pr(>Chisq)    
# timeTrack.min 46.266  1  1.032e-11 ***
r.squaredGLMM(t.m)
#           R2m       R2c
# [1,] 0.5452688 0.8711667

plot.t.m <- data.frame(t.m@frame, fitted.re = fitted(t.m))
plot.t.m2 <- data.frame(t.m2@frame, fitted.re = fitted(t.m2))

timeDEE <-ggplot(nightNets, aes(x = timeTrack.min, y = dee.kJ, color = groupID))+
  geom_point(size=3)+
  scale_color_manual(values = mycols, name = "Group")+
  geom_point(data=kunz, aes(x = mins.out, y = dee.kJday), color = "red")+
  geom_smooth(data = plot.t.m, 
              aes(x = timeTrack.min, y = fitted.re), method = "lm", se = F, color = "black")+
  theme_bw()+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(x = "time tracked per night (min)", 
       y = "daily energy\nexpenditure (kJ)", 
       subtitle =  "DEE increases with tracking time\nKunz data in red")

e.m <- lmer(Pmet.kJ~timeTrack.min+(1|batID), data=nightNets)
summary(e.m)
Anova(e.m)
# Response: Pmech.kJ
# Chisq Df Pr(>Chisq)    
# timeTrack.min 99.264  1  < 2.2e-16 ***
r.squaredGLMM(e.m)
#             R2m      R2c
# [1,] 0.7566291 0.899628
plot.e.m <- data.frame(e.m@frame, fitted.re = fitted(e.m))

timeEE <-ggplot()+
  geom_point(data = nightNets, 
             aes(x = timeTrack.min, y = Pmet.kJ, color = groupID), size=3)+
  scale_color_manual(values = mycols, name = "Group")+
  #geom_point(data=kunz, aes(x = mins.out, y = dee.kJday), color = "red")+
  geom_smooth(data = plot.e.m, 
              aes(x = timeTrack.min, y = fitted.re), method = "lm", se = F, color = "black")+
  theme_bw()+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(x = "time tracked per night (min)", 
       y = "Foraging energy\nexpenditure (kJ)", 
       subtitle =  "Foraging energy expenditure increases\nwith tracking time")



#Figure of tracking times, DEE, EE etc
pdf("./output/timeDistDEE_EE_plots.pdf", width = 8, height = 10)
plot_grid(timeDist, firstTime, lastTime, timeDEE, timeEE, ncol=2)
dev.off()


#Foraging & Patch use needs a line ####

ptchTdis.m <- lmer(nPatches~totalDistance + (1|batID), data=nightNets)
Anova(ptchTdis.m) #The positive relationship is driven by distance flown
summary(ptchTdis.m)
# esponse: nPatches
# Chisq Df Pr(>Chisq)    
# totalDistance 12.103  1  0.0005033 ***
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


# Time & Patches needs a line ####
cp.m <- lmer(nPatches~timeTrack.min+(1|batID), data=nightNets)
summary(cl.m) 
# Anova(cp.m)Response: nPatches
# Chisq Df Pr(>Chisq)   
# timeTrack.min 6.7079  1   0.009599 **
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
summary(cl.m) #no effect
Anova(cl.m)

clusPatch <- ggplot(nightNets, aes(x = nPatches, y = nFeedingClusters, color = groupID))+
  geom_point(size=3)+
  scale_color_manual(values = mycols, name = "Group")+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(x = "Number of foraging patches", 
       y = "Number of feeding clusters")
clusPat.m <- lmer(nFeedingClusters~nPatches+(1|batID), data=nightNets)
summary(cl.m) #no effect
Anova(clusPat.m)

pdf("./output/Fig S3 - PatchClusterTimeDist.pdf", width = 8, height = 8)
plot_grid(patchDist, patchTime, feedDist, clusPatch,labels = c('A', 'B', 'C', 'D'), label_size = 12,  ncol=2)
dev.off()

ptchDEE.m <- lmer(Pmet.kJ~nPatches+(1|batID), data=nightNets)
summary(cl.m)
Anova(ptchDEE.m) #no effect

patchEE <- ggplot(nightNets, aes(y = Pmet.kJ, x = nPatches, color = groupID))+
  geom_point(size=3)+
  scale_color_manual(values = mycols, name = "Group")+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(y = "foraging energy expenditure (kJ)", 
       x = "number of patches", 
       subtitle =  "Number of patches used doesn't\neffect foraging energy expenditure")


l.m <- lmer(Pmet.kJ~nFeedingClusters+(1|batID), data=nightNets)
summary(cl.m) #strong effect
Anova(cl.m)
# Response: Pmet.kJ
# Chisq Df Pr(>Chisq)    
# nFeedingClusters 31.039  1  2.529e-08 ***

plot.l.m <- data.frame(l.m@frame, fitted.re = fitted(l.m))


feedClusEE <- ggplot()+
  geom_point(data = nightNets, 
             aes(y = Pmet.kJ, x = nFeedingClusters, color = groupID), size=3)+
  scale_color_manual(values = mycols, name = "Group")+
  geom_smooth(data = plot.l.m,
              aes(x =nFeedingClusters, y = fitted.re), method = "lm", se = F, color = "black")+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(y = "foraging energy expenditure (kJ)", 
       x = "number of feeding clusters", 
       subtitle =  "No relationship between number of\nfeeding clusters & foraging energy expenditure")
c




#Activity Budget & EE ####
actBudg_w <- actBudg %>% dplyr::select(batID_day, newState, freqAct) %>%
  spread(newState, freqAct) %>% 
  dplyr::select(-`<NA>`)
nightNets <- nightNets %>% left_join(actBudg_w)
save(nightNets, file="./processed data/NightlySummaries 2020-08-10.Rdata")

eevals <- nightNets %>% dplyr::select(groupID, batID_day, Pmet.kJ)
actBudg <- actBudg %>% filter(!is.na(newState)) %>% 
  left_join(eevals)

actBudg <- actBudg %>% left_join(observed, by=c("batID" = "ID"))

actEE <- ggplot(actBudg, aes(x = freqAct, y = Pmet.kJ, color=groupID))+
  geom_point(size=3)+
  scale_color_manual(values = mycols, name = "Group")+
  theme_bw()+
  theme(legend.position = "none",
        legend.background = element_blank())+
  facet_wrap(~newState, 
             labeller = labeller(newState = behaviorStates))+
    labs(x = "proportion in activity", 
       y = "foraging energy expenditure (kJ)", 
       subtitle =  "No clear relationships between activity states\n& energy expenditure")

patchClus.m <- lmer(nFeedingClusters~nPatches+(1|batID), data=nightNets)
Anova(patchClus.m)

patchClus <- ggplot(nightNets, aes(x = nPatches, y = nFeedingClusters, color = groupID))+
  geom_point(size=4)+
  scale_color_manual(values = mycols, name = "Group")+
  theme_bw()+
  theme(legend.position = "none",
        legend.background = element_blank(),
        text = element_text(size=16),
        axis.text.x = element_text(size = 10))+
  labs(x = "number of foraging patches used", 
       y = "number of flower clusters visited", 
       subtitle =  "No relationship between number of\npatches & feeding clusters")
ggsave(patchClus, file = "./output/PatchClusters per night.pdf")


### Do close network ties influence energy expenditure? ####
socI <- lmer(Pmech.kJ~strength+(1|batID), data=nightNets) #There is a slightly positive effect of strength on Pmech.kJ (how much energy is spent outside the roost). This is driven by animals with higher strength being trak for longer.
summary(socI)
Anova(socI)
# Response: Pmech.kJ
# Chisq Df Pr(>Chisq)  
# strength 4.5103  1    0.03369 *


#Mechanical power is predicted by network strength, but metabolic power isn't. That's the one that counts. 

socI <- lmer(Pmet.kJ~strength+(1|batID), data=nightNets) #There is a slightly positive effect of strength on Pmech.kJ (how much energy is spent outside the roost). This is driven by animals with higher strength being trak for longer.
summary(socI)
Anova(socI)
# Response: Pmet.kJ
# Chisq Df Pr(>Chisq)
# strength 1.0137  1      0.314


socDEE <- lmer(dee.kJ~strength+(1|batID), data=nightNets)
Anova(socDEE)


plot.socI <- data.frame(socI@frame, fitted.re = fitted(socI))
head(plot.socI)


strEE <- ggplot()+
  geom_point(data= nightNets, aes(x = strength, y = Pmet.kJ, color = groupID), size=3)+
  scale_color_manual(values = mycols, name = "Group")+
  theme_bw()+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(x = "Network strength", 
       y = "Energy expenditure\nduring foraging (kJ)")

strDEE <- ggplot()+
  geom_point(data= nightNets, aes(x = strength, y = dee.kJ, color = groupID), size=3)+
  scale_color_manual(values = mycols, name = "Group")+
  theme_bw()+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(x = "Network strength", 
       y = "Daily energy\n expenditure (kJ)")




socItime <- lmer(Pmech.kJ~strength*timeTrack.min+(1|batID), data=nightNets)
Anova(socItime)


distI <- lmer(totalDistance~strength+(1|batID), data=nightNets)
Anova(distI)
# Response: totalDistance
# Chisq Df Pr(>Chisq)  
# strength 2.7223  1    0.09895

plot.distI <- data.frame(distI@frame, fitted.re = fitted(distI))
head(plot.distI)

strTime <- ggplot()+
  geom_point(data=nightNets, aes(y = timeTrack.min, x = strength, color = groupID), size=3)+
  scale_color_manual(values = mycols)+
  
  
  
strDist <- ggplot()+
  geom_point(data = nightNets, aes(y = totalDistance, x = strength, color = groupID), size=3)+
  scale_color_manual(values = mycols, 
                     name = "", 
                     breaks = c("blue", "brown", "yellow"),
                     labels = c("Group 1", "Group 2", "Group 3"))+
  theme(legend.position = c(0.6, 0.94), 
        legend.background = element_blank(),
        legend.title = element_blank())+
  labs(y = "Distance moved\nper night (km)", 
       x = "Network strength")


strPatch <- ggplot(nightNets, aes(y = nPatches, x = strength, color = groupID))+
  geom_point(size=3)+
  scale_color_manual(values = mycols)+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(y = "Number of patches", 
       x = "Network strength")
patchI <- lmer(nPatches~strength+(1|batID), data=nightNets)
Anova(patchI)
# Response: nPatches
# Chisq Df Pr(>Chisq)
# strength 0.0058  1     0.9395

strClust <- ggplot(nightNets, aes(y = nFeedingClusters, x = strength, color = groupID))+
  geom_point(size=3)+
  scale_color_manual(values = mycols)+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(y = "Number of feeding clusters", 
       x = "Newtork strength")
clustI <- lmer(nFeedingClusters~strength+(1|batID), data=nightNets)
Anova(clustI)
# esponse: nFeedingClusters
# Chisq Df Pr(>Chisq)
# strength 0.5218  1     0.4701

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
       x = "Newtork strength")


pdf("./output/Fig S4 - NetStrength_forFeed.pdf", width = 9, height=3)
plot_grid(strDist, strPatch, strClust, labels = c('A', 'B', 'C'), label_size = 12, ncol=3)
dev.off()

#Do networks measures influence how quickly bats hit their first patch & flower? ####

fstPtch <- ggplot()+
  geom_point(data=nightNets, aes(y = abs(trackDistPatch), x = strength, color = groupID), size=3)+
  scale_color_manual(values = mycols)+
  #geom_smooth(aes(y = fitted.re, x=strength), 
              method = "lm", data = plot.timeI, se=FALSE, color = "black")+
  theme_bw()+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(y = "time to first patch (min)", 
       x = "strength",
       subtitle = "Network strength predicts\ntime tracked")



#Activity & network measures ####
ggplot(actBudg, aes(x = freqAct , y = strength, color = groupID))+
  geom_point(size=3)+
  scale_color_manual(values = mycols, name="Group")+
  theme_bw()+
  theme(legend.position = "none",
        legend.background = element_blank())+
  facet_wrap(~newState, 
             labeller = labeller(newState = behaviorStates))+
  labs(x = "proportion in activity", 
       y = "network strength", 
       subtitle =  "Activity & network strength")


#Network strength vs individual activity measures.

actNet.r <- lmer(state1~strength + (1|batID), data = nightNets)
Anova(actNet.r)
# Response: state1
# Chisq Df Pr(>Chisq)  
# strength 4.4804  1    0.03429 *
plot.actNet.r <- data.frame(actNet.r@frame, fitted.re = fitted(actNet.r))
actNet.f <- lmer(state2~strength + (1|batID), data = nightNets)
Anova(actNet.f)
actNet.m <- lmer(state3~strength + (1|batID), data = nightNets)
Anova(actNet.m)
# Response: state3
# Chisq Df Pr(>Chisq)  
# strength 4.5993  1    0.03198 *
plot.actNet.m <- data.frame(actNet.m@frame, fitted.re = fitted(actNet.m))
actNet.c <- lmer(state4~strength + (1|batID), data = nightNets)
Anova(actNet.c)

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

pdf("./output/Fig 3 - ActBudgNetStr.pdf",  width = 8, height = 8)
plot_grid(rest, moveF, forage, commute, labels = c('A', 'B', 'C', 'D'), label_size = 12, ncol = 2)
dev.off()

pdf("./output/Fig 3 - ActBudgNetStr_wide.pdf",  width = 12, height = 4)
plot_grid(rest, moveF, forage, commute, labels = c('A', 'B', 'C', 'D'), label_size = 12, rows = 1)
dev.off()




ggplot(actBudg, aes(x = freqAct , y = strength, color = groupID))+
  geom_point(size=3)+
  scale_color_manual(values = mycols, name="Group")+
  theme_bw()+
  theme(legend.position = "none",
        legend.background = element_blank())+
  facet_wrap(~newState, 
             labeller = labeller(newState = behaviorStates))+
  labs(x = "proportion in activity", 
       y = "network strength", 
       subtitle =  "Activity & network strength")



ggplot(actBudg, aes(x = freqAct , y = centrality, color = groupID))+
  geom_point(size=3)+
  scale_color_manual(values = mycols, name="Group")+
  theme_bw()+
  theme(legend.position = "none",
        legend.background = element_blank())+
  facet_wrap(~newState, 
             labeller = labeller(newState = behaviorStates))+
  labs(x = "proportion in activity", 
       y = "network strength", 
       subtitle =  "Activity & network centrality")

ggplot(actBudg, aes(x = freqAct , y = degree, color = groupID))+
  geom_point(size=3)+
  scale_color_manual(values = mycols, name="Group")+
  theme_bw()+
  theme(legend.position = "none",
        legend.background = element_blank())+
  facet_wrap(~newState, 
             labeller = labeller(newState = behaviorStates))+
  labs(x = "proportion in activity", 
       y = "network strength", 
       subtitle =  "Activity & network degree")


#What happens close locations(less than 100m)


closeActBudg <- hast_df %>% filter(otherDist < 51) %>% 
  group_by(groupID, batID, batID_day, date, newState) %>% 
  summarise(nlocsBehav = n(),
            sumTimeBehav = sum(tlag, na.rm = TRUE)) %>% 
  mutate(freqAct = nlocsBehav / sum(nlocsBehav))

closeActBudg <- closeActBudg %>% left_join(observed, by=c("batID" = "ID"))

ggplot(closeActBudg, aes(x = freqAct , y = strength, color = groupID))+
  geom_point(size=3)+
  scale_color_manual(values = mycols, name="Group")+
  theme_bw()+
  theme(legend.position = "none",
        legend.background = element_blank())+
  facet_wrap(~newState, 
             labeller = labeller(newState = behaviorStates))+
  labs(x = "proportion in activity", 
       y = "network strength", 
       subtitle =  "Activity & network strength")

ggplot(closeActBudg, aes(x = freqAct , y = centrality, color = groupID))+
  geom_point(size=3)+
  scale_color_manual(values = mycols, name="Group")+
  theme_bw()+
  theme(legend.position = "none",
        legend.background = element_blank())+
  facet_wrap(~newState, 
             labeller = labeller(newState = behaviorStates))+
  labs(x = "proportion in activity", 
       y = "network strength", 
       subtitle =  "Activity & network centrality")

ggplot(closeActBudg, aes(x = freqAct , y = degree, color = groupID))+
  geom_point(size=3)+
  scale_color_manual(values = mycols, name="Group")+
  theme_bw()+
  theme(legend.position = "none",
        legend.background = element_blank())+
  facet_wrap(~newState, 
             labeller = labeller(newState = behaviorStates))+
  labs(x = "proportion in activity", 
       y = "network strength", 
       subtitle =  "Activity & network degree")


#Number of feeding clusters & ee?
ggplot(nightNets, aes(x = nFeedingClusters, y = Pmech.kJ, color = groupID))+
  geom_point(size=3)+
  scale_color_manual(values = mycols, name="Group")+
  theme_bw()+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(x = "number of feeding clusters", 
       y = "movement energy expenditure (kJ)",
       subtitle = "Network stregnth increases with time tracked")
m_clusEE <- lmer(Pmech.kJ~nFeedingClusters+(1|batID), data=nightNets)
Anova(m_clusEE) #no effect

ggplot(nightNets, aes(y = nFeedingClusters, x = strength, color = groupID))+
  geom_point(size=3)+
  scale_color_manual(values = mycols)+
  theme_bw()+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(y = "number of feeding clusters", x = "network strength")

ggplot(nightNets, aes(y = nFeedingClusters, x = centrality, color = groupID))+
  geom_point(size=3)+
  scale_color_manual(values = mycols)+
  theme_bw()+
  theme(legend.position = "none",
        legend.background = element_blank())+
  labs(y = "number of feeding clusters", x = "centrality")

ggplot(hast_df, aes(x=otherDist))+geom_histogram()+xlim(0, 500)+geom_vline(xintercept = 50)

