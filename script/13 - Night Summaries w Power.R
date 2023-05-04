#Create nightly summary tables.

library(tidyverse)
library(lubridate)
library(lme4)
library(car)

load("./data/11_HastatusSegStateBiodatPwr.Rdata") #hastMorph

##Night sums ####
nightSums <- hastMorph %>% group_by(groupID, batID, BatDay, batIDday) %>% 
  dplyr::summarise(minTimeUTC = min(timestamps),
            maxTimeUTC = max(timestamps),
            minTimeLocal = min(with_tz(timestamps, tz="America/Panama")),
            maxTimeLocal = max(with_tz(timestamps, tz="America/Panama")),
            timeTrack.min = round(as.numeric(as.duration(maxTimeLocal-minTimeLocal))/60, 2),
            timeLeftDay.min = 1440 - timeTrack.min,
            nlocs = n(),
            nlocs290 = length(which(otherDist < 291)),
            prop290 = nlocs290 / nlocs,
            meanTimeLag.s = mean(tlag, na.rm = TRUE),
            maxTimeLag.s = max(tlag, na.rm = TRUE),
            totalDistance = sum(stepLength, na.rm = TRUE)/1000,
            nPatches = length(unique(patchID)), 
            nClus1 = length(unique(dbClus, na.rm=TRUE)),
            nClus2 = length(unique(dbClus_feed2, na.rm=TRUE)),
            nClus3 = length(unique(dbClus_feed3, na.rm=TRUE)),
            nClus4 = length(unique(dbClus_feed4, na.rm=TRUE)),
            Pmech.kJ = sum(pwrPenny, na.rm = TRUE)/1000, 
            Pmet.kJ = sum(PmetTotal, na.rm = TRUE)/1000,
            tlagSum.min = sum(tlag, na.rm = TRUE)/60,
            bodymass = unique(mass.deploy)) %>% 
  mutate(Pmet.unTracked.kJ = (timeLeftDay.min * 23.8/60 * bodymass)/1000, 
         dee.kJ = Pmet.unTracked.kJ + Pmet.kJ, 
         nFeedingClusters = nClus2 + nClus3 + nClus4)

#Proportion of time within 290 m)
mean(nightSums$prop290, na.rm=T) #0.1110512
sd(nightSums$prop290, na.rm=T) #0.1411449
max(nightSums$prop290, na.rm = T) #0.5006399

time290 <- hastMorph %>% group_by(groupID, batID, BatDay) %>% 
  filter(otherDist < 291) %>%
  summarize(tsum.m = sum(tlag, na.rm = T)/60)
nightSums <- nightSums %>% left_join(time290)

ggplot(nightSums)+
  geom_point(aes(x = nFeedingClusters, y = tsum.m, color = groupID))
ggplot(nightSums)+
  geom_point(aes(x = dee.kJ, y = tsum.m, color = groupID))
ggplot(nightSums)+
  geom_point(aes(x = dee.kJ, y = nFeedingClusters, color = groupID))

m.290distEE <- lmer(nFeedingClusters ~ tsum.m + (1|batID:groupID), data=nightSums)
Anova(m.290distEE)

ggplot(nightSums)+
  geom_point(aes(x = timeTrack.min, y = tsum.m, color = groupID))

ggplot(nightSums)+
  geom_point(aes(x = tsum.m, y = dee.kJ, color = groupID))


#Individual summary table
individualSums <- nightSums %>% group_by(batID) %>% 
  summarize(nightsTracked = n(),
            meanTimeTracked = mean(timeTrack.min, na.rm=T),
            sdTimeTracked = sd(timeTrack.min, na.rm=T),
            meanLocs = mean(nlocs, na.rm=T),
            sdLocs = sd(nlocs, na.rm=T),
            meanDistance = round(mean(totalDistance, na.rm=T),2),
            sdDistance = round(sd(totalDistance, na.rm=T),2),
            meanPatches = mean(nPatches, na.rm=T),
            meanFeedingClusters = mean(nFeedingClusters, na.rm=T),
            sdFeedingCluster = sd(nFeedingClusters, na.rm=T),
            sdPatches = sd(nPatches, na.rm=T),
            meanDEE = mean(dee.kJ, na.rm=T),
            sdDEE = sd(dee.kJ, na.rm=T))
write.csv(individualSums, file = "./output/IndividualSums.csv", row.names = FALSE)

 #Add in activity sums for each night.
actBudg <- hastMorph %>% group_by(groupID, batID, batIDday, newState) %>% 
  dplyr::summarize(nlocs = n()) %>% 
  mutate(stateFreq = nlocs / sum(nlocs)) %>% 
  dplyr::select(-nlocs) %>% 
  pivot_wider(names_from = newState, values_from = stateFreq) %>% 
  dplyr::select(-`NA`)

actBudgTime <- hastMorph %>% 
  group_by(groupID, batID, batIDday, newState) %>% 
  dplyr::summarize(nStateObs = n(),
                   nStateTime = sum(tlag, na.rm = TRUE)) %>%
  ungroup() %>% 
  dplyr::group_by(groupID, batID, batIDday) %>% 
  mutate(dayTime = sum(nStateTime),
         stateFreqTime = nStateTime / dayTime) %>% 
  dplyr::select(-dayTime, -nStateObs, -nStateTime) %>% 
  pivot_wider(names_from = newState, values_from = stateFreqTime) %>% 
  dplyr::select(-`NA`)

nightSums <- nightSums %>% left_join(actBudgTime)

  
  #can I also add in time to first patch & first flower?

firstPtchPt <- hastMorph %>% 
  filter(!is.na(patchID)) %>% 
  group_by(groupID, batID, batIDday) %>%
  slice(which.min(timestamp)) %>% 
  dplyr::select(batIDday, firstPatch = patchID, firstPatchTimeUTC = timestamp, firstPatchDist = ptDistToCave)

firstFlwPt <- hastMorph %>% 
  group_by(groupID, batID, batIDday, patchID) %>%
  gather("clustNum", "value", dbClus_feed2:dbClus_feed4 ) %>% 
  filter(value == 1) %>% 
  slice(which.min(timestamp)) %>% 
  dplyr::select(batIDday, firstFlwPatch = patchID, firstFlwrID = value, firstFlwrTimeUTC = timestamp, firstFlwrDist = ptDistToCave)

firstPatchFlw <- left_join(firstPtchPt, firstFlwPt)
nightSums <- left_join(nightSums, firstPatchFlw)

#Add in the social network metrics to the daily summaries
load("./data/10_ObservedSocNetworkMetrics.Rdata") #observed
nightNets <- nightSums %>% left_join(observed, by = c("batID" = "ID"))

startDist <- hastMorph %>% group_by(batIDday) %>%
  filter(timestamps == min(with_tz(timestamps, tz="America/Panama"))) %>%
  dplyr::select(batIDday, firstFixDist = ptDistToCave)
lastDist <- hastMorph %>% group_by(batIDday) %>%
  filter(timestamps == max(with_tz(timestamps, tz="America/Panama"))) %>%
  dplyr::select(batIDday, lastFixDist = ptDistToCave)
firstLast <- startDist %>% left_join(lastDist) 

nightNets <- nightNets %>% left_join(firstLast)
nightNets <- nightNets %>% mutate(firstFixDist = firstFixDist/1000,
                                  lastFixDist = lastFixDist/1000,
                                  elapsedTimePatch = firstPatchTimeUTC - minTimeUTC,
                                  elapsedTimeFlw = firstFlwrTimeUTC - minTimeUTC,
                                  trackDistPatch = firstFixDist - firstPatchDist/1000,
                                  trackDistFlw = firstFixDist - firstFlwrDist/1000)
save(nightNets, file = "./data/13_NightSumValues.Rdata")

#These data are used for Figure 3, Figure 4, Figure S2
