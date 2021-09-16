#Create nightly summary tables.

library(tidyverse)
library(lubridate)
library(lme4)
library(car)

load("./data/11_HastatusSegStateBiodatPwr.Rdata") #hastMorph

##Night sums ####
nightSums <- hastMorph %>% group_by(groupID, batID, batIDday) %>% 
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
            Pmech.kJ = sum(pwrPenny, na.rm = TRUE)/1000, 
            Pmet.kJ = sum(PmetTotal, na.rm = TRUE)/1000,
            tlagSum.min = sum(tlag, na.rm = TRUE)/60,
            bodymass = unique(mass.deploy)) %>% 
  mutate(Pmet.unTracked.kJ = (timeLeftDay.min * 23.8/60 * bodymass)/1000, 
         dee.kJ = Pmet.unTracked.kJ + Pmet.kJ, 
         nFeedingClusters = nClus2 + nClus3 + nClus4)
#Add in activity sums for each night.
actBudg <- hastMorph %>% group_by(groupID, batID, batIDday, newState) %>% 
  dplyr::summarize(nlocs = n()) %>% 
  mutate(stateFreq = nlocs / sum(nlocs)) %>% 
  dplyr::select(-nlocs) %>% 
  pivot_wider(names_from = newState, values_from = stateFreq) %>% 
  dplyr::select(-`NA`)

nightSums <- nightSums %>% left_join(actBudg)

  
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


###Explore effects of time tracking ####
tDist.m <- lmer(totalDistance~timeTrack.min+(1|batID), data=nightNets)
fd.m <-  lmer(timeTrack.min~firstFixDist+(1|batID), data=nightNets)
ld.m <-  lmer(timeTrack.min~lastFixDist+(1|batID), data=nightNets)
t.m <- lmer(dee.kJ~timeTrack.min+(1|batID), data=nightNets)
e.m <- lmer(Pmet.kJ~timeTrack.min+(1|batID), data=nightNets)