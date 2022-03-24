#Patch Residency
#
#I want to know how time in a patch 1) influences the decision to stay / go and 2) is influenced by the length of time the bat is in the patch with another animal. 
library(tidyverse)
library(lubridate)
load("./data/11_HastatusSegStateBiodatPwr.Rdata") #hastMorph
load("./data/9_HastHMMdbClus_flwTrueTlagTrueGSAS.Rdata") #hast_df
load("./data/7_HastatusSegClusRelevel.Rdata") #hastP 

forage <- hastP %>% filter(behaviour=="forage")
forage$timestamp2s <- ceiling_date(forage$timestamps, unit = "2 seconds")
otherBatPatch <- forage %>% dplyr::select(timestamps, timestamp2s, batID, otherBatPatch = patchID, otherBatState = newState)
forage <- forage %>% left_join(otherBatPatch, by=c("timestamp2s" = "timestamp2s", "otherBat" = "batID")) 
batgroups <- hastP %>% group_by(batID) %>% summarize(groupID = unique(groupID))

#This will now show what patch the parter bat is in (if any), what it's behavior state is, and if it's in the same patch as the focal individual at any time. 

forage$samePatch <- NA
forageL <- split(forage, f=forage$batIDday)

forageX <- lapply(forageL, function(x){
  samePatch <- x$samePatch
  for(i in 1:length(x$samePatch)){
  samePatch[i] <- ifelse(x$patchID[i] == x$otherBatPatch[i], 1, 0)
  }
  x$samePatch <- samePatch
  return(x)
})
dat <- do.call("rbind", forageX)
forage <- dat
#NA in samePatch is when at least 1 bat *isn't* in a patch. If I'm interested in social effects of being in the same patch, these should be coded as 0
forage$samePatch <- replace_na(forage$samePatch, 0)

save(forage, file="./data/12_AllBatsSamePatch.Rdata")

#Create a variable to split everything up by if the bat is in / out of a patch continuously. This needs to be done since splitting by patch ID will pull all locations, regardless of continuity. Also calculate behavioral state changes.

forageL <- split(forage, f = forage$batIDday)

#This will take some time. Go get a coffee
forageLx <- lapply(forageL, function(x){
  x$patchStaySame <- NA
  x$stateStaySame <- NA
  #Start with 0 as a baseline. No change is 0, change is 1. This will let us count up with a cumulative sum
  x$patchStaySame[1] <- ifelse(x$patchID[1] == x$patchID[2], 0, 1)
  x$stateStaySame[1] <- ifelse(x$newState[1] == x$newState[2], 0, 1)
  for(i in 2:length(x$patchStaySame)-1){
    j = i+1
    x$patchStaySame[j] = ifelse(x$patchID[j] == x$patchID[i], 0, 1)
    x$stateStaySame[j] = ifelse(x$newState[j] == x$newState[i], 0, 1)
  }
  x$patchBout <- x$patchStaySame
  x$patchBout[!is.na(x$patchStaySame)] <- cumsum(x$patchBout[!is.na(x$patchStaySame)]) 
  x$stateBout <- x$stateStaySame
  x$stateBout[!is.na(x$stateStaySame)] <- cumsum(x$stateBout[!is.na(x$stateStaySame)])
  return(x)
})
patchStay_df <- do.call("rbind", forageLx)

save(patchStay_df, file="./data/12_AllBatsPatchStateSeq.Rdata")

### I want to summarize how much time is spent in each state per bout and how many bats are in the same patch as the focal.
### Would like to add in if there were other bats in the patch at the same time. This probably needs to come from the 7-6 patchStay that had all combinations of bats in a patch. Then reduce it down to add here.

patchBoutSums <- patchStay_df %>% group_by(batIDday, patchID, patchBout, samePatch, newState) %>%
  summarize(timeInState = sum(tlag, na.rm=T), 
            minBats = min(n.bats),
            maxBats = max(n.bats), 
            meanBats = mean(n.bats, na.rm=T)) %>% 
  mutate(timeInBout = sum(timeInState, na.rm=T)) %>% 
  mutate(stateFreq = timeInState / timeInBout)

pairDist <- hast_df %>% group_by(groupID, batID, behaviour, otherBat) %>% 
  summarize(distMean = mean(otherDist, na.rm=T),
            distMin = min(otherDist, na.rm=T), 
            distMax = max(otherDist, na.rm=T),
            distSD = sd(otherDist, na.rm=T)) 

#Bats are far apart ordered by commute_back > commute_there > forage. But there is considerable variation. Are there pairs that maybe are maybe closer slightly more?

#I now have if the patch the bat is in stays the same and the sequence of being in a patch for each individual. `patchBout` is NA when the bat is not in a patch, but numeric otherwise. `patchBout` should be the variable that each bat night is split by to test how bat behavior *in a patch* is influenced by a given variable.

#What are they doing when they are in the same patch? They are mostly resting. Resting is much more frequent when they are in the same patch.
samePatch <- patchStay_df %>%  group_by(batIDday, patchID, samePatch, newState) %>% 
  dplyr::summarise(sumTime = sum(tlag)) %>% 
  mutate(freqTime = sumTime / sum(sumTime))

m1 <- glmer(freqTime~samePatch*newState+(1|batIDday), data=samePatch, family = binomial(link=logit))
Anova(m1)
