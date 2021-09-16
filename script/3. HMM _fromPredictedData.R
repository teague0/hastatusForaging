
#HMM of full tracks using momentuHMM
library(momentuHMM)
library(tidyverse)
library(move)
library(lubridate)
library(parallel)

load("./data/2_Hastatus Crawl Fits.Rdata")

#### WARNING: Go get coffee or a snack because this will take a few minutes. If running on a slower machine, then a leisurely meal is recommended. This took approx 20 mins to run ####

#Please refer to MicClintock et al 2018 (10.1111/2041-210x.12995) and Michelot et al 2017 (10.1002/ecy.1880) for tutorials & description of moveHMM 

testmods <- mclapply(X=hastLcrw, mc.cores = detectCores()-2, FUN = function(batHMM) { 
  tryCatch({
    tmp <- prepData(batHMM)
    state4Names <- c("state1", "state2", "state3", "state4")
    # initial parameters
    step4ParMean <- c(1, 4, 7, 12)
    step4ParSD <- c(3, 3, 1, 1)
    zero4mass0 <- c(0.9, 0.5, 0, 0)
    step4Par0 <- c(step4ParMean,step4ParSD, zero4mass0)
    angle4Par0 <- c(0.5, 0.5, 0.1, 0.001)
    dist = list(step = "gamma", angle = "wrpcauchy")
    m <- fitHMM(data=tmp, 
                nbStates=4, 
                dist=dist,
                Par0=list(step=step4Par0, angle=angle4Par0),
                stateNames = state4Names)
    return(m)
  },error=function(e) NULL)
})

save(testmods, file="./data/HMM_fullNights.Rdata")

## HMM interpretation ####
#The HMMs were fit per individual bat night. Because these were fit on individual nights, the state classifications are all over the place. They need to be re-ordered (or given some sort of behavior) to know what is going on. 

#use object hast_df from script 1 or load in written csv file.
hast_df <- read.csv("./data/Phyllostomus Hastatus Processed GPS Data S1.csv")
hast_df$timestamp <- as.POSIXct(hast_df$timestamp, tz = "UTC")
hast_df$timestamps <- as.POSIXct(hast_df$timestamps, tz = "UTC")
load("./processed data/HMM_fullNights.Rdata")

#Extract table of the step values & it's state. Create a re-ordered newState to reflect low to high step lengths. 
states <- lapply(testmods, function(x){
  tryCatch({
  step_df <- data.frame(t(x$mle$step))
  step_df$stateID <- rownames(step_df)
  step_df$batIDday <- x$data$ID[1]
  o <- order(step_df$mean)
  step_df$newState <- NA
  for(i in 1:length(step_df$stateID)){
    state <- step_df$stateID[i]
    newspot <- o[i]
    step_df$newState[newspot] <- state
  }
  return(step_df)
  },error=function(e) NULL)
})
states_df <- do.call("rbind", states)

states_df %>% filter(batIDday=="X74F7D4C.1") #all of the states are less than 1 except 4 it's a small foraging area, with almost everything as state1. so recode 1 to 3, 2 to 1, 3 & 4 to NA
tooslow <- which(states_df$batIDday == "X74F7D4C.1")
states_df$newState[tooslow] <- c("state3","state1", NA, NA)

states_df %>% filter(batIDday=="X74F8E19.3")
toofast <- which(states_df$batIDday == "X74F8E19.3")
states_df$newState[toofast] <- c("state4","state3")

#Extract out the HMM dataframe and the viterbi states, and join these to the original data. This will eliminate the crawl estimated points. 
locsStates <- lapply(testmods, function(x){
  tryCatch({
  tmp <- x$data
  tmp$stateID <- viterbi(x)
  dat <- tmp %>% dplyr::select(batIDday=ID, time, stateID)
  return(dat)
  },error=function(e) NULL)
})

fullDat <- do.call("rbind", locsStates)
fullDat$batID_day <- as.character(fullDat$batIDday)
hastdfHMM <- left_join(hast_df, fullDat, by = c("batIDday", "timestamps" = "time"))
save(hastdfHMM, file="./processed data/HastatusGPS_HMMstates.Rdata")

#Add in the re-leveled states
hast_L <- split(hastdfHMM, f=hastdfHMM$batIDday)
newStates <- lapply(hast_L, function(x){
  x$states <- paste0("state",x$stateID)
  ID <- unique(x$batIDday)
  batTable <- states_df %>% filter(batIDday == ID) %>% dplyr::select(states=stateID, newState, stateMean=mean)
  dat <-left_join(x, batTable, by = "states")
})
#save(newStates, file="./data/3_HastatusGPS_HMMRelevelstates.Rdata")


