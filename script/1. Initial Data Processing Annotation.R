#This starts at the beginning and pulls in original Movebank data, drops unused locations (tags that weren't included, outlier points) & cleans up some names, etc.

#Note that depending where / how the original GPS data were downloaded from Movebank, the names of the columns will either have underscores "_" or periods/full stops ".". This is a legacy issue of the API.

library(move)
library(tidyverse)
library(lubridate)

ph <- move("./data/Greater spear-nosed bats Phyllostomus hastatus Bocas del Toro-6665207312962129521.zip") #This is a move object that includes environmental annotation from the Env-DATA service

droppedLocs <- read.csv("./data/Dropped Locations.csv") #locations that were dropped from analysis because they were the wrong tag type or were tag errors.
ph_df <- as(ph, "data.frame") #set to data frame to make this easier
`%notin%` <- function(x,y) !(x %in% y) #Define a not in function
hast_df <- ph_df[ph_df$event.id %notin% droppedLocs$x,] #There should be 237640 observations left.

behavs <- read.csv("./data/SegmentedBehavior_eventIDs.csv") #manual segmentation of major foraging vs commuting / direction
hast_df <- hast_df %>% left_join(behavs)
hast_df$height.calc <- hast_df$height.above.msl - hast_df$ASTER.ASTGTM2.Elevation

#Add local wind data from STRI (http://biogeodb.stri.si.edu/physical_monitoring/research/bocas)

wind <- read.csv(file = "./data/bocas wind data.csv")
wind$timestamp <- as.POSIXct(wind$timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "America/Panama")

wind$speed.ms <- wind$wind.speed.kmh * 1000/3600 #kmh to ms
wind.phas <- wind %>% filter(date(timestamp) > "2016-02-27" & date(timestamp) < "2016-03-12") #filter to daterange
#Calculate U & V components of wind
library(circular)
wind.phas$U.Component <-  -wind.phas$speed.ms*sin(rad(wind.phas$wind.dir)) 
wind.phas$V.Component <-  -wind.phas$speed.ms*cos(rad(wind.phas$wind.dir))
names(wind.phas) <- c("timestamp", "windDir", "windSpeed.kmh", 'windSpeed.ms', "U.Component", "V.Component")

#Bocas data were binned into 15 minute intervals but the GPS data are in 1 s. Interpolate to get weighted weather values.
t <- unique(hast_df$timestamps)
t.fl <- data.frame(floor_date(t, unit = "15 mins"))
names(t.fl) <- c("timestamp")
pre <- left_join(t.fl, wind.phas)
t.ce <- data.frame(ceiling_date(t, unit = "15 mins"))
names(t.ce) <- c("timestamp")
post <- left_join(t.ce, wind.phas)

Udif <- post$U.Component - pre$U.Component
Vdif <- post$V.Component - pre$V.Component
dtif <- as.numeric(t - pre$timestamp)
wt <- dtif / 3600
U.Component <- pre$U.Component + Udif * wt
V.Component <- pre$V.Component + Vdif * wt
timestamp <- t
winddat <- data.frame(timestamp, U.Component, V.Component)

#add in the wind components to the gps data
hast_df <- hast_df %>% left_join(winddat)

#Calculate airspeeds. The easiest way for many of these metrics is to convert back to a move object to take advantage of the functions there.

phas.ann <- move(data=hast_df,
                 x=hast_df$location.long, 
                 y=hast_df$location.lat, 
                 time=as.POSIXct(hast_df$timestamp, format="%Y-%m-%d %H:%M:%S", tz="UTC"), 
                 animal = hast_df$individual.local.identifier,
                 sensor = hast_df$sensor.type,
                 proj = CRS("+proj=longlat +ellps=WGS84"))

#Calculate variables needed for airspeed calcualtion
phas.ann$heading <- unlist(lapply(lapply(split(phas.ann), angle), "c", NA))
phas.ann$tlag <- unlist(lapply(timeLag(phas.ann, units="secs"), c, NA))
phas.ann$tagGroundSpeed.ms <- phas.ann$ground_speed
phas.ann$ground.speed <- unlist(lapply(speed(phas.ann), c, NA))
phas.ann$stepLength <- unlist(lapply(distance(phas.ann), c, NA))
phas.ann$angle <- unlist(lapply(angle(phas.ann), c, NA))
phas.ann$batID <- phas.ann@trackId

#Add in Bat Day
moveList <- lapply(split(phas.ann), function(myInd){
  datechange <- c(0, abs(diff(as.numeric(as.factor(date(myInd@timestamps-(12*60*60)))))))
  myInd$BatDay <- cumsum(datechange)+1
  return(myInd)
})
phas.ann <- moveStack(moveList, forceTz="UTC")
behaves <- phas.ann$behaviour
other <- which(behaves == "other") #Other looks like foraging to me.
phas.ann$behaviour[other] <- "forage"

#Add in the UTM XY locations, just in case they are needed
phas.ann <- spTransform(phas.ann, CRS("+proj=utm +zone=17 +datum=WGS84"))
crds <- as.data.frame(phas.ann@coords)
phas.ann$utm.x <- crds$coords.x1 #27
phas.ann$utm.y <- crds$coords.x2 #28
phas.ann <- spTransform(phas.ann, CRS("+proj=longlat +datum=WGS84")) #convert back to longlat
phas.ann$batIDday <- paste0(phas.ann$batID, ".", phas.ann$BatDay)

#Calculate airspeed
wind_support <- function(u,v,heading) {
  angle <- atan2(u,v) - heading/180*pi
  return(cos(angle) * sqrt(u*u+v*v))
}

cross_wind <- function(u,v,heading) {
  angle <- atan2(u,v) - heading/180*pi
  return(sin(angle) * sqrt(u*u+v*v))
}

airspeed <- function(x)
{
  va <- sqrt((x$ground.speed - x$ws)^2 + (x$cw)^2)
  return(va)
}

hastL <- split(phas.ann)
hasLX <- lapply(hastL, function(x){
  df <- data.frame(x)
  x$ws <- wind_support(df$U.Component, df$V.Component, df$heading)
  x$cw <- cross_wind(df$U.Component, df$V.Component, df$heading)
  x$airspeed <- airspeed(x)
  return(x)
}) 
hast_mv <- moveStack(hasLX,  forceTz="UTC") #back to a move object
hast_df <- as.data.frame(hast_mv) #and a dataframe
write.csv(hast_df, file = "./data/Phyllostomus Hastatus Processed GPS Data S1.csv", row.names = FALSE)
