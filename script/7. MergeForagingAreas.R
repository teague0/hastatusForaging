#Create consistent naming for foraging sites. This repeats some of the analysis in 4. ForagingMCPareas.R, but will add in the centroid, calculate teh distances among centroids and figure out a cutoff for merging these areas.

library(tidyverse)
library(lubridate)
library(parallel)
library(adehabitatHR)
library(geosphere)

load("./data/6_HastatusClusDistPairDup.Rdata") #The data with pairwise distances. Events duplicated.

clusFrgFlw <- split(hastPairsDupLoc, f = hastPairsDupLoc$batIDday)

#Calculate 100% MCP for each cluster 

clusMCPcent <-  mclapply(X=clusFrgFlw, mc.cores = detectCores()-2, FUN = function(x){
  tryCatch({
    z <- x %>% filter(dbClus != 0)
    ID <- unique(z$batIDday)
    clusterList <- split(z, f=z$dbClus)
    mcpList <- lapply(clusterList, function(w){
      tryCatch({
        clusName <- unique(w$dbClus)
        batIDday <- unique(w$batIDday)
        coordinates(w) <- ~utm.x+utm.y
        proj4string(w) <- CRS("+proj=utm +zone=17 +datum=WGS84")
        mcp100 <- mcp(w, percent = 100, unin="m", unout = "m2")
        MCPcentroid <- centroid(mcp100)
        MCPcentroid.x <- MCPcentroid[,1]
        MCPcentroid.y <- MCPcentroid[,2]
        ptCenter <- data.frame("ptMean.x"=mean(w@coords[,1]),"ptMean.y"=mean(w@coords[,2]))
        datOut <- cbind(data.frame(batIDday, clusName, MCPcentroid.x, MCPcentroid.y), ptCenter)
        return(datOut)
      },error=function(e) NULL)
    })
    dayClus <- do.call("rbind", mcpList)
    return(dayClus)
  }, error=function(e) NULL)
})

clusCent <- do.call("rbind", clusMCPcent)
clusCent$batIDday_clus <- paste0(clusCent$batIDday, ".", clusCent$clusName)
clusCent <- clusCent %>%  
  dplyr::mutate(centroidPtsDist=sqrt((MCPcentroid.x - ptMean.x)^2 + (MCPcentroid.y - ptMean.y)^2))

clusCent$k <- 1
centroidDist <- clusCent %>% 
  dplyr::select(batIDday_clus, MCPcentroid.x, MCPcentroid.y, k) %>% 
  dplyr::full_join(clusCent, by = "k") %>% 
  dplyr::filter(batIDday_clus.x != batIDday_clus.y) %>%
  dplyr::mutate(dist = sqrt((MCPcentroid.x.x - MCPcentroid.x.y)^2 + (MCPcentroid.y.x - MCPcentroid.y.y)^2)) %>%
  dplyr::select(-k)
centroidDist$type <- "centroid"

ptCentrDist <- clusCent %>% 
  dplyr::select(batIDday_clus, ptMean.x, ptMean.y, k) %>% 
  dplyr::full_join(clusCent, by = "k") %>% 
  dplyr::filter(batIDday_clus.x != batIDday_clus.y) %>%
  dplyr::mutate(dist = sqrt((ptMean.x.x - ptMean.x.y)^2 + (ptMean.y.x - ptMean.y.y)^2)) %>%
  dplyr::select(-k)
ptCentrDist$type <- "point"
# ggplot(ptCentrDist, aes(x=dist))+
#   geom_histogram()+
#   xlim(0,100)+
#   geom_vline(xintercept = 30)

#30 m distance looks like a good cutoff

#Recode to a systematic foraging area identity
ptCentrDist <-  ptCentrDist %>% mutate(same = ifelse(dist < 30, 0, 1))

#Order the patches by timestamps
patchInfo <- hastPairsDupLoc %>% 
  dplyr::mutate(batIDday_clus=paste0(batIDday, ".", dbClus)) %>% 
  dplyr::filter(dbClus > 0) %>%
  dplyr::group_by(batIDday_clus, dbClus) %>%
  dplyr::summarize(minTime = min(timestamps))
o <- order(patchInfo$minTime)
patchInfo$patchOrder <- o
patchInfo <- patchInfo %>% dplyr::select(batIDday_clus, patchOrder)

patches <- dplyr::left_join(ptCentrDist, patchInfo, by=c("batIDday_clus.x"="batIDday_clus"))
patches <- dplyr::left_join(patches, patchInfo, by=c("batIDday_clus.y"="batIDday_clus"))
patches <- patches %>% dplyr::arrange(patchOrder.x)
patches$patchID.x <- NA
patches$patchID.y <- NA
for(i in 1:length(patches$patchOrder.x)){
  patches$patchID.x[i] <- ifelse(patches$same[i] == 1, patches$patchOrder.x[i],
                              ifelse(patches$patchOrder.x[i] < patches$patchOrder.y[i], patches$patchOrder.x[i], patches$patchOrder.y[i]))
  patches$patchID.y[i] <- ifelse(patches$same[i] == 1, patches$patchOrder.y[i],
                              ifelse(patches$patchOrder.x[i] > patches$patchOrder.y[i], patches$patchOrder.y[i], patches$patchOrder.x[i]))
}

xvals <- patches %>%  dplyr::select(batIDday_clus = batIDday_clus.x, patchID = patchID.x, 
                                    center.x = ptMean.x.x, center.y = ptMean.y.x)
yvals <- patches %>%  dplyr::select(batIDday_clus = batIDday_clus.y, patchID = patchID.y, 
                                    center.x = ptMean.x.y, center.y = ptMean.y.y)

patchNames <- dplyr::bind_rows(xvals, yvals) %>% 
  dplyr::group_by(batIDday_clus) %>% 
  dplyr::summarize(patchID = min(patchID), 
                   patchCtr.x = mean(center.x),
                   patchCtr.y = mean(center.y)) %>% 
    arrange(patchID)
hastPairsDupLoc$batIDday_clus <- paste0(hastPairsDupLoc$batIDday, ".", hastPairsDupLoc$dbClus)
hastPairsDupLoc <- dplyr::left_join(hastPairsDupLoc, patchNames)

forage <- hastPairsDupLoc %>% filter(behaviour == "forage")
#Recalculate patch areas from the re-classified patchIDs ####
patchGPSptsL <- split(forage, f=forage$patchID) 
patchAreaL <- lapply(patchGPSptsL, function(x){
  coordinates(x) <- ~utm.x+utm.y
  proj4string(x) <- CRS("+proj=utm +zone=17 +datum=WGS84")
  mcp100 <- mcp(x, percent = 110, unin="m", unout = "m2")
  patchID <-  unique(x$patchID)
  patchArea <- mcp100@data$area
  cent <- centroid(mcp100) 
  patchCent.x <- cent[,1]
  patchCent.y <- cent[,2]
  dat <- data.frame(patchID, patchArea, patchCent.x, patchCent.y)
  return(dat)
})
patchAreas <- do.call("rbind", patchAreaL)

#Load coordinates of cave to calculate distance
coords <- cbind(-82.271541, 9.396448)
lagruta <- SpatialPoints(coords, proj4string = CRS("+proj=longlat +datum=WGS84"))
lagrutaUTM <- spTransform(lagruta, CRS("+proj=utm +zone=17 +datum=WGS84"))

patchAreas <- patchAreas %>% 
  mutate(patchDistToCave = spDists(SpatialPoints(cbind(patchCent.x, patchCent.y), proj4string = CRS("+proj=utm +zone=17 +datum=WGS84")), lagrutaUTM))

hastP <- dplyr::left_join(hastPairsDupLoc, patchAreas, by = c("patchID" = "patchID"))
save(hastP, file="./data/7_HastatusSegClusRelevel.Rdata")
