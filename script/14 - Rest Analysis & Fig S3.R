#Rest proximity & clusters to know:
# - narrow down distances between individuals when resting
# - are they resting together in the same place repeatedly, etc.
# - map of foraging & resting locations

library(tidyverse)
library(lubridate)

load("./data/12_AllBatsPatchStateSeq.Rdata") #patchStay_df
load("./data/12_AllBatsSamePatch.Rdata") #forage
load("./data/9_HastHMMdbClus_flwTrueTlagTrueGSAS.Rdata") #hast_df 

#Cluster resting sites ####
# My first step is to cluster the resting locations to find if there are sites that are more / less likely to be resting sites, and if these are then used by more than 1 bat repeatedly. Basically, are there just some places that are good for taking a break?

#Filter out the resting locations and split into a list.
resting <- hast_df %>% dplyr::filter(newState == "state1")
restL <- split(resting, f = resting$batID_day)


#Cluster the resting locations. This will be done using parallel clustering to speed it up. Regular lapply would work too. All of this code is based on the previous clustering steps for the foraging and feeding data.

library(fpc)
library(factoextra)
library(parallel)
library(cowplot)
theme_set(theme_cowplot())

clus <- mclapply(X=restL, mc.cores = detectCores()-2, FUN = function(x){
  tryCatch({
    xy <-  x %>% dplyr::select(utm.x, utm.y)
    ID <- unique(x$batIDday)
    epsVal <- 5
    ptsVal <- 15
    set.seed(123)
    db <- fpc::dbscan(xy, eps = epsVal, MinPts=ptsVal)
    x$restClus <- db$cluster
    return(x)
  },error=function(e) finally = print("No Rest"))
})
hastRest <- do.call("rbind", clus)
hastRest$batIDday_RestClus <- paste(hastRest$batIDday, hastRest$restClus, sep=".")


#Merge clusters together ####
#Now the goal is to group the clusters together. To do this I will find the centroid & weights for each resting cluster, then the find distance between the clusters, then group clusters together that are very close
library(tidyverse)
library(adehabitatHR)
library(geosphere)
library(sf)



clusMCPcent <-  mclapply(X=clus, mc.cores = detectCores()-2, FUN = function(x){
  tryCatch({
    z <- x %>% filter(restClus != 0)
    ID <- unique(z$batIDday)
    clusterList <- split(z, f=z$restClus)
    mcpList <- lapply(clusterList, function(w){
      tryCatch({
        restClusName <- unique(w$restClus)
        batIDday <- unique(w$batIDday)
        nlocs <- length(w$utm.x)
        coordinates(w) <- ~utm.x+utm.y
        proj4string(w) <- CRS("+proj=utm +zone=17 +datum=WGS84")
        mcp100 <- mcp(w, percent = 100, unin="m", unout = "m2")
        MCPcentroid <- centroid(mcp100)
        MCPcentroid.x <- MCPcentroid[,1]
        MCPcentroid.y <- MCPcentroid[,2]
        ptCenter <- data.frame("ptMean.x"=mean(w@coords[,1]),"ptMean.y"=mean(w@coords[,2]))
        datOut <- cbind(data.frame(batIDday, restClusName, MCPcentroid.x, MCPcentroid.y), ptCenter)
        return(datOut)
      },error=function(e) NULL)
    })
    dayRestClus <- do.call("rbind", mcpList)
    return(dayRestClus)
  }, error=function(e) NULL)
})

restLoc <- do.call("rbind", clusMCPcent)
restLoc$batIDday_RestClus <- paste(restLoc$batIDday, restLoc$restClusName, sep=".")
restLoc <- restLoc %>% separate(batIDday, c("batID", "batDay"))

#Now I need to group the clusters back together if they are less than a given distance apart. First find out what the distances are. From 7. MergeForagingAreas to ID & merge resting locations.

#Centroids. I don't use this, but it's here.
restLoc$k <- 1
centroidDist <- restLoc %>% 
  dplyr::select(batIDday_RestClus, MCPcentroid.x, MCPcentroid.y, k) %>% 
  dplyr::full_join(restLoc, by = "k") %>% 
  dplyr::filter(batIDday_RestClus.x != batIDday_RestClus.y) %>%
  dplyr::mutate(dist = sqrt((MCPcentroid.x.x - MCPcentroid.x.y)^2 + (MCPcentroid.y.x - MCPcentroid.y.y)^2)) %>%
  dplyr::select(-k)
centroidDist$type <- "centroid"

#Center of locations
ptCentrDist <- restLoc %>% 
  dplyr::select(batIDday_RestClus, ptMean.x, ptMean.y, k) %>% 
  dplyr::full_join(restLoc, by = "k") %>% 
  dplyr::filter(batIDday_RestClus.x != batIDday_RestClus.y) %>%
  dplyr::mutate(dist = sqrt((ptMean.x.x - ptMean.x.y)^2 + (ptMean.y.x - ptMean.y.y)^2)) %>%
  dplyr::select(-k)
ptCentrDist$type <- "point"

ggplot(ptCentrDist, aes(x=dist))+
  geom_histogram()+
  xlim(0,100)+
  geom_vline(xintercept = 25)

#25 m distance looks like a good cutoff, and it's just under 2x the mean stationary error (14.22 ± 12.05 m from Tadarida paper) in these tags

#Recode to a systematic foraging area identity
ptCentrDist <-  ptCentrDist %>% mutate(same = ifelse(dist < 25, 0, 1))

#Order the patches by timestamps
patchInfo <- hastRest %>% dplyr::mutate(batIDday_RestClus=paste(batIDday,restClus, sep=".")) %>% 
  dplyr::filter(restClus > 0) %>%
  dplyr::group_by(batIDday_RestClus, restClus) %>%
  dplyr::summarize(minTime = min(timestamps))
o <- order(patchInfo$minTime)
patchInfo$patchOrder <- o
patchInfo <- patchInfo %>% dplyr::select(batIDday_RestClus, patchOrder)

patches <- dplyr::left_join(ptCentrDist, patchInfo, by=c("batIDday_RestClus.x"="batIDday_RestClus"))
patches <- dplyr::left_join(patches, patchInfo, by=c("batIDday_RestClus.y"="batIDday_RestClus"))
patches <- patches %>% dplyr::arrange(patchOrder.x)
patches$patchID.x <- NA
patches$patchID.y <- NA
for(i in 1:length(patches$patchOrder.x)){
  patches$patchID.x[i] <- ifelse(patches$same[i] == 1, 
                                 patches$patchOrder.x[i],
                                 ifelse(patches$patchOrder.x[i] < patches$patchOrder.y[i], 
                                        patches$patchOrder.x[i],
                                        patches$patchOrder.y[i]))
  patches$patchID.y[i] <- ifelse(patches$same[i] == 1, 
                                 patches$patchOrder.y[i],
                                 ifelse(patches$patchOrder.x[i] > patches$patchOrder.y[i], 
                                        patches$patchOrder.y[i], 
                                        patches$patchOrder.x[i]))
}

xvals <- patches %>%  dplyr::select(batIDday_RestClus = batIDday_RestClus.x, patchID = patchID.x, 
                                    center.x = ptMean.x.x, center.y = ptMean.y.x)
yvals <- patches %>%  dplyr::select(batIDday_RestClus = batIDday_RestClus.y, patchID = patchID.y, 
                                    center.x = ptMean.x.y, center.y = ptMean.y.y)

patchNames <- dplyr::bind_rows(xvals, yvals) %>% 
  dplyr::group_by(batIDday_RestClus) %>% 
  dplyr::summarize(RestpatchID = min(patchID), 
                   RestPatchCtr.x = mean(center.x),
                   RestPatchCtr.y = mean(center.y)) %>% 
  arrange(RestpatchID)
hastRest <- dplyr::left_join(hastRest, patchNames, by = c("batIDday_RestClus" = "batIDday_RestClus"))

save(hastRest, file = "./data/15_Resting Site Locations.Rdata.")

#Plot the merged resting clusters ####
isla_sf <- st_read("./data/shapefiles/almiranteIsla.shp")
isla_utm <- st_transform(isla_sf, crs=32617)

#Load coordinates of cave
coords <- cbind(-82.271541, 9.396448)
lagruta <- SpatialPoints(coords, proj4string = CRS("+proj=longlat +datum=WGS84"))
lagrutaUTM <- spTransform(lagruta, CRS("+proj=utm +zone=17 +datum=WGS84"))
lagruta_sf <- st_transform(st_as_sf(lagrutaUTM))
panama_sf <- st_read("./data/Panama_Province_Boundaries-shp/Province_Boundaries.shp")
bocas_sf <- panama_sf %>% dplyr::filter(NOMBRE == "Bocas del Toro")
bats_sf <- hast_df %>%  
  st_as_sf(coords = c("utm.x", "utm.y"), crs=32617)
bats_bbox <- st_bbox(bats_sf)
bats_bboxB <- bats_bbox
bats_bboxB[1] <- bats_bbox[1] - 5000
bats_bboxB[3] <- bats_bbox[3] - 9000
bats_bboxB[2] <- bats_bbox[2] - 1000
bats_bboxB[4] <- bats_bbox[4] + 2000
bocas_crop <- st_crop(bocas_sf, bats_bboxB)

allRestsites <- hastRest %>% group_by(RestpatchID) %>% 
  dplyr::summarize(batsUsed = length(unique(batID)),
            `Number of\nbat days` = length(unique(date(timestamp))),
            x = min(RestPatchCtr.x),
            y = min(RestPatchCtr.y),
            nPoints = length(timestamp)) %>%
  filter(!is.na(RestpatchID)) %>% 
  dplyr::arrange(rev(batsUsed)) %>% 
  ggplot()+
  geom_point(aes(x = x, y = y, size = `Number of\nbat days`, fill = batsUsed), 
             alpha = 0.8, shape = 21, color = "black")+
  scale_fill_distiller(palette = "RdYlBu",
                     name = "Number bats\nusing",
                     breaks = c(1, 3, 5, 7, 9, 11))+
  scale_shape_manual(name = "Number of\nbat days")+
  geom_sf(data = bocas_crop, color = "black", fill = NA, inherit.aes = FALSE)+
  geom_sf(data = lagruta_sf, color = "black", fill = "grey", shape = 23, size = 2, inherit.aes = FALSE)+
  geom_sf(data = isla_utm, color = "black", fill = NA, inherit.aes = FALSE)+
  labs(x = "Longitude", y = "Latitude")+
  theme(legend.position = "top")

#Plot only the Changuinola island####
changuinolaRestSites <-  hastRest %>% group_by(RestpatchID) %>% 
  dplyr::summarize(batsUsed = length(unique(batID)),
                   `Number of\nbat days` = length(unique(date(timestamp))),
                   x = min(RestPatchCtr.x),
                   y = min(RestPatchCtr.y),
                   nPoints = length(timestamp)) %>%
  filter(!is.na(RestpatchID)) %>% 
  dplyr::arrange(rev(batsUsed)) %>% 
  ggplot()+
  geom_point(aes(x = x, y = y, size = `Number of\nbat days`, fill = batsUsed), 
             alpha = 0.8, shape = 21, color = "black")+
  scale_fill_distiller(palette = "RdYlBu",
                       name = "Number bats\nusing",
                       breaks = c(1, 3, 5, 7, 9, 11))+
  scale_shape_manual(name = "Number of\nbat days")+
  lims(x = c(335800.8, 338184), y = c(1040962, 1043526))+
  geom_sf(data = isla_utm, color = "black", fill = NA, inherit.aes = FALSE)+
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6))+
  labs(x = "Longitude", y = "Latitude")

#Plot only the Colón sites####
colonSites <- hastRest %>% group_by(RestpatchID) %>% 
  dplyr::summarize(batsUsed = length(unique(batID)),
                   `Number of\nbat days` = length(unique(date(timestamp))),
                   x = min(RestPatchCtr.x),
                   y = min(RestPatchCtr.y),
                   nPoints = length(timestamp)) %>%
  filter(!is.na(RestpatchID)) %>% 
  dplyr::arrange(rev(batsUsed)) %>% 
  ggplot()+
  geom_point(aes(x = x, y = y, size = `Number of\nbat days`, fill = batsUsed), 
             alpha = 0.8, shape = 21, color = "black")+
  scale_fill_distiller(palette = "RdYlBu",
                       name = "Number bats\nusing",
                       breaks = c(1, 3, 5, 7, 9, 11))+
  scale_shape_manual(name = "Number of\nbat days")+
  lims(x = c(354000, 364000), y = c(1032000, 1043600))+
  geom_sf(data = bocas_sf, color = "black", fill = NA, inherit.aes = FALSE)+
  geom_sf(data = lagruta_sf, color = "black", fill = "grey", shape = 23, size = 2, inherit.aes = FALSE)+
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6))+
  labs(x = "Longitude", y = "Latitude")


# Plot to file composite resting zoom levels ####
#Figure S3. Resting patches clustered from GPS locations. A) Resting patches are shown across the study area, scaled by the number of individual bat days that the sites were used and colored by the number of individual bats that used the patch. Resting patches are shown for B) Isla Changuinola, and C) Isla Colón.  

bottom_row <- plot_grid(changuinolaRestSites, colonSites, labels = c('B', 'C'), label_size = 18, ncol = 2) #size 18 for pdf, 18 for png

pdf("./output/Fig S3 - CompositeRestLocations.pdf", width = 10, height = 8)
plot_grid(allRestsites, bottom_row, labels = c('A', ''), 
          label_size = 18, nrow = 2, rel_heights = c(2, 1))
dev.off()


plot_grid(changuinolaRestSites, colonSites)
plot_grid(allRestsites, zoomSites, ncol=1, nrow = 2)

ggsave("./output/AllRestingSites.png", width = 12, height = 6.2)

restSites <- hastRest %>% group_by(RestpatchID) %>% 
  filter(!is.na(RestpatchID)) %>% 
  summarize(batsUsed = length(unique(batID)),
            batDays = length(unique(date(timestamp))),
            x = min(RestPatchCtr.x),
            y = min(RestPatchCtr.y),
            nPoints = length(timestamp))
# There are 43 resting locations that were identified from this analysis. One (on Isla Changuinola) was used by 11 bats over 5 different days. This site is also very close to 2-3 other positions like within another 20 m). A couple of other locations were used repeatedly by the same bat (but only by that bat).


#Impacts of social rest on movement ####
#Does resting near another bat change your behavior? Do you rest longer or shorter? Make more foraging movements away?

#This will start with a behavioral change point analysis to count up behavioral sequences.
# I think I can use the patchStay df from 9. Patch Residency Influence since it split up both the in/out of a patch and the behavior in the same loop by bat day. This only was done on the 'forage' behavior class (so no commuting included). This df has duplicates b/c all bat pairs are included for distance, patch effects, etc. I'll keep this analysis under 1000 m. 

load("./processed data/12_AllBatsPatchStateSeq.Rdata") #patchStay_DF

behavBoutTimeDist <- patchStay %>% filter(otherDist < 100) %>% 
  group_by(batID_day, stateBout, newState) %>% 
  dplyr::summarize(timeElapsed = difftime(max(timestamp), min(timestamp), units = "secs"),
            meanDistNN = mean(otherDist, na.rm = TRUE),
            minDistNN  = min(otherDist, na.rm = TRUE),
            sdDistNN = min(otherDist))
stateNo <- c("state1", "state2", "state3", "state4")
behavName <- c("rest", "slow flight" ,"move", "commute")
names(behavName) <- stateNo

behavBoutTimeDist %>%
  ggplot()+
  geom_point(aes(x = minDistNN, y = timeElapsed/60))+
  labs(x = "Distance (m) between nearest neighbors", y = "Duration (min) of bout")+
  facet_wrap(~newState, labeller = labeller(newState = behavName))

#At least with this figure, being in close proximity increases the duration of rest. Try to find a fit for rest
restProx <- behavBoutTimeDist %>% filter(newState == "state1")

fit <- nls((timeElapsed/60) ~ SSasymp(, yf, y0, log_alpha), data = restProx)

library(aomisc)
model <- drm((as.numeric(timeElapsed)/60) ~ minDistNN, fct = DRC.powerCurve(),
             data = restProx)
summary(model)
plot(model, log = "")

# predictions and confidence intervals.
demo.fits <- expand.grid(minDistNN=seq(0.1, 100, length=100))
pm <- predict(model, newdata=demo.fits, interval="confidence") 
demo.fits$time <- pm
demo.fits$pmin <- pm[,2]
demo.fits$pmax <- pm[,3]

pdf("./output/RestingBoutDuration_Distance.pdf", width = 6, height = 4)
png("./output/RestingBoutDuration_Distance.png", width = 6, height = 4, units = "in", res = 96)
ggplot(restProx)+
  geom_point(aes(x = minDistNN, y = as.numeric(timeElapsed)/60), alpha = 0.6)+
  geom_line(aes(x = minDistNN, y = time), data = demo.fits, color = "red")+
  labs(x = "Distance to nearest neighbor (m)", y = "Duration of resting bout (min)")
dev.off()

distRestDur <- ggplot(restProx)+
  geom_point(aes(x = minDistNN, y = as.numeric(timeElapsed)/60), alpha = 0.6)+
  geom_line(aes(x = minDistNN, y = time), data = demo.fits, color = "red")+
  labs(x = "Distance (m)", 
       y = "Duration (min)", 
       subtitle = "Proximity influences on resting bout duration")


top_row <- plot_grid(actBudgOut, actBudgIn, labels = c('A', 'B'), label_size = 12, ncol = 2)
bottom_row <- plot_grid(behavSameFreq, distRestDur, labels = c('C', 'D'), label_size = 12, ncol = 2)


pdf("./output/Fig 2 - Distance Effects on Behavior.pdf", width = 8, height = 8)
plot_grid(top_row, bottom_row, label_size = 12, nrow = 2)
dev.off()

png("./output/Fig 2 - Distance Effects on Behavior.png", width = 1600, height = 1600)
plot_grid(top_row, bottom_row, label_size = 12, nrow = 2)
dev.off()