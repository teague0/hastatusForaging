#Social proximity plots

#Where are bats in close social proximity?

library(tidyverse)
library(lubridate)
library(sf)
library(sp)
library(ggmap)
library(ggsn)
library(raster)
library(cowplot)
theme_set(theme_cowplot())

#Load necessary data
load("./data/12_AllBatsPatchStateSeq.Rdata") #patchStay_df
panama_sf <- st_read("./data/Panama_Province_Boundaries-shp/Province_Boundaries.shp") #CRS 32617
isla_sf <- st_read("./data/shapefiles/almiranteIsla.shp")
load("./data/6_PairwiseDistances.Rdata") #pair.dist

bocas_sf <- panama_sf %>% dplyr::filter(NOMBRE == "Bocas del Toro")#UTM zone 17 NAD84 32617
bocas_sf <- st_transform(bocas_sf, crs=32617) 
bats_sf <- hastMorph %>%  
  st_as_sf(coords = c("utm.x", "utm.y"), crs=32617)
bats_sf_ll <- hastMorph %>%  
  st_as_sf(coords = c("location.long", "location.lat"), crs=4326)

#Crop Bocas del Toro province to the study area
bats_bboxll <- st_bbox(bats_sf_ll)
bats_bbox <- st_bbox(bats_sf)
bats_bboxB <- bats_bbox
bats_bboxB[1] <- bats_bbox[1] - 5000
bats_bboxB[3] <- bats_bbox[3] + 5000
bats_bboxB[2] <- bats_bbox[2] - 2000
bats_bboxB[4] <- bats_bbox[4] + 2000
bocas_crop <- st_crop(bocas_sf, bats_bboxB) 

# Define extents and locations for background images from ggmap
e <- bbox(extent(bats_sf_ll)*1.5)
batBox <- c(left = e[1], 
            bottom = e[2] - 0.01, 
            right = e[3], 
            top = e[4])
isla_box <- extent(isla_sf)

isla_box_map <- c(left = isla_box[1], 
                  bottom = isla_box[3] - 0.005, 
                  right = isla_box[2], 
                  top = isla_box[4] - 0.005)
ctr <- c(-82.48, 9.42)

coords <- cbind(-82.271541, 9.396448)
lagruta <- SpatialPoints(coords, proj4string = CRS("+proj=longlat +datum=WGS84"))
lagrutaUTM <- spTransform(lagruta, CRS("+proj=utm +zone=17 +datum=WGS84"))
lagruta_sf <- st_transform(st_as_sf(lagrutaUTM))

closeProx <- patchStay_df %>% filter(otherDist<291)

mycols <- viridisLite::plasma(9)[c(1, 5, 8)]

bdt_290 <- ggplot()+
  geom_point(data = closeProx, aes(x = utm.x, y = utm.y, color = groupID), alpha = 0.4, size = 1)+
  geom_sf(data = bocas_crop, color = "black", fill = NA, inherit.aes = FALSE)+
  geom_sf(data = isla_sf, color = "black", fill = NA, inherit.aes = FALSE)+
  geom_sf(data = lagruta_sf, color = "black", fill = "grey", shape = 23, size = 2, inherit.aes = FALSE)+
  scale_color_manual(values = mycols, 
                     name = "", 
                     breaks = c("blue", "brown", "yellow"),
                     labels = c("Group 1", "Group 2", "Group 3"),
                     guide = guide_legend(override.aes = list(size = 3)))+
    theme(legend.position = "none", 
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6))+
  labs(x = "Longitude", y = "Latitude")


chang_290 <-ggplot()+
  geom_point(data = closeProx, aes(x = utm.x, y = utm.y, color = groupID), alpha = 0.3, size = 1)+
  geom_sf(data = bocas_crop, color = "black", fill = NA, inherit.aes = FALSE)+
  geom_sf(data = isla_sf, color = "black", fill = NA, inherit.aes = FALSE)+
  scale_color_manual(values = mycols, 
                     name = "", 
                     breaks = c("blue", "brown", "yellow"),
                     labels = c("Group 1", "Group 2", "Group 3"),
                     guide = guide_legend(override.aes = list(size = 3)))+
  lims(x = c(335800.8, 338184), y = c(1040962, 1043526))+
  theme(legend.position = "none", 
        axis.text.x = element_text(angle=45, vjust=0.5, size = 6), 
        axis.text.y = element_text(size = 6))+
  labs(x = "Longitude", y = "Latitude")
chang_290  

dens290 <- ggplot(patchStay_df)+
  geom_density(aes(x = otherDist, fill = groupID), alpha = 0.5)+
  geom_vline(xintercept = 290, lty = 2)+
  scale_fill_manual(values = mycols, 
                     name = "", 
                     breaks = c("blue", "brown", "yellow"),
                     labels = c("Group 1", "Group 2", "Group 3"),
                     guide = guide_legend(override.aes = list(size = 3)))+
  scale_x_continuous(breaks = c(1000, 5000, 10000, 15000, 20000, 25000))+
  theme(legend.position = c(0.5, 0.7), 
        axis.text.x = element_text(angle=45, vjust=0.5, size=12))+
  xlab("distance between tagged bats (m)")

#Plot of distance over time doesn't tell us anything
# library(scales)
# distTime <- patchStay_df %>% 
#   filter(!is.na(otherDist)) %>% 
#   ggplot()+
#   geom_point(aes(x = as_hms(timestamp), y = otherDist, color = groupID), alpha = 0.2)


bottom_row <- plot_grid(chang_290, dens290, labels = c('B', 'C'),
                        label_size = 16, ncol = 2) #size 18 for pdf, 22 for png

pdf("./output/Fig 2 - SocialProx290.pdf", width = 10, height = 4)
plot_grid(bdt_290, bottom_row, labels = c('A', ''), nrow = 2, label_size = 16)
dev.off()

png("./output/Fig 2 - SocialProx290.png", width = 2400, height = 1200, res = 300)
plot_grid(bdt_290, bottom_row, labels = c('A', ''), label_size = 16, nrow = 2)
dev.off()
