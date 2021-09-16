#Figure 1
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
load("./processed data/HastatusSegStateBiodatPwr_ShortNightFilt_2020-08-09.Rdata") #hast_df
panama_sf <- st_read("./QGIS/Panama_Province_Boundaries-shp/Province_Boundaries.shp")
isla_sf <- st_read("./QGIS/shapefiles/almiranteIsla.shp")
load("./processed data/PairwiseDistances.Rdata") #pair.distOG

bocas_sf <- panama_sf %>% dplyr::filter(NOMBRE == "Bocas del Toro")
bocas_sf <- st_transform(bocas_sf, crs=32617)  #26917 was used before, but that NAD83
bats_sf <- hast_df %>%  
  st_as_sf(coords = c("utm.x", "utm.y"), crs=26917)
bats_sf_ll <- hast_df %>%  
  st_as_sf(coords = c("location_long", "location_lat"), crs=4326)

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

#Note: using ggmamp and google products will require registering a personal API key.

bocas_map2 <- get_map(batBox, source="stamen", maptype = "terrain")
isla_map <- get_map(ctr, source="google", zoom = 14, maptype = "satellite")

# Create Figure 1 with all tracks map ####
mycols <- viridisLite::viridis(6)[c(3, 4, 6)]

#Manually extract the min / max extents of bocas_map2 to force scalebar into correct placement.
map_xmin <- -82.55
map_xmax <- -82.2
map_ymin <-  9.327
map_ymax <- 9.45

#Bocas Study Map ####
bocasStudyMap <- ggmap(bocas_map2)+
  ggsn::scalebar(x.min = map_xmin, x.max = map_xmax, 
                 y.min = map_ymin, y.max = map_ymax, 
                 location = "bottomleft",
                 dist = 5,
                 dist_unit = "km",
                 model = "WGS84",
                 transform = TRUE,
                 st.bottom = FALSE,
                 st.dist = 0.04,
                 st.size = 8,
                 border.size = 0.5)+
  geom_path(data = hast_df, 
            aes(x = location_long, 
                y = location_lat, 
                group = batID, 
                color = groupID), 
            alpha = 0.8, size = 1)+
  scale_color_manual(values = mycols, 
                     name = "", 
                     breaks = c("blue", "brown", "yellow"),
                     labels = c("Group 1", "Group 2", "Group 3"),
                     guide = guide_legend(override.aes = list(size = 3)))+
  theme(legend.position = c(0.85, 0.9), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 20, face = "bold"),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(colour = NA, fill = NA),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))+
  geom_rect(xmin = isla_box[1],
            xmax = isla_box[2],
            ymin = isla_box[3],
            ymax = isla_box[4], alpha = 0, size = 1, color = "black")+
  labs(x = "longitude", y = "latitude")

yellow3bats <- hast_df %>% filter(batID %in% c("X74DE4E7", "X74D932E","X74DCBCC"))

threeYlwBats <- ggmap(isla_map)+
  ggsn::scalebar(x.min = -82.50, x.max = -82.48,
                 y.min = 9.395, y.max = 9.445,
                 location = "bottomleft",
                 dist = 0.5,
                 dist_unit = "km",
                 model = "WGS84",
                 transform = TRUE,
                 st.bottom = FALSE,
                 st.size = 8,
                 st.dist = 0.04,
                 st.color = "white",
                 border.size = 0.5)+
  geom_sf(data = isla_sf, color = "white", fill = NA, inherit.aes = FALSE)+
  geom_path(data = yellow3bats,
            aes(x = location_long,
                y = location_lat,
                group = batID,
                color = batID),
            alpha = 0.8, 
            size = 2)+
  scale_color_brewer(palette = "OrRd", 
                     guide = guide_legend(override.aes = list(size = 3)))+
  theme(legend.position = c(0.05, 0.93),
        legend.title = element_blank(),
        legend.text = element_text(color = "white", size = 20, face = "bold"),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(colour = NA, fill = NA), 
        axis.text = element_text(size = 20), 
        axis.title = element_text(size = 20))+
  labs(x = "longitude", y = "latitude")

nnDist <- pair.distOG %>% filter(date(timestamps) == "2016-03-04",
                         bat1 %in% c("X74D932E", "X74DCBCC"),
                         bat2 %in% c("X74DE4E7", "X74D932E","X74DCBCC")) %>% 
  ggplot()+
  geom_line(aes(x = timestamps, y = geo.dist/1000), size = 2)+
  labs(x = "hour of night", y = "distance between bats (km)")+
  scale_x_datetime(date_labels = "%H", timezone = "America/Panama")+ 
  facet_grid(bat1~bat2)+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = rel(1.5)))+
  theme(axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)))

bottom_row <- plot_grid(threeYlwBats, nnDist, labels = c('B', 'C'),
                        label_size = 26, ncol = 2) #size 18 for pdf, 22 for png

png("./output/Fig 1 - CompositeTrackingIsla_maps_Distances.png", width = 1800, height = 1125)
plot_grid(bocasStudyMap, bottom_row, labels = c('A', ''), 
          label_size = 26, nrow = 2)
dev.off()


