##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## SeaWinds - Data processing & Analyses
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sapply(c('sf','sp','spData','tidyverse', 'data.table','rgeos','MASS','plyr', 'magrittr', 'gdistance','geosphere','raster', "doParallel", "foreach",
         'ggplot2', 'rWind','windR','adehabitatHR',"ggpubr"),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE)))

### utility functions ----

source("functions/FUNCTION_OverlapPoly.r")
source('functions/FUNCTION_QuickMap.r')
cosd <-function(x) cos(x*pi/180) #cos of an angle in degrees 
sind <-function(x) sin(x*pi/180) 

bird_path <- "data/Kittiwake_data_treated"
bird_filename <- "BLKI_tracks.rds"

map_path <- "data/baseline_data"
map_filename <- "SEAwinds_Worldmap_res-c.rds"

### Read BLKI data ----
bird_data <- readRDS(paste0(bird_path,'/', bird_filename)) %>% as.data.table

world_map <- st_as_sf(readRDS(paste0(map_path,"/",map_filename)), coords=c("lon","lat"), crs = 4326)

bird_data[, x2 := data.table::shift(x, type = 'lead'), by = ring]
bird_data[, y2 := data.table::shift(y, type = 'lead'), by = ring]
bird_data[, n_seg := length(unique(bird_data$burst[bird_data$ring==ring])), by = 1:nrow(bird_data)]
bird_data[, dist :=distGeo(c(x, y),c(x2, y2)), by = 1:nrow(bird_data)]

bird_data[, length_id := sum(bird_data$dist[bird_data$ring==ring], na.rm=T), by = 1:nrow(bird_data)]
bird_data[, length_burst := sum(bird_data$dist[bird_data$burst==burst], na.rm=T), by = 1:nrow(bird_data)]
for (col in unique(bird_data$colony))
{
  cat(paste0(col, "\n"))
  #print(summary(bird_data$length_id[bird_data$colony==col]/1000))
  #print(sd(bird_data$length_id[bird_data$colony==col], na.rm = T)/1000)
  
  #print(summary(bird_data$n_seg[bird_data$colony==col]))
  #print(sd(bird_data$n_seg[bird_data$colony==col], na.rm = T))
  
  print(summary(bird_data$length_burst[bird_data$colony==col]/1000))
  print(sd(bird_data$length_burst[bird_data$colony==col], na.rm = T)/1000)
}
