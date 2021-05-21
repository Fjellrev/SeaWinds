##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## SeaWinds - Data processing & Analyses
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sapply(c('sf','sp','spData','tidyverse', 'data.table','MASS','plyr', 'magrittr', 'gdistance','geosphere','raster', "doParallel", "foreach",
         'ggplot2', 'rWind','windR','adehabitatHR'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE)))

### utility functions ----
source("functions/FUNCTION_OverlapPoly.r")
source('functions/FUNCTION_QuickMap.r')
cosd <-function(x) cos(x*pi/180) #cos of an angle in degrees 
sind <-function(x) sin(x*pi/180) 
proj.latlon <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #classic longlat projection to use the windR function


map_path <- "data/baseline_data"
map_filename <- "SEAwinds_Worldmap_res-c.rds"
world_map <- st_as_sf(readRDS(paste0(map_path,"/",map_filename)), coords=c("lon","lat"), crs = 4326)

wind_path <- "data/ERA_Interim"
wind_filename <- "ERA_Interim_daily_sfc_09_2016_to_06_2017.rds"
wind_data <- readRDS(paste0(daily_wind_path, '/', daily_wind_filename)) %>% as.data.table

wind_data[, wdir := atan2(u, v)] #wind direction

r <- raster(xmn=-70.25, xmx=70.25, ymn=29.75, ymx=85.25, res=0.75) #Empty raster of the studied area
r.wind_var <- rasterize(wind_data[,c("x","y")], r, field = wind_data$wdir, fun = sd)
df.wind_var <- as.data.frame(r.wind_var, xy= T) %>% as.data.table

col = c('firebrick4', 'firebrick3', 'gold', 'gold', 'springgreen3', 'springgreen4')
ggplot(data=df.wind_var)+
  geom_raster(aes(x=x,y=y,fill=layer))+
  geom_sf(data=world_map)+
  coord_sf(xlim=c(-70,70),ylim=c(30,85))+
  scale_fill_gradientn(colours = col)
