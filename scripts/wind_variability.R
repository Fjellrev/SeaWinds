##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## SeaWinds - Data processing & Analyses
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sapply(c('sf','sp','spData','tidyverse', 'data.table','MASS','plyr', 'magrittr', 'gdistance','geosphere','raster', "doParallel", "foreach",
         'ggplot2', 'rWind','windR','adehabitatHR', 'circular'),
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

daily_wind_path <- "data/ERA_Interim"
daily_wind_filename <- "ERA_Interim_daily_sfc_2013_to_2018.rds"
wind_data <- readRDS(paste0(daily_wind_path, '/', daily_wind_filename)) %>% as.data.table

wind_data[, wdir := atan2(u, v)] #wind direction
wind_data[, wspeed := sqrt(u^2+v^2)]

r <- raster(xmn=-70.25, xmx=70.25, ymn=29.75, ymx=85.25, res=0.75) #Empty raster of the studied area
r.wdir_var <- rasterize(wind_data[,c("x","y")], r, field = wind_data$wdir, fun=function(x, na.rm=T){sd.circular(x)})
writeRaster(r.wdir_var,"outputs/Rasters/wdir_sd",format="raster", overwrite=TRUE)

df.wind_var <- as.data.frame(r.wdir_var, xy= T) %>% as.data.table

col = rev(c('firebrick4', 'firebrick3', 'gold', 'springgreen3', 'springgreen4'))
plot_ <- ggplot(data=df.wind_var)+
  geom_raster(aes(x=x,y=y,fill=layer))+
  geom_sf(data=world_map)+
  coord_sf(xlim=c(-70,70),ylim=c(30,85))+
  scale_fill_gradientn(colours = col)

plot_
