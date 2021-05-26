##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## SeaWinds - Data processing & Analyses
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sapply(c('sf','sp','spData','tidyverse', 'data.table','rgeos','MASS','plyr', 'magrittr', 'gdistance','geosphere','raster', "doParallel", "foreach",
         'ggplot2', 'rWind','windR','adehabitatHR'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE)))

### utility functions ----

source("functions/FUNCTION_OverlapPoly.r")
source('functions/FUNCTION_QuickMap.r')
source('functions/FUNCTION_lineDens.r')
cosd <-function(x) cos(x*pi/180) #cos of an angle in degrees 
sind <-function(x) sin(x*pi/180) 


###Get map ----

map_path <- "data/baseline_data"
map_filename <- "SEAwinds_Worldmap_res-c.rds"
world_map <- st_as_sf(readRDS(paste0(map_path,"/",map_filename)), coords=c("lon","lat"), crs = 4326)

### Read BLKI data ----

bird_path <- "data/Kittiwake_data_treated"
bird_filename <- "BLKI_tracks.rds"
bird_data <- readRDS(paste0(bird_path,'/', bird_filename)) %>% as.data.table

###Define projections ----

medlon <- median(bird_data$x, na.rm = T)
medlat <- median(bird_data$y, na.rm = T)
proj.aeqd <- paste("+proj=aeqd +lat_0=",round(medlat), " +lon_0=",round(medlon)," +units=km ", sep="")
proj.latlon <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

###Get bursts in spatialLines ----

l.bursts <- foreach(id = unique(bird_data$burst)) %do% {
  l.burst.id <- Line(bird_data[burst==id][,c("x","y")])
  Lines(list(l.burst.id), ID = id)
}
l.bursts <- SpatialLines(l.bursts)

###Get density raster ----

r <- raster(xmn=-70.25, xmx=70.25, ymn=29.75, ymx=85.25, res=0.75) #Empty raster of the studied area

r.dens <- lineDens(r, l.bursts)

df.dens.burst <- r.dens %>%as.data.frame(xy=T) #as df to plot

(ggplot()+
  geom_raster(data=df.dens.burst, aes(x=x,y=y,fill=layer))+
  geom_sf(data=world_map)+
  coord_sf(xlim=c(-70,70),ylim=c(30,85)))
