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

bird_path <- "data/Kittiwake_data_treated"
bird_filename <- "BLKI_Alkefjellet_May_tracks.rds"
data_path <- "outputs/random_tracks"
map_path <- "data/baseline_data"
map_filename <- "SEAwinds_Worldmap_res-c.rds"

### Read BLKI data ----
bird_data <- readRDS(paste0(bird_path,'/', bird_filename)) %>% as.data.table
medlon <- median(bird_data$x, na.rm = T)
medlat <- median(bird_data$y, na.rm = T)

ring<- unique(bird_data$ring) 


### Define projections ---- 
proj.aeqd <- paste("+proj=aeqd +lat_0=",round(medlat), " +lon_0=",round(medlon)," +units=m ", sep="")
proj.latlon <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #classic longlat projection to use the windR function

world_map <- st_as_sf(readRDS(paste0(map_path,"/",map_filename)), coords=c("lon","lat"), crs = 4326)
wrld <- st_transform(world_map,CRS(proj.aeqd))


#traj analysis
stat <- data.frame(track_type = c(), f_cor = c())
for (id in ring)
{
  cat(id)
  data_filename <- paste0("tracks_",id,".rds")
  traj<-readRDS(paste0(data_path, "/", data_filename)) %>% as.data.table
  traj$track_type <- "sim" 
  traj$track_type[traj$N=="gc"] <- "gc" 
  traj$track_type[traj$N=="obs"] <- "obs" 
  traj[, lon2 := data.table::shift(lon, type = 'lead'), by = N]
  traj[,lat2 := data.table::shift(lat, type = 'lead'), by = N]
  
  ws_opt<-sort((traj$mean_ws[traj$track_type!="obs"]))[95/100*nrow(traj[traj$track_type!="obs"])]
  
  opt_trajs <- st_as_sf(x = traj[track_type!="observed"&mean_ws>ws_opt][,c("lon","lat")],
                        coords = c("lon","lat"), crs = CRS(proj.latlon))
  opt_trajs <- st_transform(opt_trajs, CRS(proj.aeqd))
  opt_trajs <- as(opt_trajs, "Spatial")
  ke_opt_traj <- kernelUD(opt_trajs)
  ver_opt_traj <- getverticeshr(ke_opt_traj, 50) %>% st_as_sf%>%st_transform(. , crs=CRS(proj.latlon))
  
  f_obs <- length(which(is.land(traj[traj$track_type=="obs"]$lon,
                                traj[traj$track_type=="obs"]$lat,proj.latlon,ver_opt_traj)))/nrow(traj[traj$track_type=="obs"])
  f_gc <- length(which(is.land(traj[traj$track_type=="gc"]$lon,
                               traj[traj$track_type=="gc"]$lat,proj.latlon,ver_opt_traj)))/nrow(traj[traj$track_type=="gc"])
  
  stat <- rbind(stat, data.frame(track_type = c("obs","gc"),f_cor = c(f_obs,f_gc)))
  
  #kde <- kde2d(traj[track_type!="observed"&mean_ws>ws_opt]$lon,
  #traj[track_type!="observed"&mean_ws>ws_opt]$lat)
  #kde.df <- raster(kde)%>% as.data.frame(xy=T)
  
  #plot_ <- ggplot(data=traj[N=="obs"]) + 
  #geom_raster(data=kde.df, aes(x=x,y=y, fill=layer))+
  #geom_sf(data=ver_opt_traj,fill="green", alpha = 0, size = 2, color="green")+
  #geom_segment(data = traj[N=="obs"], size = 1.25, aes(x = lon, y = lat,xend=lon2,yend=lat2),color='black') +
  #geom_segment(data = traj[N=="gc"], size = 1.25, aes(x = lon, y = lat,xend=lon2,yend=lat2),color='red') +
  #geom_sf(data=world,fill = "grey", color = "grey") + 
  #coord_sf(xlim = c(-60, 60), ylim = c(30,85), expand = FALSE)
  #print(plot_)
  
}

ggplot(stat)+geom_boxplot(aes(x=track_type, y = f_cor), notch = T)
