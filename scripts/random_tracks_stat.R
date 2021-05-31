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
cosd <-function(x) cos(x*pi/180) #cos of an angle in degrees 
sind <-function(x) sin(x*pi/180) 

bird_path <- "data/Kittiwake_data_treated"
bird_filename <- "BLKI_tracks.rds"

data_path <- "outputs/random_tracks"

map_path <- "data/baseline_data"
map_filename <- "SEAwinds_Worldmap_res-c.rds"

### Read BLKI data ----
bird_data <- readRDS(paste0(bird_path,'/', bird_filename)) %>% as.data.table
medlon <- median(bird_data$x, na.rm = T)
medlat <- median(bird_data$y, na.rm = T)

### Define projections and get world map ---- 
proj.aeqd <- paste("+proj=aeqd +lat_0=",round(medlat), " +lon_0=",round(medlon)," +units=km ", sep="")
proj.latlon <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

world_map <- st_as_sf(readRDS(paste0(map_path,"/",map_filename)), coords=c("lon","lat"), crs = 4326)
wrld <- st_transform(world_map,CRS(proj.aeqd))

###Main script ----

#spring and autumn rasters of corridors 

r <- raster(xmn=-70.25, xmx=70.25, ymn=29.75, ymx=85.25, res=0.75) #Empty raster of the studied area
corr.r<- r
corr.r[]<-0
r.names <- c("")

#get the data filenames 
burst <- list.files(path = data_path, pattern='.rds')

#for each burst
stat <- data.frame(ring=c(), month = c(), season = c(), colony= c(), track_type = c(),
                   f_cor = c(), detour = c(), gain_corr = c())  %>%as.data.table()

for (id in burst)
{
  cat(id)
  traj<-readRDS(paste0(data_path, "/", id)) %>% as.data.table
  ring <- append(unique(substring(traj$ring,1,10)),unique(substring(traj$ring,1,11)))
  colony <- bird_data$colony[bird_data$ring%in%ring][1]
  month <- traj$migr_month[1]
  season <- "spring"
  if(as.numeric(month)>6){season <- "autumn"}
  
  traj$track_type <- "sim" 
  traj$track_type[traj$N=="gc"] <- "gc" 
  traj$track_type[traj$N=="obs"] <- "obs" #get type of track
  
  #get optimal tracks as SpatialPoints and corridor as SpatialPolygonsDF
  ws_opt <- quantile(traj$mean_ws,.95)
  opt_trajs <- SpatialPoints(traj[track_type!="observed"&mean_ws>ws_opt][,c("lon","lat")], CRS(proj.latlon))%>% 
    spTransform(., CRS(proj.aeqd))
  
  ke_opt_traj <- kernelUD(opt_trajs) #density kernel of opt trajs
  corr.poly.sp <- getverticeshr(ke_opt_traj, 50) # polygon with density above 50%
  corr.poly<- corr.poly.sp %>% st_as_sf%>%st_transform(. , crs=CRS(proj.latlon)) #into sf class latlon to plot
  
  #get length inside the corridor  
  pt_obs <- SpatialPoints(traj[track_type=="obs"][,c("lon","lat")], CRS(proj.latlon))%>%
    spTransform(., CRS(proj.aeqd))
  l_obs <- pt_obs%>%as(., "SpatialLines")
  
  pts_gc <- SpatialPoints(traj[track_type=="gc"][,c("lon","lat")], CRS(proj.latlon))%>% 
    spTransform(., CRS(proj.aeqd))
  l_gc <- pts_gc%>%as(., "SpatialLines")
  
  f_obs <- gLength(raster::intersect(l_obs,corr.poly.sp))/gLength(l_obs)
  f_gc <- gLength(raster::intersect(l_gc,corr.poly.sp))/gLength(l_gc)
  detour <- max(gDistance(pts_gc, corr.poly.sp,byid=TRUE))
  gain_corr <- ws_opt - traj$mean_ws[traj$track_type=="gc"]
  
  #save stats of this burst
  stat <- rbind(stat, data.frame(ring = traj$ring[1], month = month, season = season, colony=colony,
                                 track_type = c("obs","gc"), f_cor = c(f_obs,f_gc),
                                 detour = detour, gain_corr = gain_corr))
  
  #Save location of corridors
  corr.r.id <- rasterize(corr.poly,r)
  corr.r <- stack(corr.r, corr.r.id)
  r.names <- append(r.names,paste0(colony,"-",season))
  
  ##UNCOMMENT HERE TO PLOT THE CORRIDOR 
  #kde <- kde2d(traj[track_type!="observed"&mean_ws>ws_opt]$lon,
  #traj[track_type!="observed"&mean_ws>ws_opt]$lat)
  #kde.df <- raster(kde)%>% as.data.frame(xy=T)
  
  #traj[, lon2 := data.table::shift(lon, type = 'lead'), by = N]
  #traj[,lat2 := data.table::shift(lat, type = 'lead'), by = N]
  
  #plot_ <- ggplot(data=traj[N=="obs"]) + 
  #geom_raster(data=kde.df, aes(x=x,y=y, fill=layer))+
  #geom_sf(data=corr.poly,fill="green", alpha = 0, size = 2, color="green")+
  #geom_segment(data = traj[N=="obs"], size = 1.25, aes(x = lon, y = lat,xend=lon2,yend=lat2),color='darkgrey') +
  #geom_segment(data = traj[N=="gc"], size = 1.25, aes(x = lon, y = lat,xend=lon2,yend=lat2),color='red') +
  #geom_sf(data=world,fill = "grey", color = "grey") + 
  #coord_sf(xlim = c(-70, 70), ylim = c(30,85), expand = FALSE)
  #print(plot_)
  
}

reg_gain <- lm(stat$gain_corr~stat$detour)
stat[, corr_gain_rel := residuals(reg_gain)]

#plot stats
ggplot(stat)+geom_boxplot(aes(x=season, y = f_cor, fill = track_type),notch=TRUE)

###save wind corridors by burst and by season----
names(corr.r) <- r.names

spr.corr <- calc(corr.r[[which(str_detect(r.names, "spring"))]], fun = sum, na.rm=T)
aut.corr <- calc(corr.r[[which(str_detect(r.names, "autumn"))]], fun = sum, na.rm=T)

writeRaster(corr.r[[-1]],"outputs/Rasters/corr_by_burst",format="raster", overwrite=TRUE)
writeRaster(spr.corr,"outputs/Rasters/summ_corr_spring",format="raster", overwrite=TRUE)
writeRaster(aut.corr,"outputs/Rasters/summ_corr_autumn",format="raster", overwrite=TRUE)
