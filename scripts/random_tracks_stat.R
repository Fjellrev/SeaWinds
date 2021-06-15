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

data_path <- "outputs/random_tracks"
raster_path <- "outputs/Rasters"

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

r.var_autumn <- raster("outputs/Rasters/wdir_sd_autumn")
r.var_spring <- raster("outputs/Rasters/wdir_sd_spring")
#get the data filenames 
burst <- list.files(path = data_path, pattern='.rds')

#for each burst
stat <- data.frame(ring=c(), month = c(), season = c(), colony= c(), track_type = c(),
                   f_cor = c(), predictability = c(), predictability_null =c(), l = c())  %>%as.data.table()

for (id in burst)
{
  cat(id)

  traj<-readRDS(paste0(data_path, "/", id)) %>% as.data.table
  ring <- append(unique(substring(traj$ring,1,10)),unique(substring(traj$ring,1,11)))
  colony <- bird_data$colony[bird_data$ring%in%ring][1]
  month <- traj$migr_month[1]
  season_ <- "spring"
  if(as.numeric(month)>6){season_ <- "autumn"}
  
  traj$track_type <- "sim" 
  traj$track_type[traj$N=="gc"] <- "gc" 
  traj$track_type[traj$N=="obs"] <- "obs" #get type of track
  
  r.var <- get(paste0("r.var_",season_))
  
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
  
  wind_var_in<- extract(r.var,cellFromXY(r.var,traj[traj$track_type =="sim"&traj$mean_ws>=ws_opt][,c("lon","lat")]))%>%
    mean(.,na.rm=T)
  
  wind_var_out<- extract(r.var,cellFromXY(r.var,traj[traj$track_type =="obs"&traj$mean_ws<ws_opt][,c("lon","lat")]))%>%
    mean(.,na.rm=T)
  
  #save stats of this burst
  stat <- rbind(stat, data.frame(ring = traj$ring[1], month = month, season = season_, colony=colony,
                                 track_type = c("obs","gc"), f_cor = c(f_obs,f_gc),
                                 wind_var = c(wind_var_in, wind_var_out), l = c(gLength(l_obs),gLength(l_gc))))
  
  #Save location of corridors
  corr.r.id <- rasterize(corr.poly,r)
  corr.r <- stack(corr.r, corr.r.id)
  r.names <- append(r.names,paste0(colony,"-",season_))
  
# #UNCOMMENT HERE TO PLOT THE CORRIDOR
# kde <- kde2d(traj[track_type!="observed"&mean_ws>ws_opt]$lon,
# traj[track_type!="observed"&mean_ws>ws_opt]$lat)
# kde.df <- raster(kde)%>% as.data.frame(xy=T)
# 
#   traj[, lon2 := data.table::shift(lon, type = 'lead'), by = N]
#   traj[,lat2 := data.table::shift(lat, type = 'lead'), by = N]
#   col = c('firebrick4', 'firebrick3', 'gold', 'springgreen3', 'springgreen4')
# 
#   traj$label = "A. Construction of a wind corridor"
#   plot_corr <- ggplot() +
#   geom_segment(data = traj[N%in%unique(traj$N[traj$track_type=="sim"])[1:1000]],
#                  size = .25, aes(x = lon, y = lat,xend=lon2,yend=lat2, color= mean_ws)) +
#   geom_sf(data=corr.poly,fill="black", alpha = 0.4, size = 2, color="black")+
#   scale_color_gradientn(colors=col, name="Mean \nwind support\n(m/s)")+
#   geom_sf(data=world,fill = "grey", color = "grey") +
#   coord_sf(xlim = c(-70, 70), ylim = c(30,85), expand = FALSE)+xlab("")+ylab("")+
#   theme(strip.text = element_text(face="bold"),
#         legend.title = element_text(hjust = 0.5))+
#     theme(  axis.title.x = element_blank(),
#             axis.text.x = element_blank(),
#             axis.ticks.x = element_blank())+
#     theme(  axis.title.y = element_blank(),
#             axis.text.y = element_blank(),
#             axis.ticks.y = element_blank())+ 
#     ggtitle('B.')
#   print(plot_corr)
#   
#   plot_traj <- ggplot() +
#     geom_segment(data = traj[N=="obs"], size = 1.25, aes(x = lon, y = lat,xend=lon2,yend=lat2),color='blue', show.legend=T) +
#     geom_segment(data = traj[N=="gc"][seq(1,nrow(traj[N=="gc"]),2)], size = 1.25, aes(x = lon, y = lat,xend=lon2,yend=lat2),color='red', linetype = "dashed") +
#     geom_sf(data=world,fill = "grey", color = "grey") +
#     coord_sf(xlim = c(-70, 70), ylim = c(30,85), expand = FALSE)+xlab("")+ylab("")+
#     theme(strip.text = element_text(face="bold"),
#           legend.title = element_text(hjust = 0.5))+
#     theme(  axis.title.x = element_blank(),
#             axis.text.x = element_blank(),
#             axis.ticks.x = element_blank())+ 
#     ggtitle('A.')
#   print(plot_traj)
}

stat$n_burst <- lapply(burst,
                      FUN=function(x){rep(x,2)})%>%unlist
stat$x_lab_fig[stat$track_type=="obs"] <- "Birds' \nrealized routes"
stat$x_lab_fig[stat$track_type=="gc"] <- "Great circle lines \n(theoretical)"
#plot stats
ggpaired(stat, x= "x_lab_fig", y="f_cor", 
         color="track_type",line.color = "gray", line.size = 0.4,palette = "npg", id = "n_burst")+
  facet_wrap(.~season)+xlab("")+ylab("Proportion of the route inside the corridor")+
  theme(legend.position = "none",
        strip.text = element_text(face="bold"))

ggpaired(stat, x= "track_type", y="wind_var", id = "n_burst",
         color="track_type",line.color = "gray", line.size = 0.4,palette = "npg")+
  facet_wrap(.~season)+xlab("")+ylab("Proportion of the route in the \nlow wind variability areas")+
  theme(legend.position = "none",
        strip.text = element_text(face="bold"))


###save wind corridors by burst and by season----
names(corr.r) <- r.names

spr.corr <- calc(corr.r[[which(str_detect(r.names, "spring"))]], fun = sum, na.rm=T)
aut.corr <- calc(corr.r[[which(str_detect(r.names, "autumn"))]], fun = sum, na.rm=T)

writeRaster(corr.r[[-1]],"outputs/Rasters/corr_by_burst",format="raster", overwrite=TRUE)
writeRaster(spr.corr,"outputs/Rasters/summ_corr_spring",format="raster", overwrite=TRUE)
writeRaster(aut.corr,"outputs/Rasters/summ_corr_autumn",format="raster", overwrite=TRUE)
