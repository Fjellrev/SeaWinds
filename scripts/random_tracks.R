##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## SeaWinds - Data processing & Analyses
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sapply(c('sf','spData','tidyverse', 'data.table', 'magrittr', 'gdistance','geosphere','raster', "doParallel", "foreach",
         'ggplot2', 'rWind','windR'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE)))


### utility functions ----
source("functions/FUNCTION_OverlapPoly.r")
source('functions/FUNCTION_QuickMap.r')
cosd <-function(x) cos(x*pi/180) #cos of an angle in degrees 
sind <-function(x) sin(x*pi/180) 

bird_path <- "data/Kittiwake_data_treated"
bird_filename <- "BLKI_Alkefjellet_May_tracks.rds"
wind_path <- "data/ERA_Interim/Interannual_means" 
wind_filename <-  "ERA_Interim_interannual_monthly_mean_sfc_10_2013_to_2018.rds"
map_path <- "data/baseline_data"
map_filename <- "SEAwinds_Worldmap_res-c.rds"

### Read BLKI data ----
bird_data <- readRDS(paste0(bird_path,'/', bird_filename)) %>% as.data.table
medlon <- median(bird_data$x, na.rm = T)
medlat <- median(bird_data$y, na.rm = T)

### Define projections ---- 
proj.aeqd <- paste("+proj=aeqd +lat_0=",round(medlat), " +lon_0=",round(medlon)," +units=m ", sep="")
proj.latlon <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #classic longlat projection to use the windR function

world_map <- st_as_sf(readRDS(paste0(map_path,"/",map_filename)), coords=c("lon","lat"), crs = 4326)
wrld <- st_transform(world_map,CRS(proj.aeqd))

### Transform BLKI data ----
bird_data.sp <- st_as_sf(bird_data, coords=4:5, crs=CRS(proj.latlon))
bird_data.sp <- st_transform(bird_data.sp, CRS(proj.aeqd))
bird_data.proj <- as.data.table(bird_data.sp) 
bird_data.proj <- bird_data.proj[,-c("geometry")]
bird_data.proj[,x := unlist(map(bird_data.sp$geometry,1))]
bird_data.proj[,y := unlist(map(bird_data.sp$geometry,2))]

### Read and transform wind data ----
wind_data <- readRDS(paste0(wind_path,"/", wind_filename))%>% as.data.table
wind_data[, wdir := atan2(u, v)*180/pi, by = 1:nrow(wind_data)] #get wind direction
wind_data[, wdir := wdir+360*(wdir<0), by = 1:nrow(wind_data)] #correct wind direction to be between 0 and 360° 
wind_data[, wspeed     := sqrt(u^2 + v^2), by = 1:nrow(wind_data)] #get wind speed

### parameters of the simulation ----
n <- 5 #number of iterations to produce a track --> track with 2^n segments
N <- 10000 #number of tracks created

### Main script ----

for (id in unique(bird_data.proj$ring)){
  
  start.time <- Sys.time()
  # traj <- data.table(ring = c(), N = c(), x = c(), y = c())
  
  #first trajectory = observed track
  
  bird_data.proj %>%
    dplyr::filter(ring == id) %>%
    dplyr::mutate(N = 'obs') %>%
    dplyr::select(ring, N, x, y) %>%
    tibble() -> 
    traj_0
  
  start_pt <- traj_0[1 ,c("x","y")] 
  end_pt   <- traj_0[nrow(traj_0),c("x","y")]
  
  #N random tracks
  
  a <- sqrt((start_pt$x - end_pt$x)^2 + (start_pt$y - end_pt$y)^2)/3 #magnitude of the deviation
  
  # Registering backend for parallel computing
  n.cores <- detectCores() - detectCores() %/% 10
  cl      <- makeCluster(n.cores, sep = "")) #, type = 'FORK')
registerDoParallel(cl)
# getDoParWorkers()

traj <- foreach(k = 1:N, .errorhandling = 'pass', .packages = c("data.table", "rgdal", "sf", "sp")) %dopar% {  
  
  source("functions/FUNCTION_OverlapPoly.r")
  
  # cat(k, "\n")
  traj_x <- c(start_pt$x,end_pt$x) #x coords of the track
  traj_y <- c(start_pt$y,end_pt$y)
  
  for (j in 1:n){
    
    traj_x_ <- c(traj_x[1]) #new x coords 
    traj_y_ <- c(traj_y[1])
    
    for (i in 1:(length(traj_x)-1)) #create a new point betwen point i and i+1 
    {
      n.pts <- 100
      v  <- c(traj_x[i+1]-traj_x[i], traj_y[i+1]-traj_y[i]) 
      vn <- c(-v[2], v[1])/sqrt(v[1]^2+v[2]^2) #orthogonal vector to the segment treated 
      r  <- runif(n.pts,-1,1)
      xi <- (traj_x[i]+traj_x[i+1])/2+a/(2^(j-1))*r*vn[1] #coords of the new points
      yi <- (traj_y[i]+traj_y[i+1])/2+a/(2^(j-1))*r*vn[2]
      
      land_out <- is.land(x = xi, y = yi, prj = proj.aeqd, mask = wrld)
      
      #quick.map(xi[which(land_out == FALSE)], yi[which(land_out == FALSE)],wrld) # just checking where the new locations are
      a2 <- a
      
      while(length(which(land_out == FALSE)) < 3 & a2 < 10*a){ # check that enough points are generated above ocean
        a2 <- 1.1*a2
        r  <- runif(n.pts,-1,1)
        xi <- (traj_x[i]+traj_x[i+1])/2+a2/(2^(j-1))*r*vn[1]
        yi <- (traj_y[i]+traj_y[i+1])/2+a2/(2^(j-1))*r*vn[2]
        land_out <- is.land(x = xi, y = yi, prj = proj.aeqd, mask = wrld)
      }
      
      if(length(which(land_out == FALSE)) > 0){
        ii <- sample(which(land_out == FALSE), 1)
        xi <- xi[ii]; yi <- yi[ii]
      } else {
        ii <- sample(1:length(land_out), 1)
        xi <- xi[ii]; yi <- yi[ii]
      }
      
      traj_x_ <- append(traj_x_, c(xi,traj_x[i+1]))
      traj_y_ <- append(traj_y_, c(yi,traj_y[i+1]))
    }
    traj_x <- traj_x_
    traj_y <- traj_y_
  }
  
  traj_k <- data.table(ring=rep(id,length(traj_x)), N = formatC(rep(k, length(traj_x)), flag = '0', width = nchar(max(N))), x = traj_x, y = traj_y)
  traj_k
}

stopCluster(cl) # close connection to cluster/cores

rbindlist(traj) %>%
  bind_rows(traj_0, .) %>%
  sf::st_as_sf(. ,coords=c("x","y"), crs=CRS(proj.aeqd)) %>%
  sf::st_transform(CRS(proj.latlon)) %>%  # reproject to longlat
  dplyr::mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2]) %>%
  tibble() %>%
  dplyr::select(ring, N, lon, lat) %>%
  dplyr::filter(!N %in% unique(.$N[.$lon < -70 | .$lon > 70 | .$lat > 85 | .$lat < 30])) -> #remove above land and outside the study area
  traj


#last traj = great circle line
start_pt_latlon <- traj[traj$N == 'obs', c("lon", "lat")][1,]
end_pt_latlon   <- traj[traj$N == 'obs', c("lon", "lat")][nrow(traj[traj$N == 'obs',]),]

gctraj <- gcIntermediate(start_pt_latlon, end_pt_latlon, n = 2^n, addStartEnd = TRUE, sp = TRUE)

tibble(ring = rep(id, 2^n+2), N = rep('gc', 2^n+2)) %>%
  dplyr::mutate(lon = geom(gctraj)[,"x"], lat = geom(gctraj)[, "y"]) %>%
  dplyr::bind_rows(traj, .) ->
  traj

split(traj, as.factor(traj$N)) %>%
  purrr::map(~ mutate(., gspeed = c(geosphere::distGeo(.[1:(nrow(.)-1), c("lon","lat")], .[2:(nrow(.)), c("lon","lat")]), NA))) %>%
  purrr::map(~ mutate(., gdir = c(geosphere::bearing(.[1:(nrow(.)-1), c("lon","lat")], .[2:(nrow(.)), c("lon","lat")]), NA))) %>%
  dplyr::bind_rows() ->
  traj


# Registering backend for parallel computing
n.cores <- detectCores() - detectCores() %/% 10
cl      <- makeCluster(n.cores, outfile = paste("outputs/log_makeCluster_", format(Sys.time(), "%Y%m%d_%H%M%S"),".txt", sep = "")) #, type = 'FORK')
registerDoParallel(cl)
# getDoParWorkers()

wind_df <- foreach(i = 1:nrow(traj), .errorhandling = 'pass', .packages = c("data.table", "rgdal", "sf", "sp", "windR")) %dopar% {  
  
  uv_wind = getWind(traj$lon[i], traj$lat[i], w = wind_data, PROJ = proj.latlon)
  cbind.data.frame(u = uv_wind[[1]], v = uv_wind[[2]])
}

stopCluster(cl)

wind_df <- rbindlist(wind_df)

traj %>%
  dplyr::mutate(u_wind = as.numeric(wind_df$u),
                v_wind = as.numeric(wind_df$v),
                wdir = atan2(u_wind, v_wind)*180/pi, #wind direction
                wspeed = sqrt(u_wind^2 + v_wind^2),   #wind speed
                ws = wspeed*cosd(wdir -gdir)) %>%     #wind support
  dplyr::group_by(N) %>%
  dplyr::mutate(mean_ws = mean(ws, na.rm = T)) %>%
  dplyr::ungroup() ->
  traj

saveRDS(traj, file = paste0("outputs/observed_sim_gc_tracks/tracks_",id,"_n", n, ".rds"))

cat(id, " --- ", round(difftime(Sys.time(), start.time, units = "secs"),1), "sec\n")

}

#### Display ---- 

#col_ws = c('firebrick4', 'firebrick3', 'gold', 'gold', 'springgreen3', 'springgreen4')

#plot_ <- ggplot(data=traj[track_type!="observed"]) + 
#geom_segment(size = .75, aes(x = x, y = y,xend=x2,yend=y2,color=mean_ws)) +
#scale_colour_gradientn(colours = col_ws)+
#geom_segment(data = traj[track_type=="observed"], size = 1.25, aes(x = x, y = y,xend=x2,yend=y2),color='black') +
#geom_sf(data=world,fill = "black", color = "black") + 
#geom_point(aes(x=start_pt$x, y =start_pt$y), color = 'blue', size = 2.5) +
#geom_point(aes(x=end_pt$x, y =end_pt$y), color = 'red', size = 2.5) +
#coord_sf(xlim = c(-60, 60), ylim = c(30,85), expand = FALSE)+
#ggtitle(paste0("n = ",n, " ; ring = ",id))
#print(plot_)

#r <- raster(xmn=-70.25, xmx=70.25, ymn=29.75, ymx=85.25, res=1.5) #Empty raster of the studied area

#r.ws <- rasterize(traj[track_type!="observed"][,c("x","y")],r,field=traj[track_type!="observed"]$mean_ws,fun=mean)%>%as.data.frame(xy=T)
#plot_ <- ggplot(data=r.ws) + 
#geom_raster(aes(x = x, y = y,fill=layer)) +
#scale_fill_gradientn(colours = col_ws, name="WS" )+
#geom_segment(data = traj[track_type=="observed"], size = 1.25, aes(x = x, y = y,xend=x2,yend=y2),color='black') +
#geom_segment(data = traj[track_type=="gc"], size = 1.25, aes(x = x, y = y,xend=x2,yend=y2),color='black') +
#geom_sf(data=world,fill = "black", color = "black") + 
#coord_sf(xlim = c(-60, 60), ylim = c(30,85), expand = FALSE)+
#ggtitle(paste0("n = ",n, " ; ring = ",id))
#print(plot_)