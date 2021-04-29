##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## SeaWinds - Data processing & Analyses
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sapply(c('sf','spData','tidyverse', 'data.table', 'magrittr', 'gdistance','geosphere','raster',
         'ggplot2', 'rWind','windR'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE)))

data("world")

bird_path <- "data/Kittiwake_data_treated"
bird_filename <- "BLKI_Alkefjellet_May_tracks.rds"
wind_path <- "data/ERA-Interrim/Interanual_means" 
wind_filename <-  "ERA_Interim_interanual_monthly_mean_sfc_10_2013_to_2018.rds"

###Functions ---

is.land <- function(x,y) #return TRUE if a point is over lands
{
  pt <- expand.grid(x, y)
  pt <- st_as_sf(pt, coords=1:2, crs=CRS(proj.aeqd))
  !is.na(as.numeric(suppressMessages(st_intersects(pt, wrld))))
}

cosd <-function(x) cos(x*pi/180) #cos of an angle in degrees 
sind <-function(x) sin(x*pi/180) 

###Read BLKI data ---

bird_data <- readRDS(paste0(bird_path,'/', bird_filename)) %>% as.data.table

###Read and transform wind data ---

wind_data <- readRDS(paste0(wind_path,"/", wind_filename))%>% as.data.table

wind_data[, wdir := atan2(u, v)*180/pi, by = 1:nrow(wind_data)] #get wind direction
wind_data[, wdir := wdir+360*(wdir<0), by = 1:nrow(wind_data)] #correct wind direction to be between 0 and 360° 
wind_data[, wspeed     := sqrt(u^2 + v^2), by = 1:nrow(wind_data)] #get wind speed

r <- raster(xmn=-70.25, xmx=70.25, ymn=29.75, ymx=85.25, res=.75) #Empty raster of the studied area

r_wdir <- rasterize(wind_data[,c("x","y")], r, field=wind_data$wdir) #raster of wind direction
r_wspeed <- rasterize(wind_data[,c("x","y")], r, field=wind_data$wspeed) #raster of wind speed

wind_layer <- stack(r_wdir, r_wspeed) #raster layer
names(wind_layer) <- c("direction", "speed")
Conductance <- flow.dispersion(x=wind_layer, output="transitionLayer")

###Projections --- 

proj.latlon <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #classic longlat projection to use the windR function

medlon <- median(bird_data$x, na.rm = T)
medlat <- median(bird_data$y, na.rm = T)
proj.aeqd <- paste("+proj=aeqd +lat_0=",round(medlat), " +lon_0=",round(medlon)," +units=m ", sep="")

bird_data.sp <- st_as_sf(bird_data, coords=4:5, crs=CRS(proj.latlon))
bird_data.sp <- st_transform(bird_data.sp, CRS(proj.aeqd))
bird_data.proj <- as.data.table(bird_data.sp) 
bird_data.proj <- bird_data.proj[,-c("geometry")]
bird_data.proj[,x := unlist(map(bird_data.sp$geometry,1))]
bird_data.proj[,y := unlist(map(bird_data.sp$geometry,2))]

wrld <- st_transform(world,CRS(proj.aeqd))

###VAR --- 

#parameters of the simulation
n <- 5 #number of iterations to produce a track --> track with 2^n segments
N <- 50 #number of tracks created

###Main script ---
id <- unique(bird_data.proj$ring)[17] #temporary : will be changed into : for (id in unique(bird_data.proj$ring))
traj <- data.table(ring = c(), N=c(), track_type = c(), x = c(), y = c())

#first trajectory = observed track

traj_id <- bird_data.proj[bird_data.proj$ring==id]
start_pt <- traj_id[,c("x","y")][1] 
end_pt <- traj_id[,c("x","y")][nrow(traj_id)]
traj_0 <- data.table(ring=rep(id,nrow(traj_id)),
                     N = rep("observed",nrow(traj_id)),
                     track_type = rep("observed",nrow(traj_id)),
                     x= traj_id$x, y = traj_id$y)
traj <- rbind(traj,traj_0)

#N random tracks
rm <- c()
a <- sqrt((start_pt$x-end_pt$x)^2+(start_pt$y-end_pt$y)^2)/2 #magnitude of the deviation
for (k in seq (1:N))
{
  cat(k)
  traj_x <- c(start_pt$x,end_pt$x) #x coords of the track
  traj_y <- c(start_pt$y,end_pt$y)
  
  for (j in seq(1:n)) 
  {
    traj_x_ <- c(traj_x[1]) #new x coords 
    traj_y_ <- c(traj_y[1])
    for (i in seq(1:(length(traj_x)-1))) #create a new point betwen point i and i+1 
    {
      v <- c(traj_x[i+1]-traj_x[i], traj_y[i+1]-traj_y[i]) 
      vn <- c(-v[2], v[1])/sqrt(v[1]^2+v[2]^2) #orthogonal vector to the segment treated 
      r <- runif(1,-1,1)
      xi <- (traj_x[i]+traj_x[i+1])/2+a/(2^(j-1))*r*vn[1] #coords of the new point
      yi <- (traj_y[i]+traj_y[i+1])/2+a/(2^(j-1))*r*vn[2]
      nland=0 #number of times the point has been created above the lands (in case there is no choice)
      while(is.land(xi,yi)&nland<10)#is it above the lands
      {
        nland=nland+1
        r <- runif(1,-1,1)
        xi <- (traj_x[i]+traj_x[i+1])/2+a/(2^(j-1))*r*vn[1]
        yi <- (traj_y[i]+traj_y[i+1])/2+a/(2^(j-1))*r*vn[2]
        if (nland == 10){rm <- append(rm, k)}
      }
      traj_x_ <- append(traj_x_, c(xi,traj_x[i+1]))
      traj_y_ <- append(traj_y_, c(yi,traj_y[i+1]))
    }
    traj_x <- traj_x_
    traj_y <- traj_y_
  }
  traj_k <- data.table(ring=rep(id,length(traj_x)),N=rep(k,length(traj_x)), track_type = rep(paste0("sim"),length(traj_x)), x = traj_x, y = traj_y)
  traj <- rbind(traj,traj_k)
}

#into geographic coord
traj.sp <- st_as_sf(traj, coords=4:5, crs=CRS(proj.aeqd))
traj.sp <- st_transform(traj.sp, CRS(proj.latlon))
traj <- as.data.table(traj.sp) 
traj <- traj[,-c("geometry")]
traj[,x := unlist(map(traj.sp$geometry,1))]
traj[,y := unlist(map(traj.sp$geometry,2))]

#remove above land and outside the study area
rm_2 <- unique(traj$N[traj$x < -70|traj$x>70|traj$y>85|traj$y<30])
traj <- traj[not(traj$N%in%append(rm,rm_2))]

#last traj = great circle line
start_pt_latlon <- traj[traj$ring==id&traj$track_type=="observed"][1]
end_pt_latlon <-  traj[traj$ring==id&traj$track_type=="observed"][nrow(traj[traj$ring==id&traj$track_type=="observed"])]
traj_gc <- data.table(ring = rep(id,52),N=rep("gc",52), track_type=rep("gc",52))
gctraj <- gcIntermediate(c(start_pt_latlon$x,start_pt_latlon$y),c(end_pt_latlon$x,end_pt_latlon$y), n = 50, addStartEnd = TRUE, sp = TRUE)
traj_gc[, x :=geom(gctraj)[,"x"]]
traj_gc[, y :=geom(gctraj)[,"y"]]
traj <- rbind(traj,traj_gc)

#Analysis of tracks ---

traj[, x2 := data.table::shift(x, type = 'lead'), by = N]
traj[, y2 := data.table::shift(y, type = 'lead'), by = N]

#1. speed and direction
gdir=c()
gspeed=c()
for (i in (1:nrow(traj)))
{
  gspeed = append(gspeed,distGeo(c(traj$x[i], traj$y[i]),c(traj$x2[i], traj$y2[i])))
  gdir = append(gdir, geosphere::bearing(c(traj$x[i], traj$y[i]),c(traj$x2[i], traj$y2[i])))
  
}
traj[, gspeed := gspeed]
traj[, gdir := gdir]

#2. wind conditions
u_wind = c()
v_wind = c()
for (i in 1:nrow(traj)) #get the wind data at the given coordinates and dates
{
  cat(paste0(i," "))
  uv_wind = getWind(traj$x[i], traj$y[i], w = wind_data, PROJ = proj.latlon)
  u_wind = append(u_wind, uv_wind[1])
  v_wind = append(v_wind, uv_wind[2])
}
traj[, u_wind := as.numeric(u_wind)]
traj[, v_wind := as.numeric(v_wind)]

traj[, wdir := atan2(u_wind, v_wind)*180/pi, by = 1:nrow(traj)] #wind direction
traj[, wspeed     := sqrt(u_wind^2 + v_wind^2), by = 1:nrow(traj)] #wind speed
traj[, ws := wspeed*cosd(wdir -gdir)] #wind support
traj[, mean_ws := mean(traj$ws[traj$N==N], na.rm= T), by = 1:nrow(traj)]
#3. cost
cost <- c()
for (i in 1:nrow(traj))
{
  cat(paste0(i," "))
  if (is.na(traj$x2[i])|traj$y[i]>85|traj$y2[i]>85)
  {cost <- append(cost, NA)}
  else
  {cost <- append(cost, costDistance(Conductance, c(traj$x[i],traj$y[i]), c(traj$x2[i],traj$y2[i])))} 
}
traj$cost <- cost
traj[, cost_tot := sum(traj$cost[traj$N==N], na.rm= T), by = 1:nrow(traj)]
opt_cost <- sort(unique(traj$cost_tot))[floor(5/100*length(unique(traj$N)))]
  
#saveRDS(traj, file = paste0("outputs/observed_sim_gc_tracks/tracks_",id,".rds")


###Display --- 

col_ws = c('firebrick4', 'firebrick3', 'gold', 'gold', 'springgreen3', 'springgreen4')

plot_ <- ggplot(data=traj[track_type!="observed"]) + 
  geom_segment(size = .75, aes(x = x, y = y,xend=x2,yend=y2,color=mean_ws)) +
  scale_colour_gradientn(colours = col_ws)+
  geom_segment(data = traj[track_type=="observed"], size = 1.25, aes(x = x, y = y,xend=x2,yend=y2),color='black') +
  geom_sf(data=world,fill = "black", color = "black") + 
  geom_point(aes(x=start_pt$x, y =start_pt$y), color = 'blue', size = 2.5) +
  geom_point(aes(x=end_pt$x, y =end_pt$y), color = 'red', size = 2.5) +
  coord_sf(xlim = c(-60, 60), ylim = c(30,85), expand = FALSE)+
  ggtitle(paste0("n = ",n, " ; ring = ",id))
print(plot_)

r <- raster(xmn=-70.25, xmx=70.25, ymn=29.75, ymx=85.25, res=1) #Empty raster of the studied area

r.ws <- rasterize(traj[track_type!="observed"][,c("x","y")],r,field=traj[track_type!="observed"]$mean_ws,fun=mean)%>%as.data.frame(xy=T)
plot_ <- ggplot(data=r.ws) + 
  geom_raster(aes(x = x, y = y,fill=layer)) +
  scale_fill_gradientn(colours = col_ws, name="WS" )+
  geom_segment(data = traj[track_type=="observed"], size = 1.25, aes(x = x, y = y,xend=x2,yend=y2),color='black') +
  geom_segment(data = traj[track_type=="gc"], size = 1.25, aes(x = x, y = y,xend=x2,yend=y2),color='black') +
  geom_sf(data=world,fill = "black", color = "black") + 
  coord_sf(xlim = c(-60, 60), ylim = c(30,85), expand = FALSE)+
  ggtitle(paste0("n = ",n, " ; ring = ",id))
print(plot_)


