##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## SeaWinds - Data processing & Analyses
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sapply(c('sp','spData','raster','fields','lubridate', 'shape', 'data.table',"rnaturalearthdata","rnaturalearth", 'magrittr', "maptools", 'rWind','rworldmap', 'gdistance','geosphere', 'ggplot2'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE)))

data("world")

path_bird <- "data/Kittiwake_data_treated"


###Read BLKI data ---

birdRDS <- list.files(path=path_bird, pattern = "tracks")
bird_data <- readRDS(paste0(path_bird,'/', birdRDS)) %>% as.data.table

###Projections --- 

proj.latlon <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #classic longlat projection to use the windR function

medlon <- median(bird_data$x, na.rm = T)
medlat <- median(bird_data$y, na.rm = T)
proj.aeqd <- paste("+proj=aeqd +lat_0=",round(medlat), " +lon_0=",round(medlon)," +units=m ", sep="")

bird_data.sp <- bird_data #bird data in aeqs projection
coordinates(bird_data.sp) <- ~ x + y
proj4string(bird_data.sp) <- CRS(proj.latlon)

bird_data.sp <- spTransform(bird_data.sp, CRS(proj.aeqd))
bird_data.proj <- as.data.table(bird_data.sp) 

###Functions ---

is.land <- function(x,y) #return TRUE if a point is over lands
{
  pt <- expand.grid(x, y)
  coordinates(pt) <- ~Var1 + Var2
  proj4string(pt) <- CRS(proj.aeqd)
  pt <- spTransform(pt, CRS(proj.latlon))
  pt <- st_as_sf(pt)
  !is.na(as.numeric(suppressMessages(st_intersects(pt, world))))
}

###VAR --- 

#parameters of the simulation
n <- 5 #number of iterations to produce a track --> track with 2^n segments
N <- 5 #number of tracks created

###Main script ---
id <- unique(bird_data.proj$ring)[30] #temporary : we focus on one individual for now
traj <- data.table(ring = c(),track_type = c(), x = c(), y = c())

#first trajectory = observed track

traj_id <- bird_data.proj[bird_data.proj$ring==id]
start_pt <- traj_id[,c("x","y")][1] 
end_pt <- traj_id[,c("x","y")][nrow(traj_id)]
traj_0 <- data.table(ring=rep(id,nrow(traj_id)), track_type = rep("observed",nrow(traj_id)),
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
        if (nland == 10){rm <- append(rm, paste0("sim-",k))}
      }
      traj_x_ <- append(traj_x_, c(xi,traj_x[i+1]))
      traj_y_ <- append(traj_y_, c(yi,traj_y[i+1]))
    }
    traj_x <- traj_x_
    traj_y <- traj_y_
  }
  traj_k <- data.table(ring=rep(id,length(traj_x)),track_type = rep(paste0("sim-", k),length(traj_x)), x = traj_x, y = traj_y)
  traj <- rbind(traj,traj_k)
}

traj <- traj[not(traj$track_type%in%rm)]#remove tracks that go above the land

#into geographic coord


coordinates(traj) <- ~ x + y
proj4string(traj) <- CRS(proj.aeqd) # gives the CRS of the data

traj <- spTransform(traj, CRS(proj.latlon)) # projects the data
traj <- as.data.table(traj) # back to a dataframe

#last traj = great circle line
start_pt_latlon <- traj[traj$ring==id&traj$track_type=="observed"][1]
end_pt_latlon <-  traj[traj$ring==id&traj$track_type=="observed"][nrow(traj[traj$ring==id&traj$track_type=="observed"])]
traj_gc <- data.table(ring = rep(id,22),track_type=rep("gc",22))
gctraj <- gcIntermediate(c(start_pt_latlon$x,start_pt_latlon$y),c(end_pt_latlon$x,end_pt_latlon$y), n = 20, addStartEnd = TRUE, sp = TRUE)
traj_gc[, x :=geom(gctraj)[,"x"]]
traj_gc[, y :=geom(gctraj)[,"y"]]
traj <- rbind(traj,traj_gc)

###Display --- 
traj[, x2 := data.table::shift(x, type = 'lead'), by = track_type]
traj[, y2 := data.table::shift(y, type = 'lead'), by = track_type]

plot_ <- ggplot(data=traj[track_type!="observed"]) + 
  geom_segment(size = .5, aes(x = x, y = y,xend=x2,yend=y2,color=track_type)) +
  geom_segment(data = traj[track_type=="observed"], size = 1, aes(x = x, y = y,xend=x2,yend=y2)) +
  geom_sf(data=world,fill = "black", color = "black") + 
  geom_point(aes(x=start_pt$x, y =start_pt$y), color = 'blue', size = 2.5) +
  geom_point(aes(x=end_pt$x, y =end_pt$y), color = 'red', size = 2.5) +
  coord_sf(xlim = c(-60, 60), ylim = c(30,85), expand = FALSE)+
  ggtitle(paste0("n = ",n, " ; a = ",a))
print(plot_)

#saveRDS(traj, file = "outputs/random_tracks.rds")
