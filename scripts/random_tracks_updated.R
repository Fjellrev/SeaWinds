##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## SeaWinds - Data processing & Analyses
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sapply(c('sf','spData','tidyverse', 'data.table', 'magrittr', 'gdistance','geosphere', 'ggplot2'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE)))

source("functions/FUNCTION_OverlapPoly.r")

###Projections --- 
proj.aeqd <- paste("+proj=aeqd +lat_0=",round(medlat), " +lon_0=",round(medlon)," +units=m ", sep="")
proj.latlon <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #classic longlat projection to use the windR function

data("world")
wrld <- st_transform(world,CRS(proj.aeqd))

path_bird <- "data/Kittiwake_data_treated"


###Read BLKI data ---

birdRDS <- list.files(path=path_bird, pattern = "tracks")
bird_data <- readRDS(paste0(path_bird,'/', birdRDS)) %>% as.data.table


medlon <- median(bird_data$x, na.rm = T)
medlat <- median(bird_data$y, na.rm = T)

bird_data.sp <- st_as_sf(bird_data, coords=4:5, crs=CRS(proj.latlon))
bird_data.sp <- st_transform(bird_data.sp, CRS(proj.aeqd))
bird_data.proj <- as.data.table(bird_data.sp) 
bird_data.proj <- bird_data.proj[,-c("geometry")]
bird_data.proj[,x := unlist(map(bird_data.sp$geometry,1))]
bird_data.proj[,y := unlist(map(bird_data.sp$geometry,2))]


###VAR --- 

#parameters of the simulation
n <- 3 #number of iterations to produce a track --> track with 2^n segments
N <- 100 #number of tracks created

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
n.pts <- 50 # number of random points generated during computation of new location

for (k in 1:N)
{
  cat(k, "\n")
  traj_x <- c(start_pt$x,end_pt$x) #x coords of the track
  traj_y <- c(start_pt$y,end_pt$y)
  
  for (j in 1:n) 
  {
    cat("j = ", j)
    
    traj_x_ <- c(traj_x[1]) #new x coords 
    traj_y_ <- c(traj_y[1])
    
    for (i in 1:(length(traj_x)-1)) #create a new point betwen point i and i+1 
    {
      
      cat("i = ", i, "\n")
      
      v <- c(traj_x[i+1]-traj_x[i], traj_y[i+1]-traj_y[i]) 
      vn <- c(-v[2], v[1])/sqrt(v[1]^2+v[2]^2) #orthogonal vector to the segment treated 
      r <- runif(n.pts,-1,1)
      xi <- (traj_x[i]+traj_x[i+1])/2+a/(2^(j-1))*r*vn[1] #coords of the new points
      yi <- (traj_y[i]+traj_y[i+1])/2+a/(2^(j-1))*r*vn[2]
      # nland=0 #number of times the point has been created above the lands (in case there is no choice)
      
      land_out <- is.land(x = xi, y = yi, prj = proj.aeqd, mask = wrld)
      
      while(length(which(land_out == FALSE)) < 5) # check that enough points are generated above ocean
      {
        r <- runif(n.pts,-1,1)
        xi <- (traj_x[i]+traj_x[i+1])/2+a/(2^(j-1))*r*vn[1]
        yi <- (traj_y[i]+traj_y[i+1])/2+a/(2^(j-1))*r*vn[2]
        land_out <- is.land(x = xi, y = yi, prj = proj.aeqd, mask = wrld)
        # if (nland == 10){rm <- append(rm, paste0("sim-",k))}
      }
      
      xi <- sample(xi[which(land_out == FALSE)], 1)
      yi <- sample(yi[which(land_out == FALSE)], 1)
      
      traj_x_ <- append(traj_x_, c(xi,traj_x[i+1]))
      traj_y_ <- append(traj_y_, c(yi,traj_y[i+1]))
    }
    traj_x <- traj_x_
    traj_y <- traj_y_
  }
  
  traj_k <- data.table(ring=rep(id,length(traj_x)),track_type = rep(paste0("sim-", k),length(traj_x)), x = traj_x, y = traj_y)
  traj <- rbind(traj,traj_k)
}

# traj <- traj[not(traj$track_type%in%rm)]#remove tracks that go above the land

#into geographic coord

traj.sp <- st_as_sf(traj, coords=3:4, crs=CRS(proj.aeqd))
traj.sp <- st_transform(traj.sp, CRS(proj.latlon))
traj <- as.data.table(traj.sp) 
traj <- traj[,-c("geometry")]
traj[,x := unlist(map(traj.sp$geometry,1))]
traj[,y := unlist(map(traj.sp$geometry,2))]

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
  ggtitle(paste0("n = ",n, " ; ring = ",id))
print(plot_)

#saveRDS(traj, file = "outputs/random_tracks.rds")

test.map <- function(xi, yi){
  ggplot() +  geom_sf(data=wrld, fill = "black", color = "black") + 
    geom_point(aes(x=start_pt$x, y =start_pt$y), color = 'blue', size = 2.5) + 
    geom_point(data = end_pt, aes(x, y), colour = 'red', size = 2.5) + 
    geom_point(aes(x=xi, y =yi), color = 'green', size = 2.5)
}
