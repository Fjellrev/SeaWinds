##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## SeaWinds - Data processing & Analyses
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sapply(c('sf','spData','raster','fields','lubridate', 'shape', 'data.table',"rnaturalearthdata","rnaturalearth", 'magrittr', "maptools", 'rWind','rworldmap', 'gdistance','geosphere', 'ggplot2'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE)))

###Functions ---
data("world")
is.land <- function(x,y) #return TRUE if a point is over lands
{
  pt <- expand.grid(x, y)
  pt <- st_as_sf(pt, coords=1:2, crs = 4326)
  !is.na(as.numeric(suppressMessages(st_intersects(pt, world))))
  
}

###VAR --- 

#coord of the colony
col_x <- 19
col_y <- 80

#coord of the wintering range
wint_x <- -55
wint_y <- 60

#parameters of the simulation
n <- 3 #number of iterations to produce a track --> track with 2^n segments
a <- 15 #magnitude of the deviation
N <- 500 #number of tracks created

###Main script ---
traj <- data.table(N = c(), x = c(), y = c())

#firs traj = great circle line
traj_gc <- data.table(N=rep(0,22))
gctraj <- gcIntermediate(c(col_x,col_y),c(wint_x,wint_y), n = 20, addStartEnd = TRUE, sp = TRUE)
traj_gc[, x :=geom(gctraj)[,"x"]]
traj_gc[, y :=geom(gctraj)[,"y"]]
traj <- rbind(traj,traj_gc)

#N random tracks
rm <- c()
for (k in seq (1:N))
{
  cat(k)
  traj_x <- c(col_x,wint_x) #x coords of the track
  traj_y <- c(col_y,wint_y)
  
  for (j in seq(1:n)) 
  {
    traj_x_ <- c(traj_x[1]) #new x coords 
    traj_y_ <- c(traj_y[1])
    for (i in seq(1:(length(traj_x)-1))) #create a new point betwen point i and i+1 
    {
      v <- c(traj_x[i+1]-traj_x[i], traj_y[i+1]-traj_y[i]) 
      vn <- c(-v[2], v[1])/sqrt(v[1]^2+v[2]^2) #orthogonal vector to the segment treated 
      r <- runif(1,-1,1)
      xi <- (traj_x[i]+traj_x[i+1])/2+a/j*r*vn[1] #coords of the new point
      yi <- (traj_y[i]+traj_y[i+1])/2+a/j*r*vn[2]
      nland=0 #number of times the point has been created above the lands (in case there is no choice)
      while(is.land(xi,yi)&nland<10)#is it above the lands
      {
        nland=nland+1
        r <- runif(1,-1,1)
        xi <- (traj_x[i]+traj_x[i+1])/2+a/(j)*r*vn[1]
        yi <- (traj_y[i]+traj_y[i+1])/2+a/(j)*r*vn[2]
        if (nland == 10){rm <- append(rm, k)}
      }
      traj_x_ <- append(traj_x_, c(xi,traj_x[i+1]))
      traj_y_ <- append(traj_y_, c(yi,traj_y[i+1]))
    }
    traj_x <- traj_x_
    traj_y <- traj_y_
  }
  traj_k <- data.table(N = rep(k,length(traj_x)), x = traj_x, y = traj_y)
  traj <- rbind(traj,traj_k)
}

traj <- traj[not(traj$N%in%rm)]#remove tracks that go above the land

plot_ <- ggplot(data=traj[N]) + 
  geom_path(size = 1, aes(x = x, y = y,color=N)) +
  geom_sf(data=world,fill = "black", color = "black") + 
  geom_point(aes(x=col_x, y =col_y), color = 'blue', size = 2.5) +
  geom_point(aes(x=wint_x, y = wint_y), color = 'red', size = 2.5) +
  coord_sf(xlim = c(-60, 60), ylim = c(30,85), expand = FALSE)+
  scale_color_distiller(palette = "Set1")+ 
  ggtitle(paste0(N," simulations ; a = ",a))
print(plot_)

#saveRDS(traj, file = "outputs/random_tracks.rds")
