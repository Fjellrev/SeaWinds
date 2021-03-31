##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## SeaWinds - Data processing & Analyses
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sapply(c('raster','fields','lubridate', 'shape', 'data.table', 'magrittr', "maptools", 'rWind','rworldmap', 'gdistance', 'ggplot2'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE)))

date_departure <- as.POSIXct("2016-10-06") #Date of the start of the migration
daily_speed <- 200000 #number of metters travelled a day 
#coord of the colony
col_x <- 8
col_y <- 78

#coord of the wintering range
wint_x <- -46
wint_y <- 54

###Get wind data ---
wind_path <- "data/ERA-Interrim" 
wind_filename <-  "ERA_Interim_daily_sfc_09_2016_to_12_2016.rds"
wind_data <- readRDS(paste0(wind_path,"/", wind_filename))%>% as.data.table

###Get a conductance layer of the emerged lands ---

data(wrld_simpl)

r <- raster(xmn=-70.25, xmx=70.25, ymn=29.75, ymx=85.25, res=.75) #Empty raster of the studied area
r_land <- rasterize(wrld_simpl, r, field = 0) #raster with 0 above lands and NA above the sea
r_df <- r_land %>% as.data.frame(xy=T)
r_df[is.na(r_df)] <- 1 #turn NA into 1 above the sea
r_land <- rasterize(r_df[,c("x","y")], r, field = r_df$layer)
Conductance_land <- transition(r_land,transitionFunction=min, directions = 4) #turn this raster into a transition layer


###Main script ---
date <- date_departure
pos <- c(col_x,col_y)
traj <- data.table(x=pos[1],y=pos[2],timestamp=date)
while(distGeo(pos,c(wint_x,wint_y))>daily_speed)
{
  cat(paste0(date, " : ", distGeo(pos,c(wint_x,wint_y))))
  wind_data_<- wind_data[wind_data$datetime_ == format(date)] #select only the correct datetime
  
  wind_data_[, wdir := atan2(u, v)*180/pi, by = 1:nrow(wind_data_)] #get wind direction
  wind_data_[, wdir := wdir+360*(wdir<0), by = 1:nrow(wind_data_)] #correct wind direction to be between 0 and 360° 
  wind_data_[, wspeed     := sqrt(u^2 + v^2), by = 1:nrow(wind_data_)] #get wind speed

  r_wdir <- rasterize(wind_data_[,c("x","y")], r, field=wind_data_$wdir) #raster of wind direction
  r_wspeed <- rasterize(wind_data_[,c("x","y")], r, field=wind_data_$wspeed) #raster of wind speed
  
  wind_layer <- stack(r_wdir, r_wspeed) #raster layer
  names(wind_layer) <- c("direction", "speed")
  
  Conductance <- flow.dispersion(x=wind_layer, output="transitionLayer")
  Conductance <- Conductance*Conductance_land
  daily_sp <- shortestPath(Conductance,pos, c(wint_x, wint_y),output="SpatialLines")
  coords <- geom(daily_sp)[,c("x","y")]
  daily_sp_df <- data.table(x=coords[,"x"],y=coords[,"y"])
  daily_sp_df[, dist := distGeo(c(x,y),c(daily_sp_df$x[1],daily_sp_df$y[1])), by = 1:nrow(daily_sp_df)]
  n <- which(daily_sp_df$dist>daily_speed)[1]
  pos <- c(daily_sp_df$x[n],daily_sp_df$y[n])
  traj <- rbind(traj,data.table(x=daily_sp_df$x[1:n],y=daily_sp_df$y[1:n],timestamp=rep(date,n)))
  date <- date + days(1)
}

traj_line <-Line(cbind(traj[,c("x")], traj[,c("y")]))

ggplot(data=traj) + 
  geom_line(size = 1, aes(x = x, y = y), color = 'blue') +
  geom_sf(data=world ,fill = "black", color = "black") + 
  geom_point(aes(x=col_x, y =col_y), color = 'blue', size = 2.5) +
  geom_point(aes(x=wint_x, y = wint_y), color = 'red', size = 2.5) +
  coord_sf(xlim = c(min(bird_data_$x)-1, max(bird_data_$x)+1), ylim = c(min(bird_data_$y)-1, max(bird_data_$y)+1), expand = FALSE)
