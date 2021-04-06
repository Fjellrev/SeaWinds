##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## SeaWinds - Data processing & Analyses
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sapply(c('raster','fields','lubridate', 'shape', 'data.table', 'magrittr', "maptools", 'rWind','rworldmap', 'gdistance'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE)))

date <- "2016-10-24" #Date of the wind layer used

#coord of the colony
col_x <- 8
col_y <- 68

#coord of the wintering range
wint_x <- -37
wint_y <- 52

###Get wind data ---

wind_path <- "data/ERA-Interrim" 
wind_filename <-  "ERA_Interim_daily_sfc_09_2016_to_12_2016.rds"

wind_data <- readRDS(paste0(wind_path,"/", wind_filename))%>% as.data.table
wind_data_<- wind_data[wind_data$datetime_ == date] #select only the correct datetime

wind_data_[, wdir := atan2(u, v)*180/pi, by = 1:nrow(wind_data_)] #get wind direction
wind_data_[, wdir := wdir+360*(wdir<0), by = 1:nrow(wind_data_)] #correct wind direction to be between 0 and 360° 
wind_data_[, wspeed     := sqrt(u^2 + v^2), by = 1:nrow(wind_data_)] #get wind speed

###Transform wind data into a raster layer to calculate conductance --- 

r <- raster(xmn=-70.25, xmx=70.25, ymn=29.75, ymx=85.25, res=.75) #Empty raster of the studied area

r_wdir <- rasterize(wind_data_[,c("x","y")], r, field=wind_data_$wdir) #raster of wind direction
r_wspeed <- rasterize(wind_data_[,c("x","y")], r, field=wind_data_$wspeed) #raster of wind speed

wind_layer <- stack(r_wdir, r_wspeed) #raster layer
names(wind_layer) <- c("direction", "speed")

###Get a conductance layer of the emerged lands ---

data(wrld_simpl)
r_land <- rasterize(wrld_simpl, r, field = 0) #raster with 0 above lands and NA above the sea
r_df <- r_land %>% as.data.frame(xy=T)
r_df[is.na(r_df)] <- 1 #turn NA into 1 above the sea
r_land <- rasterize(r_df[,c("x","y")], r, field = r_df$layer)
Conductance_land <- transition(r_land,transitionFunction=min, directions = 8) #turn this raster into a transition layer

###Get the shortest path with rWind ---

Conductance <- flow.dispersion(x=wind_layer, output="transitionLayer")
Conductance <- Conductance*Conductance_land
to_wintering_range <- shortestPath(Conductance,c(col_x,col_y), c(wint_x, wint_y),output="SpatialLines")
to_colony <- shortestPath(Conductance,c(wint_x, wint_y),c(col_x,col_y), output="SpatialLines")

image.plot(wind_layer$speed, col=terrain.colors(10), zlim=c(0,0.1))
points(col_x, col_y, pch=19, cex=1.4,col="red")
points(wint_x, wint_y, pch=19, cex=1.4,col="blue")
lines(to_wintering_range, col="red", lwd=4, lty=2)
lines(to_colony, col="blue", lwd=4, lty=2)
lines(getMap(resolution="low"), lwd=4)


