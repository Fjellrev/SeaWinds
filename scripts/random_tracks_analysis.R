##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## SeaWinds - Data processing & Analyses
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sapply(c('sf','spData','raster','windR','fields','lubridate', 'shape', 'data.table',"rnaturalearthdata","rnaturalearth", 'magrittr', "maptools", 'rWind','rworldmap', 'gdistance','geosphere', 'ggplot2'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE)))

proj.latlon <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #classic longlat projection to use the windR function

cosd <-function(x) cos(x*pi/180) #cos of an angle in degrees 
sind <-function(x) sin(x*pi/180) 
###Get wind data ---

wind_path <- "data/ERA-Interrim/Interanual_means" 
wind_filename <-  "ERA_Interim_interanual_monthly_mean_sfc_10_2013_to_2018.rds"

wind_data <- readRDS(paste0(wind_path,"/", wind_filename))%>% as.data.table

wind_data[, wdir := atan2(u, v)*180/pi, by = 1:nrow(wind_data)] #get wind direction
wind_data[, wdir := wdir+360*(wdir<0), by = 1:nrow(wind_data)] #correct wind direction to be between 0 and 360° 
wind_data[, wspeed     := sqrt(u^2 + v^2), by = 1:nrow(wind_data)] #get wind speed

###Transform wind data into a raster layer to calculate conductance --- 

r <- raster(xmn=-70.25, xmx=70.25, ymn=29.75, ymx=85.25, res=.75) #Empty raster of the studied area

r_wdir <- rasterize(wind_data[,c("x","y")], r, field=wind_data$wdir) #raster of wind direction
r_wspeed <- rasterize(wind_data[,c("x","y")], r, field=wind_data$wspeed) #raster of wind speed

wind_layer <- stack(r_wdir, r_wspeed) #raster layer
names(wind_layer) <- c("direction", "speed")
Conductance <- flow.dispersion(x=wind_layer, output="transitionLayer")

###Get random tracks ---

track_path <- "outputs"
track_filename <- "random_tracks.rds"

tracks <- readRDS(paste0(track_path,"/", track_filename))%>% as.data.table
tracks[, x2 := data.table::shift(x, type = 'lead'), by = N]
tracks[, y2 := data.table::shift(y, type = 'lead'), by = N]
gspeed = c()
gdir=c()
for (i in (1:nrow(tracks)))
{
  gspeed = append(gspeed,distGeo(c(tracks$x[i], tracks$y[i]),c(tracks$x2[i], tracks$y2[i])))
  gdir = append(gdir, geosphere::bearing(c(tracks$x[i], tracks$y[i]),c(tracks$x2[i], tracks$y2[i])))
  
}
tracks[, gspeed := gspeed]
tracks[, gdir := gdir]

u_wind = c()
v_wind = c()
for (i in 1:nrow(tracks)) #get the wind data at the given coordinates and dates
{
  cat(paste0(i," "))
  uv_wind = getWind(tracks$x[i], tracks$y[i], w = wind_data, PROJ = proj.latlon)
  u_wind = append(u_wind, uv_wind[1])
  v_wind = append(v_wind, uv_wind[2])
}
tracks[, u_wind := as.numeric(u_wind)]
tracks[, v_wind := as.numeric(v_wind)]

tracks[, wdir := atan2(u_wind, v_wind)*180/pi, by = 1:nrow(tracks)] #wind direction
tracks[, wspeed     := sqrt(u_wind^2 + v_wind^2), by = 1:nrow(tracks)] #wind speed
tracks[, ws := wspeed*cosd(wdir -gdir)] #wind support

cost <- c()
for (i in 1:nrow(tracks))
{
  cat(paste0(i," "))
  if (is.na(tracks$x2[i])|tracks$y[i]>85|tracks$y2[i]>85)
  {cost <- append(cost, NA)}
  else
  {cost <- append(cost, costDistance(Conductance, c(tracks$x[i],tracks$y[i]), c(tracks$x2[i],tracks$y2[i])))} 
}
tracks$cost <- cost
tracks[, cost_tot := sum(tracks$cost[tracks$N==N], na.rm= T), by = 1:nrow(tracks)]
cost_95 <- sort(tracks$cost_tot)[5/100*nrow(tracks)]

saveRDS(tracks, file = "outputs/random_tracks_analized.rds")

world_map <- ne_countries(scale = "medium", returnclass = "sf")
plot_ <- ggplot(data=tracks[cost_tot<cost_95]) + 
  geom_path(size = 1, aes(x = x, y = y,color=N)) +
  geom_sf(data=world_map,fill = "black", color = "black") +
  coord_sf(xlim = c(-100, 100), ylim = c(0,85), expand = FALSE)+
  scale_color_distiller(palette = "Set1")
print(plot_)

