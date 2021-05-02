##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## SeaWinds - Data processing & Analyses
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sapply(c( 'raster','stringr', 'data.table', 'magrittr', 'sf', "rnaturalearthdata","rnaturalearth", 'windR', 'geosphere', 'ggplot2', 'circular'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE)))

proj.latlon <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #classic longlat projection to use the windR function

path_bird= "data/Kittiwake_data_treated"
path_wind = "data/ERA-Interrim"

cosd <-function(x) cos(x*pi/180) #cos of an angle in degrees 
sind <-function(x) sin(x*pi/180) 
mean_angle <-function(x) #calculates mean angle by decomposing it in u and v components
{
  u <- cosd(x)
  v <- sind(x)
  mean_u <- mean(u, na.rm = T)
  mean_v <- mean(v, na.rm = T)
  atan2(mean_v,mean_u)*180/pi
}

adir <- function(gspeed, gdir, wspeed, wdir) #get air direction
{
  u = wspeed * sind(wdir)
  v= wspeed * cosd(wdir)
  
  x_g = gspeed * sind(gdir)
  y_g = gspeed * cosd(gdir)
  
  x_a = x_g-u
  y_a = y_g-v
  
  atan2(x_a,y_a)*180/pi
}

###GET DATA ----
month = "2016-11" #month studied in this example
birdRDS = list.files(path=path_bird, pattern = "5col")
bird_data = readRDS(paste0(path_bird,'/', birdRDS)) %>% as.data.table
bird_data[, timestamp := as.POSIXct(bird_data$timestamp)] 
bird_data_ = bird_data#[as.logical(str_count(bird_data$timestamp, pattern = month))] #reduction of the dataset to the month and segment studied
setnames(bird_data_, c("lon","lat"),c("x", "y"))

windRDS = list.files(path=path_wind, pattern= 'ERA_Interim_daily_sfc_09_2016_to_12_2016.RDS')
wind_data = readRDS(paste0(path_wind,"/", windRDS)) %>% as.data.table
wind_data_ = wind_data#[as.logical(str_count(wind_data$datetime_, pattern = month))]  #reduction of the dataset to the month and segment studied

###GET MIGRATORY SEGMENTS ---- based on the residency time
bird_data_[, migrx := data.table::shift(x, n = 6L, type = 'lead'), by = ring] #get coordinates 3 days later
bird_data_[, migry := data.table::shift(y, n = 6L, type = 'lead'), by = ring]

dist_migr=c()
for (i in (1:nrow(bird_data_)))
{
dist_migr = append(dist_migr,distGeo(c(bird_data_$x[i], bird_data_$y[i]),c(bird_data_$migrx[i], bird_data_$migry[i])))
}
bird_data_[, migr_rt := as.numeric(dist_migr>3e+05)] #migratory segments = move more than 300km in 3 days
  
###GET BIRD SPEED AND DIRECTION ----
bird_data_[, timestamp2 := data.table::shift(timestamp, type = 'lead'), by = ring]
bird_data_[, x2 := data.table::shift(x, type = 'lead'), by = ring]
bird_data_[, y2 := data.table::shift(y, type = 'lead'), by = ring]
bird_data_[, dt := as.numeric(difftime(timestamp2, timestamp, units = 'sec'))]

gspeed = c()
gdir=c()
wdir=c()
for (i in (1:nrow(bird_data_)))
{
  gspeed = append(gspeed,distGeo(c(bird_data_$x[i], bird_data_$y[i]),c(bird_data_$x2[i], bird_data_$y2[i]))/bird_data_$dt[i])
  gdir = append(gdir, geosphere::bearing(c(bird_data_$x[i], bird_data_$y[i]),c(bird_data_$x2[i], bird_data_$y2[i])))
  
}
bird_data_[, gspeed := gspeed/(1-std_conductivity)]
bird_data_[, gdir := gdir]

#GET WIND ----
w_datetime=unique(wind_data_$datetime_) #dates with available wind data
bird_data_[, wind_time := closestDatetime(timestamp, w_datetime), by = 1:nrow(bird_data_)] #closest date with available wind data for each bird observation
u_wind = c()
v_wind = c()
for (i in 1:nrow(bird_data_)) #get the wind data at the given coordinates and dates
{
  cat(paste0(i," "))
  uv_wind = getWind(bird_data_$x[i], bird_data_$y[i], w = wind_data_[wind_data_$datetime_ == bird_data_$wind_time[i]], PROJ = proj.latlon)
  u_wind = append(u_wind, uv_wind[1])
  v_wind = append(v_wind, uv_wind[2])
}
bird_data_[, u_wind := as.numeric(u_wind)]
bird_data_[, v_wind := as.numeric(v_wind)]

#WIND SUPPORT, CROSS WIND, AIR SPEED ----

bird_data_[, wdir := atan2(u_wind, v_wind)*180/pi, by = 1:nrow(bird_data_)] #wind direction
bird_data_[, wspeed     := sqrt(u_wind^2 + v_wind^2), by = 1:nrow(bird_data_)] #wind speed
bird_data_[, ws := wspeed*cosd(wdir -gdir)] #wind support
bird_data_[, cw := bird_data_$wspeed*sind(bird_data_$wdir-bird_data_$gdir)] # Crosswind > 0 towards the right
bird_data_[,adir := adir(bird_data_$gspeed, bird_data_$gdir, bird_data_$wspeed, bird_data_$wdir)] #air speed
bird_data_[,aspeed := sqrt(bird_data_$gspeed^2+bird_data_$wspeed^2-2*bird_data_$gspeed*bird_data_$wspeed*cosd(bird_data_$wdir-bird_data_$gdir))]

###SAVING ---- 
saveRDS(bird_data_, file = "outputs/daily_tracks_analysed.rds")
