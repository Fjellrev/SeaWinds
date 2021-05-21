##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## SeaWinds - Data processing & Analyses
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sapply(c('sf','sp','spData','tidyverse', 'data.table','MASS','plyr', 'magrittr', 'gdistance','geosphere','raster', "doParallel", "foreach",
         'ggplot2', 'rWind','windR','adehabitatHR'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE)))

### utility functions ----
source("functions/FUNCTION_OverlapPoly.r")
source('functions/FUNCTION_QuickMap.r')
cosd <-function(x) cos(x*pi/180) #cos of an angle in degrees 
sind <-function(x) sin(x*pi/180) 
proj.latlon <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #classic longlat projection to use the windR function

bird_path <- "data/Kittiwake_data_treated"
bird_filename <- "BLKI_tracks.rds"

map_path <- "data/baseline_data"
map_filename <- "SEAwinds_Worldmap_res-c.rds"

daily_wind_path <- "data/ERA_Interim"
daily_wind_filename <- "ERA_Interim_daily_sfc_09_2016_to_06_2017.rds"
monthly_wind_path <-"data/ERA_Interim/Interannual_means" 

### Read BLKI data ----
bird_data <- readRDS(paste0(bird_path,'/', bird_filename)) %>% as.data.table
bird_data[, x2 := data.table::shift(x, type = 'lead'), by = ring]
bird_data[, y2 := data.table::shift(y, type = 'lead'), by = ring]
bird_data[, gdir := geosphere::bearing(c(x,y),c(x2,y2)), by = 1:nrow(bird_data)]

### GET REAL WIND ----
wind_data <- readRDS(paste0(daily_wind_path, '/', daily_wind_filename)) %>% as.data.table
w_datetime <- unique(wind_data$datetime_) #dates with available wind data
bird_data[, wind_time := closestDatetime(timestamp, w_datetime), by = 1:nrow(bird_data)] #closest date with available wind data for each bird observation
u_wind <- c()
v_wind <- c()
for (i in 1:nrow(bird_data)) #get the wind data at the given coordinates and dates
{
  cat(paste0(i," "))
  uv_wind <- getWind(bird_data$x[i], bird_data$y[i], w = wind_data[wind_data$datetime_ == bird_data$wind_time[i]], PROJ = proj.latlon)
  u_wind <- append(u_wind, uv_wind[1])
  v_wind <- append(v_wind, uv_wind[2])
}
bird_data[, u_wind := as.numeric(u_wind)]
bird_data[, v_wind := as.numeric(v_wind)]

bird_data[, wdir := atan2(u_wind, v_wind)*180/pi, by = 1:nrow(bird_data)] #wind direction
bird_data[, ws := cosd(wdir - gdir)] #wind support
### GET MEAN MONTHLY WINDs ----
bird_data[, u_monthly_wind := 0]
bird_data[, v_monthly_wind := 0]
months <- unique(bird_data$migr_month)
for (month in months){
  monthly_wind_filename <-  paste0("ERA_Interim_interannual_monthly_mean_sfc_",month,"_2013_to_2018.rds")
  wind_data <- readRDS(paste0(monthly_wind_path,"/", monthly_wind_filename))%>% as.data.table
  
  u_monthly_wind <- c()
  v_monthly_wind <- c()
  for (i in 1:nrow(bird_data[bird_data$migr_month==month])) #get the wind data at the given coordinates and dates
  {
    cat(paste0(i," "))
    uv_wind <- getWind(bird_data$x[i], bird_data$y[i], w = wind_data, PROJ = proj.latlon)
    u_monthly_wind <- append(u_monthly_wind, uv_wind[1])
    v_monthly_wind <- append(v_monthly_wind, uv_wind[2])
  }
  bird_data$u_monthly_wind[bird_data$migr_month==month] <- as.numeric(u_monthly_wind)
  bird_data$v_monthly_wind[bird_data$migr_month==month] <- as.numeric(v_monthly_wind)
  
}
bird_data[, monthly_wdir := atan2(u_monthly_wind, v_monthly_wind)*180/pi, by = 1:nrow(bird_data)] #wind direction
bird_data[, monthly_ws := cosd(monthly_wdir - gdir)] #wind support

