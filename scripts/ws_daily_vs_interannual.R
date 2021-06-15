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

raster_path <- "outputs/Rasters"

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
bird_data[, wspeed := sqrt(u_wind^2+v_wind^2)]
bird_data[, ws := wspeed*cosd(wdir - gdir)] #wind support

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
bird_data[, monthly_wspeed := sqrt(u_monthly_wind^2+v_monthly_wind^2)]
bird_data[, monthly_ws := monthly_wspeed*cosd(monthly_wdir - gdir)] #wind support

bird_data[, ws_by_burst := mean(bird_data$ws[bird_data$burst==burst], na.rm=T), by= 1:nrow(bird_data)]
bird_data[, ws_monthly_by_burst := mean(bird_data$monthly_ws[bird_data$burst==burst], na.rm=T), by= 1:nrow(bird_data)]

bird_data$season <- "spring"
bird_data$season[as.numeric(bird_data$migr_month)>6] <- "autumn"

r.var_spring <- raster(paste0(raster_path,"/wdir_sd_spring"))
r.var_autumn <- raster(paste0(raster_path,"/wdir_sd_spring"))

bird_data$wind_var <- NA
bird_data$wind_var[bird_data$season =="spring"] <- extract(r.var_spring,
                                                 cellFromXY(r.var_spring,bird_data[season =="spring"][,c("x","y")]))
bird_data$wind_var[bird_data$season =="autumn"] <- extract(r.var_autumn,
                                                           cellFromXY(r.var_autumn,bird_data[season =="autumn"][,c("x","y")]))
bird_data[, wind_var_track := mean(bird_data$wind_var[bird_data$burst==burst], na.rm=T), by= 1:nrow(bird_data)]
bird_data[, d_ws := mean(bird_data$ws[bird_data$burst==burst]-bird_data$monthly_ws[bird_data$burst==burst], na.rm=T), by= 1:nrow(bird_data)]

#lmm.dws_windvar <- lmer(ws-monthly_ws ~ wind_var + (1|colony/ring) + (1|migr_month), data = bird_data)

#stat_ws_d <- data.frame(ring=bird_data$ring,colony=bird_data$colony, season = bird_data$season,
                      #wind_dataset = "daily", ws = bird_data$ws)

#stat_ws_i <- data.frame(ring=bird_data$ring,colony=bird_data$colony, season = bird_data$season,
                        #wind_dataset = "interanual", ws = bird_data$monthly_ws)

#stat <- rbind(stat_ws_d,stat_ws_i)
#stat$season <- "spring"
#stat$season[as.numeric(stat$month)>6] <- "autumn"


burst_data <- data.table(burst = unique(bird_data$burst))
burst_data[,colony := bird_data$colony[bird_data$burst==burst][1], by = 1:nrow(burst_data)]
burst_data[,season := bird_data$season[bird_data$burst==burst][1], by = 1:nrow(burst_data)]
burst_data[,dataset := "Daily wind data"]
burst_data[,mean_ws := mean(bird_data$ws[bird_data$burst == burst], na.rm=T), by = 1:nrow(burst_data)]

burst_data_monthly <- burst_data[,-c("mean_ws")]
burst_data_monthly$dataset <- "Average monthly wind data"
burst_data_monthly[,mean_ws := mean(bird_data$monthly_ws[bird_data$burst == burst], na.rm=T),
                   by = 1:nrow(burst_data)]
burst_data <- rbind(burst_data, burst_data_monthly)

ggpaired(burst_data, x= "dataset", y="mean_ws", 
         color="dataset",line.color = "gray", line.size = 0.4,palette = "npg", id = "burst")+
  facet_wrap(.~season)+xlab("")+ylab(bquote("Wind Support (" ~m.s^-1~")"))+
  theme(legend.position = "none",
        strip.text = element_text(face="bold"))


