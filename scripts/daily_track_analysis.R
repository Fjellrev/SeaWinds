##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## SeaWinds - Data processing & Analyses
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sapply(c( 'raster','stringr', 'data.table', 'magrittr', 'sf', "rnaturalearthdata","rnaturalearth", 'windR', 'geosphere', 'ggplot2', 'circular'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE)))

proj.latlon <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #classic longlat projection to use the windR function

cosd <-function(x) cos(x*pi/180) #degrees 
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
migr = c(0,1)     # type of segment studied (0 = stationnary, 1 = migratory)
birdRDS = list.files(path="data/Kittiwake_data", pattern = "Alkefjellet_cond_nov2016")
bird_data = readRDS(paste0("data/Kittiwake_data/", birdRDS)) %>% as.data.table
bird_data[, timestamp := as.POSIXct(bird_data$timestamp)] 
bird_data_ = bird_data[as.logical(str_count(bird_data$timestamp, pattern = month))&bird_data$migr10 %in% migr] #reduction of the dataset to the month and segment studied
bird_data_ = bird_data_[abs(as.numeric(format(bird_data_$timestamp, "%H"))-12)<6] #keep one daily position

windRDS = list.files(path="data/ASCAT", pattern= 'ASCAT_daily_201611_wind.RDS')
wind_data = readRDS(paste0("data/ASCAT/", windRDS)) %>% as.data.table
wind_data_ = wind_data[as.logical(str_count(wind_data$datetime_, pattern = month))]  #reduction of the dataset to the month and segment studied

#GET BIRD SPEED AND DIRECTION ----
setnames(bird_data_, c("lon","lat"),c("x", "y"))
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

###ANALYSIS OF INDIVIDUALS TRACKS ----
ws_track <- c() #mean wind support on the individual tracks
for (i in unique(bird_data_$ring))
{
  ws_track <- append(ws_track, sum(bird_data_$ws[bird_data_$ring==i], na.rm = T))
}
track <- data.table(ring = unique(bird_data_$ring), ws = ws_track)

###PLOTS###
ggplot(data=bird_data_)+
  geom_histogram(aes(x=ws, color = migr10==1), fill = 'white')+
  geom_vline(data = bird_data_[migr10==1],aes(xintercept=mean(ws, na.rm = T)),
             linetype="dashed", color = 'blue')+
  geom_vline(data = bird_data_[migr10==0],aes(xintercept=mean(ws, na.rm = T)),
             linetype="dashed", color = 'red')

ggplot(data=track)+
  geom_histogram(aes(x=ws))


###DISPLAY###
world <- ne_countries(scale = "medium", returnclass = "sf")
a=4 #Arrow size
ggplot(data = bird_data_[bird_data_$std_conductivity<0.75]) + 
  geom_segment(aes(x = x, y = y, xend = x2, yend = y2, color=ws),arrow = arrow(length = unit(.1, "cm"))) +
  scale_color_gradient("wind support", low = "red", high = "green") +
  geom_sf(data=world ,fill = "black", color = "black") + 
  coord_sf(xlim = c(min(bird_data_$x)-1, max(bird_data_$x)+1), ylim = c(min(bird_data_$y)-1, max(bird_data_$y)+1), expand = FALSE)


i = 110

print(bird_data_[which(not(is.na(bird_data_$gdir)))[i]])

ggplot(data=bird_data_[which(not(is.na(bird_data_$gdir)))[i]]) + 
  geom_segment(aes(x = 0, y = 0, xend = gspeed*sind(gdir), yend = gspeed*cosd(gdir)), color = 'red', size = 1, arrow = arrow(length = unit(1, "cm"))) + #Ground 
  geom_segment(aes(x = 0, y = 0, xend = wspeed*sind(wdir), yend = wspeed*cosd(wdir)), color = 'blue', size = 1, arrow = arrow(length = unit(1, "cm"))) + #Wind 
  geom_segment(aes(x = 0, y = 0, xend = ws*sind(gdir), yend = ws*cosd(gdir)), color = 'black', size = 1, arrow = arrow(length = unit(1, "cm"))) + #wind support 
  geom_segment(aes(x = 0, y = 0, xend = cw*sind(gdir+90), yend = cw*cosd(gdir+90)), color = 'black', size = 1, arrow = arrow(length = unit(1, "cm"))) + #crosswind
  geom_segment(aes(x=0, y =0, xend =  aspeed*sind(adir), yend = aspeed*cosd(adir)), color = 'green', size = 1, arrow = arrow(length = unit(1, "cm")))+ #air
  coord_fixed()
