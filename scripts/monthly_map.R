
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## SeaWinds - Data processing & Analyses
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sapply(c( 'raster','stringr', 'data.table', 'magrittr', 'sf', "rnaturalearthdata","rnaturalearth", 'windR', 'geosphere', 'ggplot2', 'circular'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE)))

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
  birdRDS = list.files(path="data/Kittiwake_data", pattern = "Alkefjellet_nov2016")
  bird_data = readRDS(paste0("data/Kittiwake_data/", birdRDS)) %>% as.data.table
  bird_data_ = bird_data[as.logical(str_count(bird_data$timestamp, pattern = month))&bird_data$migr10 %in% migr] #data used
  
  windRDS = list.files(path="data/ASCAT", pattern= 'ASCAT_monthly_201309_201912_wind.RDS')
  wind_data = readRDS(paste0("data/ASCAT/", windRDS)) %>% as.data.table
  wind_data_ = wind_data[as.logical(str_count(wind_data$datetime_, pattern = month))] #data used

#GET SPEED AND DIRECTION ----
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
  for (i in (1:nrow(wind_data_)))
  {
    wdir = append(wdir, windR::bearing(0,0,wind_data_$u[i],wind_data_$v[i]))
  }
  bird_data_[, gspeed := gspeed]
  bird_data_[, gdir := gdir]
  wind_data_[,wspeed := sqrt(wind_data_$u^2+wind_data_$v^2)]
  wind_data_[,wdir := wdir*180/pi]
  
###GET MIGRATORY SEGMENTS ----
  #bird_data_[, migrx := data.table::shift(x, n = 3L, type = 'lead'), by = ring]
  #bird_data_[, migry := data.table::shift(y, n = 3L, type = 'lead'), by = ring]
  
  #dist_migr=c()
  #for (i in (1:nrow(bird_data_)))
  #{
    #dist_migr = append(dist_migr,distGeo(c(bird_data_$x[i], bird_data_$y[i]),c(bird_data_$migrx[i], bird_data_$migry[i])))
  #}
  #bird_data_[, migr10 := as.numeric(dist_migr>4e+05)]
  
###CONVERT TO RASTER ----
  r <- raster(xmn=(min(bird_data_$x)-1), ymn=(min(bird_data_$y)-1), xmx=(max(bird_data_$x)+1), ymx=(max(bird_data_$y)+1), res=2)
  r.gdir <- rasterize(bird_data_[,c("x","y")], r, field = bird_data_$gdir, fun = function(x, na.rm=T) mean_angle(x))
  r.gspeed <- rasterize(bird_data_[,c("x","y")], r, field = bird_data_$gspeed, fun = mean)
  r.n <- rasterize(bird_data_[,c("x","y")], r, field = bird_data_$x2, fun = 'count')
  r.wdir <- rasterize(wind_data_[,c("x","y")],r, field = wind_data_$wdir, fun = function(x, na.rm=T) mean_angle(x))
  r.wspeed <- rasterize(wind_data_[,c("x","y")],r, field = wind_data_$wspeed,fun = mean)

###BACK TO DATA FRAME (pxdata = data for each pixel) ----
  pxdata <- merge(as.data.frame(r.gdir, xy = T), as.data.frame(r.gspeed, xy = T), by = c("x", "y"))
  pxdata <- merge(pxdata, as.data.frame(r.n, xy = T), by = c("x","y"))
  pxdatawind <- merge(as.data.frame(r.wdir, xy = T), as.data.frame(r.wspeed, xy = 'T'), by = c("x","y"))
  pxdata <- merge(pxdata, pxdatawind, by = c("x","y"))
  setnames(pxdata, c("layer.x.x","layer.y.x","layer", "layer.x.y", "layer.y.y"), c("gdir","gspeed","n","wdir","wspeed"))
  pxdata <- pxdata %>% as.data.table

###DATA ANALYSIS ----
  pxdata[, adir := adir(pxdata$gspeed, pxdata$gdir, pxdata$wspeed, pxdata$wdir)]
  pxdata[, ws := pxdata$wspeed*cosd(pxdata$wdir-pxdata$gdir)]
  pxdata[, cw := pxdata$wspeed*sind(pxdata$wdir-pxdata$gdir)] # CW > 0 towards the right
  pxdata[,aspeed := sqrt((pxdata$gspeed-pxdata$ws)^2+pxdata$cw^2)]
  saveRDS(pxdata, file = "outputs/monthly_pixel_data.rds")

  ggplot(data = pxdata)+ #histogram of wind support 
    geom_histogram(aes(x=ws))+
    xlab("Wind support")
  
  ggplot(data = pxdata)+ #Air speed as a function of wind support
    geom_point(aes(x=ws, y = aspeed, size = gspeed))+
    geom_smooth(aes(x=ws, y = aspeed), method = "lm")+
    xlab("Wind support")+ylab("air speed")
  
  ggplot(data = pxdata)+ #Ground speed as a function of wind support
    geom_point(aes(x=ws, y = gspeed))+
    geom_smooth(aes(x=ws, y = gspeed), method = "lm")+
    xlab("Wind support")+ylab("ground speed")
  
  
###DISPLAY###
world <- ne_countries(scale = "medium", returnclass = "sf")
a=4 #Arrow size
plottraj =  ggplot(data = bird_data_) + 
  geom_segment(size = 1, aes(x = x, y = y, xend = x2, yend = y2, color=migr10),arrow = arrow(length = unit(.15, "cm"))) +
  geom_point(aes(x,y),size =.7) +
  geom_sf(data=world ,fill = "black", color = "black") + 
  coord_sf(xlim = c(min(bird_data_$x)-1, max(bird_data_$x)+1), ylim = c(min(bird_data_$y)-1, max(bird_data_$y)+1), expand = FALSE)

plotbird =  ggplot(data = pxdata[not(is.na(pxdata$n))]) + 
  geom_segment(size = 1, aes(x = x, y = y, xend = x+gspeed/a*sind(gdir), yend = y+gspeed/a*cosd(gdir), color=ws),arrow = arrow(length = unit(.15, "cm"))) +
  scale_color_gradient("wind support", low = "red", high = "green") +
  geom_sf(data=world ,fill = "black", color = "black") + 
  coord_sf(xlim = c(min(bird_data_$x)-1, max(bird_data_$x)+1), ylim = c(min(bird_data_$y)-1, max(bird_data_$y)+1), expand = FALSE)
plotwind = ggplot(data = pxdata) +
  geom_segment(size = 1, aes(x = x, y = y, xend = x+wspeed/a*sind(wdir), yend = y+wspeed/a*cosd(wdir)),arrow = arrow(length = unit(.15, "cm"))) +
  geom_sf(data=world ,fill = "black", color = "black") + 
  coord_sf(xlim = c(min(bird_data_$x)-1, max(bird_data_$x)+1), ylim = c(min(bird_data_$y)-1, max(bird_data_$y)+1), expand = FALSE)


