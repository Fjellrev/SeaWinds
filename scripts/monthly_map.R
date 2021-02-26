
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## SeaWinds - Data processing & Analyses
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sapply(c('rgdal', 'raster', 'data.table', 'magrittr', 'sp', 'rgeos', 'raster', 'foreach',
         'ncdf4', 'gdalUtilities',"rnaturalearthdata","rnaturalearth", 'windR', 'geosphere', 'ggplot2'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE)))
cosd <-function(x) cos(x*pi/180) #degrees 
sind <-function(x) sin(x*pi/180) # radians


###GET DATA ----
  birdRDS = list.files(path="Kittiwake_data", pattern = "Alkefjellet_nov2016")
  bird_data = readRDS(paste0("Kittiwake_data/", birdRDS)) %>% as.data.table
  bird_data_ = bird_data["2016-11" %in% bird_data$timestamp]
  
  windRDS = list.files(path="ASCAT", pattern= '.RDS')
  wind_data = readRDS(paste0("ASCAT/", windRDS)) %>% as.data.table
  wind_data_ = wind_data[wind_data$datetime_=="2016-11-16 CET" ] #dataset used in this example (11-2016)

#GET SPEED AND DIRECTION ----
  setnames(bird_data, c("lon","lat"),c("x", "y"))
  bird_data[, timestamp2 := data.table::shift(timestamp, type = 'lead'), by = ring]
  bird_data[, x2 := data.table::shift(x, type = 'lead'), by = ring]
  bird_data[, y2 := data.table::shift(y, type = 'lead'), by = ring]
  bird_data[, dt := as.numeric(difftime(timestamp2, timestamp, units = 'sec'))]
  
  gspeed = c()
  gdir=c()
  wdir=c()
  for (i in (1:nrow(bird_data)))
  {
    gspeed = append(gspeed,distGeo(c(bird_data$x[i], bird_data$y[i]),c(bird_data$x2[i], bird_data$y2[i]))/bird_data$dt[i])
    gdir = append(gdir, geosphere::bearing(c(bird_data$x[i], bird_data$y[i]),c(bird_data$x2[i], bird_data$y2[i])))
    
  }
  for (i in (1:nrow(wind_data_)))
  {
    wdir = append(wdir, windR::bearing(0,0,wind_data_$u[i],wind_data_$v[i]))
  }
  bird_data[, gspeed := gspeed]
  bird_data[, gdir := gdir]
  wind_data_[,wspeed := sqrt(wind_data_$u^2+wind_data_$v^2)]
  wind_data_[,wdir := wdir*180/pi]

###CONVERT TO RASTER ----
  r <- raster(xmn=(min(bird_data$x)-1), ymn=(min(bird_data$y)-1), xmx=(max(bird_data$x)+1), ymx=(max(bird_data$y)+1), res=2)
  r.gdir <- rasterize(bird_data[,c("x","y")], r, field = bird_data$gdir, fun = mean, na.rm = T)
  r.gspeed <- rasterize(bird_data[,c("x","y")], r, field = bird_data$gspeed, fun = mean, na.rm = T)
  r.n <- rasterize(bird_data[,c("x","y")], r, field = bird_data$x2, fun = 'count', na.rm = T)
  r.wdir <- rasterize(wind_data_[,c("x","y")],r, field = wind_data_$wdir,fun = mean)
  r.wspeed <- rasterize(wind_data_[,c("x","y")],r, field = wind_data_$wspeed,fun = mean)

###BACK TO DATA FRAME (pxdata = data for each pixel) ----
  pxdata <- merge(as.data.frame(r.gdir, xy = T), as.data.frame(r.gspeed, xy = T), by = c("x", "y"))
  pxdata <- merge(pxdata, as.data.frame(r.n, xy = T), by = c("x","y"))
  pxdatawind <- merge(as.data.frame(r.wdir, xy = T), as.data.frame(r.wspeed, xy = 'T'), by = c("x","y"))
  pxdata <- merge(pxdata, pxdatawind, by = c("x","y"))
  setnames(pxdata, c("layer.x.x","layer.y.x","layer", "layer.x.y", "layer.y.y"), c("gdir","gspeed","n","wdir","wspeed"))
  pxdata <- pxdata %>% as.data.table


###DISPLAY###
world <- ne_countries(scale = "medium", returnclass = "sf")
a=4 #Arrow size
plottraj =  ggplot(data = bird_data) + 
  geom_segment(aes(x = x, y = y, xend = x2, yend = y2, color=migr10),arrow = arrow(length = unit(.1, "cm"))) +
  geom_point(aes(x,y),size =.7) +
  geom_sf(data=world ,fill = "black", color = "black") + 
  coord_sf(xlim = c(min(bird_data$x)-1, max(bird_data$x)+1), ylim = c(min(bird_data$y)-1, max(bird_data$y)+1), expand = FALSE)

plotbird =  ggplot(data = pxdata[not(is.na(pxdata$n))]) + 
  geom_segment(aes(x = x, y = y, xend = x+gspeed/a*sind(gdir), yend = y+gspeed/a*cosd(gdir), color=n),arrow = arrow(length = unit(.1, "cm"))) +
  scale_color_gradient("n", low = "grey", high = "red") +
  geom_sf(data=world ,fill = "black", color = "black") + 
  coord_sf(xlim = c(min(bird_data$x)-1, max(bird_data$x)+1), ylim = c(min(bird_data$y)-1, max(bird_data$y)+1), expand = FALSE)
plotwind = ggplot(data = pxdata) +
  geom_segment(aes(x = x, y = y, xend = x+wspeed/a*sind(wdir), yend = y+wspeed/a*cosd(wdir)),arrow = arrow(length = unit(.1, "cm"))) +
  geom_sf(data=world ,fill = "black", color = "black") + 
  coord_sf(xlim = c(min(bird_data$x)-1, max(bird_data$x)+1), ylim = c(min(bird_data$y)-1, max(bird_data$y)+1), expand = FALSE)


