##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## SeaWinds - Testing how many na are present in a wind dataset
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sapply(c( 'raster','stringr', 'data.table', 'magrittr', 'sf', "rnaturalearthdata","rnaturalearth", 'windR', 'geosphere', 'ggplot2', 'circular'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE)))

cosd <-function(x) cos(x*pi/180) #cos of an angle in degrees
sind <- function(x) sin(x*pi/180)

path_wind<- "data/ERA-Interrim"

RDS_sfc <- list.files(path=path_wind, pattern = "sfc")
RDS_100 <- list.files(path=path_wind, pattern = "100m")

wind_sfc <- readRDS(paste0(path_wind,"/",RDS_sfc)) %>% as.data.table
wind_100 <- readRDS(paste0(path_wind,"/",RDS_100)) %>% as.data.table

###Display an example of the wind conditions---

world <- ne_countries(scale = "medium", returnclass = "sf")
a <- 0.1 #size of the arrows
day <- 5 #day we want to see the wind conditions

#ggplot(data = wind_sfc[datetime_==unique(wind_sfc$datetime_)[day]])+
  #geom_segment(size = 1, aes(x = x, y = y, xend = x+a*u, yend = y+a*v),arrow = arrow(length = unit(.09, "cm"))) +
  #geom_sf(data=world ,fill = "black", color = "black") + 
  #coord_sf(xlim = c(min(wind_sfc$x)-1, max(wind_sfc$x)+1), ylim = c(min(wind_sfc$y)-1, max(wind_sfc$y)+1), expand = FALSE)

#ggplot(data = wind_100[datetime_==unique(wind_100$datetime_)[day]])+
  #geom_segment(size = 1, aes(x = x, y = y, xend = x+a*u, yend = y+a*v),arrow = arrow(length = unit(.09, "cm"))) +
  #geom_sf(data=world ,fill = "black", color = "black") + 
  #coord_sf(xlim = c(min(wind_sfc$x)-1, max(wind_sfc$x)+1), ylim = c(min(wind_sfc$y)-1, max(wind_sfc$y)+1), expand = FALSE)

###Compare the speed and direction---

wind_sfc[, wspeed     := sqrt(u^2 + v^2), by = 1:nrow(wind_sfc)] #wind speed
wind_100[, wspeed     := sqrt(u^2 + v^2), by = 1:nrow(wind_100)] #wind speed

wind_sfc[, wdir := atan2(u, v)*180/pi, by = 1:nrow(wind_sfc)] #wind direction
wind_100[, wdir := atan2(u, v)*180/pi, by = 1:nrow(wind_100)] #wind direction

wind_sfc[, dspeed := wspeed-wind_100$wspeed] #absolute difference in speed
wind_sfc[, ddir := asin(sind(wdir-wind_100$wdir))*180/pi] #angle between the directions

###Statistics---

ggplot(wind_sfc)+geom_histogram(aes(x=ddir)) #direction
summary(wind_sfc$ddir)
t.test(wind_sfc$ddir, mu = 0)

ggplot(wind_sfc)+geom_histogram(aes(x=dspeed)) #speed
summary(wind_sfc$dspeed)
t.test(wind_sfc$dspeed, mu = 0)
