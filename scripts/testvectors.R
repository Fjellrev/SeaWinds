##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## SeaWinds - Testing the calculated vector
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sapply(c( 'raster','stringr', 'data.table', 'magrittr', 'sf', "rnaturalearthdata","rnaturalearth", 'windR', 'geosphere', 'ggplot2'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE)))

pxdata = readRDS("outputs/monthly_pixel_data.rds")
i = 11

print(pxdata[which(not(is.na(pxdata$gdir)))[i]])

ggplot(data=pxdata[which(not(is.na(pxdata$gdir)))[i]]) + 
  geom_segment(aes(x = 0, y = 0, xend = gspeed*sind(gdir), yend = gspeed*cosd(gdir)), color = 'red', size = 1, arrow = arrow(length = unit(1, "cm"))) + #Ground 
  geom_segment(aes(x = 0, y = 0, xend = wspeed*sind(wdir), yend = wspeed*cosd(wdir)), color = 'blue', size = 1, arrow = arrow(length = unit(1, "cm"))) + #Wind 
  geom_segment(aes(x = 0, y = 0, xend = ws*sind(gdir), yend = ws*cosd(gdir)), color = 'black', size = 1, arrow = arrow(length = unit(1, "cm"))) + #wind support 
  geom_segment(aes(x = 0, y = 0, xend = cw*sind(gdir+90), yend = cw*cosd(gdir+90)), color = 'black', size = 1, arrow = arrow(length = unit(1, "cm"))) + #crosswind
  geom_segment(aes(x=wspeed*sind(wdir), y =wspeed*cosd(wdir), xend = wspeed*sind(wdir) + aspeed*sind(adir), yend = wspeed*cosd(wdir) + aspeed*cosd(adir)), color = 'green', size = 1, arrow = arrow(length = unit(1, "cm"))) #air
