##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## SeaWinds - Data processing & Analyses
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sapply(c('sf','sp','spData','tidyverse', 'data.table','MASS','plyr', 'magrittr', 'gdistance','geosphere','raster', "doParallel", "foreach",
         'ggplot2', 'rWind','windR','adehabitatHR', 'circular', 'maptools'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE)))

### utility functions ----
source("functions/FUNCTION_OverlapPoly.r")
source('functions/FUNCTION_QuickMap.r')
cosd <-function(x) cos(x*pi/180) #cos of an angle in degrees 
sind <-function(x) sin(x*pi/180) 
proj.latlon <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #classic longlat projection to use the windR function

  
map_path <- "data/baseline_data"
map_filename <- "SEAwinds_Worldmap_res-c.rds"
world_map <- st_as_sf(readRDS(paste0(map_path,"/",map_filename)), coords=c("lon","lat"), crs = 4326)
data(wrld_simpl)

wind_path <- "data/ERA_Interim/monthly_raw_data"
wind_filename <- list.files(path=wind_path)
wind_filename <- wind_filename[grep(".RDS", wind_filename)]
wind_data <- do.call('rbind', lapply(paste0(wind_path,"/", wind_filename), readRDS)) %>% as.data.table

wind_data[, wdir := atan2(u, v)] #wind direction
wind_data[, wspeed := sqrt(u^2+v^2)]
wind_data$season <- "autumn"
wind_data$season[as.numeric(format(wind_data$datetime_, "%m"))<6] <- "spring"
seasons <- unique(wind_data$season)

r <- raster(xmn=-70.25, xmx=70.25, ymn=29.75, ymx=85.25, res=0.75) #Empty raster of the studied area

r.land <- rasterize(wrld_simpl, r, field = 0) #raster with 0 above lands and NA above the sea
df.land <- r.land %>% as.data.frame(xy=T)
df.land[is.na(df.land)] <- 1 #turn NA into 1 above the sea
df.land$layer[df.land$layer==0] <- NA #turn 0 into na above the lands
r.land <- rasterize(df.land[df.land$layer!=0,][,c("x","y")], r, field = df.land$layer[df.land$layer!=0], background=NA,)

for (season_ in seasons)
{
cat(season_)
wind_data_ <- wind_data[season==season_]
r.wdir_var <- r.land*rasterize(wind_data_[,c("x","y")], r, field = wind_data_$wdir, fun=function(x, na.rm=T){sd.circular(x)})
writeRaster(r.wdir_var,paste0("outputs/Rasters/wdir_sd_", season_),format="raster", overwrite=TRUE)
}

r.wdir_var_spr <- raster("outputs/Rasters/wdir_sd_spring")
r.wdir_var_aut <- raster("outputs/Rasters/wdir_sd_autumn")

df.wind_var_aut <- as.data.frame(r.wdir_var_aut, xy= T) %>% as.data.table
df.wind_var_aut$season <- "Autumn"
df.wind_var_spr <- as.data.frame(r.wdir_var_spr, xy= T) %>% as.data.table
df.wind_var_spr$season <- "Spring"
df.wind_var <- rbind(df.wind_var_spr, df.wind_var_aut)
df.wind_var_min <- df.wind_var
df.wind_var_min$layer[df.wind_var$layer>quantile(df.wind_var$layer,.05, na.rm=T)]<-NA

r.corridors_aut <- raster("outputs/Rasters/summ_corr_autumn")
r.corridors_aut <- r.corridors_aut/cellStats(r.corridors_aut, max)
r.corridors_spr <- raster("outputs/Rasters/summ_corr_spring")
r.corridors_spr <- r.corridors_spr/cellStats(r.corridors_spr, max)
poly.corr_aut <- rasterToPolygons(clump(r.corridors_aut > 0.75), dissolve = TRUE)%>%st_as_sf
poly.corr_spr <- rasterToPolygons(clump(r.corridors_spr > 0.55), dissolve = TRUE)%>%st_as_sf
poly.corr_aut$season <-"Autumn"
poly.corr_spr$season <-"Spring"

col = rev(c('firebrick4', 'firebrick3', 'gold', 'springgreen3', 'springgreen4'))
plot_ <- ggplot(data=df.wind_var)+
  geom_raster(aes(x=x,y=y,fill=layer))+
  geom_sf(data=world_map)+
  geom_sf(data=poly.corr_aut, color="black", size = .75, fill=NA)+
  geom_sf(data=poly.corr_spr, color="black", size = .75,fill=NA)+
  coord_sf(xlim=c(-70,70),ylim=c(30,85))+
  scale_fill_gradientn(colours = col, name = "sd (rad)")+
  facet_grid(.~season)+  
  xlab("")+ylab("")+
  theme(strip.text = element_text(face="bold"))

plot_


saveRDS(df.wind_var_min, file = "outputs/Rasters/wdir_sd_min.rds")
