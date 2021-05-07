##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## SeaWinds - Data processing & Analyses
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sapply(c('raster','stringr', 'data.table', 'magrittr', 'sf', "rnaturalearthdata","rnaturalearth", 'rWind', 'ggplot2'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE)))

path_wind <- "data/ERA-Interrim/"
path_save <- "data/ERA-Interrim/wind_by_traj/"

path_bird <- "data/Kittiwake_data_treated/"
bird_filename <- "BLKI_tracks.rds"

windRDS <- list.files(path=path_wind, pattern= 'ERA_Interim_daily_sfc_2013_to_2018')
windRDS <- windRDS[grep(".RDS", windRDS)]
#wind_data <- readRDS(paste0(path_wind,"/", windRDS)) %>% as.data.table

bird_data <- readRDS(paste0(path_bird, bird_filename)) %>% as.data.table
bird_data[, timestamp := as.POSIXct(bird_data$timestamp)] 

traj_list <- unique(bird_data$journey)

years <- seq(2013,2018)

for (id in traj_list[155])
{
  cat(id)
  n_pos <- nrow(bird_data[bird_data$journey==id])
  day_1 <- bird_data$timestamp[bird_data$journey==id][1]
  day_n <- bird_data$timestamp[bird_data$journey==id][n_pos]
  days <- format(seq(day_1,day_n, by="days"), "%m-%d")
  wind_mean <- data.frame(x=c(),y=c(),u=c(),v=c())
  
  for (year in years)
  {
    dates <- paste0(year,"-",days)
    wind_data_ <- wind_data[FALSE]
    for (date in dates)
    {
      wind_data_ <- rbind(wind_data_,wind_data[wind_data$datetime_==dates])
    }
    
    r <- raster(xmn=-70, xmx=70, ymn=30, ymx=85, res=.75) #Empty raster of the studied area
    
    u_mean <- rasterize(wind_data_[,c("x","y")], r, field = wind_data_$u, fun = mean)%>% as.data.frame(xy=T)%>% data.table
    v_mean <- rasterize(wind_data_[,c("x","y")], r, field = wind_data_$v, fun = mean)%>% as.data.frame(xy=T)%>% data.table
    setnames(u_mean, c('x', 'y','u'))
    setnames(v_mean, c('x', 'y','v'))
    wind_mean <- rbind(wind_mean,merge(u_mean, v_mean, by = c("x","y")))
    
  }
  
  u_mean <- rasterize(wind_mean[,c("x","y")], r, field = wind_mean$u, fun = mean)%>% as.data.frame(xy=T)%>% data.table
  v_mean <- rasterize(wind_mean[,c("x","y")], r, field = wind_mean$v, fun = mean)%>% as.data.frame(xy=T)%>% data.table
  setnames(u_mean, c('x', 'y','u'))
  setnames(v_mean, c('x', 'y','v'))

  wind_mean <- merge(u_mean, v_mean, by = c("x","y"))
  wind_mean[, date := rep(paste0(format(day_1, "%m-%d"), " to ",format(day_n, "%m-%d")), nrow(wind_mean))]
  
  ggplot(wind_mean)+
    geom_segment(aes(x=x,y=y,xend=x+a*u,yend=y+a*v), arrow=arrow(length=unit(0.10,"cm")))+
    geom_sf(data=world)+
    coord_sf(xlim=c(-60,60),ylim=c(30,85))
  
#saveRDS(wind_mean, paste0(path_save, 'ERA_Interim_',id,'_interanual_mean_sfc_2013_to_2018.RDS'), compress='xz')
}
