##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## SeaWinds - Data processing & Analyses
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sapply(c('raster','stringr', 'data.table', 'magrittr', 'sf', "rnaturalearthdata","rnaturalearth", 'rWind', 'ggplot2'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE)))

path_wind <- "data/ERA-Interrim/monthly_raw_data"
path_save <- "data/ERA-Interrim/Interanual_means/"

windRDS <- list.files(path=path_wind, pattern= 'ERA_Interim_monthly_sfc_10_11_2013_to_2018')
windRDS <- windRDS[grep(".RDS", windRDS)]
wind_data <- readRDS(paste0(path_wind,"/", windRDS)) %>% as.data.table

months <- unique(format(wind_data$datetime_, "%m"))

r <- raster(xmn=-70, xmx=70, ymn=30, ymx=85, res=.75) #Empty raster of the studied area

for (month in months)
{
  wind_data_ <- wind_data[as.logical(str_count(wind_data$datetime_, pattern = month))]
  u_mean <- rasterize(wind_data_[,c("x","y")], r, field = wind_data_$u, fun = mean)%>% as.data.frame(xy=T)%>% data.table
  v_mean <- rasterize(wind_data_[,c("x","y")], r, field = wind_data_$v, fun = mean)%>% as.data.frame(xy=T)%>% data.table
  setnames(u_mean, c('x', 'y','u'))
  setnames(v_mean, c('x', 'y','v'))
  wind_mean <- merge(u_mean, v_mean, by = c("x","y"))
  wind_mean[, month := rep(month, nrow(wind_mean))]
  saveRDS(wind_mean, paste0(path_save, 'ERA_Interim_interanual_monthly_mean_sfc_',month,'_2013_to_2018.RDS'), compress='xz')

}
