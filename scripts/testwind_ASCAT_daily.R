##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## SeaWinds - Testing wind data
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sapply(c('rgdal', 'raster', 'data.table', 'magrittr', 'sp', 'rgeos', 'raster', 'foreach', 'ncdf4', 'gdalUtilities'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE) ) )

netcdf_ = list.files(path="data/ASCAT", pattern = "ASCAT_daily")

u = brick(paste0("data/ASCAT/", netcdf_), varname = "x_wind")
v = brick(paste0("data/ASCAT/", netcdf_), varname = "y_wind")

daily_na = c()
for (i in (1:nlayers(u)))
{
  daily_na = append(daily_na, length(which(is.na(u[[i]][]))))
}

study_area = length(u[[1]][])-min(daily_na) #size of the largest zone that can be studied in the sample
daily_na = (daily_na - min(daily_na))/study_area #proportion of the study area that is lost
print(daily_na)

