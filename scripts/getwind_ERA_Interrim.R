rm(list = ls())
setwd("Era_Interrim")

sapply(c('rgdal', 'raster', 'data.table', 'magrittr', 'sp', 'rgeos', 'raster', 'foreach', 'ncdf4', 'gdalUtilities'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE) ) )

netcdf_ = list.files(pattern = "ERA_Interrim_monthly")

u = brick(netcdf_, varname = "u10")
v = brick(netcdf_, varname = "v10")

dd = data.table(rastNam = names(u))
dd[, datetime_ := substring(rastNam, 2, 11) %>% as.POSIXct(., format = '%Y.%m.%d') ]

write.table(dd,'ERA_Interrim_names.txt', sep = ';', row.names = FALSE)

writeRaster(u, filename = 'ERA_Interim_monthly_201201_201908_u_wind.tif', overwrite=TRUE)
writeRaster(v, filename = 'ERA_Interim_monthly_201201_201908_v_wind.tif', overwrite=TRUE)

names(u) = dd$rastNam
names(v) = dd$rastNam

o = foreach(i = 1:nlayers(u) ) %do% {
  # u wind component
  xu = as(u[[i]] , "SpatialPixelsDataFrame")
  xu = xu %>% as.data.frame %>% data.table
  xu[, datetime_ := dd$datetime_[i]]
  setnames(      xu, c('u', 'x', 'y', 'datetime_'))
  setcolorder(   xu, c('x', 'y', 'datetime_', 'u'))
  xu[, xy        := paste0(x, '_',y)]
  
  # v wind component
  xv = as(v[[i]] , "SpatialPixelsDataFrame")
  xv = xv %>% as.data.frame %>% data.table
  xv[, datetime_ := dd$datetime_[i]]
  setnames(      xv, c('v', 'x', 'y', 'datetime_'))
  setcolorder(   xv, c('x', 'y', 'datetime_', 'v'))
  xv[, xy        := paste0(x, '_',y)]
  
  # merge u and v
  xuv = merge(xu, xv[, c('xy', 'v', 'datetime_'), with = FALSE], by = c('xy', 'datetime_') )
  
  # get rid of useless columns
  xuv[, ':=' (xy = NULL)]
  cat(i)
  xuv
}
w = rbindlist(o)

saveRDS(w, 'ERA_Interim_monthly_201201_201908_wind.RDS', compress='xz')
