rm(list = ls())
setwd("ASCAT")

sapply(c('rgdal', 'raster', 'data.table', 'magrittr', 'sp', 'rgeos', 'raster', 'foreach', 'ncdf4', 'gdalUtilities'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE) ) )

netcdf_ = list.files(pattern = "ASCAT_monthly")

u = brick(netcdf_, varname = "x_wind")
v = brick(netcdf_, varname = "y_wind")

dd = data.table(rastNam = names(u))
dd[, datetime_ := substring(rastNam, 2, 11) %>% as.POSIXct(., format = '%Y.%m.%d') ]

write.table(dd,'ASCAT_names.txt', sep = ';', row.names = FALSE)

writeRaster(u, filename = 'ASCAT_monthly_201309_201912_u_wind.tif', overwrite=TRUE)
writeRaster(v, filename = 'ASCAT_monthly_201309_201912_v_wind.tif', overwrite=TRUE)

names(u) = dd$rastNam
names(v) = dd$rastNam

o = foreach(i = 1:nlayers(u) ) %do% {
  if(i!=9){ #No data on layer 9, i need to find a condition more general to specify if no data
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
}}
w = rbindlist(o)

saveRDS(w, 'ASCAT_monthly_201309_201912_wind.RDS', compress='xz')
