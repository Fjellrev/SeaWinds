
dirWind <- "data/ERA-Interrim/"

sapply(c('rgdal', 'raster', 'data.table', 'magrittr', 'sp', 'rgeos', 'raster', 'foreach', 'ncdf4', 'gdalUtilities'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE) ) )

netcdf_ <- list.files(dirWind, pattern = "ERA_Interim_daily_sfc_09_2016_to_12_2016")
netcdf_ <- netcdf_[grep(".nc", netcdf_)]

u <- brick(paste0(dirWind, netcdf_), varname = "u10")
v <- brick(paste0(dirWind, netcdf_), varname = "v10")

dd <- data.table(rastNam = names(u))
dd[, datetime_ := substring(rastNam, 2, 11) %>% as.POSIXct(., format = '%Y.%m.%d') ]

names(u) <- dd$rastNam
names(v) <- dd$rastNam

o <- foreach(i = 1:nlayers(u) ) %do% {
  
  if(length(which(is.na(u[[i]][]))) != length(u[[i]][])){ # check if no data
    
    # u wind component
    xu <- as(u[[i]] , "SpatialPixelsDataFrame")
    xu <- xu %>% as.data.frame %>% data.table
    xu[, datetime_ := dd$datetime_[i]]
    setnames(xu, c('u', 'x', 'y', 'datetime_'))
    setcolorder(xu, c('x', 'y', 'datetime_', 'u'))
    xu[, xy        := paste0(x, '_',y)]
    
    # v wind component
    xv <- as(v[[i]] , "SpatialPixelsDataFrame")
    xv <- xv %>% as.data.frame %>% data.table
    xv[, datetime_ := dd$datetime_[i]]
    setnames(      xv, c('v', 'x', 'y', 'datetime_'))
    setcolorder(   xv, c('x', 'y', 'datetime_', 'v'))
    xv[, xy        := paste0(x, '_',y)]
    
    # merge u and v
    xuv <- merge(xu, xv[, c('xy', 'v', 'datetime_'), with = FALSE], by = c('xy', 'datetime_') )
    
    # get rid of useless columns
    xuv[, ':=' (xy = NULL)]
    cat(i)
    xuv
  }
}

w <- rbindlist(o)

saveRDS(w, paste0(dirWind, 'ERA_Interim_daily_sfc_09_2016_to_12_2016.RDS'), compress='xz')

