
## Creates a raster of the density of spatial lines
## Returns a density raster

#' @param r empty raster of the studied area
#' @param l.bursts Spatial Lines object 

lineDens <- function(r, l.bursts) 
{
  r.dens <- r
  rpoly <- rasterToPolygons(r, na.rm=T) #Polygons with raster celles
  density <- foreach(i = 1:ncell(r.dens))%do%
  {
    cat(paste0(i,"\n"))
    poly_i <- Polygon(rpoly@polygons[[i]]@Polygons[[1]]@coords)
    poly_i <- Polygons(list(poly_i),1)%>%list%>%SpatialPolygons() #poly of cell i
    length(raster::intersect(l.bursts, poly_i)) #number of lines in the cell
  }
  r.dens[] <- as.numeric(density)
  r.dens
}