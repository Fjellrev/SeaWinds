
## Check whether points overlap background polygons
## Returns TRUE/FALSE

#' @param x x-coordinates
#' @param y y-coordinates
#' @param prj projection (as proj.4 string)
#' @param mask underlying polygon (e.g. background map with landmasses) 

is.land <- function(x, y, prj, mask) #return TRUE if a point is over lands
{
  pt <- cbind.data.frame(x, y)
  pt <- st_as_sf(pt, coords=1:2, crs=CRS(prj))
  !is.na(as.numeric(suppressMessages(st_intersects(pt, mask))))
}