
# quick & dirty function to map points against world map background



quick.map <- function(xi, yi, wrld){
  ggplot() +  geom_sf(data=wrld, fill = "black", color = "black") + 
    geom_point(aes(x=start_pt$x, y =start_pt$y), color = 'blue', size = 2.5) + 
    geom_point(data = end_pt, aes(x, y), colour = 'red', size = 2.5) + 
    geom_point(aes(x=xi, y =yi), color = 'green', size = 2.5)
}