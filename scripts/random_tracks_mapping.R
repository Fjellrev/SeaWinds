

myfiles <- list.files("outputs/observed_sim_gc_tracks", full.names = TRUE, pattern = "_n5.rds")
col_ws = c('firebrick4', 'firebrick3', 'gold', 'gold', 'springgreen3', 'springgreen4')

for(i in 1:length(myfiles)){

    traj <- readRDS(myfiles[i])
    
    plot <- ggplot(data=traj) +
    geom_path(size = .75, aes(x = lon, y = lat,color=mean_ws, group = N)) +
    scale_colour_gradientn(colours = col_ws)+
    geom_path(data = traj[traj$N=="obs",], size = 1.25, aes(x = lon, y = lat),color='black') +
    geom_path(data = traj[traj$N=="gc",], size = 1.25, aes(x = lon, y = lat),color='blue') +
    geom_sf(data=world_map,fill = "black", color = "black") +
    coord_sf(xlim = c(-70, 40), ylim = c(40,85), expand = FALSE)
      
    ggsave(plot = plot, filename = paste0(gsub(".rds", "", myfiles[i]), ".png"))
    
}



