##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## SeaWinds - Activity data
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sapply(c( 'raster','stringr', 'data.table', 'magrittr', 'sf', "rnaturalearthdata","rnaturalearth", 'windR', 'geosphere', 'ggplot2', 'circular'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE)))
path_act = "data/activity"
path_bird = "data/Kittiwake_data"
  
###GET BIRD DATA ----
birdRDS <- list.files(path=path_bird, pattern = "Alkefjellet_nov2016")
bird_data <- readRDS(paste0(path_bird,"/", birdRDS)) %>% as.data.table
bird_data[, timestamp := as.POSIXct(bird_data$timestamp)]

###ADD THE ACTIVITY DATA ----
bird_act <- data.table(ring=c(),timestamp=c(),std_cond=c())

ring <- unique(bird_data$ring)
for (id in ring)
{
actRDS <- list.files(path=path_act, pattern = id)
act_data <- readRDS(paste0(path_act,"/",actRDS))
timestamp <- unique(format(act_data$timestamp, "%D"))
mean_act <- c()
for (d in timestamp)
{
mean_act <- append(mean_act, mean(act_data$std_conductivity[format(act_data$timestamp, "%D") == d], na.rm = T))
}
bird_act <- rbind(bird_act, data.table(ring=rep(id, length(timestamp)),timestamp=timestamp, std_cond=mean_act))
}

bird_data[,std_conductivity := bird_act$std_cond[ring==bird_act$ring&format(timestamp, "%D")==bird_act$timestamp], by = 1:nrow(bird_data)]

saveRDS(bird_data, paste0(path_bird,"/","BLKI_Alkefjellet_cond_nov2016.RDS"))