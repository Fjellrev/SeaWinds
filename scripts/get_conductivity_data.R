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
cat(id)
actRDS <- list.files(path=path_act, pattern = id)
act_data <- readRDS(paste0(path_act,"/",actRDS)) %>% as.data.table
act_data[, timestamp := as.POSIXct(act_data$timestamp)]
timestamp <- (unique(format(bird_data$timestamp[bird_data$ring==id])))
tmean_act <- c()
for (t in timestamp)
{
dt<-difftime(format(act_data$timestamp),t,units="hours")
mean_act <- append(mean_act, mean(act_data$std_conductivity[dt<12&dt>=0], na.rm = T))
}
bird_act <- rbind(bird_act, data.table(ring=rep(id, length(timestamp)),timestamp=timestamp, std_cond=mean_act))
}

bird_data[,std_conductivity := bird_act$std_cond[ring==bird_act$ring&timestamp==bird_act$timestamp], by = 1:nrow(bird_data)]

saveRDS(bird_data, paste0(path_bird,"/","BLKI_Alkefjellet_cond_nov2016.RDS"))