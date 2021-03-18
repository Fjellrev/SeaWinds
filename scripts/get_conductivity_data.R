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
birdRDS <- list.files(path=path_bird, pattern = "5col")
bird_data <- readRDS(paste0(path_bird,"/", birdRDS)) %>% as.data.table
bird_data[, timestamp := as.POSIXct(bird_data$timestamp)]
bird_data = bird_data[abs(as.numeric(format(bird_data$timestamp, "%H"))-12)>6] #keep only the midnight positions
###ADD THE ACTIVITY DATA ----
bird_act <- data.table(ring=c(),timestamp=c(),std_cond=c())

ring <- unique(bird_data$ring)
for (id in ring)
{
cat(id)
actRDS <- list.files(path=path_act, pattern = id)
act_data <- readRDS(paste0(path_act,"/",actRDS)) %>% as.data.table
act_data[, timestamp := as.POSIXct(act_data$timestamp)]
id_timestamp <- as.POSIXct(bird_data$timestamp[bird_data$ring==id])
mean_act <- c()
for (t in format(id_timestamp, "%D"))
{
mean_act <- append(mean_act, mean(act_data$std_conductivity[format(act_data$timestamp, "%D") ==t], na.rm = T))
}
bird_act <- rbind(bird_act, data.table(ring=rep(id, length(id_timestamp)),timestamp=id_timestamp, std_cond=mean_act))
}

bird_data[,std_conductivity := bird_act$std_cond]

saveRDS(bird_data, paste0("data/Kittiwake_data_treated","/","BLKI_5col_cond_sept_dec2016.RDS"))
