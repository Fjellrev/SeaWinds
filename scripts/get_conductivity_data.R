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
bird_data <- bird_data[not(abs(as.numeric(format(bird_data$timestamp, "%H"))-12)<6)] #keep only the midnight positions


###ADD THE ACTIVITY DATA ----

bird_act <- data.table(ring=c(),timestamp=c(),std_cond=c())

ring <- unique(bird_data$ring)
for (id in ring)
{
cat(id)
actRDS <- list.files(path=path_act, pattern = id)
act_data <- readRDS(paste0(path_act,"/",actRDS)) %>% as.data.table #activity for the individual id
act_data[, timestamp := as.POSIXct(act_data$timestamp)] 
id_timestamp <- as.POSIXct(bird_data$timestamp[bird_data$ring==id]) #list of the GLS timestamp for the individual id
id_act <- c() 

for (i in 1:(length(id_timestamp)-1)) #for each timestamp get the conductivity between t and t+1
{
t <- id_timestamp[i]
t_dt <- id_timestamp[i+1]
act <- mean(act_data$std_conductivity[act_data$timestamp>t&act_data$timestamp<t_dt], na.rm = T)#mean activity between t and t+1 
id_act <- append(id_act, act)
}

bird_act <- rbind(bird_act, data.table(ring=rep(id, length(id_timestamp)),timestamp=id_timestamp, std_cond=append(id_act,NA)))
}

bird_data[,std_conductivity := bird_act$std_cond]

saveRDS(bird_data, paste0("data/Kittiwake_data_treated","/","BLKI_5col_cond_sept_dec2016.RDS"))
