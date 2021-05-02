##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## SeaWinds - Data processing & Analyses
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sapply(c( 'raster','stringr', 'data.table', 'magrittr', 'sf', "rnaturalearthdata","rnaturalearth", 'windR', 'geosphere', 'ggplot2', 'circular'),
       function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE)))

proj.latlon <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #classic longlat projection to use the windR function

path_bird <- "data/Kittiwake_data_treated"
path_wind <- "data/ERA-Interrim"
world <- ne_countries(scale = "medium", returnclass = "sf")

colonies <- c("Isle of May", "Alkefjellet") #studied colonies
min_length <- 2e06
min_npos <- 5

###GET DATA ----

birdRDS <- list.files(path=path_bird, pattern = "5col")
bird_data <- readRDS(paste0(path_bird,'/', birdRDS)) %>% as.data.table
bird_data[, timestamp := as.POSIXct(bird_data$timestamp)] 
setnames(bird_data, c("lon","lat"),c("x", "y"))
bird_data_ <- bird_data[bird_data$colony%in%colonies]
bird_data_ <- bird_data_[bird_data_$migr10==1]

###PREPARE DATASET ----

rm <- c() #tracks to be removed
for (id in unique(bird_data_$ring))
{
  track_id <- bird_data_[bird_data_$ring==id]
  npos <- nrow(track_id)
  length_track_id <- distGeo(c(track_id$x[1],track_id$y[1]),c(track_id$x[npos],track_id$y[npos]))
  if (npos<min_npos|length_track_id<min_length)
  {
    rm <- append(rm,id)
  }
}

bird_data_ <- bird_data_[not(bird_data_$ring%in%rm)]

saveRDS(bird_data_,"data/Kittiwake_data_treated/BLKI_Alkefjellet_May_tracks.rds")
