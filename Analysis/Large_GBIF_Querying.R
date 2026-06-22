library(tidyverse)
library(gbifdb)
library(rgbif) #Only used to find taxonkey
library(CoordinateCleaner)

gbif <- gbif_local(dir='/home/shared/occurrence/2024-04-01') #Tell R where to find our local GBIF copy
#I apologize to the reproducibility gods for absolute paths

matches <- read.csv(file="backbone_matches.csv")

nmlist <- data.frame(Name=NA, ID=0, NROW=0, time=Sys.time()) #Start an empty df to keep track of query progress

keyquery <- NULL
for(name in matches$verbatim_name){
  filname <- name %>% tolower() %>% gsub(pattern=" ", replacement="_",.) %>%
      paste("data/GBIF/occs/", ., ".csv", sep="") #Create a filename
    if(file.exists(filname)==TRUE){ #If the occ file already exists locally and is filled, skip and move to the next species
      if(file.info(filname)$size > 150){
      next
      }
    }
  key <- matches[matches$verbatim_name==name,] #Save the key to search by
  keyquery <- rbind(keyquery, data.frame(name=name, key=key$speciesKey, class=key$class))
}

###TEMP: Let's do non-birds first
#keyquery <- dplyr::filter(keyquery, class != "Aves")
#Each query seems to take ~8min, pretty much reguardless of size. So to speed things up, we're going to call our gbifdb query in chunks and then split the outputs up afterwards
#idex <- c(0, seq(from=1, to=nrow(keyquery), by=5)[-1], nrow(keyquery))
idex <- c(seq(1, nrow(keyquery)))

#keyquery <- gbif |> #adding in counts of bird species records
#  filter(taxonkey %in% keyquery$key) |> 
#  count(taxonkey) |>
#  collect() |>
#  dplyr::rename("key"=taxonkey) %>%
#  left_join(keyquery, ., by="key")


lat_breaks <- c(-90, 0, 90)
lon_breaks <- c(-180, -90, 0, 90, 180)

panels <- expand.grid(
  lat_i = seq_len(length(lat_breaks) - 1),
  lon_i = seq_len(length(lon_breaks) - 1)
)

panels$lat_min <- lat_breaks[panels$lat_i]
panels$lat_max <- lat_breaks[panels$lat_i + 1]
panels$lon_min <- lon_breaks[panels$lon_i]
panels$lon_max <- lon_breaks[panels$lon_i + 1]

#keyquery <- keyquery |>
#  dplyr::arrange(n) #start with smallest first and go up from there
#for(j in idex[-length(idex)]){
for (p in seq_len(nrow(panels))) {
  
  lat_min <- panels$lat_min[p]
  lat_max <- panels$lat_max[p]
  lon_min <- panels$lon_min[p]
  lon_max <- panels$lon_max[p]
  
  for(j in idex){
    #keystemp <- keyquery$key[(j+1):(j+5)]
    keystemp <- keyquery$key[j]
    Sys.time()
    
    spdat <- gbif |> 
      filter(taxonkey %in% keystemp) |> 
      filter(
        !is.na(decimallatitude),
        !is.na(decimallongitude),
        decimallatitude >= lat_min,
        decimallatitude <  lat_max,
        decimallongitude >= lon_min,
        decimallongitude <  lon_max
      ) |>
      collect()
    
    if (nrow(spdat) == 0) {
      print(paste("no records for key", j, "panel", p))
      rm(spdat)
      gc()
      next
    }
      spdat <- dplyr::filter(spdat, is.na(decimallongitude)==F) #remove recs without latlong
      
      spdat <- spdat %>% #Do some data cleaning
          dplyr::filter(occurrencestatus  == "PRESENT") %>%
          dplyr::filter(!basisofrecord %in% c("FOSSIL_SPECIMEN")) %>% #Remove fossils
          CoordinateCleaner::cc_cen( #Remove records within 1km buffer of country centroids
          lon = "decimallongitude", 
          lat = "decimallatitude", 
          buffer = 1000, # radius of circle around centroid to look for centroids
          value = "clean",
          test="both")  %>% 
          cc_sea( #Remove oceanic records
        lon = "decimallongitude",
        lat = "decimallatitude"
        )# %>%
        #dplyr::filter(., establishmentMeans != "introduced") #Remove known introduced  spp
      
      occs <- dplyr::select(spdat, genus, species,scientificname, decimallatitude, decimallongitude, elevation, day, month, year, taxonkey) %>% unique() 
      #subset across species and save into seperate files
      filname <- keyquery$name[j] %>% tolower() %>% gsub(pattern=" ", replacement="_",.) %>%
        paste("data/GBIF/occs/", ., "_panel", p, ".csv", sep="") #Create a filename
      write.csv(occs, file=filname, row.names = FALSE) 
      print(paste("done with", j, sep=" "))
      rm(occs)
      gc()
  }
}
