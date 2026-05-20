library(tidyverse)
library(gbifdb)
library(CoordinateCleaner)

load(file="missingtoQuery.RDA") #Load in our total species list
gbif <- gbif_local(dir='/home/shared/occurrence/2024-04-01') #Tell R where to find our local GBIF copy
#I apologize to the reproducibility gods for absolute paths

  spdat <- gbif |> 
    filter(species %in% unlist(missing)) |> 
    collect() #Grab all occurence records for the species of interest and call is spdat
  
  spdat <- dplyr::filter(spdat, is.na(decimallongitude)==F)
  
  spdat <- spdat %>% #Do some data cleaning
    dplyr::filter(occurrencestatus  == "PRESENT") %>%
    dplyr::filter(!basisofrecord %in% c("FOSSIL_SPECIMEN")) %>% #Remove fossils
    CoordinateCleaner::cc_cen( #Remove records within 1k buffer of country centroids
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
  
  missing_occs <- dplyr::select(spdat, genus, species,scientificname, decimallatitude, decimallongitude, elevation, day, month, year, taxonkey) %>% unique()
  write.csv(missing_occs, file="missingOccs.csv", row.names = FALSE)