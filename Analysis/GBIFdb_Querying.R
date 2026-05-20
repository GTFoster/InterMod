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

#keyquery <- keyquery |>
#  dplyr::arrange(n) #start with smallest first and go up from there

#for(j in idex[-length(idex)]){
for(j in idex){
  #keystemp <- keyquery$key[(j+1):(j+5)]
  keystemp <- keyquery$key[j]
  Sys.time()
    spdat <- gbif |> 
      #filter(taxonkey %in% key) |> 
      filter(taxonkey %in% keystemp) |> 
      collect() #Grab all occurence records for the species of interest and call is spdat, using our backbone key to search
    Sys.time()
    outrun <- spdat
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
    
    #for(i in (j+1):(j+5)){
    for(i in j){
      filname <- keyquery$name[i] %>% tolower() %>% gsub(pattern=" ", replacement="_",.) %>%
      paste("data/GBIF/occs/", ., ".csv", sep="") #Create a filename
      tempocc <- dplyr::filter(occs, taxonkey==keyquery$key[i])
      write.csv(tempocc, file=filname, row.names = FALSE) 
      IDkey <- matches$verbatim_index[matches$verbatim_name==keyquery$name[i]]
      if(length(IDkey)==0){
        IDkey <- NA
      }
      tempcheck <- data.frame(Name=keyquery$name[i], ID=IDkey, NROW=try(nrow(tempocc)), time=Sys.time())
    }
    print(paste("done with", j, sep=" "))
    nmlist <- rbind(nmlist, tempcheck)
    rm(occs, tempocc)
    gc()
    write.csv(nmlist, file="Progresscheck_birds_big.csv")
}
