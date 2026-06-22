options(java.parameters = "-Xmx20g")
library(rJava)
library(tidyverse)
library(terra)
library(dismo)
library(predicts)

setwd("~/Documents/InterMod/Analysis")

PCA <- terra::rast("worldclim/Clima_PCA_1to8.tif")

trainMaxent <- function(Occ, seed=1, split=0.75, flnm){
  Occ <- as.data.frame(Occ)
  colnames(Occ) <- tolower(colnames(Occ))
  Occ <- dplyr::select(Occ, species, decimallongitude, decimallatitude) %>% unique() #remove redundant records
  Occ['cells']<-terra::cellFromXY(object= PCA[[1]], xy = Occ[,2:3]) #Grab the occurrence cells
  
  Occ <- Occ %>% dplyr::group_by(cells) %>%
    slice_sample() #Grab one unique coordinate from each cell
  
  if(nrow(Occ) > 25000){ #Sample back down to 25k if you're above it
    Occ <- Occ %>% ungroup() %>% sample_n(25000)
  }
  
  #creating a 200km2 buffer to select the psudo-absence points
  OccTerra <- terra::vect(x=Occ[,2:3], crs = "+proj=longlat +datum=WGS84") #Create SpatVector object
  OccBuffer <- terra::buffer(OccTerra, width=200000) #Create 200km buffer
  
  OccBuffer<- terra::aggregate(OccBuffer) #Merge buffer objects into
  
  OccBufferRas <-terra::rasterize(x = OccBuffer, y=PCA) #snap to Clima grid
  OccBufferRas <-terra::mask(x = OccBufferRas, mask = PCA[[1]]) #Mask to remove oceanic points within the buffer.
  
  ##Selecting the potential pseudo absence sites
  AllPseudo <- terra::as.points(OccBufferRas)
  AllPseudo <- terra::extract(x=OccBufferRas[[1]], AllPseudo[,1:2],cell=T)
  AllPseudo <- AllPseudo[-(which(AllPseudo$cell%in%Occ$cells)),] #remove true presence from potential psuedo points
  
  #Actually sample psuedoabsences
  set.seed(seed)
  PseudoA <- sample(x = AllPseudo$cell, size = nrow(Occ)) #Sample a number of psuedoA equal to # of occurances
  PseudoA <- as.data.frame(terra::xyFromCell(object = PCA[[1]], cell = PseudoA)) #Make a df
  
  #Now, train a Maxent model
  Occurrences<-cbind(Occ[,2:3], PseudoA[,1:2]) #bind into 4 cols; left is real, right is associated psuedo
  colnames(Occurrences)<-c('Long','Lat','Long_Psuedo','Lat_Psuedo')
  
  id.training<-sample(1:nrow(Occurrences), round(split*nrow(Occurrences),0)) #create test-train split (split=percent train)
  training <- dismo::prepareData(x=PCA, p=Occurrences[id.training,c("Long", "Lat")], b=Occurrences[id.training,c("Long_Psuedo", "Lat_Psuedo")]) #Use dismo prepare data to set up df
  testing <-dismo::prepareData(x=PCA, p=Occurrences[-id.training,c("Long", "Lat")], b=Occurrences[-id.training,c("Long_Psuedo", "Lat_Psuedo")])
  training<-na.omit(training)
  testing<-na.omit(testing)
  
  training <- dplyr::select(training, -ID)
  ##Maxent
  Sys.setenv(NOAWT=TRUE)
  Maxent.Model <- predicts::MaxEnt(x = as.data.frame(training[,-1]), p = training[,1]) #Train Maxent Model
  Maxent.eval <- dismo::evaluate(p = testing[testing[,"pb"]==1,-1], a = testing[testing[,"pb"]==0,-1], model = Maxent.Model)
  
  terra::predict(
    PCA,
    Maxent.Model,
    filename = flnm,
    overwrite = TRUE,
    ncores=6
  )
  output <- list(MaxentModel=Maxent.Model, eval=Maxent.eval, flnm=flnm)
  return(output)
}

#Ok, let's apply this across a bunch of species

files <- list.files("data/GBIF/occs")
bigfiles <- files[grepl("panel", files)==TRUE]

dfmap <- data.frame(bigfile=bigfiles, sp=sub("^(([^_]*_){1}[^_]*).*", "\\1", bigfiles))

Occ <- NULL
for(i in 1:length(unique(dfmap$sp))){
  map <- dplyr::filter(dfmap, sp==unique(dfmap$sp[i]))
  tifname <- paste("suitabilityOutputs/", map$sp[1], ".tif", sep="")
  for(j in 1:nrow(map)){
    Occ <- rbind(Occ, read.csv(file=paste("data/GBIF/occs/", map$bigfile[j], sep="")))
    print(j)
}
  output <- trainMaxent(Occ, flnm=tifname)
  sdmname <- gsub(".csv", ".RDA", nm) %>% paste("SdmOut/",., sep="")
  save(output, file=sdmname) 
  print(paste(i, "complete of ", length(files)))
}
