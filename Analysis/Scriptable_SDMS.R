library(dplyr)
library(tidyr)
library(geodata)
occs <- read.csv(file="data/trimmed_occurances.csv")[,-1]
load(file="data/TestPoints.Rda")

Clima <- geodata::worldclim_global("worldclim", var="bio", res=2.5)

sp <- unique(c(testpoints$Animal_Species, testpoints$Plant_Species))
spfiles <- tolower(sp) %>% gsub(" ", "_", .) %>% paste(., ".csv", sep="")
#table(spfiles %in% occs$file)
#table(sp %in% occs$species)

#missing <- sort(spfiles[spfiles %in% occs$file==F])

psuedos <- NULL
splist <- spfiles[spfiles %in% occs$file]
for(i in 99:length(splist)){
  temp <- occs[occs$file==splist[i],]
  if(nrow(temp)>10000){ #If too many points, sample down to 10k
    temp <- temp[sample(1:nrow(temp), 10000, replace=F),]
  }
  
  temp['cells']<-terra::cellFromXY(object= Clima[[1]], xy = temp[,5:4]) 
  
  presences <- terra::vect(as.matrix(temp[,5:4]), crs=crs(Clima))
  v_buf <- terra::buffer(presences, width=200000)
  if(i==98){ #Need to expand bounds for endemic puerto rican bird
  v_buf <- terra::buffer(presences, width=400000)
  }
  
  bg <- Clima[[1]]
  values(bg)[!is.na(values(bg))] <- 1
  
  region_buf <- terra::mask(bg, v_buf)
  
  
  # Exclude presence locations:
  sp_cells <- terra::extract(region_buf, presences, cells=T)$cell
  region_buf_exclp <- region_buf
  values(region_buf_exclp)[sp_cells] <- NA

  bg_rand_buf <- terra::spatSample(region_buf_exclp, nrow(temp), "random", na.rm=T, as.points=TRUE, exhaustive=T)
  temp['psudeocells'] <- terra::extract(region_buf, bg_rand_buf, cells=T)$cell
  

  ##plot(Clima[[1]])
  ##points(terra::xyFromCell(Clima[[1]], temp$psudeocells))
  
  
  ##plot(bg_rand_buf)
  #creating a 200km2 buffer to select the pseudo-absence points
  
  ##v_buf <- terra::buffer(pop_imag, width=200000)
  ##terra::vect(temp, coords = c("decimalLongitude",'decimalLatitude'), 
              ##crs=sp::CRS("EPSG:4326"))
  
  ##OccSf <- sf::st_as_sf(x=temp, coords = c("decimalLongitude",'decimalLatitude'), 
  ##crs=sp::CRS("EPSG:4326")) #crs=sp::CRS("+proj=longlat +ellps=WGS84
  
##OccBuffer <-sf::st_buffer(OccSf, dist=2) #Create buffer objects
##OccBuffer <-sf::st_union(OccBuffer) #Merge buffer objects into one
  #plot(OccBuffer)
  
  ##Masking the buffer so that it does not occur in the ocean
##OccBuffer2 <-as(OccBuffer, 'Spatial')
  #terra::plot(OccBuffer2)
  ##OccBuffer2 <- terra::vect(OccBuffer)
  #OccBufferRas <-terra::rasterize(x = OccBuffer2, y=Clima)
  #OccBufferRas<-raster::mask(x = OccBufferRas, mask = Clima[[1]])
  #plot(OccBufferRas)
  
  ##Selecting the pseudo absences. Since I'm doing 1:1 Presence:Psuedoabsence, I'm saving the cell number in the df as a seperate column. 
  ##AllPseudo <- terra::as.points(OccBuffer2)
  ##AllPseudo <- raster::extract(x= Clima[[1]], terra::as.points(OccBuffer2),cell=T)
  ##AllPseudo2 <-  AllPseudo[-which(AllPseudo$cells %in% temp$cells),]

  ##temp$psuedocell <- sample(AllPseudo$cell, size=nrow(temp), replace=F)
  psuedos <- rbind(psuedos, temp)
  if(i%%10==0){
    print(i)
    write.csv(psuedos, file="data/psuedos.csv")
  }
}

#Now we have all of our data, let's perform our PCA. 

