
# Download NLCD
# C. Lane Scher
# clanescher@gmail.com


devtools::install_github("ropensci/FedData")

library(tidyverse)
library(sf)

source("nlcdFunctions.r")

# load counts
load("DATA/birds/counts.rdata")

# load state codes
stateCodes <- read.csv("DATA/birds/BBS/bbsStateCodes.csv")
colnames(stateCodes)[1] <- "StateNum"



## get locations ----

# select routes that we have counts for during this time period
useRoutes <- counts %>%
  ungroup() %>%
  mutate(id = gsub("-.*", "", idYear)) %>%
  select(id) %>%
  distinct()

# read in route locations and select the ones we want
locs <- read.csv("DATA/birds/BBS/routes.csv") %>%
  filter(BCR %in% c(12, 22, 23),
         CountryNum == 840) %>%
  mutate(id = paste0(CountryNum, "_", StateNum, "_", Route)) %>%
  filter(id %in% useRoutes$id) %>%
  inner_join(stateCodes, by = "StateNum")

# downloaded from https://www.birdscanada.org/bird-science/nabci-bird-conservation-regions/
bcr <- sf::st_read("DATA/birds/BCR/bcr_terrestrial_shape/BCR_Terrestrial_master.shp")

bcrUse <- bcr %>%
  filter(BCR %in% c(12, 22, 23),
         COUNTRY == "USA")


## download ----

yearsNLCD <- c("2016")
redoNLCD <- F

nlcdAll <- c()
for (y in 1:length(yearsNLCD)) {
  year <- yearsNLCD[y]
  
  for (s in 1:nrow(bcrUse)) {
    st <- bcrUse$BCR[s]
    
    cat("\n")
    
    print(Sys.time())
    
    print(paste("Shape ", s, " out of ", nrow(bcrUse)))
    print(paste0("BCR ", bcrUse$BCR[s]))
    print(as.character(bcrUse$PROVINCE_S[s]))
    
    
    utmZone <- 16
    crs1 <- paste0("+proj=utm +zone=", utmZone," +ellps=WGS84 +datum=WGS84 +units=m +no_defs ")
    
    shp <- bcrUse[s,]
    # st_crs(shp)
    # shp <- st_transform(shp, st_crs(crs1))
    
    done <- list.files("DATA/birds/nlcdShapes/")
    thisShape <- paste0(s, "-",
                        bcrUse$BCR[s], "-",
                        bcrUse$PROVINCE_S[s])
    
    
    # all files for this shape
    num <- grep(thisShape, done)
    
    # raster for this shape
    num1 <- grep(paste0("raw-", thisShape), done)
    
    # processed nlcd for this shape
    num2 <- setdiff(num, num1)
    
    # if redoNLCD is T, reset num2 to null
    if (redoNLCD == T) {
      num2 <- NULL
      print("Recalculating NLCD")
    }
    
    if (s == 1) {num <- NULL}
    
    # all done
    if(length(num1) >= 1 &
       length(num2) >= 1) {
      print("Loading nlcd")
      
      toLoad <- done[num2]
      
      for (n in 1:length(toLoad)) {
        load(paste0("DATA/birds/nlcdShapes/", toLoad[n]))
        nlcdAll <- rbind(nlcdAll, nlcdShape1)
        
      }
      
      print(paste0(nrow(nlcdAll), " out of ", nrow(locs), " done"))
      
      #next
    }
    
    
    
    
    # raster(s) downloaded, just calculate nlcd
    if(length(num1) >= 1 &
       length(num2) == 0) {
      print("Loading raster")
      
      toLoad <- done[num1]
      
      for (n in 1:length(toLoad)) { # for each raster, load and calculate nlcd
        load(paste0("DATA/birds/nlcdShapes/", toLoad))
        
        # assign crs to rts
        locs1 <- locs %>%
          filter(BCR == bcrUse$BCR[s],
                 state == bcrUse$PROVINCE_S[s]) %>%
          select(id, Longitude, Latitude) %>%
          st_as_sf(coords = c("Longitude", "Latitude"),
                   crs = crs1) %>%
          st_buffer(.1) # approx 10km
        
        #locs2 <- st_intersection(locs1, shp)
        if(nrow(locs1) == 0) next
        
        # assign crs to raster
        raster::projection(tmp2) <- crs1
        
        print("Calculating NLCD")
        
        nlcdShape <- c()
        pb <- txtProgressBar(min = 0, max = nrow(locs1), style = 3)
        
        for (r in 1:nrow(locs1)) {
          
          loc <- locs1[r,]
          
          
          # get nlcd for route
          nlcdRoutes <- freqNLCD(loc, tmp2, buffer = 0, max = T)
          if(length(nlcdRoutes) == 0) next
          nlcdRoutes$id <- loc$id
          
          
          nlcdShape <- rbind(nlcdShape, nlcdRoutes)
          
          setTxtProgressBar(pb, r)
          
        }
        cat("\n")
        
        nlcdShape1 <- nlcdShape[which(rowSums(nlcdShape[,1:15]) > 0),]
        save(nlcdShape1, file = paste0("DATA/birds/nlcdShapes/", s, "-",
                                       bcrUse$BCR[s], "-",
                                       bcrUse$PROVINCE_S[s], "-", n, ".rdata"))
        
        nlcdAll <- rbind(nlcdAll, nlcdShape1)
        
        print(paste0(nrow(nlcdAll), " out of ", nrow(locs), " done"))
      }
      
    }
    
    
    
    
    # need to download raster and calculate nlcd
    if(length(num) == 0) {
      
      if (bcrUse$Shape_Area[[s]] <= 5) {
        spl <- list(shp) 
      }
      
      if (bcrUse$Shape_Area[[s]] > 5 &
          bcrUse$Shape_Area[[s]] <= 9) {
        print("Too big!! Splitting!!")
        
        spl <- splitPoly(shp)
      } 
      
      if (bcrUse$Shape_Area[[s]] > 9) {
        print("Way too big!! Splitting twice!!")
        
        spl1 <- splitPoly(shp)
        spl2 <- splitPoly(spl1[[1]])
        spl3 <- splitPoly(spl1[[2]])
        
        spl <- c(spl2, spl3)
        
      }
      
      for (c in 1:length(spl)) {
        shp1 <- spl[[c]]
        
        print(paste0("Downloading raster approx half of size ", bcrUse$Shape_Area[[s]]))
        tmp1 <- FedData::get_nlcd(template = shp1,
                                  label = st,
                                  year = year,
                                  dataset = "landcover",
                                  force.redo = T,
                                  extraction.dir = paste0("DATA/birds/stateNLCD/",
                                                          thisShape))
        
        
        print("Projecting raster")
        tmp2 <- raster::projectRaster(tmp1, crs="+proj=longlat +datum=WGS84",
                                      method = "ngb")
        save(tmp2, file = paste0("DATA/birds/nlcdShapes/raw-",
                                 thisShape, "-", c, ".rdata"))
        
        # assign crs to rts
        locs1 <- locs %>%
          filter(BCR == bcrUse$BCR[s],
                 state == bcrUse$PROVINCE_S[s]) %>%
          select(id, Longitude, Latitude) %>%
          st_as_sf(coords = c("Longitude", "Latitude"),
                   crs = crs1) %>%
          st_buffer(.1) # approx 10km
        
        #locs2 <- st_intersection(locs1, shp)
        if(nrow(locs1) == 0) next
        
        # assign crs to raster
        raster::projection(tmp2) <- crs1
        
        print("Calculating NLCD")
        
        nlcdShape <- c()
        pb <- txtProgressBar(min = 0, max = nrow(locs1), style = 3)
        
        for (r in 1:nrow(locs1)) {
          
          loc <- locs1[r,]
          
          
          # get nlcd for route
          nlcdRoutes <- freqNLCD(loc, tmp2, buffer = 0, max = T)
          if(length(nlcdRoutes) == 0) next
          nlcdRoutes$id <- loc$id
          
          
          nlcdShape <- rbind(nlcdShape, nlcdRoutes)
          
          setTxtProgressBar(pb, r)
          
        }
        cat("\n")
        
        nlcdShape1 <- nlcdShape[which(rowSums(nlcdShape[,1:15]) > 0),]
        save(nlcdShape1, file = paste0("DATA/birds/nlcdShapes/", s, "-",
                                       bcrUse$BCR[s], "-",
                                       bcrUse$PROVINCE_S[s], "-", c, ".rdata"))
        
        nlcdAll <- rbind(nlcdAll, nlcdShape1)
        
        print(paste0(nrow(nlcdAll), " out of ", nrow(locs), " done"))
      }
      
      
    }
    
   
  }
  
  

}

nlcdAll$lc <- NULL

nlcd <- nlcdAll %>%
  group_by(id) %>%
  summarize_all(sum) #%>%

nlcd <- nlcd %>%
  mutate(across(2:16, ~ .x/rowSums(nlcd[,2:16])))
nlcd$lc <- colnames(nlcd)[max.col(nlcd[,2:16], ties.method = "random")]

save(nlcd, file = "DATA/birds/nlcd.rdata")


