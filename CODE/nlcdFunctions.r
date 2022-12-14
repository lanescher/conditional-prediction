
# NLCD functions!!!!
#### original

freqNLCD <- function(points, raster, buffer, max = F) {
  #polyBuffered <- rgeos::gBuffer(points, width = buffer, byid = T)
  #tmp1 <- raster::extract(raster, polyBuffered)
  
  tmp1 <- raster::extract(raster, points)
  
  ind <- lapply(tmp1, function(x) {which(length(x) == 0)})
  tmp1[ind == 1] <- -99
  
  #freq <- lapply(tmp1, function(x) {table(x)/length(x)})
  freq <- lapply(tmp1, function(x) {table(x)})
  
  if(length(freq[[1]]) != 0) {
    allFreq <- as.data.frame(dplyr::bind_rows(freq))
    
    lcNums <- c("11", "21", "22", "23", "24", "31", "41", "42", "43",
                "52", "71", "81", "82", "90", "95", "-99")
    lcNumsMissing <- lcNums[-which(lcNums %in% colnames(allFreq))]
    if (length(lcNumsMissing > 0)) {
      for (i in 1:length(lcNumsMissing)){
        allFreq[lcNumsMissing[i]] <- 0
      }
    }
    
    
    allFreq <- allFreq[,c("11", "21", "22", "23", "24", "31", "41", "42", "43",
                          "52", "71", "81", "82", "90", "95", "-99")]
    
    allFreq[is.na(allFreq)] <- 0
    
    allFreq[which(allFreq$`-99` == 1),] <- NA
    allFreq$`-99` <- NULL
    
    if (max == T) {
      allFreq$lc <- colnames(allFreq)[max.col(allFreq, ties.method = "random")]
    }
    
    return(allFreq)
  }
  
  
}



splitPoly <- function(polygon) {
  xmin <- st_bbox(polygon)[1]
  xmax <- st_bbox(polygon)[3]
  
  ymin <- st_bbox(polygon)[2]
  ymax <- st_bbox(polygon)[4]
  
  xrange <- xmax - xmin
  yrange <- ymax - ymin
  
  if(xrange > yrange) {
    split <- xmin + 0.5*xrange
    
    bbox1 <- c(xmin-0.001, ymin, split, ymax)
    names(bbox1) <- c("xmin", "ymin", "xmax", "ymax")
    bbox1 <- st_bbox(bbox1)
    
    bbox2 <- c(split, ymin, xmax, ymax)
    names(bbox2) <- c("xmin", "ymin", "xmax", "ymax")
    bbox2 <- st_bbox(bbox2)
    
    a <- st_crop(polygon, bbox1)
    b <- st_crop(polygon, bbox2)
    splits <- list(a, b)
  }
  if(xrange < yrange) {
    split <- ymin + 0.5*yrange
    
    bbox1 <- c(xmin, ymin, xmax, split)
    names(bbox1) <- c("xmin", "ymin", "xmax", "ymax")
    bbox1 <- st_bbox(bbox1)
    
    bbox2 <- c(xmin, split, xmax, ymax)
    names(bbox2) <- c("xmin", "ymin", "xmax", "ymax")
    bbox2 <- st_bbox(bbox2)
    
    a <- st_crop(polygon, bbox1)
    b <- st_crop(polygon, bbox2)
    splits <- list(a, b)
  }
  return(splits)
}
