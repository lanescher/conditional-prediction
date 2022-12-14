

# Get BBS data for bird application
# C. Lane Scher
# clanescher@gmail.com

years <- 2010:2020
library(tidyverse)


#### read in raw BBS data ####
dpath <- 'X:/conditionalPrediction/DATA/birds/BBS/'
files <- list.files(dpath, pattern = "fifty")

data <- read.csv(paste0(dpath, "/", files[1]))
colnames(data) <- tolower(colnames(data))

for (j in 2:length(files)) {
  print(j)
  tmp <- read.csv(paste0(dpath, "/", files[j]))
  colnames(tmp) <- tolower(colnames(tmp))
  data <- rbind(data, tmp)
}

weather <- read.csv("../DATA/birds/BBS/weather.csv") %>%
  filter(RunType == 1) %>%
  mutate(idYear = paste0(CountryNum, "_", StateNum, "_", Route, "-", Year))

#### read in AOU codes ####
spList <- read.fwf(file = "X:/conditionalPrediction/DATA/birds/BBS/SpeciesList.txt",
                   widths = c(7, 5, 50, 50, 50, 
                              50, 50, 50, 50),
                   skip = 11,
                   col.names = c("Seq", "AOU", "English",
                                 "French", "Spanish",
                                 "Order", "Family", "Genus", "Species"),
                   strip.white = T)

spList$gs <- paste(spList$Genus, spList$Species, sep = " ")
spList <- spList[,c("AOU", "gs")]
colnames(spList)[1] <- "aou"

spUse <- c("Molothrus ater",
           "Setophaga petechia", "Melospiza melodia",
           "Vireo olivaceus", "Spizella passerina",
           "Sayornis phoebe", "Pipilo erythrophthalmus",
           "Seiurus aurocapilla", "Geothlypis trichas",
           "Setophaga ruticilla", "Passerina cyanea",
           "Icteria virens", "Agelaius phoeniceus",
           "Oporornis formosus", "Empidonax trailli",
           "Vireo bellii", "Vireo flavifrons",
           "Spizella pusilla", "Setophaga kirtlandii")

# format for gjam
dat <- data %>%
  mutate(count = rowSums(data[,8:57]),
         idYear = paste0(countrynum, "_", statenum, "_", route, "-", year)) %>%
  filter(year %in% years) %>%
  inner_join(spList, by = "aou") %>%
  select(idYear, gs, count) %>%
  filter(gs %in% spUse) %>%
  group_by(idYear, gs) %>%
  summarize(count = sum(count)) %>%
  pivot_wider(names_from = gs, values_from = count) %>%
  filter(idYear %in% weather$idYear)
dat[is.na(dat)] <- 0

colnames(dat) <- c("idYear",
                   "red-winged blackbird", "common yellowthroat",
                   "song sparrow", "american yellow warbler",
                   "chipping sparrow", "brown-headed cowbird",
                   "red-eyed vireo", "american redstart",
                   "eastern phoebe", "ovenbird",
                   "yellow-breasted chat", "yellow-throated vireo",
                   "field sparrow", "indigo bunting",
                   "eastern towhee", "bell's vireo",
                   "kirtland's warbler")
counts <- dat

save(counts, file = "DATA/birds/counts.rdata")
