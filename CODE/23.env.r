

# Download environmental covariates
# C. Lane Scher
# clanescher@gmail.com

library(tidyverse)

load("DATA/birds/counts.rdata")
load("DATA/birds/nlcd.rdata")


## get locations ----
useRoutes <- counts %>%
  ungroup() %>%
  mutate(id = gsub("-.*", "", idYear)) %>%
  select(id) %>%
  distinct()

locs <- read.csv("DATA/birds/BBS/routes.csv") %>%
  filter(BCR %in% c(12, 22, 23)) %>%
  mutate(id = paste0(CountryNum, "_", StateNum, "_", Route)) %>%
  filter(id %in% useRoutes$id) %>%
  select(id, Latitude, Longitude)


## get climate ----

# download
done <- list.files("DATA/birds/")

if ("env.rdata" %in% done) {
  load("DATA/birds/env.rdata")
} else {
  write.csv(locs, file = "../DATA/birds/locs.csv", row.names = F)  
  
  yearStart <- 2015
  yearEnd <- 2020
  
  tmp <- daymetr::download_daymet_batch(file_location = "../DATA/birds/locs.csv",
                                        start = yearStart,
                                        end = yearEnd,
                                        internal = T,
                                        path = "../DATA/birds/envRAW/")
  save(tmp, file = "../DATA/birds/env.rdata")
}


# process 
env <- c()
for (l in 1:length(tmp)) {
  data <- tmp[[l]]$data%>%
    mutate(#date = lubridate::as_date(paste(data$yday, data$year, sep = "/"),
           #                          format = "%j/%Y"),
           season = case_when(yday <= 82 ~ "winter",
                              yday > 82 & 
                                yday <= 174 ~ "spring",
                              yday > 174 &
                                yday <= 266 ~ "summer",
                              yday > 266 &
                                yday <= 357 ~ "fall",
                              yday > 357 &
                                yday <= 365 ~ "winter"))
  data$year1 <- data$year
  data$year1[data$yday > 213] <- data$year[data$yday > 213] + 1
  
  env1 <- data %>%
    group_by(year1) %>%
    summarize(prcp = sum(prcp..mm.day.),
           tmin = mean(tmin..deg.c.)) %>%
    ungroup() %>%
    mutate(prcpMean = mean(prcp),
           tminMean = mean(tmin),
           prcpAnom = prcp - prcpMean,
           tminAnom = tmin - tminMean,
           id = tmp[[l]]$site)
  
  env <- bind_rows(env, env1)
  
}


# get elevation ----

if ("elev.rdata" %in% done) {
  load("DATA/birds/elev.rdata")
} else {
  locs1 <- locs %>%
    sf::st_as_sf(coords = c("Longitude", "Latitude"),
                 crs = "+proj=longlat +datum=WGS84")
  
  
  tmp <- elevatr::get_elev_point(locs1,
                                 prj="+proj=longlat +datum=WGS84")
  
  elevation <- tmp %>%
    as.data.frame() %>%
    select(id, elevation)
  
  save(elevation, file = "DATA/birds/elev.rdata")
}





## combine and save ----

covar <- env %>%
  mutate(idYear = paste0(id, "-", year1)) %>%
  inner_join(locs, by = "id") %>%
  inner_join(elevation, by = "id") %>%
  inner_join(nlcd, by = "id") %>%
  mutate(water = `11`,
         developed = `21` + `22` + `23` + `24`,
         barren = `31`,
         forest = `41` + `42` + `43`,
         shrub = `52`,
         herbaceous = `71`,
         planted = `81` + `82`,
         wetlands = `95`)


save(covar, file = "DATA/birds/covar.rdata")
