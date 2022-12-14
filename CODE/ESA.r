

#### Figures for ESA

source("conditiionalPredictionExamples.R")
library(gjam)
library(tidyverse)
load("X:/conditionalPrediction/OUT/gjamOutput.rdata")


#### study area map ----

load("../DATA/covar.rdata")
covar <- covar[which(covar$idYear %in% out$inputs$xdata$idYear),]

bcr <- sf::st_read("../DATA/BCR/bcr_terrestrial_shape/BCR_Terrestrial_master.shp")
bcr <- bcr %>%
  filter(BCR %in% c(12, 22, 23),
         COUNTRY == "USA")

bcr <- bcr %>%
  mutate(col = case_when(BCR %in% c(12, 22, 23) ~ "yes",
                         T ~ "no")) %>%
  filter(COUNTRY == "USA")

map1 <- data.frame(lat = covar$Latitude,
                   lon = covar$Longitude,
                   true = out$inputs$y[,"brownheadedcowbird"],
                   year = out$inputs$xdata$year1) %>%
  #filter(year == 2019) %>%
  pivot_longer(cols = c(true))


bcr1 <- bcr %>%
  group_by(PROVINCE_S) %>%
  summarize(geometry = sf::st_union(geometry)) 
bcr2 <- bcr %>%
  filter(BCR %in% c(12, 22, 23)) %>%
  #group_by(BCR) %>%
  summarize(geometry = sf::st_union(geometry))

ggplot(map1) +
  geom_sf(data = bcr1, fill = "darkolivegreen4") +
  geom_sf(data = bcr2, fill = NA, size = 1) +
  geom_point(aes(x = lon, y = lat)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title = element_blank())
  


#### conditional prediction ----

# predict  traditional and conditional

true <- as.data.frame(out$inputs$y)
trad_pred <- as.data.frame(out$prediction$ypredMu)
cond_pred <- data.frame(obs = 1:nrow(true))

trad_rmspe <- as.data.frame(out$fit$rmspeBySpec)
cond_rmspe <- data.frame(sp = colnames(trad_pred),
                         rmspe = NA)
                         
ydata <- out$inputs$y
for (s in 1:ncol(out$inputs$y)) {
  
  condCols <- 1:17
  condCols <- condCols[-s]
  pCols <- c(1:ncol(ydata))[-condCols]
  
  new <- list(ydataCond = ydata[,condCols], nsim = 200)
  cond_pred2 <- gjamPredict(output = out, newdata = new)

  cond_pred[,1+s] <- cond_pred2$sdList$yMu[,s]
  

}



colnames(cond_pred) <- c("obs", colnames(trad_pred))
trad_pred$obs <- cond_pred$obs
true$obs <- cond_pred$obs

cond_pred_long <- cond_pred %>%
  pivot_longer(!obs) %>%
  mutate(type = "cond")
trad_pred_long <- trad_pred %>%
  pivot_longer(!obs) %>%
  mutate(type = "trad")
true_long <- true %>%
  pivot_longer(!obs) %>%
  mutate(type = "true")

all <- bind_rows(cond_pred_long, trad_pred_long, true_long) %>%
  pivot_wider(names_from = type, values_from = value) %>%
  filter(name %in% c("chippingsparrow", "indigobunting", "redeyedvireo", "songsparrow")) %>%
  pivot_longer(cols = c(trad, cond), names_to = "pred")

brks <- seq(-1, max(all$value)+1, by = 10)
all <- all %>%
  mutate(bins = cut(true,
                    breaks = brks)) %>%
  mutate(bin = round(true, digits = -1))

all$pred <- factor(all$pred, levels = c("trad", "cond"))
levels(all$pred) <- c("Traditional", "Conditional")
all$name <- factor(all$name, levels = c("chippingsparrow", "indigobunting", "songsparrow", "redeyedvireo"))
levels(all$name) <- c("Chipping sparrow", "Indigo Bunting", "Song Sparrow", "Red-eyed Vireo")

ggplot(all) +
  geom_hline(yintercept = 0) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed", color = "red", size = 1) +
  geom_boxplot(aes(x = bin, y = value, group = bin), fill = "gray80") +
  facet_grid(rows = vars(pred), cols = vars(name), scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  labs(x = "Observed", y = "Predicted") +
  coord_cartesian(ylim = c(0, 80), xlim = c(-3, 83))

ggsave("../OUT/ESA/predObs.jpg", height = 5, width = 10)


output <- out

ydata <- output$inputs$y
cnames <- colnames(ydata)
rmspe <- numeric(0)
c <- data.frame(cond = cnames)

# goes through each species and predicts conditionally on all others
for(i in 1:dim(ydata)[2]){ 
  
  yc <- ydata[ , !(colnames(ydata) %in% cnames[i]), drop = F] # Loop through each species for conditional predictions
  newdata <- list(ydataCond = yc,  nsim=200)
  
  comp <- conditionalComparison( output, newdata )
  ivals <- comp[[1]]
  rmspe <- rbind(rmspe, signif(ivals, 3) )
  
  c1 <- as.data.frame(comp[[2]]) %>%
    rownames_to_column(var = "cond")
  colnames(c1)[2] <- cnames[i]
  
  c <- c %>%
    full_join(c1, by = "cond")
  
  
  #print(rmspe)
}

rmspe <- as.data.frame(rmspe)
rmspe$abundance <- 100-rmspe$percZero

ggplot(rmspe) +
  geom_point(aes(x = percZero, y = `Perc Diff`))

ggplot(rmspe) +
  geom_point(aes(x = abundance, y = (`Frac u > c`)*100)) +
  geom_smooth(aes(x = abundance, y = (`Frac u > c`)*100), method = "lm", se = F) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "% of observations", y = "Fraction of observations where prediction is improved by conditioning")
cor(rmspe$abundance, rmspe$`Frac u > c`)
ggsave(file = "../OUT/ESA/abundanceImprovement.jpg", height = 4, width = 4)


##### prediction map
bcr <- sf::st_read("../DATA/BCR/bcr_terrestrial_shape/BCR_Terrestrial_master.shp")
bcr <- bcr %>%
  filter(BCR %in% c(12, 22, 23),
         COUNTRY == "USA")

bcr <- bcr %>%
  mutate(col = case_when(BCR %in% c(12, 22, 23) ~ "yes",
                         T ~ "no")) %>%
  filter(COUNTRY == "USA")
bcr1 <- bcr %>%
  group_by(PROVINCE_S) %>%
  summarize(geometry = sf::st_union(geometry)) 
bcr2 <- bcr %>%
  filter(BCR %in% c(12, 22, 23)) %>%
  #group_by(BCR) %>%
  summarize(geometry = sf::st_union(geometry))

loc <- covar %>%
  ungroup() %>%
  select(Latitude, Longitude, year1) %>%
  mutate(obs = 1:nrow(out$inputs$xdata))
  
all <- bind_rows(cond_pred_long, trad_pred_long, true_long) %>%
  pivot_wider(names_from = type, values_from = value) %>%
  filter(name %in% c("chippingsparrow")) %>%
  inner_join(loc, by = "obs") %>%
  filter(year1 == 2019) %>%
  pivot_longer(cols = c(trad, cond, true), names_to = "pred")

all$pred <- factor(all$pred, levels = c("true", "trad", "cond"))
levels(all$pred) <- c("True", 'Traditional', "Conditional")

cols <- c("aquamarine1", 
          "blue3",
          "magenta2", "firebrick1")

cols <- c("blue3", 
          "aquamarine1",
          "magenta2", "firebrick1")

cols <- c("blue3", 
          "purple",
          "orange", "brown")


box <- data.frame(lat = c(38, 38, 42.7, 42.7),
                  lon = c(-91.5, -87.4, -87.4, -91.5)) %>%
  sf::st_as_sf(coords = c("lon", "lat"),
               crs = sf::st_crs(bcr)) %>%
  summarise(geometry = sf::st_combine(geometry)) %>%
  sf::st_cast("POLYGON")

ggplot(all) +
  geom_sf(data = bcr1, fill = "gray95") +
  geom_sf(data = bcr2, fill = NA, size = 1) +
  geom_point(aes(x = Longitude, y = Latitude, color = sqrt(value))) +
  facet_wrap(~pred) +
  scale_color_gradientn(colors = cols,
                        values = c(0, 0.45, 0.6, 1)) +
  coord_sf() +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "", y = "", color = "Abundance")
ggsave(file = "../OUT/ESA/map.jpg", height = 6, width = 10)

ggplot(all) +
  geom_sf(data = bcr1, fill = "gray95") +
  geom_sf(data = bcr2, fill = NA, size = 1) +
  geom_point(aes(x = Longitude, y = Latitude, color = sqrt(value))) +
  geom_sf(data = box, color = "black", fill = NA,
          size = 1) +
  facet_wrap(~pred) +
  scale_color_gradientn(colors = cols,
                        values = c(0, 0.45, 0.6, 1)) +
  coord_sf() +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "", y = "", color = "Abundance")
ggsave(file = "../OUT/ESA/map-box.jpg", height = 6, width = 10)

ggplot(all) +
  geom_sf(data = bcr1, fill = "gray95") +
  geom_sf(data = bcr2, fill = NA, size = 1) +
  geom_point(aes(x = Longitude, y = Latitude, color = sqrt(value)), size = 5) +
  facet_wrap(~pred) +
  scale_color_gradientn(colors = cols,
                        values = c(0, 0.45, 0.6, 1)) +
  coord_sf(ylim = c(38, 42.7), xlim = c(-91.3, -87.5)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "", y = "", color = "Abundance")


ggsave(file = "../OUT/ESA/map-zoom.jpg", height = 6, width = 10)




#### make plot for A
# predict cowbird abundance conditionally on all other species
ydata <- out$inputs$y
conOn <- colnames(ydata)
conOn <- conOn[-grep('brownheadedcowbird', conOn)]
tmp <- gjamConditionalParameters(out, conditionOn = conOn)

A <- tmp$Atab %>%
  mutate(species = conOn)

ggplot(A) +
  geom_abline(intercept = 0, slope = 0) +
  geom_point(aes(x = reorder(species, Estimate), y = Estimate, color = sig95)) +
  geom_segment(aes(x = reorder(species, Estimate), xend = reorder(species, Estimate),
                   y = CI_025, yend = CI_975, color = sig95)) +
  scale_color_manual(values = c("gray70", "darkblue"), labels = c("Not significant",
                                                                  "Significant"),
                     name = "") +
  theme_bw() +
  theme(axis.text = element_text(angle = 45,
                                 hjust = 1,
                                 vjust = 0.5)) +
  labs(x = "", y = "")



ggsave(file = "../OUT/ESA/Cmatrix.jpg", height = 3.5, width = 8)








# gjamconditionalparameters
ydata <- out$inputs$y

yc <- ydata[ , !colnames(ydata) %in% c("brownheadedcowbird")] # condition on cowbird
newdata <- list(ydataCond = yc,  nsim=200)

comp <- conditionalComparison( out, newdata )



c1 <- as.data.frame(comp[[2]]) 
c1$species <- rownames(c1)
colnames(c1)[1] <- "brownheadedcowbird"

ggplot(c1) +
  geom_hline(yintercept = 0, color = "gray20") +
  geom_point(aes(x = reorder(species, brownheadedcowbird), 
                 y = brownheadedcowbird)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5)) +
  labs(x = "Species", y = "Association with cowbird")



#### this code should do the same thing but with uncertainty
yc <- ydata[ , !colnames(ydata) %in% c("brownheadedcowbird")] # condition on cowbird
newdata <- list(ydataCond = yc,  nsim=200)


ydataCond <- newdata$ydataCond
ni <- match(colnames(ydataCond), colnames(ydata))   # condition on these columns
ci <- c(1:ncol(ydata))[-ni]


sigma <- output$parameters$sigMu
sigmaSe <- output$parameters$sigSe
sigmaStd <- sigmaSe * sqrt(nrow(output$inputs$y))

nsim <- 20000
C <- matrix(nrow = 16, ncol = nsim)
for (n in 1:nsim) {
  sigma1 <- matrix(mapply(rnorm, n = 1, mean = sigma, sd = sigmaStd), 
                   ncol = 17, nrow = 17)

  C[,n]  <- solve(sigma1[ni,ni])%*%sigma1[ni,ci] 
}

# mean C
cmean <- rowMeans(C)
c975 <- apply(C, 1, quantile, 0.975)
c025 <- apply(C, 1, quantile, 0.025)
sp <- colnames(ydata)[ni]

ggplot() +
  geom_hline(yintercept = 0, color = "gray20") +
  geom_point(aes(x = reorder(sp, cmean), y = cmean)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5)) +
  geom_segment(aes(x = reorder(sp, cmean), xend = reorder(sp, cmean),
                   y = c025, yend = c975))






