
library(gjam)
library(tidyverse)


# 1. management insight from conditioning on cowbirds ----

load("OUT/birds-gjamOutput.rdata")
output <- out
ydata <- output$inputs$y

# condition cowbirds on hosts
conOn <- colnames(ydata)
conOn <- conOn[-grep('brownheadedcowbird', conOn)]

tmp <- gjamConditionalParameters(out, conditionOn = conOn, nsim = 20000)

A <- tmp$Atab %>%
  mutate(species = conOn,
         sigSE = case_when((Estimate + SE)*(Estimate - SE) > 0 ~ "*",
                           T ~ ""),
         type = "hostToCow")

A <- A %>%
  mutate(Species = case_when(species == "commonyellowthroat" ~ "Common Yellowthroat",
                             species == "songsparrow" ~ "Song Sparrow",
                             species == "chippingsparrow" ~ "Chipping Sparrow",
                             species == "americanyellowwarbler" ~ "American Yellow Warbler",
                             species == "redwingedblackbird" ~ "Red-winged Blackbird",
                             species == "redeyedvireo" ~ "Red-eyed Vireo",
                             species == "americanredstart" ~ "American Redstart",
                             species == "easternphoebe" ~ "Eastern Phoebe",
                             species == "ovenbird" ~ "Ovenbird",
                             species == "yellowbreastedchat" ~ "Yellow-breasted Chat",
                             species == "yellowthroatedvireo" ~ "Yellow-throated Vireo",
                             species == "fieldsparrow" ~ "Field Sparrow",
                             species == "easterntowhee" ~ "Eastern Towhee",
                             species == "indigobunting" ~ "Indigo Bunting",
                             species == "bellsvireo" ~ "Bell's Vireo",
                             species == "kirtlandswarbler" ~ "Kirtland's Warbler"))





# condition hosts on cowbird
conOn <- "brownheadedcowbird"

tmp <- gjamConditionalParameters(out, conditionOn = conOn, nsim = 20000)

sp <- colnames(ydata)
sp <- sp[-grep('brownheadedcowbird', sp)]

A1 <- tmp$Atab %>%
  mutate(species = sp,
         sigSE = case_when((Estimate + SE)*(Estimate - SE) > 0 ~ "*",
                           T ~ ""),
         type = "cowToHost") 

A1 <- A1 %>%
  mutate(Species = case_when(species == "commonyellowthroat" ~ "Common Yellowthroat",
                             species == "songsparrow" ~ "Song Sparrow",
                             species == "chippingsparrow" ~ "Chipping Sparrow",
                             species == "americanyellowwarbler" ~ "American Yellow Warbler",
                             species == "redwingedblackbird" ~ "Red-winged Blackbird",
                             species == "redeyedvireo" ~ "Red-eyed Vireo",
                             species == "americanredstart" ~ "American Redstart",
                             species == "easternphoebe" ~ "Eastern Phoebe",
                             species == "ovenbird" ~ "Ovenbird",
                             species == "yellowbreastedchat" ~ "Yellow-breasted Chat",
                             species == "yellowthroatedvireo" ~ "Yellow-throated Vireo",
                             species == "fieldsparrow" ~ "Field Sparrow",
                             species == "easterntowhee" ~ "Eastern Towhee",
                             species == "indigobunting" ~ "Indigo Bunting",
                             species == "bellsvireo" ~ "Bell's Vireo",
                             species == "kirtlandswarbler" ~ "Kirtland's Warbler"))


# now plot them, ordered by host to cowbird
ord <- A$Species[order(A$Estimate)]

A1$Species <- factor(A1$Species, levels = ord)
A$Species <- factor(A$Species, levels = ord)


figA <- ggplot(A) +
  geom_abline(intercept = 0, slope = 0) +
  geom_point(aes(x = Species, y = Estimate, color = sigSE)) +
  geom_segment(aes(x = Species, xend = Species,
                   y = Estimate+SE, yend = Estimate-SE, color = sigSE)) +
  scale_color_manual(values = c("gray70", "darkblue"), labels = c("Not significant",
                                                                  "Significant"),
                     name = "") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 9)) +
  labs(x = "", y = "Species \ncoefficient") +
  coord_cartesian(ylim = c(-6.9, 3.3))



figB <- ggplot(A1) +
  geom_abline(intercept = 0, slope = 0) +
  geom_point(aes(x = Species, y = Estimate, color = sigSE)) +
  geom_segment(aes(x = Species, xend = Species,
                   y = Estimate+SE, yend = Estimate-SE, color = sigSE)) +
  scale_color_manual(values = c("gray70", "darkblue"), labels = c("Not significant",
                                                                  "Significant"),
                     name = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5,
                                 size = 9),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 9)) +
  labs(x = "", y = "Species \ncoefficient") +
  coord_cartesian(ylim = c(-6.9, 3.3))


# combine and plot

ggpubr::ggarrange(figA, figB, nrow = 2,
                  labels = "auto",
                  common.legend = T, legend = "none",
                  heights = c(1, 1.8),
                  label.x = 0, label.y = 1.04,
                  font.label = list(size = 10))

ggsave(file = "../OUT/figures/Amatrix.png", height = 120, width = 80,
       units = "mm", dpi = 600)


## now just do a

ggplot(A) +
  geom_abline(intercept = 0, slope = 0) +
  geom_point(aes(x = Species, y = Estimate, color = sigSE)) +
  geom_segment(aes(x = Species, xend = Species,
                   y = Estimate+SE, yend = Estimate-SE, color = sigSE)) +
  scale_color_manual(values = c("gray70", "darkblue"), labels = c("Not significant",
                                                                  "Significant"),
                     name = "",
                     guide = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   size = 7),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7), 
        legend.text=element_text(size = 7)) +
  labs(x = "", y = "Species \ncoefficient") +
  coord_cartesian(ylim = c(-6.9, 3.3))

ggsave(file = "../OUT/figures/Amatrix.png", height = 60, width = 80,
       units = "mm", dpi = 600)



### now see if there are more or less cowbirds at focal species sites

#conditionOn <- "kirtlandswarbler"
#conditionOnData <- select(as.data.frame(ydata), all_of(conditionOn))
conditionOnData <- ydata[,-which(colnames(ydata) == "brownheadedcowbird")]

newdata <- list(ydataCond = conditionOnData,  nsim = 500)
cond_pred <- gjamPredict(output = out, newdata = newdata)
cond_pred1 <- as.data.frame(cond_pred$sdList$yMu)


cowKirtPres <- cond_pred1 %>%
  filter(kirtlandswarbler > 0) %>%
  mutate(group = "Kirtland's Warbler \npresent")
mean(cowKirtPres$brownheadedcowbird)

cowBellPres <- cond_pred1 %>%
  filter(bellsvireo > 0) %>%
  mutate(group = "Bell's Vireo \npresent")
mean(cowBellPres$brownheadedcowbird)

cowAll <- cond_pred1 %>%
  mutate(group = "All")
mean(cowAll$brownheadedcowbird)

all <- bind_rows(cowKirtPres, cowBellPres) %>%
  bind_rows(cowAll)

trad_pred <- cond_pred1 %>%
  mutate(`brownheadedcowbird` = out$prediction$ypredMu[,"brownheadedcowbird"])

cowKirtPres1 <- trad_pred %>%
  filter(kirtlandswarbler > 0) %>%
  mutate(group = "Kirtland's Wabler \npresent")
mean(cowKirtPres1$brownheadedcowbird)

cowBellPres1 <- trad_pred %>%
  filter(bellsvireo > 0) %>%
  mutate(group = "Bell's Vireo \npresent")
mean(cowBellPres1$brownheadedcowbird)

cowAll1 <- trad_pred %>%
  mutate(group = "All")
mean(cowAll1$brownheadedcowbird)

all1 <- bind_rows(cowKirtPres1, cowBellPres1) %>%
  bind_rows(cowAll1)

a <- ggplot(all1) +
  geom_boxplot(aes(x = group, y = brownheadedcowbird, fill = group)) +
  scale_fill_manual(values = c("gray", "brown", "yellow")) +
  theme_bw() +
  labs(x = "Observations", y = "Brown-headed Cowbird \nabundance", fill = "") +
  coord_cartesian(ylim = c(0, 45))

b <- ggplot(all) +
  geom_boxplot(aes(x = group, y = brownheadedcowbird, fill = group)) +
  scale_fill_manual(values = c("gray", "brown", "yellow")) +
  theme_bw() +
  labs(x = "Observations", y = "", fill = "") +
  coord_cartesian(ylim = c(0, 45))


ggpubr::ggarrange(a, b, nrow = 1, legend = "none",
                  labels = "auto", align = "hv")

ggsave(file = "../OUT/figures/predictedAbundance.jpg",
       height = 5, width = 10)



# 2. which/how many species are best to condition on ----
load("../OUT/birds-gjamOutput.rdata")
output <- out
ydata <- output$inputs$y
xdata <- output$inputs$xdata


species <- colnames(ydata)

trad <- out$prediction$ypredMu

# keep lists of condition improvement, who was conditioned on, 
# and who was predicted

# this takes a while to run, can load output below.
condImprovement <- list()
condOnList <- list()
predictSpList <- list()

for (s in 1:length(species)) {

  print(paste0("predicting species ", s, ", ", species[s]))
  pred <- species[s]
  sp <- species[-s]
  
  
  for (i in 1:length(sp)) {
    print(paste0("conditioning on ", i, " species"))
    
    tmp <- as.data.frame(combn(sp, i))
    
    # if there are more than 10 possible combinations, pick 10
    if (ncol(tmp) > 10) {
      tmp <- tmp[,sample(1:ncol(tmp), 10)]
    }
    
    # now go through those 10 (or fewer) combinations
    for (j in 1:ncol(tmp)) {
      print(paste0("combination # ", j))
      
      conditionOn <- tmp[,j]
      conditionOnData <- select(as.data.frame(ydata), all_of(conditionOn))
      
      
      newdata <- list(ydataCond = conditionOnData,  nsim = 500)
      cond_pred <- gjamPredict(output = out, newdata = newdata)
      cond <- cond_pred$sdList$yMu[,s]
      
      # compare conditional and unconditional
      #comp <- conditionalComparison( output, newdata )
      #comp <- comp[[1]]
      
      # get spe
      coni <- (ydata[,pred] - cond )^2
      unci <- (ydata[,pred] - trad[,s] )^2
      
      # how many spe are greater for traditional prediction
      umc  <- unci - coni
      uno  <- length(which(umc > 0))/length(umc)  # fraction larger error in uncon
      
      # get rmspe
      rmspecond <- sqrt(mean(coni))
      rmspetrad <- sqrt(mean(unci))
      
      
      # add to compare dataframe
      comp <- data.frame("Unc RMSPE" = rmspetrad, 
                         "Con RMSPE" = rmspecond, 
                         "Perc Diff" = (rmspetrad-rmspecond)/rmspetrad * 100, 
                         "Frac u > c" = uno)
      colnames(comp) <- c("Unc RMSPE", "Con RMSPE", "Perc Diff", "Frac u > c")
      
      
      
      condImprovement[[length(condImprovement) + 1]] <- comp
      condOnList[[length(condOnList) + 1]] <- conditionOn
      predictSpList[[length(predictSpList) + 1]] <- species[s]
      
    }
    
    
  }
}


#save(condImprovement, condOnList, predictSpList, file = "../OUT/birds-conditionOn.rdata")
load("../OUT/birds-conditionOn.rdata")


all <- c()
for (i in 1:length(condImprovement)) {
  thisOne <- as.data.frame(condImprovement[[i]])
  
  pred <- predictSpList[[i]]
  
  #thisOne <- thisOne[which(rownames(thisOne) %in% pred),]
  
  keep <- data.frame(species = pred,
                     fracImproved = (thisOne[,c("Frac u > c")])*100,
                     improvedBy = thisOne[,c("Perc Diff")],
                     predOn = length(condOnList[[i]]))
  
  all <- bind_rows(all, keep)
  
}


all_summary <- all %>%
  group_by(species, predOn) %>%
  summarize(fracMean = mean(fracImproved),
            fracSe = sd(fracImproved)/sqrt(10),
            byMean = mean(improvedBy),
            bySe = sd(improvedBy)/sqrt(10)) %>%
  mutate(fracHi = fracMean + fracSe,
         fracLo = fracMean - fracSe,
         byHi = byMean + bySe,
         byLo = byMean - bySe)

all <- all_summary %>%
  pivot_longer(cols = c("fracMean", "fracHi", "fracLo",
                                          "byMean", "byHi", "byLo")) %>%
  mutate(stat = gsub("frac|by", "", name),
         cat = gsub("Mean|Hi|Lo", "", name)) %>%
  select(species, predOn, value, stat, cat) %>%
  pivot_wider(names_from = stat, values_from = value)


percZero <- ydata
percZero[percZero > 0] <- 1
percZero <- colMeans(percZero)

percZero <- data.frame(species = names(percZero),
                       percZero = percZero)

all <- inner_join(all, percZero, by = c("species"))



rmspeImprovedBy <- all %>%
  filter(cat == "by")
a <- ggplot(rmspeImprovedBy) +
  geom_hline(yintercept = 0) +
  geom_ribbon(aes(x = predOn, ymin = Lo, ymax  = Hi, fill = percZero, group = percZero), alpha = 0.4) +
  geom_line(aes(x = predOn, y = Mean, color = percZero, group = percZero)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid.minor.x = element_blank(),
        #strip.placement = "outside",
        # legend.key.height = unit(.7, "cm"),
        # legend.key.width = unit(.6, "cm"),
        # legend.title = element_text(size = 9,
        #                             hjust = 0.5,
        #                             vjust = 1),
        legend.key.height = unit(.5, "cm"),
        legend.key.width = unit(1.1, "cm"),
        legend.title = element_text(size = 9,
                                    hjust = 0.5),
        legend.direction = "horizontal",
        legend.text= element_text(size = 9),
        # legend.title = element_text(size = 9),
        # legend.text= element_text(size = 9),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        plot.margin = margin(.5,.5,0,.5, "cm")) +
  labs(x = "Number of species conditioned on", y = "Percent improvement \nin RMSPE", 
       color = "Percent non-zero observations") +
  scale_x_continuous(breaks = c(1, 5, 10, 15), labels = c(1, 5, 10, 15)) +
  scale_color_gradientn(values = c(0, 1),
                        colors = c('blue', "orange"), 
                        breaks = c(0, 0.25, 0.5, 0.75, 1),
                        labels = c(0, 25, 50, 75, 100),
                        guide = guide_colorbar(title.position = "top")) +
  scale_fill_gradientn(values = c(0, 1),
                       colors = c('blue', "orange"), 
                       breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c(0, 25, 50, 75, 100),
                       guide = "none") +
  coord_cartesian(clip = "off")
# percPredsImproved <- all %>%
#   filter(cat == "frac")
# b <- ggplot(percPredsImproved) +
#   geom_hline(yintercept = 50) +
#   geom_ribbon(aes(x = predOn, ymin = Lo, ymax  = Hi, fill = percZero, group = percZero), alpha = 0.4) +
#   geom_line(aes(x = predOn, y = Mean, color = percZero, group = percZero)) +
#   theme_bw() +
#   theme(strip.background = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         #strip.placement = "outside",
#         # legend.key.height = unit(.7, "cm"),
#         # legend.key.width = unit(.6, "cm"),
#         # legend.title = element_text(size = 9,
#         #                             hjust = 0.5,
#         #                             vjust = 1),
#         legend.key.height = unit(.5, "cm"),
#         legend.key.width = unit(1.1, "cm"),
#         legend.title = element_text(size = 9,
#                                     hjust = 0.5),
#         legend.direction = "horizontal",
#         legend.text= element_text(size = 9),
#         # legend.title = element_text(size = 9),
#         # legend.text= element_text(size = 9),
#         axis.text = element_text(size = 9),
#         axis.title = element_text(size = 9),
#         plot.margin = margin(.5,.5,0,.5, "cm")) +
#   labs(x = "Number of species conditioned on", y = "Percent of \npredictions improved", 
#        color = "Percent non-zero observations") +
#   scale_x_continuous(breaks = c(1, 5, 10, 15), labels = c(1, 5, 10, 15)) +
#   scale_color_gradientn(values = c(0, 1),
#                         colors = c('blue', "orange"), 
#                         breaks = c(0, 0.25, 0.5, 0.75, 1),
#                         labels = c(0, 25, 50, 75, 100),
#                         guide = guide_colorbar(title.position = "top")) +
#   scale_fill_gradientn(values = c(0, 1),
#                        colors = c('blue', "orange"), 
#                        breaks = c(0, 0.25, 0.5, 0.75, 1),
#                        labels = c(0, 25, 50, 75, 100),
#                        guide = "none") +
#   coord_cartesian(clip = "off")


ggpubr::ggarrange(a, ncol = 1,
                  common.legend = T,
                  legend = 'top', 
                  #labels = 'auto',
                  label.x = 0, label.y = 1.01,
                  font.label = list(size = 10))


ggsave(file = "../OUT/figures/conditionOn.png", height = 90, width = 100, 
       units = "mm", dpi = 600)





### see which species decline as more incidental species are added

numSp <- all %>%
  filter(predOn %in% c(1, 16),
         cat == "frac") %>%
  select(species, predOn, Mean) %>%
  pivot_wider(values_from = "Mean",
              names_from = "predOn") %>%
  mutate(diff = `16` - `1`) %>%
  inner_join(percZero, by = "species")
# if diff is negative, predictions get worse as numSp increases
print(numSp)

plot(numSp$percZero, numSp$diff)



### see what's going on with rare species ----
RMSE_func <- function(preds, actual){
  return(sqrt(mean((actual - preds)^2)))
}


library(gjam)
library(tidyverse)
library(patchwork)

#predSp <- c('kirtlandswarbler', 'redwingedblackbird')

load("X:/conditionalPrediction/OUT/birds-gjamOutput.rdata")
ydata <- out$inputs$y

percZero <- ydata
percZero[percZero > 0] <- 1
percZero <- colMeans(percZero)


percZero <- data.frame(species = names(percZero),
                       percZero = percZero)


predSp <- colnames(ydata)



predsAll <- c()
for (s in 1:length(predSp)) {
  pred <- predSp[s]
  sp <- colnames(ydata)
  sp <- sp[-grep(pred, sp)]


  # pick species combos to condition on
  condOn1 <- combn(sp, 1)
  condOn5 <- combn(sp, 5)


  for (i in 1:2) {
    
    if (i == 1) {
      # condition on all species
      ydataCond <- select(as.data.frame(ydata), !all_of(pred))
      newdata <- list(ydataCond = ydataCond, nsim = 300)
      c      <- gjamPredict(out, newdata = newdata)
      c <- c$sdList$yMu[,pred]
      
      df <- data.frame(species = pred, 
                       obs = 1:length(c),
                       combo = i,
                       count = out$inputs$y[,pred],
                       name = 'condall',
                       value = c)
      
      predsAll <- bind_rows(predsAll, df)
    }
    
    # condition on 5 species
    condOn <- condOn5[,i]
    ydataCond <- select(as.data.frame(ydata), all_of(condOn))
    newdata <- list(ydataCond = ydataCond, nsim = 300)
    c5      <- gjamPredict(out, newdata = newdata)
    c5 <- c5$sdList$yMu[,pred]

    df <- data.frame(species = pred, 
                     obs = 1:length(c5),
                     combo = i,
                     count = out$inputs$y[,pred],
                     name = 'cond5',
                     value = c5)
    
    predsAll <- bind_rows(predsAll, df)

    # condition on one species
    condOn <- condOn1[,i]
    ydataCond <- select(as.data.frame(ydata), all_of(condOn))
    newdata <- list(ydataCond = ydataCond, nsim = 300)
    c1      <- gjamPredict(out, newdata = newdata)
    c1 <- c1$sdList$yMu[,pred]
    
    df <- data.frame(species = pred, 
                     obs = 1:length(c1),
                     combo = i,
                     count = out$inputs$y[,pred],
                     name = 'cond1',
                     value = c1)
    
    predsAll <- bind_rows(predsAll, df)

    # traditionanl prediction
    t <- out$prediction$ypredMu[,pred]

    df <- data.frame(species = pred, 
                     obs = 1:length(t),
                     combo = i,
                     count = out$inputs$y[,pred],
                     name = 'trad',
                     value = t)
    
    predsAll <- bind_rows(predsAll, df)
    # all1 <- data.frame(species = pred,
    #                    obs = 1:nrow(out$inputs$y),
    #                    combo = i,
    #                    count = out$inputs$y[,pred],
    #                    trad = t,
    #                    cond1 = c1,
    #                    cond5 = c5,
    #                    condall = c) %>%
    #   pivot_longer(cols = c("trad", "cond1", "cond5", "condall"))

    #predsAll <- bind_rows(predsAll,  all1)
  }

}

predsAll <- predsAll %>%
  mutate(Name = case_when(name == "trad" ~ "Traditional",
                          name == "cond1" ~ "Conditioned on 1 species",
                          name == "cond5" ~ "Conditioned on 5 species",
                          name == "condall" ~ "Conditioned on all other species"))
predsAll$Name <- factor(predsAll$Name, levels = c("Traditional", "Conditioned on 1 species",
                                        'Conditioned on 5 species', "Conditioned on all other species"))

all1 <- predsAll %>%
  pivot_wider(values_from = value, names_from = name,
              id_cols = c("species", "count", "combo", "obs")) %>%
  mutate(count1 = count) %>%
  pivot_longer(cols = c(count1, trad, cond1, cond5, condall)) %>%
  inner_join(percZero, by = "species")


# get % obs improved for obs = 0 and obs > 0
percImpr <- c()
for (s in 1:length(predSp)) {
  #sp <- colnames(ydata)[s]
  sp <- predSp[s]
  
  tmp <- all1 %>%
    filter(species == sp) %>%
    mutate(diff = (value - count)^2) %>%
    select(species, combo, obs, name, diff) %>%
    pivot_wider(names_from = "name",
                values_from = "diff") %>%
    pivot_longer(cols = c(cond1, cond5, condall)) %>%
    mutate(impr = trad - value,
           better = case_when(impr > 0 ~ 1,
                              T ~ 0)) %>%
    group_by(name, combo) %>%
    summarize(percImpr = mean(better)) %>%
    filter((name == "condall" & combo > 1) == F) %>%
    group_by(name) %>%
    summarize(percImprSe = sd(percImpr)/sqrt(5),
              percImpr = mean(percImpr)) %>%
    mutate(group = 'All',
           species = sp)
  
  percImpr <- bind_rows(percImpr, tmp)
  
  tmp <- all1 %>%
    filter(species == sp,
           count == 0) %>%
    mutate(diff = (value - count)^2) %>%
    select(species, combo, obs, name, diff) %>%
    pivot_wider(names_from = "name",
                values_from = "diff") %>%
    pivot_longer(cols = c(cond1, cond5, condall)) %>%
    mutate(impr = trad - value,
           better = case_when(impr > 0 ~ 1,
                              T ~ 0)) %>%
    group_by(name, combo) %>%
    summarize(percImpr = mean(better)) %>%
    filter((name == "condall" & combo > 1) == F) %>%
    group_by(name) %>%
    summarize(percImprSe = sd(percImpr)/sqrt(5),
              percImpr = mean(percImpr)) %>%
    mutate(group = 'Count = 0',
           species = sp)
  
  percImpr <- bind_rows(percImpr, tmp)
  
  tmp <- all1 %>%
    filter(species == sp,
           count > 0) %>%
    mutate(diff = (value - count)^2) %>%
    select(species, combo, obs, name, diff) %>%
    pivot_wider(names_from = "name",
                values_from = "diff") %>%
    pivot_longer(cols = c(cond1, cond5, condall)) %>%
    mutate(impr = trad - value,
           better = case_when(impr > 0 ~ 1,
                              T ~ 0)) %>%
    group_by(name, combo) %>%
    summarize(percImpr = mean(better)) %>%
    filter((name == "condall" & combo > 1) == F) %>%
    group_by(name) %>%
    summarize(percImprSe = sd(percImpr)/sqrt(5),
              percImpr = mean(percImpr)) %>%
    mutate(group = 'Count > 0',
           species = sp)
  
  percImpr <- bind_rows(percImpr, tmp)
  
}


# get rmspes for when obs = 0 and obs > 0
rmspes <- c()
for (s in 1:length(predSp)) {
  #sp <- colnames(ydata)[s]
  sp <- predSp[s]
  
  tmp <- all1 %>%
    filter(species == sp,
           count == 0,
           name != "count1") %>%
    group_by(name, combo) %>%
    summarize(rmspe = RMSE_func(value, count)) %>%
    pivot_wider(values_from = rmspe, names_from = name) %>%
    pivot_longer(!c("trad", "combo")) %>%
    mutate(percDiff = (trad-value)/trad * 100) %>%
    filter((name == "condall" & combo > 1) == F) %>%
    group_by(name) %>%
    summarize(percDiffSe = sd(percDiff)/sqrt(5),
              percDiff = mean(percDiff)) %>%
    mutate(group = 'Count = 0',
           species = sp)

  rmspes <- bind_rows(rmspes, tmp)
  
  tmp <- all1 %>%
    filter(species == sp,
           count > 0,
           name != "count1") %>%
    group_by(name, combo) %>%
    summarize(rmspe = RMSE_func(value, count)) %>%
    pivot_wider(values_from = rmspe, names_from = name) %>%
    pivot_longer(!c("trad", "combo")) %>%
    mutate(percDiff = (trad-value)/trad * 100) %>%
    filter((name == "condall" & combo > 1) == F) %>%
    group_by(name) %>%
    summarize(percDiffSe = sd(percDiff)/sqrt(5),
              percDiff = mean(percDiff)) %>%
    mutate(group = 'Count > 0',
           species = sp)
  
  rmspes <- bind_rows(rmspes, tmp)
  
  tmp <- all1 %>%
    filter(species == sp,
           name != "count1")  %>%
    group_by(name, combo) %>%
    summarize(rmspe = RMSE_func(value, count)) %>%
    pivot_wider(values_from = rmspe, names_from = name) %>%
    pivot_longer(!c("trad", "combo")) %>%
    mutate(percDiff = (trad-value)/trad * 100) %>%
    filter((name == "condall" & combo > 1) == F) %>%
    group_by(name) %>%
    summarize(percDiffSe = sd(percDiff)/sqrt(5),
              percDiff = mean(percDiff)) %>%
    mutate(group = 'All',
           species = sp)
  
  rmspes <- bind_rows(rmspes, tmp)
  

}
#colnames(rmspes) <- c("species", "trad", "cond1", "cond5", 'condall', "group")
# 
# percZero <- ydata
# percZero[percZero > 0] <- 1
# percZero <- colMeans(percZero)
# 
# 
# percZero <- data.frame(species = names(percZero),
#                        percZero = percZero)
# 

rmspes1 <- rmspes %>%
  inner_join(percZero, by = "species") %>%
  mutate(percZero = percZero*100)

a <- ggplot(filter(rmspes1, group == "All")) +
  geom_hline(yintercept = 0) +
  geom_ribbon(aes(x = factor(name, levels = c('cond1', "cond5", "condall")),
                  ymax = percDiff + percDiffSe, 
                  ymin = percDiff - percDiffSe, group = species,
                  fill = percZero), alpha = 0.4) +
  geom_line(aes(x = factor(name, levels = c('cond1', "cond5", "condall")),
                y = percDiff, color = percZero, group = species)) +
  #facet_wrap(~group, scales = "free") +
  scale_x_discrete(labels = c("1", "5", "16")) +
  scale_color_gradientn(values = c(0, 1),
                        colors = c('blue', "orange")) +
  scale_fill_gradientn(values = c(0, 1),
                       colors = c('blue', "orange"), guide = "none") +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  labs(x = "", y = "Percent improvement \nin RMSPE", 
       color = "Percent \nnon-zero \nobservations",
       title = "All") +
  coord_cartesian(ylim = c(-27, 50))

b <- ggplot(filter(rmspes1, group == "Count = 0")) +
  geom_hline(yintercept = 0) +
  geom_ribbon(aes(x = factor(name, levels = c('cond1', "cond5", "condall")),
                  ymax = percDiff + percDiffSe, 
                  ymin = percDiff - percDiffSe, group = species,
                  fill = percZero), alpha = 0.4) +
  geom_line(aes(x = factor(name, levels = c('cond1', "cond5", "condall")),
                y = percDiff, color = percZero, group = species)) +
  #facet_wrap(~group, scales = "free") +
  scale_x_discrete(labels = c("1", "5", "16")) +
  scale_color_gradientn(values = c(0, 1),
                        colors = c('blue', "orange")) +
  scale_fill_gradientn(values = c(0, 1),
                       colors = c('blue', "orange"), guide = "none") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "", y = "", 
       color = "Percent \nnon-zero \nobservations",
       title = "Focal species absent") +
  coord_cartesian(ylim = c(-27, 50))

c <- ggplot(filter(rmspes1, group == "Count > 0")) +
  geom_hline(yintercept = 0) +
  geom_ribbon(aes(x = factor(name, levels = c('cond1', "cond5", "condall")),
                  ymax = percDiff + percDiffSe, 
                  ymin = percDiff - percDiffSe, group = species,
                  fill = percZero), alpha = 0.4) +
  geom_line(aes(x = factor(name, levels = c('cond1', "cond5", "condall")),
                y = percDiff, color = percZero, group = species)) +
  #facet_wrap(~group, scales = "free") +
  scale_x_discrete(labels = c("1", "5", "16")) +
  scale_color_gradientn(values = c(0, 1),
                        colors = c('blue', "orange")) +
  scale_fill_gradientn(values = c(0, 1),
                       colors = c('blue', "orange"), guide = "none") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "", y = "", 
       color = "Percent \nnon-zero \nobservations",
       title = "Focal species present") +
  coord_cartesian(ylim = c(-27, 50))


percImpr1 <- percImpr %>%
  inner_join(percZero, by = "species") %>%
  mutate(percZero = percZero*100)
d <- ggplot(filter(percImpr1, group == "All")) +
  geom_hline(yintercept = 50) +
  geom_ribbon(aes(x = factor(name, levels = c('cond1', "cond5", "condall")),
                  ymax = (percImpr + percImprSe)*100, 
                  ymin = (percImpr - percImprSe)*100, group = species,
                  fill = percZero), alpha = 0.4) +
  # geom_point(aes(x = factor(name, levels = c('cond1', "cond5", "condall")),
  #                y = percImpr, color = percZero)) +
  geom_line(aes(x = factor(name, levels = c('cond1', "cond5", "condall")),
                y = percImpr*100, color = percZero, group = species)) +
  #facet_wrap(~group, scales = "free") +
  scale_x_discrete(labels = c("1", "5", "16")) +
  scale_color_gradientn(values = c(0, 1),
                        colors = c('blue', "orange")) +
  scale_fill_gradientn(values = c(0, 1),
                        colors = c('blue', "orange"), guide = "none") +
  theme_bw() +
  theme() +
  labs(x = "", y = "Percent of \npredictions improved", 
       color = "Percent \nnon-zero \nobservations",
       title = "All") +
  coord_cartesian(ylim = c(0, 100))

e <- ggplot(filter(percImpr1, group == "Count = 0")) +
  geom_hline(yintercept = 50) +
  geom_ribbon(aes(x = factor(name, levels = c('cond1', "cond5", "condall")),
                  ymax = (percImpr + percImprSe)*100, 
                  ymin = (percImpr - percImprSe)*100, group = species,
                  fill = percZero), alpha = 0.4) +
  # geom_point(aes(x = factor(name, levels = c('cond1', "cond5", "condall")),
  #                y = percImpr, color = percZero)) +
  geom_line(aes(x = factor(name, levels = c('cond1', "cond5", "condall")),
                y = percImpr*100, color = percZero, group = species)) +
  #facet_wrap(~group, scales = "free") +
  scale_x_discrete(labels = c("1", "5", "16")) +
  scale_color_gradientn(values = c(0, 1),
                        colors = c('blue', "orange")) +
  scale_fill_gradientn(values = c(0, 1),
                       colors = c('blue', "orange"), guide = "none") +
  theme_bw() +
  theme(axis.text.y = element_blank()) +
  labs(x = "Number of species conditioned on", y = "", 
       color = "Percent \nnon-zero \nobservations",
       title = "Focal species absent") +
  coord_cartesian(ylim = c(0, 100))

f <- ggplot(filter(percImpr1, group == "Count > 0")) +
  geom_hline(yintercept = 50) +
  geom_ribbon(aes(x = factor(name, levels = c('cond1', "cond5", "condall")),
                  ymax = (percImpr + percImprSe)*100, 
                  ymin = (percImpr - percImprSe)*100, group = species,
                  fill = percZero), alpha = 0.4) +
  # geom_point(aes(x = factor(name, levels = c('cond1', "cond5", "condall")),
  #                y = percImpr, color = percZero)) +
  geom_line(aes(x = factor(name, levels = c('cond1', "cond5", "condall")),
                y = percImpr*100, color = percZero, group = species)) +
  #facet_wrap(~group, scales = "free") +
  scale_x_discrete(labels = c("1", "5", "16")) +
  scale_color_gradientn(values = c(0, 1),
                        colors = c('blue', "orange")) +
  scale_fill_gradientn(values = c(0, 1),
                       colors = c('blue', "orange"), guide = "none") +
  theme_bw() +
  theme(axis.text.y = element_blank()) +
  labs(x = "", y = "", 
       color = "Percent \nnon-zero \nobservations",
       title = "Focal species present") +
  coord_cartesian(ylim = c(0, 100))

ggpubr::ggarrange(a, b, c, d, e, f, labels = "auto", common.legend = T, legend = "right")



ggsave(file = "../OUT/figures/metrics-split.jpg", height = 8, width = 12)
