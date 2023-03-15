

# Run simulations for conditional prediction paper
# C. Lane Scher
# clanescher@gmail.com


# set up packages and functions
library(gjam)
library(tidyverse)

.binaryScore <- function(p, x){
  
  #brier and logarithmic score, prediction prob p, event x = 0 or 1
  
  a <- mean((x - p)^2)
  b <- -mean( x*log(p) + (1 - x)*log(1 - p))
  
  list(brierScore = a, logScore = b)
}



compareEachSpecies <- function(out, nsim = 2000) {
  compare <- c()
  ydata <- out$inputs$y
  
  for (s in 1:ncol(out$inputs$y)){
    
    true <- ydata[,s]
    trad <- out$prediction$ypredMu[,s]
    
    
    # conditional prediction
    
    conditionOn <- colnames(ydata)[-s]
    conditionOnData <- select(as.data.frame(ydata), all_of(conditionOn))
    
    newdata <- list(ydataCond = conditionOnData,  nsim=nsim)
    cond_pred <- gjamPredict(output = out, newdata = newdata)
    cond <- cond_pred$sdList$yMu[,s]
    
    # get spe
    con <- (ydata[,s] - cond)^2
    unc <- (ydata[,s] - trad)^2
    
    # how many spe are greater for traditional prediction
    umc  <- unc - con
    uno  <- length(which(umc > 0))/length(umc)  # fraction larger error in uncon
    
    # get rmspe
    rmspecond <- sqrt(mean(con))
    rmspetrad <- sqrt(mean(unc))
    
    if(unique(out$modelList$typeNames) == "PA") {
      tt <- .binaryScore(cond, true)
      rmspecond <- tt$brierScore

      tt <- .binaryScore(trad, true)
      rmspetrad <- tt$brierScore

    }
    
    
    
    
    # add to compare dataframe
    comp <- data.frame(species = colnames(ydata)[s],
                       "Unc RMSPE" = rmspetrad, 
                       "Con RMSPE" = rmspecond, 
                       "Perc Diff" = (rmspetrad-rmspecond)/rmspetrad * 100, 
                       "Frac u > c" = uno)
    colnames(comp) <- c("species", "Unc RMSPE", "Con RMSPE", "Perc Diff", "Frac u > c")
    
    compare <- bind_rows(compare, comp)
  }
  
  return(compare)
  
} 


myrmvnorm <- function (nn, mu, sigma){
  
  # nn - no. samples from one mu vector or nrow(mu) for matrix
  
  if(!is.matrix(mu)) mu <- matrix(mu,1)
  if(length(mu) == 1)mu <- matrix(mu,1,nrow(sigma))
  if(ncol(mu) == 1)  mu <- t(mu)
  
  m <- ncol(sigma)
  
  if(ncol(mu) != m)stop('dimension mismatch mu, sigma')
  
  if(nn > 1 & nrow(mu) == 1)mu <- matrix(mu,nn,m,byrow=T)
  
  if(nn != nrow(mu))stop('sample size does not match mu')
  
  testv <- try(svd(sigma),T)
  
  if( inherits(testv,'try-error') ){
    ev <- eigen(sigma, symmetric = TRUE)
    testv <- t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))
  } else {
    testv <- t(testv$v %*% (t(testv$u) * sqrt(testv$d)))
  }
  
  p <- matrix(rnorm(nn * m), nn) %*% testv
  p + mu
}

tnormMVNmatrix <- function(avec, muvec, smat, 
                           lo=matrix(-1000,nrow(muvec),ncol(muvec)), 
                           hi=matrix(1000,nrow(muvec),ncol(muvec)),
                           whichSample = c(1:nrow(smat)) ){
  
  #lo, hi must be same dimensions as muvec,avec
  
  lo[lo < -1000] <- -1000
  hi[hi > 1000]  <- 1000
  
  if(max(whichSample) > length(muvec))
    stop('whichSample outside length(muvec)')
  
  r <- avec
  a <- trMVNmatrixRcpp(avec, muvec, smat, lo, hi, whichSample, 
                       idxALL = c(0:(nrow(smat)-1)) )  
  r[,whichSample] <- a[,whichSample]
  r
}


# set up variables for gjam runs
ng <- 5000
burnin <- ng/2

n <- 1000
S <- 10
Q <- 5

reps <- 10


# 1. Simulate the data ----

data <- list()
for (i in 1:reps) {
  
  d <- list()
  
  # simulate data with all covariance
  f1 <- gjamSimData(n = n, 
                   S = S, 
                   Q = Q,
                   typeNames = 'CON')
  d[[1]] <- f1
  
  # simulate data with one covariance
  f2 <- gjamSimData(n = n, 
                   S = S, 
                   Q = Q,
                   typeNames = 'CON')
  
  # modify so that only two species have residual covariance
  beta  <- f2$trueValues$beta
  stmp  <- diag( 1, S )
  stmp[2,3] <- stmp[3,2] <- -.7
  
  xdata     <- as.matrix( f2$xdata )
  mu    <- xdata%*%beta
  ydata     <- mu + myrmvnorm(nrow(xdata), rep(0,S), stmp)
  
  f2$xdata <- xdata
  f2$ydata <- ydata
  
  d[[2]] <- f2
  
  data[[i]] <- d
  
}


# 2. Models ----

MVNall <- c()
MVN2 <- c()

CAall <- c()
CA2 <- c()

DAall <- c()
DA2 <- c()

PAall <- c()
PA2 <- c()

for (i in 1:reps) {
  
  rep <- data[[i]]
  
  ## MVN ----
  
  ### all covariances ----
  f <- rep[[1]]
  
  formula <- f$formula
  xdata   <- f$xdata
  ydata   <- f$ydata
  
  # fit model
  modelList <- list(ng = ng, burnin = burnin, typeNames = "CON")
  out <- gjam( formula = formula, xdata, ydata, modelList )
  
  # compare
  #compare <- conditionalComparisonEach(output = out)
  compare <- compareEachSpecies(out = out)
  
  
  MVNall <- bind_rows(MVNall, compare)

  ### one covariance ----
  f <- rep[[2]]
  
  formula <- f$formula
  xdata   <- f$xdata
  ydata   <- f$ydata
  
  # fit model
  modelList <- list(ng = ng, burnin = burnin, typeNames = "CON")
  out <- gjam( formula = formula, xdata, ydata, modelList )
  
  # compare
  #compare <- conditionalComparisonEach(output = out)
  compare <- compareEachSpecies(out = out)
  
  
  
  
  MVN2 <- bind_rows(MVN2, compare)
  
  ## CA ----
  
  ### all covariances ----
  f <- rep[[1]]
  
  formula <- f$formula
  xdata   <- f$xdata
  ydata   <- f$ydata
  
  # make it CA
  ydata[ydata < 0] <- 0
  
  # fit model
  modelList <- list(ng = ng, burnin = burnin, typeNames = "CA")
  out <- gjam( formula = formula, xdata, ydata, modelList )
  
  # compare
  #compare <- conditionalComparisonEach(output = out)
  compare <- compareEachSpecies(out = out)
  
  CAall <- bind_rows(CAall, compare)
  
  
  ### one covariance ----
  
  f <- rep[[2]]
  
  formula <- f$formula
  xdata   <- f$xdata
  ydata   <- f$ydata
  
  # make it CA
  ydata[ydata < 0] <- 0
  
  # fit model
  modelList <- list(ng = ng, burnin = burnin, typeNames = "CA")
  out <- gjam( formula = formula, xdata, ydata, modelList )
  
  # compare
  #compare <- conditionalComparisonEach(output = out)
  compare <- compareEachSpecies(out = out)
  
  CA2 <- bind_rows(CA2, compare)
  
  ## DA ----
  
  ### all covariances ----
  f <- rep[[1]]
  
  formula <- f$formula
  xdata   <- f$xdata
  ydata   <- f$ydata
  
  # make it DA
  ydata[ydata < 0] <- 0
  ydata <- round(ydata)
  
  # fit model
  modelList <- list(ng = ng, burnin = burnin, typeNames = "DA")
  out <- gjam( formula = formula, xdata, ydata, modelList )
  
  # compare
  compare <- compareEachSpecies(out = out)
  
  DAall <- bind_rows(DAall, compare)
  
  
  ### one covariance ----
  
  f <- rep[[2]]
  
  formula <- f$formula
  xdata   <- f$xdata
  ydata   <- f$ydata
  
  # make it DA
  ydata[ydata < 0] <- 0
  ydata <- round(ydata)
  
  # fit model
  modelList <- list(ng = ng, burnin = burnin, typeNames = "DA")
  out <- gjam( formula = formula, xdata, ydata, modelList )
  
  # compare
  compare <- compareEachSpecies(out = out)
  
  DA2 <- bind_rows(DA2, compare)
  
  
  
  ## PA ----
  
  ### all covariances ----
  f <- rep[[1]]
  
  formula <- f$formula
  xdata   <- f$xdata
  ydata   <- f$ydata
  
  # make it PA
  ydata[ydata < 0] <- 0
  ydata[ydata > 0] <- 1
  
  # fit model
  modelList <- list(ng = ng, burnin = burnin, typeNames = "PA")
  out <- gjam( formula = formula, xdata, ydata, modelList )
  
  # compare
  #compare <- conditionalComparisonEach(output = out)
  compare <- compareEachSpecies(out = out)
  
  PAall <- bind_rows(PAall, compare)
  
  
  ### one covariance ----
  
  f <- rep[[2]]
  
  formula <- f$formula
  xdata   <- f$xdata
  ydata   <- f$ydata
  
  # make it PA
  ydata[ydata < 0] <- 0
  ydata[ydata > 0] <- 1
  
  # fit model
  modelList <- list(ng = ng, burnin = burnin, typeNames = "PA")
  out <- gjam( formula = formula, xdata, ydata, modelList )
  
  # compare
  #compare <- conditionalComparisonEach(output = out)
  compare <- compareEachSpecies(out = out)
  
  PA2 <- bind_rows(PA2, compare)
  
}

# add data type and covariance to each model
MVNall$type <- "CON"
MVNall$cov <- "all"

MVN2$type <- "CON"
MVN2$cov <- "2"

CAall$type <- "CA"
CAall$cov <- "all"

CA2$type <- "CA"
CA2$cov <- "2"

DAall$type <- "DA"
DAall$cov <- "all"

DA2$type <- "DA"
DA2$cov <- "2"

PAall$type <- "PA"
PAall$cov <- "all"

PA2$type <- "PA"
PA2$cov <- "2"


# 3. Plot ----
all <- bind_rows(MVNall, MVN2, 
                 CAall, CA2,
                 DAall, DA2,
                 PAall, PA2) %>%
  mutate(covSp = case_when(species %in% c("S2", "S3") &
                             cov == "2" ~ 1,
                           cov == "all" ~ 9,
                           T ~ 0),
         fracuc = `Frac u > c`*100,
         type1 = case_when(type == "CON" ~ "Continuous",
                           type == "CA" ~ "Continuous abundance",
                           type == "DA" ~ "Discrete abundance",
                           type == "PA" ~ "Presence-absence"))


all$type1 <- factor(all$type1, levels = c("Continuous", "Continuous abundance", 
                                         "Discrete abundance", "Presence-absence"))

#save(all, file = "../OUT/simulation.rdata")
load("OUT/simulation.rdata")


all$type1 <- factor(all$type, levels = c("CON", "CA", "DA", "PA"))


cond.labs1 <- c("a", "b", "c")
names(cond.labs1) <- c(0, 1, 9)

cols <- c("hotpink4", "dodgerblue3", "forestgreen", "orange2")
cols1 <- colorspace::lighten(cols)
cols2 <- colorspace::lighten(cols, 0.6)


base <- ggplot() +
  scale_color_manual(values = cols1) +
  scale_fill_manual(values = cols2) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        plot.title = element_text(size = 9,
                                  hjust = 0.5))

top <- base  +
  coord_cartesian(ylim = c(-150, 100)) +
  labs(x = "", y = "Percent improvement \nin RMSPE", fill = "", color = "") +
  geom_hline(yintercept = 0, color = "black") +
  theme(plot.margin = margin(15, 0.1, 0, 0.1),
        axis.text.x = element_blank())

bottom <- base  +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "", y = "Percent of \npredictions improved", fill = "", color = "") +
  geom_hline(yintercept = 50, color = "black") +
  theme(plot.margin = margin(0, 0.1, 0, 0.1)) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100))


a <- top +
  geom_boxplot(data = filter(all, covSp == 9),
               aes(x = type1, fill = type1, color = type1,
                   y = `Perc Diff`)) +
  labs(title = "Residual covariance with \n9 species")

b <- top +
  geom_boxplot(data = filter(all, covSp == 1),
               aes(x = type1, fill = type1, color = type1,
                   y = `Perc Diff`)) +
  labs(y = "", title = "Residual covariance with \n1 species")

c <- top +
  geom_boxplot(data = filter(all, covSp == 0),
               aes(x = type1, fill = type1, color = type1,
                   y = `Perc Diff`)) +
  labs(y = "", title = "Residual covariance with \n0 species")

d <- bottom +
  geom_boxplot(data = filter(all, covSp == 9),
               aes(x = type1, fill = type1, color = type1,
                   y = (`Frac u > c`)*100))

e <- bottom +
  geom_boxplot(data = filter(all, covSp == 1),
               aes(x = type1, fill = type1, color = type1,
                   y = (`Frac u > c`)*100)) +
  labs(y = "", x = "Data type")

f <- bottom +
  geom_boxplot(data = filter(all, covSp == 0),
               aes(x = type1, fill = type1, color = type1,
                   y = (`Frac u > c`)*100)) +
  labs(y = "")

ggpubr::ggarrange(a, b, c, d, e, f,
          labels = "auto",
          label.x = c(0.13, 0.09, 0.09, 0.13, 0.09, 0.09), 
          label.y = c(0.85, 0.85, 0.85, 1.1, 1.1, 1.1),
          font.label = list(size = 10),
          heights = c(1.2, 1))

ggsave("OUT/figures/simulationSummary.png", 
       height = 100, width = 180, units = "mm", dpi = 600)



