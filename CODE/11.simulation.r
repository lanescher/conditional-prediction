

# Run simulations for conditional prediction paper
# C. Lane Scher
# clanescher@gmail.com


# set up packages and functions
library(gjam)
library(tidyverse)
#Rcpp::sourceCpp( 'R:/clark/clark.unix/GJAM/makeGJAMcurrent/RcppFunctions/cppFns.cpp' )

.binaryScore <- function(p, x){
  
  #brier and logarithmic score, prediction prob p, event x = 0 or 1
  
  a <- mean((x - p)^2)
  b <- -mean( x*log(p) + (1 - x)*log(1 - p))
  
  list(brierScore = a, logScore = b)
}


conditionalComparison <- function( output, newdata ){
  
  # compare conditional prediction with unconditional prediction in output
  
  ydataCond <- newdata$ydataCond
  colnames(ydataCond) <- .cleanNames( colnames(ydataCond) )
  
  y           <- as.matrix( output$inputs$y )
  colnames(y) <- .cleanNames( colnames(y) )
  
  typeNames <- output$modelList$typeNames
  zeroCols  <- which(typeNames %in% c('DA','CA','OC','CC'))
  
  x         <- output$inputs$xUnstand
  S         <- ncol(y)
  n         <- nrow(x)
  if( is.null(ydataCond) )stop(' newdata must include ydataCond ' )
  if( is.null(colnames(ydataCond)) )stop( 'ydataCond must have colnames matching ydata' )
  
  y0 <- y
  y0[ y0 > 0 ] <- 1
  percZero <- signif(1 - apply(y0, 2, sum)/nrow(y0),3)*100
  
  ni <- match(colnames(ydataCond), colnames(y))   # condition on these columns
  ci <- c(1:S)[-ni]
  if( length(ni) < ncol(ydataCond) )stop( 'ydataCond must have colnames matching ydata' )
  
  edata <- output$inputs$effMat
  if( is.null(edata) )edata <- newdata$effort$values
  if( is.null(edata) )edata <- y*0 + 1
  wUn <- output$prediction$ypredMu/edata
  w   <- y/edata
  
  preds  <- gjamPredict(output, newdata = newdata, FULL=T) 
  condy  <- preds$sdList$yMu                                # on data scale
  condw  <- condy/edata                                     # latent scale
  
  coni <- (w[,ci] - condw[,ci] )^2
  unci <- ( w[,ci] - wUn[,ci] )^2
  umc  <- unci - coni
  uno  <- length(which(umc > 0))/length(umc)  # fraction larger error in uncon
  
  uno <- umc
  uno[ uno > 0 ] <- 1
  uno[ uno <= 0 ] <- 0
  
  
  if( length(ci) == 1 ){
    con <- mean( coni  )
    unc <- mean( unci )
    uno <- sum(uno)/length(uno)
  }else{ 
    con <- diag( var( w[,ci] - condw[,ci] ) )
    unc <- diag( var( w[,ci] - wUn[,ci] ) )
    uno <- apply(uno, 2, sum)/nrow(umc)
  }
  
  wi <- which(typeNames[ci] == 'PA')
  if( length(wi) > 0){
    tt <- .binaryScore(condw[,ci[wi]], y[,ci[wi]])
    con <- tt$brierScore
    #  ccLog   <- tt$logScore
    
    tt <- .binaryScore(wUn[,ci[wi]], y[,ci[wi]])
    unc <- tt$brierScore
    #  uuLog   <- tt$logScore
  }
  
  
  sigma <- output$parameters$sigMu
  beta  <- output$parameters$betaMu
  
  C  <- solve(sigma[ni,ni])%*%sigma[ni,ci]      # check for correlation matrix
  mu <- x%*%beta
  
  if(length(zeroCols) > 0){
    mz <- mu[,zeroCols]
    mz[ mz < 0 ] <- 0
    mu[,zeroCols] <- mz
  }
  Rs <- w[,ci,drop=F] - mu[,ci,drop=F]
  Rn <- w[,ni,drop=F] - mu[,ni,drop=F]
  D  <- Rn%*%C
  # z  <- -(2*t(Rs) + t(D))%*%D
  
  z <- ( 2*t(Rs) - t(C)%*%t(Rn) )%*%Rn%*%C
  colnames(z) <- rownames(z)
  
  vmu <- var(mu[,ci])
  vd  <- var(D)
  cf  <- vd/(vmu + vd)      #variance fraction from conditional
  
  fe <- output$fit$fractionExplained
  ynames <- colnames(fe)
  ynames <- ynames[ ynames %in% colnames(condy)[ci] ]
  ynames <- ynames[ !ynames == 'other' ]
  
  if( length(wi) == 0 ){   # not PA data
    
    if( length(ci) == 1 ){
      
      imat <- cbind( unc, con, 
                     (unc - con)/unc*100, uno,
                     fe[1,ynames], diag(cf), diag(z)[ynames]/n )
    }else{
      imat <- cbind( unc[ynames], con[ynames], 
                     (unc - con)[ynames]/unc[ynames]*100, uno[ynames],
                     fe[1,ynames], diag(cf)[ynames], diag(z)[ynames]/n )
    }
    colnames(imat) <- c("Unc MSPE", "Con MSPE", "Perc Diff", "Frac u > c",
                        "Frac Expl", 
                        "Frac Con", "CIS")
  }else{  # PA data
    
    if( length(ci) == 1 ){
      
      imat <- cbind( unc, con, 
                     (unc - con)/unc*100, uno,
                     fe[1,ynames], diag(cf) )
    }else{
      imat <- cbind( unc[ynames], con[ynames], (unc - con)[ynames]/unc[ynames]*100, 
                     uno[ynames],
                     fe[1,ynames], diag(cf)[ynames] )
    }
    colnames(imat) <- c("Unc Brier", "Con Brier", "Perc Diff", "Frac u > c",
                        "Frac Expl", 
                        "Frac Con")
  }
  if(typeNames[1] != 'CON'){
    percZero <- percZero[ynames]
    imat <- cbind(imat, percZero)
  }
  if( nrow(imat) == 1)rownames(imat) <- colnames(y)[ci]
  out1 <- list(imat, C)
  out1
}

# this function applies conditionalComparison() iteratively
# to each species in a model, conditioning on all other species
conditionalComparisonEach <- function(output) {
  sp <- colnames(out$inputs$y)
  
  tmp <- c()
  for (s in 1:length(sp)) {
    print(paste0("Species ", s, " of ", S))
    
    # get species to predict and condition on
    pred <- sp[s]
    condOn <- sp[-s]
    
    # set up data
    condOnData <- ydata[,condOn]
    newdata <- list(ydataCond = condOnData,  nsim=200)
    
    # compare conditional and unconditional
    comp <- conditionalComparison(out, newdata)
    comp <- comp[[1]]
    comp <- as.data.frame(comp)
    
    comp$species <- sp[s]
    #comp$out <- paste0("out-", o)
    comp$run <- i

    tmp <- bind_rows(tmp, comp)
  }
  return(tmp)
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

.cleanNames <- function(xx){
 
  xx <- gsub('-','',xx)
  xx <- gsub('_','',xx)
  xx <- gsub(' ','',xx)
  xx <- gsub("'",'',xx)
  
  xx
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
  compare <- conditionalComparisonEach(output = out)
  
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
  compare <- conditionalComparisonEach(output = out)
  
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
  compare <- conditionalComparisonEach(output = out)
  
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
  compare <- conditionalComparisonEach(output = out)
  
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
  compare <- conditionalComparisonEach(output = out)
  
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
  compare <- conditionalComparisonEach(output = out)
  
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
  compare <- conditionalComparisonEach(output = out)
  
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
  compare <- conditionalComparisonEach(output = out)
  
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
load("../OUT/simulation.rdata")

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
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 8,
                                  hjust = 0.5))

top <- base  +
  coord_cartesian(ylim = c(-25, 100)) +
  labs(x = "", y = "Percent improvement \nin RMSPE", fill = "", color = "") +
  geom_hline(yintercept = 0, color = "black") +
  theme(plot.margin = margin(15, 0.1, 0, 0.1),
        axis.text.x = element_blank())

bottom <- base  +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "", y = "Percent of \nobservations improved", fill = "", color = "") +
  geom_hline(yintercept = 50, color = "black") +
  theme(plot.margin = margin(0, 0.1, 0, 0.1))


a <- top +
  geom_boxplot(data = filter(all, covSp == 9),
               aes(x = type, fill = type, color = type,
                   y = `Perc Diff`)) +
  labs(title = "Residual covariance with 9 species")

b <- top +
  geom_boxplot(data = filter(all, covSp == 1),
               aes(x = type, fill = type, color = type,
                   y = `Perc Diff`)) +
  labs(y = "", title = "Residual covariance with 1 species")

c <- top +
  geom_boxplot(data = filter(all, covSp == 0),
               aes(x = type, fill = type, color = type,
                   y = `Perc Diff`)) +
  labs(y = "", title = "Residual covariance with 0 species")

d <- bottom +
  geom_boxplot(data = filter(all, covSp == 9),
               aes(x = type, fill = type, color = type,
                   y = (`Frac u > c`)*100))

e <- bottom +
  geom_boxplot(data = filter(all, covSp == 1),
               aes(x = type, fill = type, color = type,
                   y = (`Frac u > c`)*100)) +
  labs(y = "")

f <- bottom +
  geom_boxplot(data = filter(all, covSp == 0),
               aes(x = type, fill = type, color = type,
                   y = (`Frac u > c`)*100)) +
  labs(y = "")

ggpubr::ggarrange(a, b, c, d, e, f,
          labels = "auto",
          label.x = 0.05, label.y = .95,
          font.label = list(size = 10),
          align = "hv")

ggsave("../OUT/figures/simulationSummary.png", 
       height = 100, width = 180, units = "mm", dpi = 600)



