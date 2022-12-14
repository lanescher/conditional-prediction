

# Run simulations for conditional prediction paper
# C. Lane Scher
# clanescher@gmail.com


# set up packages and functions
library(gjam)
library(tidyverse)
Rcpp::sourceCpp( 'R:/clark/clark.unix/GJAM/makeGJAMcurrent/RcppFunctions/cppFns.cpp' )

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
  
  # xx <- .replaceString(xx,'-','')
  # xx <- .replaceString(xx,'_','')
  # xx <- .replaceString(xx,' ','')
  # xx <- .replaceString(xx,"'",'')
  
  xx <- gsub('-','',xx)
  xx <- gsub('_','',xx)
  xx <- gsub(' ','',xx)
  xx <- gsub("'",'',xx)
  
  xx
}

#source("R:/clark/clark.unix/GJAM/makeGJAMcurrent/Rfunctions/gjamHfunctions.r")

# set up variables for gjam runs
ng <- 500
burnin <- ng/2

n <- 100
S <- 5
Q <- 3

reps <- 3


# Simulate the data ----

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


# Models ----

MVNall <- c()
MVN2 <- c()

for (i in 1:reps) {
  
  rep <- data[[i]]
  
  ## MVN ----
  
  # all covariances
  f <- rep[[1]]
  
  formula <- f$formula
  xdata   <- f$xdata
  ydata   <- f$ydata
  
  # fit model
  modelList <- list(ng = ng, burnin = burnin, typeNames = f$typeNames )
  out <- gjam( formula = formula, xdata, ydata, modelList )
  
  # compare
  compare <- conditionalComparisonEach(output = out)
  
  MVNall <- bind_rows(MVNall, compare)

  # one covariance
  f <- rep[[2]]
  
  formula <- f$formula
  xdata   <- f$xdata
  ydata   <- f$ydata
  
  # fit model
  modelList <- list(ng = ng, burnin = burnin, typeNames = f$typeNames )
  out <- gjam( formula = formula, xdata, ydata, modelList )
  
  # compare
  compare <- conditionalComparisonEach(output = out)
  
  MVN2 <- bind_rows(MVN2, compare)
  
  
  ## CA ----
  
  
  
  
}
MVNall$type <- "CON"
MVNall$cov <- "all"

## 2 species covariance ----

for (i in 1:reps) {
  
  # simulate data
  f <- gjamSimData(n = n, 
                   S = S, 
                   Q = Q,
                   typeNames = 'CON')
  
  # modify so that only two species have residual covariance
  formula <- f$formula
  beta  <- f$trueValues$beta
  stmp  <- diag( 1, S )
  stmp[2,3] <- stmp[3,2] <- -.7
  
  xdata     <- as.matrix( f$xdata )
  mu    <- xdata%*%beta
  ydata     <- mu + myrmvnorm(nrow(xdata), rep(0,S), stmp)
  
  # fit model
  modelList <- list(ng = ng, burnin = burnin, typeNames = f$typeNames )
  out <- gjam(formula = formula, xdata, ydata, modelList )
  
  # compare
  compare <- conditionalComparisonEach(output = out)
  
  MVN2 <- bind_rows(MVN2, compare)
}
MVN2$type <- "CON"
MVN2$cov <- "2"



# 2. CA ----

## all covariance ----

CAall <- c()
for (i in 1:reps) {
  
  # simulate data
  f <- gjamSimData(n = n, 
                   S = S, 
                   Q = Q,
                   typeNames = 'CA')
  
  formula <- f$formula
  xdata   <- f$xdata
  ydata   <- f$ydata
  
  # fit model
  modelList <- list(ng = ng, burnin = burnin, typeNames = f$typeNames )
  out <- gjam( formula = formula, xdata, ydata, modelList )
  
  # compare
  compare <- conditionalComparisonEach(output = out)
  
  CAall <- bind_rows(CAall, compare)
}
CAall$type <- "CA"
CAall$cov <- "all"


## 2 species covariance ----

CA2 <- c()
for (i in 1:reps) {
  
  # simulate data
  f <- gjamSimData(n = n, 
                   S = S, 
                   Q = Q,
                   typeNames = 'CA')
  
  # modify so that only two species have residual covariance
  formula <- f$formula
  beta  <- f$trueValues$beta
  stmp  <- diag( 1, S )
  stmp[2,3] <- stmp[3,2] <- -.7
  
  xdata     <- as.matrix( f$xdata )
  mu    <- xdata%*%beta
  ydata     <- mu + myrmvnorm(nrow(xdata), rep(0,S), stmp)
  ydata[ydata < 0] <- 0
  
  # fit model
  modelList <- list(ng = ng, burnin = burnin, typeNames = f$typeNames )
  out <- gjam(formula = formula, xdata, ydata, modelList )
  
  # compare
  compare <- conditionalComparisonEach(output = out)
  
  CA2 <- bind_rows(CA2, compare)
}
CA2$type <- "CA"
CA2$cov <- "2"


# 3. DA ----

## all covariance ----

DAall <- c()
for (i in 1:reps) {
  
  # simulate data
  f <- gjamSimData(n = n, 
                   S = S, 
                   Q = Q,
                   typeNames = 'DA')
  
  formula <- f$formula
  xdata   <- f$xdata
  ydata   <- f$ydata
  
  # fit model
  modelList <- list(ng = ng, burnin = burnin, typeNames = f$typeNames )
  out <- gjam( formula = formula, xdata, ydata, modelList )
  
  # compare
  compare <- conditionalComparisonEach(output = out)
  
  DAall <- bind_rows(DAall, compare)
}
DAall$type <- "DA"
DAall$cov <- "all"


## 2 species covariance ----

DA2 <- c()
for (i in 1:reps) {
  
  # simulate data
  f <- gjamSimData(n = n, 
                   S = S, 
                   Q = Q,
                   typeNames = 'DA')
  
  # modify so that only two species have residual covariance
  formula <- f$formula
  beta  <- f$trueValues$beta
  stmp  <- diag( 1, S )
  stmp[2,3] <- stmp[3,2] <- -.7
  
  xdata     <- as.matrix( f$xdata )
  mu    <- xdata%*%beta
  ydata     <- mu + myrmvnorm(nrow(xdata), rep(0,S), stmp)
  ydata     <- round(ydata)
  ydata[ydata < 0 ] <- 0
  
  # fit model
  modelList <- list(ng = ng, burnin = burnin, typeNames = f$typeNames )
  out <- gjam(formula = formula, xdata, ydata, modelList )
  
  # compare
  compare <- conditionalComparisonEach(output = out)
  
  DA2 <- bind_rows(DA2, compare)
}
DA2$type <- "DA"
DA2$cov <- "2"



# 4. PA ----

## all covariance ----

PAall <- c()
for (i in 1:reps) {
  
  # simulate data
  f <- gjamSimData(n = n, 
                   S = S, 
                   Q = Q,
                   typeNames = 'PA')
  
  formula <- f$formula
  xdata   <- f$xdata
  ydata   <- f$ydata
  
  # fit model
  modelList <- list(ng = ng, burnin = burnin, typeNames = f$typeNames )
  out <- gjam( formula = formula, xdata, ydata, modelList )
  
  # compare
  compare <- conditionalComparisonEach(output = out)
  
  PAall <- bind_rows(PAall, compare)
}
PAall$type <- "PA"
PAall$cov <- "all"


## 2 species covariance ----

PA2 <- c()
for (i in 1:reps) {
  
  # simulate data
  f <- gjamSimData(n = n, 
                   S = S, 
                   Q = Q,
                   typeNames = 'PA')
  
  # modify so that only two species have residual covariance
  formula <- f$formula
  beta  <- f$trueValues$beta
  stmp  <- diag( 1, S )
  stmp[2,3] <- stmp[3,2] <- -.7
  
  xdata     <- as.matrix( f$xdata )
  mu    <- xdata%*%beta
  w     <- mu + myrmvnorm(nrow(xdata), rep(0,S), stmp)
  w     <- tnormMVNmatrix(avec = w, muvec = mu,
                          smat = stmp, lo = w*0 - 3,
                          hi = w*0 + 3)
  ydata     <- round( pnorm(w) )
  
  ydata     <- round(ydata)
  ydata[ydata < 0 ] <- 0
  
  # fit model
  modelList <- list(ng = ng, burnin = burnin, typeNames = f$typeNames )
  out <- gjam(formula = formula, xdata, ydata, modelList )
  
  # compare
  compare <- conditionalComparisonEach(output = out)
  
  PA2 <- bind_rows(PA2, compare)
}
PA2$type <- "PA"
PA2$cov <- "2"

# ISSUE: 
# Error in trMVNmatrixRcpp(avec, muvec, smat, lo, hi, whichSample, idxALL = c(0:(nrow(smat) -  : 
# could not find function "trMVNmatrixRcpp"




# 5. Plot ----
all <- bind_rows(MVNall, MVN2, 
                 CAall, CA2,
                 DAall, DA2) %>%
  mutate(covSp = case_when(species %in% c("S2", "S3") &
                             cov == "2" ~ 1,
                           cov == "all" ~ 9,
                           T ~ 0),
         fracuc = `Frac u > c`*100,
         type1 = case_when(type == "CON" ~ "Continuous",
                           type == "CA" ~ "Continuous abundance",
                           type == "DA" ~ "Discrete abundance"))


all$type <- factor(all$type1, levels = c("Continuous", "Continuous abundance", "Discrete abundance"))

ggplot(all) +
  geom_hline(yintercept = 0, color = "black") +
  geom_boxplot(aes(x = as.factor(covSp), y = `Perc Diff`), fill = "gray") +
  facet_grid(cols = vars(type)) +
  coord_cartesian(ylim = c(0, 100)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank()) +
  labs(x = "", y = "Percent improvement in RMSPE")

ggplot(all) +
  geom_hline(yintercept = 50, color = "black") +
  geom_boxplot(aes(x = as.factor(covSp), y = fracuc), fill = "gray") +
  facet_grid(cols = vars(type)) +
  coord_cartesian(ylim = c(0, 100)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank()) +
  labs(x = "", y = "Percent of observations improved")

ggplot(all) +
  geom_point(aes(x = percZero, y = `Perc Diff`, color = as.factor(covSp))) +
  facet_grid(cols = vars(type))


ggplot(all) +
  geom_histogram(aes(x = percZero)) +
  facet_wrap(~type)

### ISSUE: different 0s with different simulated data







# 
# 
# 
# ############# Presence-absence
# f   <- gjamSimData( n = n, S = S, typeNames = 'PA' )
# 
# modelList  <- list( ng = 5000, burnin = 500, typeNames = f$typeNames )
# output <- gjam( formula = f$formula, xdata = f$xdata, ydata = f$ydata,
#                 modelList = modelList )
# outPA <- output
# 
# ############## PA, one covariance
# 
# beta  <- f$trueValues$beta
# #sigma <- f$trueValues$sigma
# stmp  <- diag( 1, S )
# stmp[2,3] <- stmp[3,2] <- -.7
# x     <- as.matrix( f$xdata )
# mu    <- x%*%beta
# w     <- mu + myrmvnorm(nrow(x), rep(0,S), stmp)
# w     <- tnormMVNmatrix(avec = w, muvec = mu,
#                          smat = stmp, lo = w*0 - 3,
#                          hi = w*0 + 3)
# y     <- round( pnorm(w) )
# 
# modelList  <- list( ng = 5000, burnin = 500, typeNames = f$typeNames )
# output <- gjam( formula = f$formula, xdata = x, ydata = y, modelList = modelList )
# outPA1 <- output
# 
# 
# ##### 1b. conditional prediction
# 
# # for each version of out, conditionally predict each species
# all <- c()
# for (o in 1:length(outs)) {
#   out <- outs[[o]]
#   print(paste0("Out ", o, " of ", length(outs)))
#   
#   ydata <- out$inputs$y
#   sp <- colnames(ydata)
#   
#   for (s in 1:length(sp)) {
#     print(paste0("Species ", s, " of ", length(sp)))
#     
#     for (i in 1:10) {
#       print(paste0("Run ", i, " of 10"))
#       pred <- sp[s]
#       condOn <- sp[-s]
#       
#       condOnData <- ydata[,condOn]
#       newdata <- list(ydataCond = condOnData,  nsim=200)
#       
#       # compare conditional and unconditional
#       comp <- conditionalComparison( out, newdata )
#       comp <- comp[[1]]
#       comp <- as.data.frame(comp)
#       
#       comp$species <- sp[s]
#       comp$out <- paste0("out-", o)
#       comp$run <- i
#       
#       all <- bind_rows(all, comp)
#     }
#     
#   }
#   
# }
# 
# all1 <- all %>%
#   group_by(species, out) %>%
#   summarize(percDiffMean = mean(`Perc Diff`),
#             percDiffSe = sd(`Perc Diff`)/sqrt(n()),
#             percImprMean = mean(`Frac u > c`*100),
#             percImprSe = sd(`Frac u > c`*100)/sqrt(n()))
