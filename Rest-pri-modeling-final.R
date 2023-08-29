# remove objects from workspace
rm(list=ls())

# load libraries
library(tidyverse)
library(LatticeKrig)
library(geoR)
library(raster)
library(sf)
library(stargazer)
library(geobr)
select <- dplyr::select


# load environment from script: "Rest-pri-helper-final.R"
load('./Data/Restoration_models_2021.Rdata')

## create main df
fp <- forest_points

x <- t(as.data.frame(lapply(fp$geom, as.numeric)))
fp <- fp %>% select(IDF, Forest_cover_1985, #z
                    Forest_cover_2012,
                    Percent_forest_1985,
                    Percent_forest_2012,
                    Type, #z
                    name_muni, #z
                    Size, #z
                    #Gini, #z
                    Distance_roads, #z
                    Avg_precip, #z
                    #Transition_1985, #y
                    Restoration_1985, #y
                    Restoration_2012,
                    #Deforestation_1985, #y
                    Percent_restoration_1985, #y
                    Percent_restoration_2012
)
fp <- fp %>% mutate(x1 = x[,1], x2 = x[,2])
fp <- distinct(fp, x1, x2, .keep_all=TRUE) %>% na.omit()

fpp <- fp %>% mutate(Percent_restoration_1985 = Percent_restoration_1985,
                     Size = Size,
                     Restoration_1985 = Restoration_1985)

fp <- fp %>% mutate(Forest_cover_1985 = Forest_cover_1985/900,
                    Forest_cover_2012 = Forest_cover_2012/900,
                    Size = Size,
                    Restoration_1985 = Restoration_1985/900,
                    Restoration_2012 = Restoration_2012/900,
                    Distance_roads = Distance_roads/30
                    )


## prepare subset data
fpsub <- forest_points_subset
x <- t(as.data.frame(lapply(fpsub$geom, as.numeric)))
fpsub <- fpsub %>% select(Forest_cover_1985, #z
                    Forest_cover_2012,
                    Percent_forest_1985,
                    Percent_forest_2012,
                    Type, #z
                    name_muni, #z
                    Size, #z
                    Distance_roads, #z
                    Avg_precip, #z
                    Restoration_1985, #y
                    Restoration_2012,
                    Percent_restoration_1985, #y
                    Percent_restoration_2012
)
fpsub <- fpsub %>% mutate(x1 = x[,1], x2 = x[,2])
fpsub <- distinct(fpsub, x1, x2, .keep_all=TRUE) %>% na.omit()
fpsub <- fpsub %>% mutate(Forest_cover_1985 = Forest_cover_1985/900,
                    Forest_cover_2012 = Forest_cover_2012/900,
                    Size = Size,
                    Restoration_1985 = Restoration_1985/900,
                    Restoration_2012 = Restoration_2012/900,
                    Distance_roads = Distance_roads/1000
)

resps <- c('Restoration_1985',
           'Percent_restoration_1985',
           'Restoration_1985',
           'Restoration_2012',
           'Restoration_1985'
)

fc <- c('Forest_cover_1985', 'Percent_forest_1985', 'Percent_forest_1985', 'Forest_cover_2012', 'Forest_cover_1985')

models <- list()
for(i in 1:length(resps)){
  #i <- 1
  # create matrices and train test split
  if(i==5){
    fp <- fpsub
  }
  x <- fp %>% select(x1, x2)
  st_geometry(x) <- NULL
  names(x) <- NULL
  x <- data.matrix(x)
  y <- fp %>% .[[resps[i]]]
  z <- fp %>% select(fc[i], #Forest_cover_1985,
                     Size,
                     #Gini,
                     Distance_roads,
                     Avg_precip)# %>%
  # mutate(Type = fct_lump(as.factor(Type), n=1))
  z$geometry <- NULL
  z <- model.matrix( ~ . - 1, data=z)
  #z <- z[,-4]
  dimnames(z) <- NULL
  #summary(factor(z[,4]))
  z <- data.matrix(z)
  
  set.seed(12)
  train.ind <- sample(1:nrow(fp), size = floor(0.9 * nrow(fp)) )
  x.train <- x[train.ind,]
  y.train <- y[train.ind]
  z.train <- z[train.ind,]
  x.test <- x[-train.ind,]
  y.test <- y[-train.ind]
  z.test <- z[-train.ind,]
  
  ### fit linear model
  lm.dat <- data.frame(lon = x.train[,1],
                       lat = x.train[,2],
                       y = y.train,
                       ForestCover = z.train[,1],
                       Size = z.train[,2],
                       #Gini = z.train[,3],
                       Distance_roads = z.train[,3],
                       Avg_precip = z.train[,4])
  m1 <- lm(y ~ . , data=data.frame(lm.dat))
  summary(m1)
  
  to.pred <- data.frame(lon = x.test[,1],
                        lat = x.test[,2],
                        ForestCover = z.test[,1],
                        Size = z.test[,2],
                        #Gini = z.test[,3],
                        Distance_roads = z.test[,3],
                        Avg_precip = z.test[,4])
  m1.pred.dist <- predict(m1, to.pred, se.fit=TRUE)
  y.hat.m1 <- m1.pred.dist$fit
  y.hat.m1.var <- m1.pred.dist$se.fit^2
  
  
  ### fit lattice krig spatial model
  a.w.init <- 5
  LKinfo<- LKrigSetup( x.train,  nlevel=3, nu=0.5,
                       NC=6, NC.buffer=2, a.wght=a.w.init, # 4.001, 4.01, 4.1, 4.3, 4.5, 5
                       fixedFunctionArgs = list(m = 2) )
  m2 <- LatticeKrig(x=x.train, y=y.train, Z=z.train,
                    findAwght = TRUE,
                    LKinfo=LKinfo,
                    verbose=TRUE)
  
  y.hat.m2 <- predict(m2, xnew = x.test, Znew = z.test)
  y.hat.m2.var <- predictSE(m2, xnew = x.test, Znew = z.test)^2
  
  
  # report parameters
  betas <- m2$d.coef #* 107639
  process.var <- m2$MLE$summary[5]
  nugget.var <- m2$MLE$summary[4]^2
  lambda <- m2$MLE$summary[6]
  if(!is.na(m2$MLE$summary['a.wght.MLE'])){
    a.weight <- m2$MLE$summary['a.wght.MLE']
    sar.range <- 1/sqrt(m2$MLE$summary['a.wght.MLE'] - 4)
  }else{
    a.weight <- a.w
    sar.range <- 1/sqrt(a.w - 4)
  }
  #sar.range <- 1/sqrt(m2$MLE$summary[7] - 4)
  pars <-c(process.var, nugget.var, lambda, sar.range, a.weight)
  format(signif(pars, digits=3), scientific=FALSE)
  m2$eff.df
  
  # wald tests, lrt test
  # test fixed effects H_0 : beta_i = 0
  cov.m <- m2$rho.MLE * (m2$Omega)
  t.stats <- m2$d.coef / sqrt(diag(cov.m))
  p.vals <- 2*pt( -abs(t.stats), df=m2$eff.df) # or
  p.vals <- 2*pnorm( -abs(t.stats))
  round(p.vals, digits=5)
  p.vals < 0.05
  p.stars <- gtools::stars.pval(p.vals)
  data.frame( round(cbind(m2$d.coef, sqrt(diag(cov.m)), t.stats, p.vals ), digits=4), p.stars)
  # likelihood ratio test H_0 : rho=0
  lrt.stat <- -2*( c(logLik(m1)) - m2$lnProfileLike )
  p.val <- pchisq( lrt.stat, df=1, lower.tail = FALSE)
  round(p.val, digits=5)
  str(m2$MLE$lnLike.eval)
  
  
  # RSS, r^2
  m1.train.rmse <- sqrt(mean(m1$residuals^2))
  m2.train.rmse <- sqrt(mean(m2$residuals^2))
  
  # CV RMSE, how many true values fall in simulated CI using lm and LK
  m1.test.rmse <- sqrt(mean((y.test - y.hat.m1)^2))
  m2.test.rmse <- sqrt(mean((y.test - y.hat.m2)^2))
  # CRPS
  crps.normal <- function(mu, sigma2, x){
    ### mu, sigma2, and x are vectors of equal length
    ### mu and sigma2 are the mean and variance of the normal predictive distribution
    ### x is the actual realization that came to fruition
    s <- sqrt(sigma2)
    Z <- (x-mu)/s
    score <- -s*( (1/sqrt(pi)) - 2*dnorm(Z) - Z*( 2*pnorm(Z)-1 ) )
    return(score)
  }
  m1.crps <- crps.normal(mu=y.hat.m1, sigma2=y.hat.m1.var, x=y.test)
  m2.crps <- crps.normal(mu=y.hat.m2, sigma2=y.hat.m2.var, x=y.test)
  m1.crps.mean <- mean(m1.crps)
  m2.crps.mean <- mean(m2.crps)
  
  summary.stats <- matrix( c(m1.train.rmse, m2.train.rmse,
                             m1.test.rmse, m2.test.rmse,
                             m1.crps.mean, m2.crps.mean), nrow=3, ncol=2,
                           dimnames=list(c('Train RMSE', 'Test RMSE', 'Test CRPS'),
                                         c('LM', 'LK')))
  
  # spatial model
  linregxy <- predict( m2, drop.Z = TRUE, just.fixed=TRUE)
  linregxyz <- predict( m2, drop.Z = FALSE, just.fixed=TRUE)
  spautoxy <- predict( m2, drop.Z = TRUE, just.fixed=FALSE)
  spautoxyz <- predict( m2, drop.Z = FALSE, just.fixed=FALSE)
  
  fdf <- data.frame(lon = x.train[,1],
                    lat = x.train[,2],
                    y = y.train,
                    ForestCover = z.train[,1],
                    Size = z.train[,2],
                    #Gini = z.train[,3],
                    Distance_roads = z.train[,3],
                    Avg_precip  = z.train[,4],
                    linregxy = linregxy,
                    linregxyz = linregxyz,
                    spautoxy = spautoxy,
                    spautoxyz = spautoxyz,
                    fitted = m2$fitted.values,
                    resids = m2$residuals,
                    lin.mod.preds = predict(m1),
                    row.names = NULL)
  
  ggplot(fdf) + geom_histogram(aes(resids), bins=100) + lims(x=c(-30,30))
  ggplot(fdf) + geom_point(aes(fitted, resids))
  ggplot(fdf) + geom_point(aes(lon, lat, color=resids)) + scale_colour_gradientn(colors=tim.colors(),limits = c(-30,30))
  ggplot(fdf) + geom_point(aes(lon, lat, color=fitted)) + scale_colour_gradientn(colors=tim.colors(),limits = c(-30,30))
  
  
  
  s1 <- cbind(betas, p.vals)
  
  models[[i]] <- list(linmod = m1, lkmod = m2,
                      data = fdf,
                      summstats = summary.stats,
                      tbl=s1, vars=paste(fc[i], resps[i], sep=' - '))
  
}





sslist <- list()
for(i in 1:length(models)){
  sslist[[i]] <- models[[i]]$summstats
}



i <- 1
s <- models[[i]]$lkmod
prs <- s$parameters
mles <- s$MLE$summary
tbl <- list()
tbl[[i]] <- round(models[[i]]$tbl, digits=2)
for(i in 2:5){
  s <- models[[i]]$lkmod
  prs <- cbind(prs, s$parameters[,2])
  mles <- rbind(mles, s$MLE$summary)
  tbl[[i]] <- round(models[[i]]$tbl, digits=2)
}


### View parameter values

i <- 1
m1 <- models[[i]]$linmod
fdf <- models[[i]]$data
m2 <- models[[i]]$lkmod

# set i <- 1 to get model 1
  
# report parameters
betas <- m2$d.coef #* 107639
process.var <- m2$MLE$summary[5]
nugget.var <- m2$MLE$summary[4]^2
lambda <- m2$MLE$summary[6]
if(!is.na(m2$MLE$summary['a.wght.MLE'])){
  a.weight <- m2$MLE$summary['a.wght.MLE']
  sar.range <- 1/sqrt(m2$MLE$summary['a.wght.MLE'] - 4)
}else{
  a.weight <- a.w
  sar.range <- 1/sqrt(a.w - 4)
}
#sar.range <- 1/sqrt(m2$MLE$summary[7] - 4)
pars <-c(process.var, nugget.var, lambda, sar.range, a.weight)
format(signif(pars, digits=3), scientific=FALSE)
m2$eff.df

# wald tests, lrt test
# test fixed effects H_0 : beta_i = 0
cov.m <- m2$rho.MLE * (m2$Omega)
t.stats <- m2$d.coef / sqrt(diag(cov.m))
p.vals <- 2*pt( -abs(t.stats), df=m2$eff.df) # or
p.vals <- 2*pnorm( -abs(t.stats))
round(p.vals, digits=5)
p.vals < 0.05
p.stars <- gtools::stars.pval(p.vals)
data.frame( round(cbind(m2$d.coef, sqrt(diag(cov.m)), t.stats, p.vals ), digits=4), p.stars)
# likelihood ratio test H_0 : rho=0
lrt.stat <- -2*( c(logLik(m1)) - m2$lnProfileLike )
p.val <- pchisq( lrt.stat, df=1, lower.tail = FALSE)
round(p.val, digits=5)
str(m2$MLE$lnLike.eval)



