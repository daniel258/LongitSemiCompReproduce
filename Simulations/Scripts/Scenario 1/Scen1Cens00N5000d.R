######################################################################################################
### Simulation Script for the paper: 
# Modeling semi-competing risks data as a longitudinal bivariate process by Nevo et al.
######################################################################################################
# This script runs 1000 iterations of simulating the dataset and estimation using different model specifications.
# Scenario: (I) The null scenario
# Censoring: 0%
# Sample size: 5000

rm(list = ls())
library(LongitSemiComp)
set.seed(15124)

# Number of intervals
times <- seq(1,14,1)
## True parameter values
alpha.nt <- LongitSemiComp:::logit(times*0.005  + 0.005*(times-2)^2 - (0.0002*(times + 1)^3) + 0.005)
alpha.t <- LongitSemiComp:::logit(times*0.0075  + 0.001*(times^2)  + 0.03)
alpha.or <- 0.9 + 0.175*times - 0.02*times^2 
alpha.or[times >= 13] <- 0

beta.nt <- log(c(1, 1))
beta.t <- log(c(1, 1))
beta.or <- log(c(1, 1))
beta.y <- log(1)
true.params <- c(beta.y, alpha.nt, alpha.t, alpha.or, beta.nt, beta.or, beta.t)
p <- 2
n.sample <- 5000
n.sim <- 100
n.params.BS5 <- 1 + 3*5 + 3*p
n.params.BS10 <- 1 + 3*10 + 3*p
n.params.param <- 1 + 3*14 + 3*p
pars.init.BS5 <- rep(-0.1, n.params.BS5)
pars.init.BS10 <- rep(-0.1, n.params.BS10)
pars.init.param <- rep(-0.1, n.params.param)
lambdas <- c(0, 0.1, 1, 5, 10, 25)
n.lambdas <- length(lambdas)
est.BS5 <- se.BS5 <- array(dim = c(n.sim,n.lambdas, n.params.BS5))
est.BS10 <- se.BS10 <- array(dim = c(n.sim,n.lambdas, n.params.BS10))
est.param <- se.param <- matrix(nr = n.sim, nc =  n.params.param)
conver.BS5 <- err.BS5 <- df.BS5 <- aic.BS5  <- matrix(nr = n.sim, nc = n.lambdas)
conver.BS10 <- err.BS10 <- df.BS10 <- aic.BS10  <- matrix(nr = n.sim, nc = n.lambdas)
conver.param <- err.param <- df.param <- aic.param  <- cens.rate <- vector(length = n.sim)
status.mat <- matrix(nr = n.sim, nc = 4) # 4 "types" of statuses; see FindStatusTimeDep.R

for (i in 1:n.sim)
{
  a <- proc.time()
  cat("i =", i, "\n")
  my.data <- SimLongitDataTimeDep(n.sample, times = times,  beta.y,  
                                  alpha.nt, alpha.t, alpha.or, 
                                  beta.nt, beta.t, beta.or, cens.poten.rate = 0)
  df.data <- my.data$df.return
  cens.rate[i] <- mean(my.data$cens)
  status.mat[i, ] <- LongitSemiComp:::FindStatusTimeDep(ID = df.data$ID, YNT = df.data$YNT, YT = df.data$YT)
  for (j in 1:n.lambdas)
  {
    my.lambda = lambdas[j]
    res.BS5 <- tryCatch(LongitSCtimeDep(times = times, 
                                       formula.NT = YNT ~ X.1 + X.2, 
                                       formula.T = YT ~ X.1 + X.2, 
                                       formula.OR = ~ X.1 + X.2, 
                                       data = df.data, epsOR = 10^(-10),
                                       knots = 5, lambda = my.lambda, 
                                       init = pars.init.BS5, 
                                       maxit.optim = 50000),
                         error=function(e) {e})
    if(inherits(res.BS5, "error")){
      err.BS5[i, j] <- 1
    } else {
      est.BS5[i, j, ] <- res.BS5$est
      se.BS5[i, j, ] <- res.BS5$se.rob
      conver.BS5[i, j] <- res.BS5$optim.conv
      df.BS5[i, j] <- res.BS5$df
      aic.BS5[i, j] <- res.BS5$aic
    }
    res.BS10 <- tryCatch(LongitSCtimeDep(times = times, 
                                        formula.NT = YNT ~ X.1 + X.2, 
                                        formula.T = YT ~ X.1 + X.2, 
                                        formula.OR = ~ X.1 + X.2, 
                                        data = df.data, epsOR = 10^(-10),
                                        knots = 10, lambda = my.lambda, 
                                        init = pars.init.BS10, 
                                        maxit.optim = 50000),
                        error=function(e) {e})
    if(inherits(res.BS10, "error")){
      err.BS10[i, j] <- 1
    } else {
      est.BS10[i, j, ] <- res.BS10$est
      se.BS10[i, j, ] <- res.BS10$se.rob
      conver.BS10[i, j] <- res.BS10$optim.conv
      df.BS10[i, j] <- res.BS10$df
      aic.BS10[i, j] <- res.BS10$aic
    }}
  res.param <- tryCatch(LongitSCparamTimeDep(times = times, 
                                             formula.NT = YNT ~ X.1 + X.2, 
                                             formula.T = YT ~ X.1 + X.2, 
                                             formula.OR = ~ X.1 + X.2, 
                                             data = df.data, epsOR = 10^(-10),
                                             init = pars.init.param, 
                                             maxit.optim = 50000),
                        error=function(e) {e})
  if(inherits(res.param, "error")){
    err.param[i] <- 1
  } else {
    est.param[i, ] <- res.param$est
    se.param[i, ] <- res.param$se.rob
    conver.param[i] <- res.param$optim.conv
    df.param[i] <- res.param$df
    aic.param[i] <- res.param$aic
  }
  print(proc.time() - a)
  gc()
  }

packageVersion("LongitSemiComp")

setwd("/home/dn84/LongitSemiComp/Results")
save.image("Scen1Cens00N5000d.RData")
