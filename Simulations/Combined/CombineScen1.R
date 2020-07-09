########################################################################################
# Reproduce simulation results from
# "Modeling semi-competing risks data as a longitudinal bivariate process"
# by Nevo et al.
# This script combines all simulation results for scenario 1, and output a result file
# ResScen1.csv
########################################################################################
rm(list = ls())
library(reshape2)
library(xtable)
library(ggplot2)
library(dplyr)
library(Daniel)
# Number of columns:
# censoring rate + sample size + knots + lambda + beta.y + 15*3 time.dep intercepts + 2*3 betas + 
# corresponding 1 + 45 + 6 beta SEs
# Number of rows:
# 1000 iterations for each:  3 sample sizes, 3 censoring rates. For the B-spline estimators, 6 lambda values, 
scen1.all.results.J5 <- matrix(nr = 1000*6*3*3, nc = 1 + 1 + 1 + 1 + 1 + 14*3 + 6 + 1 + 6)
scen1.all.results.J10 <- matrix(nr = 1000*6*3*3, nc = 1 + 1 + 1 + 1 + 1 + 14*3 + 6 + 1 + 6)
scen1.all.results.unstruct <- matrix(nr = 1000*3*3, nc = 1+ 1 + 1 + 1 + 1 + 14*3 + 6 + 1 + 6) 

setwd("Simulations/Results/Scenario 1/")

my.files.Scen1 <- list.files()

free <- free.unstruct <- 1
for(ii in 1:length(my.files.Scen1))
{
  load(paste0(my.files.Scen1[ii]))
  # The same as in LongitSCtimedep, the following lines constructs the baseline time-varying functions
  smooth.aux5 <- mgcv::smooth.construct.ps.smooth.spec(mgcv::s(times,bs="ps", k = 5), data = list(times = times),
                                                       knots = list(times = c(min(times),max(times))))
  smooth.aux10 <- mgcv::smooth.construct.ps.smooth.spec(mgcv::s(times,bs="ps", k = 10), data = list(times = times),
                                                        knots = list(times = c(min(times),max(times))))
  S.penal5 <- smooth.aux5$S[[1]]
  Bsplines5 <- smooth.aux5$X
  S.penal10 <- smooth.aux10$S[[1]]
  Bsplines10 <- smooth.aux10$X
  Q5 <- ncol(Bsplines5)
  Q10 <- ncol(Bsplines10)
  ###
  for (jj in 1:n.lambdas)
  {
    n.iter <- dim(est.BS5)[1]
    start <- free
    end <- free + n.iter - 1
    scen1.all.results.J5[start:end, 1] <- ifelse(substr(paste0(my.files.Scen1[ii]), 10, 11)=="00", 0, # "No Censoring"
                                           ifelse(substr(paste0(my.files.Scen1[ii]), 10, 11)=="20", 1, # "Medium Censoring",
                                                  ifelse(substr(paste0(my.files.Scen1[ii]), 10, 11)=="30", 
                                                         2, stop("Nooo!")))) #"High Censoring"
    scen1.all.results.J5[start:end, 2] <- n.sample  
    scen1.all.results.J5[start:end, 3] <- 5
    scen1.all.results.J5[start:end, 4] <- lambdas[jj]
    scen1.all.results.J5[start:end, 5] <- est.BS5[, jj, 1] # beta.y.est
    scen1.all.results.J5[start:end, 6] <- est.BS5[, jj, 1 + 3*Q5 + 1] # betaNT1.est
    scen1.all.results.J5[start:end, 7] <- est.BS5[, jj, 1 + 3*Q5 + 2] # betaNT2.est
    scen1.all.results.J5[start:end, 8] <- est.BS5[, jj, 1 + 3*Q5 + 3] # betaT1.est
    scen1.all.results.J5[start:end, 9] <- est.BS5[, jj, 1 + 3*Q5 + 4] # betaT2.est
    scen1.all.results.J5[start:end, 10] <- est.BS5[, jj, 1 + 3*Q5 + 5] # betaOR1.est
    scen1.all.results.J5[start:end, 11] <- est.BS5[, jj, 1 + 3*Q5 + 6] # betaOR2.est
    
    scen1.all.results.J5[start:end, 12] <- se.BS5[, jj, 1] # beta.y.est.se
    scen1.all.results.J5[start:end, 13] <- se.BS5[, jj, 1 + 3*Q5 + 1] # betaNT1.est.se
    scen1.all.results.J5[start:end, 14] <- se.BS5[, jj, 1 + 3*Q5 + 2] # betaNT2.est.se
    scen1.all.results.J5[start:end, 15] <- se.BS5[, jj, 1 + 3*Q5 + 3] # betaT1.est.se
    scen1.all.results.J5[start:end, 16] <- se.BS5[, jj, 1 + 3*Q5 + 4] # betaT2.est.se
    scen1.all.results.J5[start:end, 17] <- se.BS5[, jj, 1 + 3*Q5 + 5] # betaOR1.est.se
    scen1.all.results.J5[start:end, 18] <- se.BS5[, jj, 1 + 3*Q5 + 6] # betaOR2.est.se
    
    scen1.all.results.J10[start:end, 1] <- scen1.all.results.J5[start:end, 1]
    scen1.all.results.J10[start:end, 2] <- n.sample
    scen1.all.results.J10[start:end, 3] <- 10
    scen1.all.results.J10[start:end, 4] <- lambdas[jj]
    scen1.all.results.J10[start:end, 5] <- est.BS10[, jj, 1] # beta.y.est
    scen1.all.results.J10[start:end, 6] <- est.BS10[, jj, 1 + 3*Q10 + 1] # betaNT1.est
    scen1.all.results.J10[start:end, 7] <- est.BS10[, jj, 1 + 3*Q10 + 2] # betaNT2.est
    scen1.all.results.J10[start:end, 8] <- est.BS10[, jj, 1 + 3*Q10 + 3] # betaT1.est
    scen1.all.results.J10[start:end, 9] <- est.BS10[, jj, 1 + 3*Q10 + 4] # betaT2.est
    scen1.all.results.J10[start:end, 10] <- est.BS10[, jj, 1 + 3*Q10 + 5] # betaOR1.est
    scen1.all.results.J10[start:end, 11] <- est.BS10[, jj, 1 + 3*Q10 + 6] # betaOR2.est
    
    scen1.all.results.J10[start:end, 12] <- se.BS10[, jj, 1] # beta.y.est.se
    scen1.all.results.J10[start:end, 13] <- se.BS10[, jj, 1 + 3*Q10 + 1] # betaNT1.est.se
    scen1.all.results.J10[start:end, 14] <- se.BS10[, jj, 1 + 3*Q10 + 2] # betaNT2.est.se
    scen1.all.results.J10[start:end, 15] <- se.BS10[, jj, 1 + 3*Q10 + 3] # betaT1.est.se
    scen1.all.results.J10[start:end, 16] <- se.BS10[, jj, 1 + 3*Q10 + 4] # betaT2.est.se
    scen1.all.results.J10[start:end, 17] <- se.BS10[, jj, 1 + 3*Q10 + 5] # betaOR1.est.se
    scen1.all.results.J10[start:end, 18] <- se.BS10[, jj, 1 + 3*Q10 + 6] # betaOR2.est.se
    
    #  Keeping the baseline varying functions
    for (kk in 1:n.iter)
    {
      # K = 5
      scen1.all.results.J5[start + kk - 1 , (18 + 1):(18 + 14)] <-
        expit(Bsplines5%*%est.BS5[kk, jj, 2:(Q5 + 1)])  
      scen1.all.results.J5[start + kk - 1, (18 + 14 + 1):(18 + 14*2)] <-
        expit(Bsplines5%*%est.BS5[kk, jj, (1 + Q5 + 1):(1 + 2*Q5)])
      scen1.all.results.J5[start + kk - 1, (18 + 14*2 + 1):(18 + 14*3)] <-
        exp(Bsplines5%*%est.BS5[kk, jj, (1 + 2*Q5 + 1):(1 + 3*Q5)])
      #K = 10
      scen1.all.results.J10[start + kk - 1 , (18 + 1):(18 + 14)] <-
        expit(Bsplines10%*%est.BS10[kk, jj, 2:(Q10+1)])  
      scen1.all.results.J10[start + kk - 1, (18 + 14 + 1):(18 + 14*2)] <-
        expit(Bsplines10%*%est.BS10[kk, jj, (1 + Q10 + 1):(1 + 2*Q10)])
      scen1.all.results.J10[start + kk - 1, (18 + 14*2 + 1):(18 + 14*3)] <-
        exp(Bsplines10%*%est.BS10[kk, jj, (1 + 2*Q10 + 1):(1 + 3*Q10)])
    }
    free <- end + 1
  }
  start.unstruct <- free.unstruct
  end.unstruct <- free.unstruct + n.iter - 1
  ### Now unstruct estimators
  scen1.all.results.unstruct[start.unstruct:end.unstruct, 1] <-  ifelse(substr(paste0(my.files.Scen1[ii]), 10, 11)=="00", 0, # "No Censoring"
                                                         ifelse(substr(paste0(my.files.Scen1[ii]), 10, 11)=="20", 1, # "Medium Censoring",
                                                                ifelse(substr(paste0(my.files.Scen1[ii]), 10, 11)=="30", 2,stop("Nooo!")))) #"High Censoring"
  scen1.all.results.unstruct[start.unstruct:end.unstruct, 2] <- n.sample
  scen1.all.results.unstruct[start.unstruct:end.unstruct, 3] <- 0
  scen1.all.results.unstruct[start.unstruct:end.unstruct, 4] <- 0
  scen1.all.results.unstruct[start.unstruct:end.unstruct, 5] <- est.param[, 1 + 3*14] # beta.y.est
  scen1.all.results.unstruct[start.unstruct:end.unstruct, 6] <- est.param[, 2 + 3*14] # betaNT1.est
  scen1.all.results.unstruct[start.unstruct:end.unstruct, 7] <- est.param[, 3 + 3*14] # betaNT2.est
  scen1.all.results.unstruct[start.unstruct:end.unstruct, 8] <- est.param[, 4  + 3*14] # betaT1.est
  scen1.all.results.unstruct[start.unstruct:end.unstruct, 9] <- est.param[, 5 + 3*14] # betaT2.est
  scen1.all.results.unstruct[start.unstruct:end.unstruct, 10] <- est.param[, 6 + 3*14] # betaOR1.est
  scen1.all.results.unstruct[start.unstruct:end.unstruct, 11] <- est.param[, 7 + 3*14] # betaOR2.est
  
  scen1.all.results.unstruct[start.unstruct:end.unstruct, 12] <- se.param[, 1] # beta.y.est.se
  scen1.all.results.unstruct[start.unstruct:end.unstruct, 13] <- se.param[, 2 + 3*14] # betaNT1.est.se
  scen1.all.results.unstruct[start.unstruct:end.unstruct, 14] <- se.param[, 3 + 3*14] # betaNT2.est.se
  scen1.all.results.unstruct[start.unstruct:end.unstruct, 15] <- se.param[, 4 + 3*14] # betaT1.est.se
  scen1.all.results.unstruct[start.unstruct:end.unstruct, 16] <- se.param[, 5 + 3*14] # betaT2.est.se
  scen1.all.results.unstruct[start.unstruct:end.unstruct, 17] <- se.param[, 6 + 3*14] # betaOR1.est.se
  scen1.all.results.unstruct[start.unstruct:end.unstruct, 18] <- se.param[, 7 + 3*14] # betaOR2.est.se   

  scen1.all.results.unstruct[start.unstruct:end.unstruct, (18 + 1):(18 + 14)] <-  expit(est.param[, 2:(2 + 13)])  
  scen1.all.results.unstruct[start.unstruct:end.unstruct, (18 + 14 + 1):(18 + 14*2)] <-  expit(est.param[, (2 + 14):(2 + 14 + 13)])
  scen1.all.results.unstruct[start.unstruct:end.unstruct, (18 + 14*2 + 1):(18 + 14*3)] <- exp(est.param[, (2 + 2*14):(2 + 2*14 + 13)])
  free.unstruct <- end.unstruct + 1
  }

colnames(scen1.all.results.J5) <- colnames(scen1.all.results.J10) <- colnames(scen1.all.results.unstruct) <- 
  c("censoring","n.sample", "n.knots", "lambda", "beta.y.est", "betaNT1.est" ,"betaNT2.est", 
    "betaT1.est", "betaT2.est", "betaOR1.est", "betaOR2.est", "beta.y.se", "betaNT1.se", 
    "betaNT2.se", "betaT1.se", "betaT2.se", "betaOR1.se", "betaOR2.se", 
    paste0("base.prob.NT", 1:14), paste0("base.prob.T", 1:14), paste0("base.OR", 1:14))#,


df <- rbind(scen1.all.results.J5, scen1.all.results.J10, scen1.all.results.unstruct) %>% as.data.frame 
df <- df %>% mutate(Method = ifelse(n.knots==10, "B-splines, 10 knots \n ", 
                                    ifelse(n.knots==5, "B-splines, 5 knots \n ", 
                                           "Non B-spline")),
                    Penal = ifelse(lambda > 0, "Penalized (lambda=5)", "Non-Penalized" ))

setwd("../../Combined/")
write.csv(df, "ResScen1.csv")
