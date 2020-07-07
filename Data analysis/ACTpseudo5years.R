########################################################################################################
###### Analysis of the ACT psuedo data using longitudinal bivariate modeling
# This script mimics the main analysis presented in the paper:
# "Modeling semi-competing risks data as a longitudinal bivariate process"
# By Nevo et al.
# for interval length of 5 years
# NOTE: This psuedo data set sole use it to demonstrate how was the ACT data analyzed
# The results obtained from analyzing this dataset are not the same as those reported in the paper
# for the ACT data
########################################################################################################
set.seed(17)
rm(list = ls())

library(xtable)
library(dplyr)
library(ggplot2)
library(LongitSemiComp)

setwd("Data analysis/")

ACTpseudo <- read.table("ACTpseudo.csv")

########################################################################################################
##### Change the data to longitudinal data with time-frame of 5 years #####
inter.vec <- seq(65, 100, 5)
longit.data <- TimesToLongit(T1 = ACTpseudo$T1, T2 = ACTpseudo$T2, 
                             delta1 = ACTpseudo$delta1, delta2 = ACTpseudo$delta2,
                            inter.vec = inter.vec, TruncData = T, TruncTime = ACTpseudo$EnrollAge)

# ########################################################################################################
########################################################################################################
######  Bsplines Preperations ##### 
times <- seq(70, 100, 5)
smooth.aux <- mgcv::smooth.construct.ps.smooth.spec(mgcv::s(times,bs="ps", k = 5), 
                                                    data = list(times = times),
                                                    knots = list(times = c(min(times),max(times))))
Bsplines <- smooth.aux$X
n.knots <- 5
########################################################################################################
##### Covariate and model definitions ####
formula.NT <- ~ RaceWhite + CollGrad  + MaritalStatus + Depression + Female*APOE
formula.T <- ~  RaceWhite + CollGrad  + MaritalStatus + Depression + Female*APOE
formula.OR <- ~ Female*APOE
formula.inter <- ~ Female*APOE

#### Analyses under various penalty values #####
res0 <- LongitSemiComp:::LongitSC(longit.data = longit.data, formula.NT = formula.NT, formula.T = formula.T, 
                                  formula.OR = formula.OR, formula.inter = formula.inter, data = ACTpseudo,  
                                  lambda = 0, knots = n.knots)
# Output
res0

res0.1 <- LongitSemiComp:::LongitSC(longit.data = longit.data, formula.NT =  formula.NT, formula.T = formula.T, 
                                    formula.OR = formula.OR, formula.inter = formula.inter, data = ACTpseudo,  
                                    lambda = 0.1, knots = n.knots)

res0.5 <- LongitSemiComp:::LongitSC(longit.data = longit.data, formula.NT =  formula.NT, formula.T = formula.T, 
                                    formula.OR = formula.OR, formula.inter = formula.inter, data = ACTpseudo,  
                                    lambda = 0.5, knots = n.knots)

res1 <- LongitSemiComp:::LongitSC(longit.data = longit.data, formula.NT =  formula.NT, formula.T = formula.T, 
                                  formula.OR = formula.OR, formula.inter = formula.inter, data = ACTpseudo,  
                                  lambda = 1, knots = n.knots)

res2.5 <- LongitSemiComp:::LongitSC(longit.data = longit.data, formula.NT =  formula.NT, formula.T = formula.T, 
                                    formula.OR = formula.OR, formula.inter = formula.inter, data = ACTpseudo,  
                                    lambda = 2.5, knots = n.knots)

res5 <- LongitSemiComp:::LongitSC(longit.data = longit.data, formula.NT =  formula.NT, formula.T = formula.T, 
                                  formula.OR = formula.OR, formula.inter = formula.inter, data = ACTpseudo,  
                                  lambda = 5, knots = n.knots)

###### Analysis in the unstrucutred method 
res.unstruct <- LongitSCparam(longit.data = longit.data, formula.NT =  formula.NT, formula.T = formula.T, 
                              formula.OR = formula.OR, formula.inter = formula.inter, data = ACTpseudo)


##############################################################################################
##### AIC Table (analogous to Table A.17 in the SM) 
all.aic <- c(res0$aic, res0.1$aic, res0.5$aic, res1$aic, res2.5$aic, res5$aic)
A <- rbind(c(0, 0.1, 0.5, 1, 2.5, 5), all.aic) # aic is logL - 2#params
rownames(A) <- c("lambda", "AIC")
xtable(A)
##############################################################################################
 
####### Analyze results ####### 
# The LongitSC function returns point estimates and CIs, including for the time-dependent 
# intercept functions \alpha. Here we extract them for ease of plotting and 
##############################################################################################

##############################################################################################
##### Point-estimates for the time-dependent functions #####
# Non-terminal event
alpha.nt0 <- res0$time.int.NT
alpha.nt0.1 <- res0.1$time.int.NT
alpha.nt0.5 <- res0.5$time.int.NT
alpha.nt1 <- res1$time.int.NT
alpha.nt2.5 <- res2.5$time.int.NT
alpha.nt5 <- res5$time.int.NT
alpha.nt.unstruct <- res.unstruct$time.int.NT

# Terminal event
alpha.t0 <- res0$time.int.T
alpha.t0.1 <- res0.1$time.int.T
alpha.t0.5 <- res0.5$time.int.T
alpha.t1 <- res1$time.int.T
alpha.t2.5 <- res2.5$time.int.T
alpha.t5 <- res5$time.int.T
alpha.t.unstruct <- res.unstruct$time.int.T

# Odds ratios for the baseline levels
alpha.or0 <- res0$time.int.OR
alpha.or0.1 <- res0.1$time.int.OR
alpha.or0.5 <- res0.5$time.int.OR
alpha.or1 <- res1$time.int.OR
alpha.or2.5 <- res2.5$time.int.OR
alpha.or5 <- res5$time.int.OR
alpha.or.unstruct <- res.unstruct$time.int.OR

##############################################################################################
####### Confidence intervals for the time-varying functions
# Non-terminal event
ci.L.alpha.nt0 <- res0$ci.L.alpha.NT
ci.L.alpha.nt0.1 <- res0.1$ci.L.alpha.NT
ci.L.alpha.nt0.5 <- res0.5$ci.L.alpha.NT
ci.L.alpha.nt1 <- res1$ci.L.alpha.NT
ci.L.alpha.nt2.5 <- res2.5$ci.L.alpha.NT
ci.L.alpha.nt5 <-  res5$ci.L.alpha.NT
ci.L.alpha.nt.unstruct <-  res.unstruct$ci.L.alpha.NT
ci.H.alpha.nt0 <- res0$ci.H.alpha.NT
ci.H.alpha.nt0.1 <- res0.1$ci.H.alpha.NT
ci.H.alpha.nt0.5 <- res0.5$ci.H.alpha.NT
ci.H.alpha.nt1 <- res1$ci.H.alpha.NT
ci.H.alpha.nt2.5 <- res2.5$ci.H.alpha.NT
ci.H.alpha.nt5 <-  res5$ci.H.alpha.NT
ci.H.alpha.nt.unstruct <-  res.unstruct$ci.H.alpha.NT

# Terminal event
ci.L.alpha.t0 <- res0$ci.L.alpha.T
ci.L.alpha.t0.1 <- res0.1$ci.L.alpha.T
ci.L.alpha.t0.5 <- res0.5$ci.L.alpha.T
ci.L.alpha.t1 <- res1$ci.L.alpha.T
ci.L.alpha.t2.5 <- res2.5$ci.L.alpha.T
ci.L.alpha.t5 <-  res5$ci.L.alpha.T
ci.L.alpha.t.unstruct <-  res.unstruct$ci.L.alpha.T
ci.H.alpha.t0 <- res0$ci.H.alpha.T
ci.H.alpha.t0.1 <- res0.1$ci.H.alpha.T
ci.H.alpha.t0.5 <- res0.5$ci.H.alpha.T
ci.H.alpha.t1 <- res1$ci.H.alpha.T
ci.H.alpha.t2.5 <- res2.5$ci.H.alpha.T
ci.H.alpha.t5 <-  res5$ci.H.alpha.T
ci.H.alpha.t.unstruct <-  res.unstruct$ci.H.alpha.T

# Odds ratio 
ci.L.alpha.or0 <- res0$ci.L.alpha.OR
ci.L.alpha.or0.1 <- res0.1$ci.L.alpha.OR
ci.L.alpha.or0.5 <- res0.5$ci.L.alpha.OR
ci.L.alpha.or1 <- res1$ci.L.alpha.OR
ci.L.alpha.or2.5 <- res2.5$ci.L.alpha.OR
ci.L.alpha.or5 <-  res5$ci.L.alpha.OR
ci.L.alpha.or.unstruct <-  res.unstruct$ci.L.alpha.OR
ci.H.alpha.or0 <- res0$ci.H.alpha.OR
ci.H.alpha.or0.1 <- res0.1$ci.H.alpha.OR
ci.H.alpha.or0.5 <- res0.5$ci.H.alpha.OR
ci.H.alpha.or1 <- res1$ci.H.alpha.OR
ci.H.alpha.or2.5 <- res2.5$ci.H.alpha.OR
ci.H.alpha.or5 <-  res5$ci.H.alpha.OR
ci.H.alpha.or.unstruct <-  res.unstruct$ci.H.alpha.OR

##############################################################################################
#### Point estimates for the betas
# Non-terminal
beta.nt0 <- res0$coef.NT
beta.nt0.1 <- res0.1$coef.NT
beta.nt0.5 <- res0.5$coef.NT
beta.nt1 <- res1$coef.NT
beta.nt2.5 <- res2.5$coef.NT
beta.nt5 <- res5$coef.NT
beta.nt.unstruct <- res.unstruct$coef.NT
# Terminal
beta.t0 <- res0$coef.T
beta.t0.1 <- res0.1$coef.T
beta.t0.5 <- res0.5$coef.T
beta.t1 <- res1$coef.T
beta.t2.5 <- res2.5$coef.T
beta.t5 <- res5$coef.T
beta.t.unstruct <- res.unstruct$coef.T
# Odds ratio
beta.or0 <- res0$coef.OR
beta.or0.1 <- res0.1$coef.OR
beta.or0.5 <- res0.5$coef.OR
beta.or1 <- res1$coef.OR
beta.or2.5 <- res2.5$coef.OR
beta.or5 <- res5$coef.OR
beta.or.unstruct <- res.unstruct$coef.OR
# Interactions
beta.inter0 <- res0$coef.inter
beta.inter0.1 <- res0.1$coef.inter
beta.inter0.5 <- res0.5$coef.inter
beta.inter1 <- res1$coef.inter
beta.inter2.5 <- res2.5$coef.inter
beta.inter5 <- res5$coef.inter
beta.inter.unstruct <- res.unstruct$coef.inter

#### SE estimates for the (non-exponentiated) betas
# Non-terminal
se.beta.nt0 <- res0$se.coef.rob.NT
se.beta.nt0.1 <- res0.1$se.coef.rob.NT
se.beta.nt0.5 <- res0.5$se.coef.rob.NT
se.beta.nt1 <- res1$se.coef.rob.NT
se.beta.nt2.5 <- res2.5$se.coef.rob.NT
se.beta.nt5 <- res5$se.coef.rob.NT
se.beta.nt.struct <- res.unstruct$se.coef.NT

# Terminal
se.beta.t0 <- res0$se.coef.rob.T
se.beta.t0.1 <- res0.1$se.coef.rob.T
se.beta.t0.5 <- res0.5$se.coef.rob.T
se.beta.t1 <- res1$se.coef.rob.T
se.beta.t2.5 <- res2.5$se.coef.rob.T
se.beta.t5 <- res5$se.coef.rob.T
se.beta.t.struct <- res.unstruct$se.coef.T

# Odds ratio
se.beta.or0 <- res0$se.coef.rob.OR
se.beta.or0.1 <- res0.1$se.coef.rob.OR
se.beta.or0.5 <- res0.5$se.coef.rob.OR
se.beta.or1 <- res1$se.coef.rob.OR
se.beta.or2.5 <- res2.5$se.coef.rob.OR
se.beta.or5 <- res5$se.coef.rob.OR
se.beta.or.struct <- res.unstruct$se.coef.OR

# Interactions
se.beta.inter0 <- res0$se.coef.rob.inter
se.beta.inter0.1 <- res0.1$se.coef.rob.inter
se.beta.inter0.5 <- res0.5$se.coef.rob.inter
se.beta.inter1 <- res1$se.coef.rob.inter
se.beta.inter2.5 <- res2.5$se.coef.rob.inter
se.beta.inter5 <- res5$se.coef.rob.inter
se.beta.inter.struct <- res.unstruct$se.coef.inter
#### Arrange in data frames ####

df.alphas <- data.frame(alphas = c(alpha.nt0, alpha.nt0.1,  alpha.nt0.5, alpha.nt1, alpha.nt2.5, alpha.nt5, alpha.nt.unstruct,
                                   alpha.t0, alpha.t0.1,  alpha.t0.5, alpha.t1, alpha.t2.5, alpha.t5, alpha.t.unstruct,
                                   alpha.or0, alpha.or0.1, alpha.or0.5, alpha.or1, alpha.or2.5, alpha.or5, alpha.or.unstruct),
                        ci.L = c(ci.L.alpha.nt0, ci.L.alpha.nt0.1,  ci.L.alpha.nt0.5, ci.L.alpha.nt1, 
                                 ci.L.alpha.nt2.5, ci.L.alpha.nt5, ci.L.alpha.nt.unstruct,
                                 ci.L.alpha.t0, ci.L.alpha.t0.1,  ci.L.alpha.t0.5, ci.L.alpha.t1, 
                                 ci.L.alpha.t2.5, ci.L.alpha.t5, ci.L.alpha.t.unstruct,
                                 ci.L.alpha.or0, ci.L.alpha.or0.1, ci.L.alpha.or0.5, ci.L.alpha.or1, 
                                 ci.L.alpha.or2.5, ci.L.alpha.or5, ci.L.alpha.or.unstruct),
                        ci.U = c(ci.H.alpha.nt0, ci.H.alpha.nt0.1,  ci.H.alpha.nt0.5, ci.H.alpha.nt1, 
                                 ci.H.alpha.nt2.5, ci.H.alpha.nt5, ci.H.alpha.nt.unstruct,
                                 ci.H.alpha.t0, ci.H.alpha.t0.1,  ci.H.alpha.t0.5, ci.H.alpha.t1, 
                                 ci.H.alpha.t2.5, ci.H.alpha.t5,ci.H.alpha.t.unstruct,
                                 ci.H.alpha.or0, ci.H.alpha.or0.1, ci.H.alpha.or0.5, ci.H.alpha.or1, 
                                 ci.H.alpha.or2.5, ci.H.alpha.or5, ci.H.alpha.or.unstruct),
                        times = rep(times, 3*7),
                        type = rep(c("AD", "Death", "OR"), each = 7*length(times)),
                        lambda = rep(rep(c(0, 0.1,  0.5, 1, 2.5, 5, "Unstructured"), each = length(times)), 3))

# Make a factor version of lambda
df.alphas$lambda.fact <- factor(df.alphas$lambda, levels = c("Unstructured", 0, 0.1, 0.5, 1, 2.5, 5),
                                         labels = c("Unstructured", paste("lambda==", c(0, 0.1, 0.5, 1, 2.5, 5))))
### Add a column with the middle of the intervals as the time instead of their end-point
df.alphas$times.mid <- df.alphas$times - 2.5/2

 
# Arrange covariate coefficients in a data.frame
pNT <- length(names(res0$coef.NT))
pT <-  length(names(res0$coef.T))
pOR <- length(names(res0$coef.OR))
pInter <- length(names(res0$coef.inter))
df.betas <- data.frame(betas = c(beta.nt0, beta.nt0.1, beta.nt0.5, beta.nt1, beta.nt2.5, beta.nt5, beta.nt.unstruct,
                                 beta.t0, beta.t0.1, beta.t0.5, beta.t1, beta.t2.5, beta.t5, beta.t.unstruct,
                                 beta.or0, beta.or0.1, beta.or0.5, beta.or1, beta.or2.5, beta.or5, beta.or.unstruct,
                                 beta.inter0, beta.inter0.1, beta.inter0.5, beta.inter1, beta.inter2.5, beta.inter5, 
                                 beta.inter.unstruct),
                       se = c(se.beta.nt0, se.beta.nt0.1, se.beta.nt0.5, se.beta.nt1, se.beta.nt2.5, se.beta.nt5, 
                              se.beta.nt.struct,
                              se.beta.t0, se.beta.t0.1, se.beta.t0.5, se.beta.t1, se.beta.t2.5, se.beta.t5,
                              se.beta.t.struct,
                              se.beta.or0, se.beta.or0.1, se.beta.or0.5, se.beta.or1, se.beta.or2.5, se.beta.or5,
                              se.beta.or.struct,
                              se.beta.inter0, se.beta.inter0.1, se.beta.inter0.5, se.beta.inter1, se.beta.inter2.5, se.beta.inter5,
                              se.beta.inter.struct),
                       var = c(rep(names(beta.nt0), 7), rep(names(beta.t0), 7), 
                               rep(names(beta.or0), 7), rep(paste0(names(beta.inter0)), 7)),
                       type = rep(c("Alzheimer", "Death", "OR", "Death"), times = 7*c(pNT, pT, pOR, pInter)),
                       lambda = c(rep(c(0, 0.1,  0.5, 1, 2.5, 5, "Unstructured"), each = pNT), 
                                  rep(c(0, 0.1,  0.5, 1, 2.5, 5, "Unstructured"), each = pT),
                                  rep(c(0, 0.1,  0.5, 1, 2.5, 5, "Unstructured"), each = pOR),
                                  rep(c(0, 0.1,  0.5, 1, 2.5, 5, "Unstructured"), each = pInter)))
########################################################################################################################
########## Figures analogous of the figures in the paper ##########
########################################################################################################################
#### Truncate/censor large values for plotting
#df.alphas$ci.U[df.alphas$type != "OR" & df.alphas$ci.U > 0.6] <- 0.6
df.alphas$ci.U[df.alphas$type == "OR" & df.alphas$ci.U > 3] <- 3
 

### Analogous to Figure 3 of the paper ####
my_breaks <- function(x) { if (max(x) < 1) seq(0, 0.3, 0.1) else seq(0, 3, 1) }
ggplot(filter(df.alphas, lambda %in% c("Unstructured",0,  2.5)) , aes(x = times.mid, y = alphas)) + 
  theme_bw() + facet_grid(type~lambda.fact, labeller  = label_parsed, scales = "free") + geom_point(size = 2) +
  labs(x = "Age", y = "Probability/OR conditioned on no events") +# ylim(c(0, 0.6)) + 
  theme(text = element_text(size = 23), panel.spacing = unit(1.25, "lines")) +
  geom_errorbar(aes(ymin = ci.L, ymax = ci.U), alpha = 0.4)  + 
  geom_point(data = data.frame(times.mid = 80, lambda = 0, alphas = 0.3), aes(x = times.mid, y = alphas),
             color = "white") +
  scale_y_continuous(breaks = my_breaks)


### And for the supplmentary, for all lambdas (analogous to Figure A.8)
my_breaks <- function(x) { if (max(x) < 1) seq(0, 0.3, 0.1) else seq(0, 3, 1) }
ggplot(df.alphas , aes(x = times.mid, y = alphas)) + 
  theme_bw() + facet_grid(type~lambda.fact, labeller  = label_parsed, scales = "free") + geom_point(size = 2) +
  labs(x = "Age", y = "Probability/OR conditioned on no events") +# ylim(c(0, 0.6)) + 
  theme(text = element_text(size = 20), panel.spacing = unit(1.25, "lines")) +
  geom_errorbar(aes(ymin = ci.L, ymax = ci.U), alpha = 0.4)  + 
  geom_point(data = data.frame(times.mid = 80, lambda = 0, alphas = 0.3), aes(x = times.mid, y = alphas),
             color = "white") +
  scale_y_continuous(breaks = my_breaks)


########################################################################################################################
##########  Tables for the paper ##########
########################################################################################################################

######################################################################################################################## 
df.betas <- df.betas %>% mutate(ci = paste0("(",round(exp((betas) - 1.96*se), 2), ", ", 
                                            round(exp((betas) + 1.96*se), 2), ")"),
                                Exp.est = round(exp(betas),2))



###### Tables of the coefficents beta for the paper and SM
betas.table <- df.betas
betas.table$betas <- betas.table$betas %>% round(2)
betas.table$se <- betas.table$se %>% round(3)

betas.out <- betas.table[c("lambda","type", "var", "Exp.est", "ci")]# %>% 

main.vars <- c("Female", "APOE", "Female:APOE", "NT-event", "NT-Female:gender", "NT-event:APOE", "NT-event:Female:APOE")
# Analogous to Table 2 (2nd, 4th and 6th columns)
betas.out %>% filter(lambda %in% c(2.5, "Unstructured") & (var %in% main.vars)) %>% arrange(lambda,type)
# Analogous to Table A.18
betas.out %>% filter(lambda %in% c(2.5, "Unstructured")) %>% arrange(lambda,type)



