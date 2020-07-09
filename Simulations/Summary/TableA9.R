########################################################################################
# Table A.9
########################################################################################
rm(list = ls())
library(xtable)
library(dplyr)

df <- read.csv(file = "Simulations/Combined/ResScen1.csv",header = T)

## True parameter values
beta.nt <- log(c(1, 1))
beta.t <- log(c(1, 1))
beta.or <- log(c(1, 1))
beta.y <- log(1)

my.table <- df  %>% mutate(n.sample = factor(n.sample), lambda = factor(lambda)) %>% 
  filter(lambda %in% c(1) & n.knots %in% c(10)) %>%  # & censoring %in% c(2)) %>% 
  select(n.sample, censoring,n.knots,  lambda, beta.y.est, beta.y.se, 
         betaNT1.est, betaNT1.se, betaNT2.est, betaNT2.se,
         betaT1.est, betaT1.se, betaT2.est, betaT2.se,
         betaOR1.est, betaOR1.se, betaOR2.est, betaOR2.se) %>% 
  group_by(n.sample, censoring) %>%
  summarise(bias.beta.y = mean(beta.y.est - beta.y) %>% round(3),
            est.se.beta.y = mean(beta.y.se) %>% round(3),
            emp.se.beta.y = sd(beta.y.est) %>% round(3),
            ci.beta.y.in = mean((beta.y.est - qnorm(0.975)*beta.y.se)< beta.y & 
                                 (beta.y.est + qnorm(0.975)*beta.y.se) > beta.y) %>% round(3),
            bias.NT1 = mean(betaNT1.est - beta.nt[1]) %>% round(3),
            est.se.NT1 = mean(betaNT1.se) %>% round(3),
            emp.se.NT1 = sd(betaNT1.est) %>% round(3),
            ci.betaNT1.in = mean((betaNT1.est - qnorm(0.975)*betaNT1.se)< beta.nt[1] & 
                                   (betaNT1.est + qnorm(0.975)*betaNT1.se) > beta.nt[1]) %>% round(3),
             bias.NT2 = mean(betaNT2.est - beta.nt[2]) %>% round(3),
             est.se.NT2 = mean(betaNT2.se) %>% round(3),
             emp.se.NT2 = sd(betaNT2.est) %>% round(3),
             ci.betaNT2.in = mean((betaNT2.est - qnorm(0.975)*betaNT2.se)< beta.nt[2] & 
                                   (betaNT2.est + qnorm(0.975)*betaNT2.se) > beta.nt[2]) %>% round(3),
            bias.T1 = mean(betaT1.est- beta.t[1]) %>% round(3),
            est.se.T1 = mean(betaT1.se) %>% round(3),
            emp.se.T1 = sd(betaT1.est) %>% round(3),
            ci.betaT1.in = mean((betaT1.est - qnorm(0.975)*betaT1.se)< beta.t[1] & 
                                  (betaT1.est + qnorm(0.975)*betaT1.se) > beta.t[1]) %>% round(3),
            bias.T2 = mean(betaT2.est - beta.t[2]) %>% round(3),
            est.se.T2 = mean(betaT2.se) %>% round(3),
            emp.se.T2 = sd(betaT2.est) %>% round(3),
            ci.betaT2.in = mean((betaT2.est - qnorm(0.975)*betaT2.se)< beta.t[2] & 
                                   (betaT2.est + qnorm(0.975)*betaT2.se) > beta.t[2]) %>% round(3),
            bias.OR1 = mean(betaOR1.est - beta.or[1]) %>% round(3),
            est.se.OR1 = mean(betaOR1.se) %>% round(3),
            emp.se.OR1 = sd(betaOR1.est) %>% round(3),
            ci.betaOR1.in = mean((betaOR1.est - qnorm(0.975)*betaOR1.se)< beta.or[1] & 
                                   (betaOR1.est + qnorm(0.975)*betaOR1.se) > beta.or[1]) %>% round(3),
            bias.OR2 = mean(betaOR2.est - beta.or[2]) %>% round(3),
            est.se.OR2 = mean(betaOR2.se) %>% round(3),
            emp.se.OR2 = sd(betaOR2.est) %>% round(3),
            ci.betaOR2.in = mean((betaOR2.est - qnorm(0.975)*betaOR2.se)< beta.or[2] & 
                                   (betaOR2.est + qnorm(0.975)*betaOR2.se) > beta.or[2]) %>% round(3)#,
  ) 
t(my.table) %>% xtable

my.table.beta.y <- t(my.table)[c(1:2, 3:6),]
my.table.nt1 <- t(my.table)[c(1:2, 7:10),]
my.table.nt2 <- t(my.table)[c(1:2, 11:14),]
my.table.t1 <- t(my.table)[c(1:2, 15:18),]
my.table.t2 <- t(my.table)[c(1:2, 19:22),]
my.table.or1 <- t(my.table)[c(1:2, 23:26),]
my.table.or2 <- t(my.table)[c(1:2, 27:30),]

# Table A.9
my.table.nt1 %>% xtable
my.table.nt2 %>% xtable()
my.table.t1 %>% xtable()
my.table.t2 %>% xtable()

