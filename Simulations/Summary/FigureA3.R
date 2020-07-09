########################################################################################
# Figure A.3 all three panels
########################################################################################
rm(list = ls())
library(reshape2)
library(xtable)
library(ggplot2)
library(dplyr)

#### For plots, thanks to https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

##### Load data
df <- read.csv(file = "Simulations/Combined/ResScen2.csv",header = T)
df <- df %>% filter(censoring==2) # keep only the high censoring rate

###### Summarize function estimation #####
df.baselines <- df %>% dplyr::select(n.sample, Method, Penal, n.knots, lambda, starts_with("base."))
df.baselines.means <- df.baselines %>% group_by(n.sample, Method, Penal, n.knots, lambda) %>% summarise_all(mean, na.rm = T) 
df.long <- melt(df.baselines.means, id.vars=c("n.sample", "Method", "Penal", "n.knots", "lambda")) %>% na.omit()
df.long$time <- rep(rep(65 + 2.5*(1:14), each = 39),3) 
df.long$type <- rep(c("Non-terminal baseline \n probability", "Terminal baseline \n probability", 
                      "Baseline OR"), each = nrow(df.long)/3)

#### True values ####
# Number of intervals
times <- seq(1, 14, 1)
alpha.nt <- LongitSemiComp:::logit(times*0.005  + 0.005*(times-2)^2 - (0.0002*(times + 1)^3) + 0.005)
alpha.t <- LongitSemiComp:::logit(times*0.0075  + 0.001*(times^2)  + 0.03)
alpha.or <- 0.9 + 0.175*times - 0.02*times^2 
alpha.or[times >= 13] <- 0


df.true.nt <- data.frame(time = 65 + 2.5*(1:14), value = expit(alpha.nt))
df.true.t <- data.frame(time = 65 + 2.5*(1:14), value = expit(alpha.t))
df.true.or <- data.frame(time = 65 + 2.5*(1:14), value = exp(alpha.or))
df.true <- cbind(rbind(df.true.nt, df.true.t, df.true.or), Method = rep("True baseline", 14*3), Penal = "Non-Penalized" ,
                 type = rep(c("Non-terminal baseline \n probability","Terminal baseline \n probability","Baseline OR"),each = 14))
df.true.for.plot <- rbind(df.true, df.true, df.true)
df.true.for.plot$n.sample <- rep(c(500, 1000, 5000), each = nrow(df.true))

###### Compare true baseline curves and estimated based using 10 knots #####
df.est <- df.long  
df.est <- df.est %>%  dplyr::select(time,Method, Penal, n.knots, n.sample, lambda, type, value)
df.true.for.plot$n.knots <- 999
df.true.for.plot$lambda <- 999

df.plot.baseline <- rbind(df.true.for.plot, df.est)
df.plot.baseline$N <- factor(paste0("N = ", df.plot.baseline$n.sample), levels = c("N = 500", "N = 1000", "N = 5000"))
df.plot.baseline <- df.plot.baseline %>% 
  mutate(Group = ifelse(Method %in% c("True baseline", "Non B-spline"), paste0(Method),
                        paste0(Method, ", ", Penal)))

df.plot.baseline$type.f = factor(df.plot.baseline$type, levels=c("Non-terminal baseline \n probability",
                                                                 "Terminal baseline \n probability",
                                                                 "Baseline OR"))

df.plot.baseline$n.sample.f <- factor(df.plot.baseline$n.sample, levels = c(500, 1000, 5000), 
                                      label = paste0("N = ", c(500, 1000, 5000)))
df.plot.baseline$n.knots.f <- factor(df.plot.baseline$n.knots, levels = sort(unique(df.plot.baseline$n.knots)), 
                                      label = paste0(sort(unique(df.plot.baseline$n.knots)), " Knots"))
df.plot.baseline$lambda.f <- factor(df.plot.baseline$lambda, levels = sort(unique(df.plot.baseline$lambda)), 
                                     label = paste0("lambda = ", sort(unique(df.plot.baseline$lambda))))

#get 6 ggplot colors
my.colors = gg_color_hue(6)


## Panel (a)
df.plot.baseline %>% filter(type.f=="Non-terminal baseline \n probability") %>% 
  filter(substr(Method,1,1) %in% c("B") & n.knots > 0) %>% 
  ggplot(aes(x = time, y = value, color = lambda.f, shape = lambda.f)) + 
  facet_grid(n.sample.f~n.knots.f, scales = "free") + 
  theme_bw() + geom_point(size = 2, shape = 19) + 
  ylab("Non-terminal baseline probability") + xlab("Age")  +
  theme(text = element_text(size = 24),  
        legend.title=element_blank(),
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 24),
        legend.key.size =  unit(1, "cm"),
        legend.key.height = unit(0.5,"cm")) +
  geom_point(data = filter(df.true,type=="Non-terminal baseline \n probability"), 
            aes(x = time,y = value, color = "True baseline"), size = 3,  shape = 19) + 
  scale_color_manual(values = c("True baseline" = "black", "lambda = 0" = my.colors[1],
                                "lambda = 0.1" = my.colors[2], "lambda = 1" = my.colors[3],
                                "lambda = 5" = my.colors[4], "lambda = 10"=  my.colors[5], 
                                "lambda = 25" = my.colors[6]),
                     breaks=c("True baseline","lambda = 0", "lambda = 0.1", "lambda = 1", "lambda = 5",
                                "lambda = 10", "lambda = 25"))
## Panel (b)
df.plot.baseline %>% filter(type.f=="Terminal baseline \n probability") %>% 
  filter(substr(Method,1,1) %in% c("B") & n.knots > 0) %>% 
  ggplot(aes(x = time, y = value, color = lambda.f, shape = lambda.f)) + 
  facet_grid(n.sample.f~n.knots.f, scales = "free") + 
  theme_bw() + geom_point(size = 2, shape = 19) + 
  ylab("Terminal baseline probability") + xlab("Age")  +
  theme(text = element_text(size = 24),  
        legend.title=element_blank(),
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 24),
        legend.key.size =  unit(1, "cm"),
        legend.key.height = unit(0.5,"cm")) +
  geom_point(data=filter(df.true,type=="Terminal baseline \n probability"), 
             aes(x = time,y = value, color = "True baseline"), size = 3 , shape = 19) + 
  scale_color_manual(values = c("True baseline" = "black", "lambda = 0" = my.colors[1],
                                "lambda = 0.1" = my.colors[2], "lambda = 1" = my.colors[3],
                                "lambda = 5" = my.colors[4], "lambda = 10"=  my.colors[5],
                                "lambda = 25" = my.colors[6]),
                     breaks=c("True baseline","lambda = 0", "lambda = 0.1", "lambda = 1", "lambda = 5",
                              "lambda = 10", "lambda = 25")) + 
  scale_shape_manual(breaks=c("True baseline", "lambda = 0", "lambda = 0.1", "lambda = 1", "lambda = 5",
                              "lambda = 10", "lambda = 25"), 
                     values = c(19, rep(19, 6)))
## Panel (c)
df.plot.baseline %>% filter(type.f=="Baseline OR" & value < 5) %>% 
  filter(substr(Method, 1, 1) %in% c("B") & n.knots > 0) %>% 
  ggplot(aes(x = time, y = value, color = lambda.f, shape = lambda.f)) + 
  facet_grid(n.sample.f~n.knots.f, scales = "free") + 
  theme_bw() + geom_point(size = 2, shape = 19) + 
  ylab("OR") + xlab("Age")  +
  theme(text = element_text(size = 24),  
        legend.title=element_blank(),
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 24),
        legend.key.size =  unit(1, "cm"),
        legend.key.height = unit(0.5,"cm")) +
  geom_point(data = filter(df.true,type=="Baseline OR"), 
             aes(x = time,y = value, color = "True baseline"), size = 3, shape = 19) + 
  scale_color_manual(values = c("True baseline" = "black", "lambda = 0" = my.colors[1],
                                "lambda = 0.1" = my.colors[2], "lambda = 1" = my.colors[3],
                                "lambda = 5" = my.colors[4], "lambda = 10"=  my.colors[5],
                                "lambda = 25" = my.colors[6]),
                     breaks=c("True baseline", "lambda = 0", "lambda = 0.1", "lambda = 1", "lambda = 5",
                              "lambda = 10", "lambda = 25")) 
