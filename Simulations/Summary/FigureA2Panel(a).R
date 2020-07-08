########################################################################################
# Figure A.2 panel (b)
########################################################################################
rm(list = ls())

read.csv2(file = "Simulations/Summary/ResScen1.csv")
df <- df %>% filter(censoring==2) # keep only the high censoring rate

###### Summarize function estimation #####
df.baselines <- df %>% dplyr::select(n.sample, Method, Penal, n.knots, lambda, starts_with("base."))
df.baselines.means <- df.baselines %>% group_by(n.sample, Method, Penal, n.knots, lambda) %>% summarise_all(mean, na.rm = T) 
df.long <- melt(df.baselines.means, id.vars=c("n.sample", "Method", "Penal", "n.knots", "lambda")) %>% na.omit()
df.long$time <- rep(rep(65 + 2.5*(1:14), each = 39),3) 
df.long$type <- rep(c("Non-terminal baseline \n probability", "Terminal baseline \n probability", 
                      "Baseline OR"), each = nrow(df.long)/3)
df.true.nt <- data.frame(time = 65 + 2.5*(1:14), value = expit(alpha.nt))
df.true.t <- data.frame(time = 65 + 2.5*(1:14), value = expit(alpha.t))
df.true.or <- data.frame(time = 65 + 2.5*(1:14), value = exp(alpha.or))
df.true <- cbind(rbind(df.true.nt, df.true.t, df.true.or), Method = rep("True baseline", 14*3), Penal = "Non-Penalized" ,
                 type = rep(c("Non-terminal baseline \n probability","Terminal baseline \n probability","Baseline OR"),each = 14))

df.true.for.plot <- rbind(df.true, df.true, df.true)
df.true.for.plot$n.sample <- rep(c(500, 1000, 5000), each = nrow(df.true))
###### Compare true baseline curves and estimated based using 10 knots #####
df.est <- df.long   %>% filter(lambda %in% c(0,5))
df.est <- df.est %>%  dplyr::select(time, Method, Penal, n.sample, type, value)
df.plot.baseline <- rbind(df.true.for.plot, df.est)
df.plot.baseline$N <- factor(paste0("N = ", df.plot.baseline$n.sample), levels = c("N = 500", "N = 1000", "N = 5000"))
df.plot.baseline <- df.plot.baseline %>% 
  mutate(Group = ifelse(Method %in% c("True baseline", "Non B-spline"), paste0(Method),
                        paste0(Method, ", ", Penal)))

df.plot.baseline$type.f = factor(df.plot.baseline$type, levels=c("Non-terminal baseline \n probability",
                                                                 "Terminal baseline \n probability",
                                                                 "Baseline OR"))

#################### Figure A.2, Panel (b)####################
df.plot.baseline %>% filter(value < 5) %>% 
  ggplot(aes(x = time, y = value, color = Group, shape = Group, size = Group)) + 
  facet_grid(type.f ~ N, scales = "free") + 
  theme_bw() + geom_point()  +  #scale_colour_grey(start = 0, end = .75)  + 
  ylab("Probability/OR") + xlab("Age")  +
  theme(text = element_text(size = 20),  
        legend.title=element_blank(),
        legend.text = element_text(size = 16),
        strip.text.x = element_text(size = 18),
        legend.key.size =  unit(1.3, "cm"),
        legend.key.height = unit(1,"cm")) + 
  #  legend.position="bottom") +
  guides(color = guide_legend(override.aes = list(size = 1.5))) +
  scale_colour_manual(name = "",
                      labels = c("B-splines, 10 knots, \n Non-Penalized", "B-splines, 10 knots, \n Penalized (lambda=5)",
                                 "B-splines, 5 knots, \n Non-Penalized", "B-splines, 5 knots, \n Penalized (lambda=5)",
                                 "Non B-spline", "True baseline"),
                      values = c("grey60","grey60","grey40", "grey40", "grey75", "black")) +   
  scale_shape_manual(name = "",
                     labels = c("B-splines, 10 knots, \n Non-Penalized", "B-splines, 10 knots, \n Penalized (lambda=5)",
                                "B-splines, 5 knots, \n Non-Penalized", "B-splines, 5 knots, \n Penalized (lambda=5)",
                                "Non B-spline", "True baseline"),
                     values = c(0, 3, 0, 3, 16, 16)) + 
  scale_size_manual(name = "",
                    labels = c("B-splines, 10 knots, \n Non-Penalized", "B-splines, 10 knots, \n Penalized (lambda=5)",
                               "B-splines, 5 knots, \n Non-Penalized", "B-splines, 5 knots, \n Penalized (lambda=5)",
                               "Non B-spline", "True baseline"),
                    values = c(1.3, 1.3, 1.3, 1.3, 1.5, 2.5)) + 
  guides(colour = guide_legend(override.aes = list(size = 4)))



