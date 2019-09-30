
############################################
#### code for Brambilla et. al analysis ####
############################################

## upload libraries
library(brms)
library(rstan)
library(ggplot2)
library(tidyverse)

## read data
df.tot<-read.csv("recruitment_data_resubmit.csv", head = TRUE)

str(df.tot)
df.tot$srug <- scale(df.tot$sr) ## mean center surface rugosity 
df.tot$rack <- as.factor(df.tot$rack)   # factorize rack

#### fitting models ####

# Rstan settings

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

### status, rugosity, depth ### 
set.seed(1234)
mod_stat_rug_depth <- brm(s.pres ~ alive.dead + srug + site + (1|rack),
                 data= df.tot, family = bernoulli(),
                 control = list(adapt_delta = 0.99),
                 iter = 20000, warmup = 1000, thin = 10)

summary(mod_stat_rug_depth)
plot(mod_stat_rug_depth) # chains check
#plot(marginal_effects(mod_stat_rug_depth)) # visualize effects

### status, depth ###
set.seed(1234)
mod_stat_depth <- brm(s.pres ~ alive.dead + site + (1|rack),
                          data= df.tot, family = bernoulli(),
                          control = list(adapt_delta = 0.99),
                          iter = 20000, warmup = 1000, thin = 10)

summary(mod_stat_depth)
plot(mod_stat_depth) # chains check
#plot(marginal_effects(mod_stat_depth)) # visualize effects

### rugosity, depth ###
set.seed(1234)
mod_rug_depth <- brm(s.pres ~ srug + site + (1|rack),
                      data= df.tot, family = bernoulli(),
                      control = list(adapt_delta = 0.99),
                      iter = 20000, warmup = 1000, thin = 10)

summary(mod_rug_depth)
plot(mod_rug_depth) # chains check
#plot(marginal_effects(mod_rug_depth)) # visualize effects

### rugosity, status ###
set.seed(1234)
mod_rug_stat <- brm(s.pres ~ srug + alive.dead + (1|rack),
                     data= df.tot, family = bernoulli(),
                     control = list(adapt_delta = 0.99),
                     iter = 20000, warmup = 1000, thin = 10)

summary(mod_rug_stat)
plot(mod_rug_stat) # chains check
#plot(marginal_effects(mod_rug_stat)) # visualize effect

### rugosity ###
set.seed(1234)
mod_rug <- brm(s.pres ~ srug + (1|rack),
                    data= df.tot, family = bernoulli(),
                    control = list(adapt_delta = 0.99),
                    iter = 20000, warmup = 1000, thin = 10)

summary(mod_rug)
plot(mod_rug) # chains check
#plot(marginal_effects(mod_rug)) # visualize effects

### status ###
set.seed(1234)
mod_stat <- brm(s.pres ~ alive.dead + (1|rack),
               data= df.tot, family = bernoulli(),
               control = list(adapt_delta = 0.99),
               iter = 20000, warmup = 1000, thin = 10)

summary(mod_stat)
plot(mod_stat) # chains check
#plot(marginal_effects(mod_stat), ask = FALSE) # visualize effects

### depth ###
set.seed(1234)
mod_depth <- brm(s.pres ~ site + (1|rack),
                data= df.tot, family = bernoulli(),
                control = list(adapt_delta = 0.99),
                iter = 20000, warmup = 1000, thin = 10)

summary(mod_depth)
plot(mod_depth) # chains check
#plot(marginal_effects(mod_depth), ask = FALSE) # visualize effects

### summary table ###
model.results <- data.frame(rbind(summary(mod_stat_rug_depth)$fixed,
summary(mod_stat_depth)$fixed,
summary(mod_rug_depth)$fixed,
summary(mod_rug_stat)$fixed,
summary(mod_rug)$fixed,
summary(mod_stat)$fixed,
summary(mod_depth)$fixed))

model.results <- data.frame(rbind(summary(mod_stat_rug_depth)$fixed,
                                  summary(mod_stat_depth)$fixed,
                                  summary(mod_rug_depth)$fixed,
                                  summary(mod_rug_stat)$fixed,
                                  summary(mod_rug)$fixed,
                                  summary(mod_stat)$fixed,
                                  summary(mod_depth)$fixed))

# compute and compare  WAIC and LOOic

WAIC(mod_stat_rug_depth, mod_stat_depth, mod_rug_depth, mod_rug_stat,mod_rug,
     mod_stat, mod_depth) # model 5 is the best fit

l1 <- loo(mod_stat_rug_depth)
l2 <- loo(mod_stat_depth)
l3 <- loo(mod_rug_depth)
l4 <- loo(mod_rug_stat)
l5 <- loo(mod_rug)
l6 <- loo(mod_stat)
l7 <- loo(mod_depth)

loo_compare(l1,l2,l3,l4,l5,l6,l7) # model 5 is the best fit

# create table
models <- c(rep("presence ~ alive.dead + srug + site + (1 | rackID)", 4),
            rep("presence ~ alive.dead + site + (1 | rackID)", 3),
            rep("presence ~ srug + site + (1 | rackID)", 3),
            rep("presence ~ srug + nubbin.status + (1 | rackID)", 3),
            rep("presence ~ srug + (1 | rackID)", 2),
            rep("presence ~ alive.dead + (1 | rackID)", 2),
            rep("presence ~ site + (1 | rackID)", 2))
            
model.results <- cbind(model = models, variable = rownames(model.results), model.results)

# write.csv(model.results,"results_summary.csv", row.names = FALSE, col.names = TRUE)

#### fig.3 ####

# add tile category based on combination of presence/absence of settlers 
# and status of coral nubbin
df.tot$rec.graph[df.tot$s.pres == 1 & df.tot$alive.dead == "alive"] <- "a"
df.tot$rec.graph[df.tot$s.pres == 1 & df.tot$alive.dead == "dead"] <- "b"
df.tot$rec.graph[df.tot$s.pres == 0 & df.tot$alive.dead == "alive"] <- "c"
df.tot$rec.graph[df.tot$s.pres == 0 & df.tot$alive.dead == "dead"] <-  "d"

# add sum of the tiles in each category and round sr at 3 digits
df.tot <- df.tot %>% group_by(rack) %>% group_by(alive.dead) %>% 
  mutate(spat.per = sum(s.pres), srr = round(sr, 3))

df.tot <- arrange(df.tot,sr,site)

# plot pie chart
ggplot(df.tot, aes(x=1, y= spat.per, fill=factor(rec.graph))) +
  geom_bar(stat="identity", position = position_fill()) +
  scale_fill_manual(values = viridis(9)[c(1,4,7,8)],
                    labels = c("  tiles with adult coral and settlers", 
                               "  tiles with dead coral and settlers", 
                               "  tiles with adult coral and no settlers", 
                               "  tiles with dead coral and no settlers"))+
  facet_wrap (site~srr, nrow = 2) +
  coord_polar(theta = "y", direction = -1)+
  ggtitle(" Tiles proportion in each rack\n")+
  theme_minimal()+
  theme(axis.title = element_blank(), axis.text = element_blank(), 
        title =  element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA,), # bg of the plot
        legend.title=element_blank(),
        legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", color = NA)) # get rid of legend panel


#### fig.4 ####

## fig 4a - marginal effect
m <- plot(marginal_effects(mod_rug, col = c("black","red")))
fig4a <- m[[1]] +
    geom_count(inherit.aes = FALSE,data = df.tot, aes(x = srug, y = s.pres), fill = "black", alpha = 0.2) +
    theme_minimal()+ ylim(c(0,1))+
    scale_color_grey() +
    scale_fill_grey() +
    labs(x = "\n standardized rugosity", y = "presence probability \n") +
    #geom_point(aes(x = srug, y = tot.ab, color = alive.dead), inherit.aes = FALSE, size = 2, alpha = 0.1,
    #           position = position_jitter(width = 0.01)) +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 13), 
          legend.margin = margin(3,3,3,3), legend.title = element_text(size = 13),
          legend.text = element_text(size = 12), legend.position ="top", 
          panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
          legend.box.background = element_rect(fill = "transparent", color = NA)) # get rid of legend panel

# fig 4b - posterior distribution of the slopes of mod_rug
posterior_slope <- as.array(mod_rug)
dim(posterior_slope)
dimnames(posterior_slope) 
p <-plot(posterior_slope, pars = "b_srug")

library(bayesplot)
d <- mcmc_areas(
  posterior_slope, 
  pars = "b_srug",
  prob = .95, # 80% intervals
  prob_outer = 1, # 99%
  point_est = "mean"
) 

d + 
  scale_color_grey()+
  scale_fill_grey()+
  scale_y_discrete(breaks="b_srug",
    labels="surface \n rugosity")+
  labs(x = "\n slope", y = " frequency\n")+
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 13), 
        legend.text = element_text(size = 12),legend.margin = margin(3,3,3,3), legend.title = element_text(size = 13),
        legend.text = element_text(size = 12), legend.position ="top", 
        panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", color = NA)) # get rid of legend panel

# fig 4c - marginal effect of status of the nubbin in mod_stat

s <- plot(marginal_effects(mod_stat))

fig4c <- s[[1]] +
    #geom_count(inherit.aes = FALSE,data = df.tot, aes(x = srug, y = s.pres), fill = "black", alpha = 0.2) +
    #theme_minimal()+ ylim(c(0,1))+
    #scale_color_grey() +
    #scale_fill_grey() +
    labs(x = "\n adult coral status", y = "presence probability \n") +
    #geom_point(aes(x = srug, y = tot.ab, color = alive.dead), inherit.aes = FALSE, size = 2, alpha = 0.1,
    #           position = position_jitter(width = 0.01)) +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 13), 
          legend.margin = margin(3,3,3,3), legend.title = element_text(size = 13),
          legend.text = element_text(size = 12), legend.position =c(0.8,0.6), 
          panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
          legend.box.background = element_rect(fill = "transparent", color = NA),
          axis.ticks = element_blank())) # get rid of legend panel)) # get rid of legend panel