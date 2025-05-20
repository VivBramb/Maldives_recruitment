
############################################
#### code for Brambilla et. al analysis ####
############################################

## upload libraries
library(brms)
library(rstan)
library(ggplot2)

## read data
df.tot<-read.csv("recruitment_data_resubmit.csv", head = TRUE)

str(df.tot)
df.tot$srug <- scale(df.tot$sr) ## mean center surface rugosity 
df.tot$rack <- as.factor(df.tot$rack)   # factorize rack

#### fitting models ####

# Rstan settings

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

### status, rugosity, depth #### 
set.seed(1234)
mod_stat_rug_depth <- brm(s.pres ~ alive.dead + srug + site + (1|rack),
                          data= df.tot, family = bernoulli(),
                          control = list(adapt_delta = 0.99),
                          iter = 20000, warmup = 1000, thin = 10)

summary(mod_stat_rug_depth)
plot(mod_stat_rug_depth) # chains check
#plot(marginal_effects(mod_stat_rug_depth)) # visualize effects

### status, depth ####
set.seed(1234)
mod_stat_depth <- brm(s.pres ~ alive.dead + site + (1|rack),
                      data= df.tot, family = bernoulli(),
                      control = list(adapt_delta = 0.99),
                      iter = 20000, warmup = 1000, thin = 10)

summary(mod_stat_depth)
plot(mod_stat_depth) # chains check
#plot(marginal_effects(mod_stat_depth)) # visualize effects

### rugosity, depth ####
set.seed(1234)
mod_rug_depth <- brm(s.pres ~ srug + site + (1|rack),
                     data= df.tot, family = bernoulli(),
                     control = list(adapt_delta = 0.99),
                     iter = 20000, warmup = 1000, thin = 10)

summary(mod_rug_depth)
plot(mod_rug_depth) # chains check
#plot(marginal_effects(mod_rug_depth)) # visualize effects

### rugosity, status ####
set.seed(1234)
mod_rug_stat <- brm(s.pres ~ srug + alive.dead + (1|rack),
                    data= df.tot, family = bernoulli(),
                    control = list(adapt_delta = 0.99),
                    iter = 20000, warmup = 1000, thin = 10)

summary(mod_rug_stat)
plot(mod_rug_stat) # chains check
#plot(marginal_effects(mod_rug_stat)) # visualize effect

### rugosity ####
set.seed(1234)
mod_rug <- brm(s.pres ~ srug + (1|rack),
               data= df.tot, family = bernoulli(),
               control = list(adapt_delta = 0.99),
               iter = 20000, warmup = 1000, thin = 10)

summary(mod_rug)
plot(mod_rug) # chains check
#plot(marginal_effects(mod_rug)) # visualize effects

### status ####
set.seed(1234)
mod_stat <- brm(s.pres ~ alive.dead + (1|rack),
                data= df.tot, family = bernoulli(),
                control = list(adapt_delta = 0.99),
                iter = 20000, warmup = 1000, thin = 10)

summary(mod_stat)
plot(mod_stat) # chains check
#plot(marginal_effects(mod_stat), ask = FALSE) # visualize effects

### depth ####
set.seed(1234)
mod_depth <- brm(s.pres ~ site + (1|rack),
                 data= df.tot, family = bernoulli(),
                 control = list(adapt_delta = 0.99),
                 iter = 20000, warmup = 1000, thin = 10)

summary(mod_depth)
plot(mod_depth) # chains check
#plot(marginal_effects(mod_depth), ask = FALSE) # visualize effects

### summary table SM####
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

#### fig.2 ####
ggplot2::theme_set(theme_bw())
theme_set(theme_tidybayes() + panel_border() + background_grid())
theme_set(theme_tidybayes())
m <- plot(marginal_effects(mod_rug))
(mod1 <- m[[1]] + 
    theme_minimal()+ ylim(c(0,1))+
    scale_color_grey() +
    scale_fill_grey() +
    geom_count(inherit.aes = FALSE,data = df.tot, aes(x = srug, y = s.pres), fill = "black", alpha = 0.2) +
    labs(x = "\n standardized rugosity", y = "presence probability \n") +
    # geom_point(data = df.tot,aes(x = alive.dead, y = tot.ab, color = alive.dead), inherit.aes = FALSE, size = 2, alpha = 0.1,
    #            position = position_jitter(width = 0.01)) +
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 13), 
          legend.margin = margin(3,3,3,3), legend.title = element_text(size = 13),
          legend.text = element_text(size = 12), legend.position = c(0.8,0.65), 
          panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
          legend.box.background = element_rect(fill = "transparent", color = NA))) # get rid of legend panel
#
s <- plot(marginal_effects(mod_stat))

(mod2<- s[[1]] +
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



# effect size
posterior_slope_r <- as.array(mod_rug)
dim(posterior_slope_r)
dimnames(posterior_slope_r) 
p <-plot(posterior_slope_r, pars = "b_srug")

p
library(bayesplot)

color_scheme_set("darkgray")
r <- mcmc_areas(
  posterior_slope_r, 
  pars = c("b_srug"),
  prob = .95, # 80% intervals
  prob_outer = 1, # 99%
  point_est = "mean"
) 
r1 <- r + 
  #scale_y_discrete(breaks="b_srug",
  #                 labels="surface \n rugosity")+
  labs(x = "\n slope (SR)", y = " frequency distribution")+
  #xlim(c(-0.9,1.9))+
  vline_0(linetype = "dotted")+
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 13), 
        #legend.margin = margin(3,3,3,3), legend.title = element_text(size = 13),
        #legend.text = element_text(size = 12), legend.position ="top", 
        panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", color = NA),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) # get rid of legend panel



r1
posterior_slope_s <- as.array(mod_stat)
dim(posterior_slope_s)
dimnames(posterior_slope_s) 
p <-plot(posterior_slope_s, pars = "b_srug")


d <- mcmc_areas(
  posterior_slope_s, 
  pars = "b_alive.deaddead",
  prob = .95, # 80% intervals
  prob_outer = 1, # 99%
  point_est = "mean"
) 

d1 <-d + 
  scale_color_grey()+
  vline_0(linetype = "dotted")+
  xlim(c(-2,2))+
  labs(x = "\n slope", y = " frequency\n")+
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 13), 
        #legend.margin = margin(3,3,3,3), legend.title = element_text(size = 13),
        #legend.text = element_text(size = 12), legend.position ="top", 
        panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", color = NA)) # get rid of legend panel
d1


D <- ggplotGrob(d1)
plot(mod1) + annotation_custom(
     grob = D,
     xmin = -1,
     xmax = 0,
     ymin = 0.5,
     ymax = 0.875
     )

ggarrange(plot(mod1) + annotation_custom(
  grob = D,
  xmin = -1,
  xmax = 0,
  ymin = 0.5,
  ymax = 0.875
))


library(scatterpie)

df <- as.data.frame(df.tot %>% group_by(srr, site, rack, rec.graph) %>% count()) 
df <- df%>% pivot_wider(names_from = rec.graph, values_from = n)
df$b[is.na(df$b)]<- 0
df$depth <- 0
df$depth[df$site == "deep"]<- -0.08

ggplot() +
  geom_scatterpie(data=df, aes(x=srr, y= depth, group = rack), cols = names(df[4:7])) + coord_equal()+
  scale_fill_manual(values = viridis(9)[c(1,4,7,8)],
                    labels = c("  tiles with adult coral and settlers",
                               "  tiles with dead coral and settlers",
                               "  tiles with adult coral and no settlers",
                               "  tiles with dead coral and no settlers"))+
  # facet_wrap (site~srr, nrow = 2) +
  # coord_polar(theta = "y", direction = -1)+
  # ggtitle(" Tiles proportion in each rack\n")+
  theme_bw()+ 
  scale_x_continuous("surface rugosity",trans='log') +
  scale_y_continuous( "depth", breaks = c(0,-0.1),
                      label = c("6m", "18m"))+
  theme(legend.position = "bottom",
        legend.direction="vertical",
        text = element_text(size = 12),
        legend.text = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA,), # bg of the plot
        legend.title=element_blank(),
        legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", color = NA)) # get rid of legend panel
