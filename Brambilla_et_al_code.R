#### code for Brambilla et. al analysis ####

## read data
df.tot <- read.csv("recruitment_data.csv")
df.tot$srug <- scale(df.tot$sr)  # scale surface rugosity

## upload libraries
library(brms)
library(rstan)
library(ggplot2)

## Rstan settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#### model fitting ####

## model 1: juvenile presence
set.seed(1234)
modbrmsjp_sr <- brm(presence ~ srug + (1|rack),
                    data= df.tot, family = bernoulli(),
                    control = list(adapt_delta = 0.99),
                    iter = 20000, warmup = 1000, thin = 10)
summary(modbrmsjp_sr)
summary(modbrmsjp_sr)[17]
plot(marginal_effects(modbrmsjp_sr))

### model 2: settlers presence
set.seed(1234)

modbrmssp_sr <- brm(rec.presence ~ srug + (1|rack),
                    data= df.tot, family = bernoulli(),
                    control = list(adapt_delta = 0.99),
                    iter = 20000, warmup = 1000, thin = 10)
summary(modbrmssp_sr)
summary(modbrmssp_sr)[17]
plot(marginal_effects(modbrmssp_sr))

#### Fig 2 ####
## model 1
m1 <- plot(marginal_effects(modbrmsjp_sr, spaghetti = TRUE, col = c("black","red")))
(mod1 <- m1[[1]] +
    geom_count(inherit.aes = FALSE,data = df.tot, aes(x = srug, y = presence), alpha = 0.2) +
    theme_minimal()+ ylim(c(0,1))+
    scale_color_manual(values="black") +
    scale_fill_manual(values="black") +
    labs(x = "\n standardized rugosity", y = "presence probability \n") +
    #geom_point(aes(x = srug, y = tot.ab, color = alive.dead), inherit.aes = FALSE, size = 2, alpha = 0.1,
    #           position = position_jitter(width = 0.01)) +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 11), 
          legend.margin = margin(3,3,3,3), legend.title = element_text(size = 11),
          legend.text = element_text(size = 10), legend.position ="top", 
          panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
          legend.box.background = element_rect(fill = "transparent", color = NA))) # get rid of legend panel

## model 2
m2 <- plot(marginal_effects(modbrmssp_sr, spaghetti = TRUE, 
                             spaghetti_args	= list(inherit.aes = FALSE, aes(colour = "black", fill = "black"))))

(mod2 <- m2$`srug` + ylim(c(0,1))+
    geom_count(inherit.aes = FALSE, data = df.tot, aes(x = srug, y = rec.presence), alpha = 0.2)+
    scale_color_manual(values=c("#3CBB75FF", "#404788FF"))+
    scale_fill_manual(values=c("#3CBB75FF", "#404788FF"))+
    theme_minimal() +
    labs(x = "\n standardized rugosity", y = "presence probability \n") +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 11), 
          legend.margin = margin(3,3,3,3), legend.title = element_text(size = 11),
          legend.text = element_text(size = 10),legend.position = "top",
          panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
          legend.box.background = element_rect(fill = "transparent", color = NA))) # get rid of legend panel

