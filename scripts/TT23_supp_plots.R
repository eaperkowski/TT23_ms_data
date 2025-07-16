##############################################################################
## Prep  
##############################################################################
## Load libraries
library(tidyverse)
library(ggpubr)
library(lme4)
library(car)
library(emmeans)
library(multcomp)
library(scales)

## Read compiled data file
df <- read.csv("../data/TT23_data.csv") %>%
  mutate(canopy_plot = ifelse(canopy == "open",
                              "open (April-May)",
                              "closed (June)"))
head(df)

# How many measurements per ID?
n_measurements <- df %>%
  group_by(id, spp, plot, subplot, gm.trt) %>%
  summarize(n_measurements = length(id)) %>%
  ungroup() %>%
  dplyr::select(id, spp, n_meas = n_measurements)
head(n_measurements)

# Join n measurements into photo traits
df2 <- df %>%
  left_join(n_measurements, by = "id") %>%
  filter(n_meas > 1) %>%
  dplyr::select(id, spp = spp.x, plot:canopy, canopy_plot,
                anet:inorg_n_ppm, n_meas) %>%
  mutate(gm.trt = factor(gm.trt, levels = c("weeded", "ambient")),
         canopy = factor(canopy, levels = c("open", "closed")),
         canopy_plot = factor(canopy_plot, 
                              levels = c("open (April-May)",
                                         "closed (June)")))

# Turn off digit rounding in emmean args
emm_options(opt.digits = FALSE)

## Read and subset soil dataset
df.soil <- df %>%
  distinct(composite, plot, canopy, .keep_all = TRUE) %>%
  mutate(gm.trt = factor(gm.trt, levels = c("weeded", "ambient")),
         canopy = factor(canopy, levels = c("open", "closed")),
         canopy_plot = factor(canopy_plot, 
                              levels = c("open (April-May)",
                                         "closed (June)")),
         np.ratio = inorg_n_ppm/phosphate_ppm)

## Read daily soil moisture dataset
df.sm <- read.csv("../data/TT23_tomst_probe_sm_daily.csv")

## Remove outliers
df2$anet[43] <- NA
df2$gsw[43] <- NA
df2$l[68] <- NA
df2$vcmax25[c(67, 96)] <- NA
df2$jmax25[c(67, 96)] <- NA
df2$jmax.vcmax[c(97)] <- NA

## Create models for soil data
nitrate <- lmer(
  nitrate_ppm ~ gm.trt * canopy_plot + (1 | plot), data = df.soil)

ammonium <- lmer(
  sqrt(ammonium_ppm) ~ gm.trt * canopy_plot + (1 | plot), data = df.soil)

phosphate <- lmer(
  phosphate_ppm ~ gm.trt * canopy_plot + (1 | plot), data = df.soil)

plant_availableN <- lmer(
  log(inorg_n_ppm) ~ gm.trt * canopy_plot + (1 | plot), data = df.soil)

n_to_p_ratio <- lmer(
  log(np.ratio) ~ gm.trt * canopy_plot + (1 | plot), data = df.soil)

sm_model <- lmer(daily_vwc ~ gm.trt * doy + (1 | plot),
                 data = df.sm)

## Create models for photosynthesis data
anet.tri <- lmer(log(anet) ~ gm.trt * canopy_plot + (1 | plot) + (1 | id), 
                 data = subset(df2, spp == "Tri"))

gsw.tri <- lmer(log(gsw) ~ gm.trt * canopy_plot + (1 | plot) + (1 | id), 
                data = subset(df2, spp == "Tri"))

l.tri <- lmer(l ~ gm.trt * canopy_plot + (1 | plot) + (1 | id), 
              data = subset(df2, spp == "Tri"))

vcmax25.tri <- lmer(log(vcmax25) ~ gm.trt * canopy_plot + (1 | plot) + (1 | id), 
                    data = subset(df2, spp == "Tri"))

jmax25.tri <- lmer(log(jmax25) ~ gm.trt * canopy_plot + (1 | plot) + (1 | id), 
                   data = subset(df2, spp == "Tri"))

jmax25_vcmax25.tri <- lmer(jmax.vcmax ~ gm.trt * canopy_plot + (1 | plot) + (1 | id), 
                           data = subset(df2, spp == "Tri"))

spad.tri <- lmer(SPAD ~ gm.trt * canopy_plot + (1 | plot) + (1 | id), 
                 data = subset(df2, spp == "Tri"))

anet.mai <- lmer(anet ~ gm.trt * canopy_plot + (1 | plot) + (1 | id), 
                 data = subset(df2, spp == "Mai"))

gsw.mai <- lmer(gsw ~ gm.trt * canopy_plot + (1 | plot) + (1 | id), 
                data = subset(df2, spp == "Mai"))

l.mai <- lmer(l ~ gm.trt * canopy_plot + (1 | plot) + (1 | id), 
              data = subset(df2, spp == "Mai"))

vcmax25.mai <- lmer(log(vcmax25) ~ gm.trt * canopy_plot + (1 | plot) + (1 | id), 
                    data = subset(df2, spp == "Mai"))

jmax25.mai <- lmer(log(jmax25) ~ gm.trt * canopy_plot + (1 | plot) + (1 | id), 
                   data = subset(df2, spp == "Mai"))

jmax25_vcmax25.mai <- lmer(jmax.vcmax ~ gm.trt * canopy_plot + (1 | plot) + (1 | id), 
                           data = subset(df2, spp == "Mai"))

spad.mai <- lmer(SPAD ~ gm.trt * canopy_plot + (1 | plot) + (1 | id), 
                 data = subset(df2, spp == "Mai"))

## Add code for facet labels
facet.labs <- c("Trillium spp.", "M. racemosum")
names(facet.labs) <- c("Tri", "Mai")

## Color palettes
gm.colors <- c("#00B2BE", "#F1B700")


##############################################################################
## Net photosynthesis - Trillium (gm.trt)
##############################################################################

anet_tri_gm_results <- cld(emmeans(anet.tri, pairwise~gm.trt, type = "response"), 
                        Letters = LETTERS, reversed = TRUE) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

anet_tri_gm_plot <- ggplot(data = subset(df2, spp == "Tri"),
                        aes(x = gm.trt, y = anet, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = anet_tri_gm_results, 
            aes(y = 18, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors) +
  scale_x_discrete(labels = label_wrap(10)) +
  scale_y_continuous(limits = c(0, 18), breaks = seq(0, 18, 6)) +
  labs(x = "Tree canopy status",
       y = expression(bold(italic("A")["net"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 18),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(color = "white"))
anet_tri_gm_plot


##############################################################################
## Net photosynthesis - Mai (gm.trt)
##############################################################################

anet_mai_gm_results <- cld(emmeans(anet.mai, pairwise~gm.trt, type = "response"), 
                           Letters = LETTERS, reversed = TRUE) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

anet_mai_gm_plot <- ggplot(data = subset(df2, spp == "Mai"),
                           aes(x = gm.trt, y = anet, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = anet_mai_gm_results, 
            aes(y = 18, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors) +
  scale_x_discrete(labels = label_wrap(10)) +
  scale_y_continuous(limits = c(0, 18), breaks = seq(0, 18, 6)) +
  labs(x = "Tree canopy status",
       y = expression(bold(italic("A")["net"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 18),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(color = "white"))
anet_mai_gm_plot

##############################################################################
## Net photosynthesis - Trillium (gm.trt)
##############################################################################

anet_tri_gm_results <- cld(emmeans(anet.tri, pairwise~gm.trt, type = "response"), 
                           Letters = LETTERS, reversed = TRUE) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

anet_tri_gm_plot <- ggplot(data = subset(df2, spp == "Tri"),
                           aes(x = gm.trt, y = anet, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = anet_tri_gm_results, 
            aes(y = 18, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors) +
  scale_x_discrete(labels = label_wrap(10)) +
  scale_y_continuous(limits = c(0, 18), breaks = seq(0, 18, 6)) +
  labs(x = "Tree canopy status",
       y = expression(bold(italic("A")["net"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 18),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(color = "white"))
anet_tri_gm_plot


##############################################################################
## Stomatal conductance - Mai (gm.trt)
##############################################################################

gsw_mai_gm_results <- cld(emmeans(gsw.mai, pairwise~gm.trt, type = "response"), 
                           Letters = LETTERS, reversed = TRUE) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

gsw_mai_gm_plot <- ggplot(data = subset(df2, spp == "Mai"),
                           aes(x = gm.trt, y = gsw, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = gsw_mai_gm_results, 
            aes(y = 0.3, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors) +
  scale_x_discrete(labels = label_wrap(10)) +
  scale_y_continuous(limits = c(0, 0.3), breaks = seq(0, 0.3, 0.1)) +
  labs(x = "Tree canopy status",
       y = expression(bold(italic("g")["sw"]*" (mol m"^"-2"*" s"^"-1"*")")),
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 18),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(color = "white"))
gsw_mai_gm_plot

