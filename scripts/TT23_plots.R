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
library(gghalves)

## Read compiled data file
df <- read.csv("../data/TT23_data.csv") %>%
  mutate(gm.trt = factor(gm.trt, levels = c( "weeded", "ambient")),
         canopy = factor(canopy, levels = c("open", "closed")),
         spp = factor(spp, levels = c("Tri", "Mai", "Ari"))) %>%
  unite(col = "gm.canopy", gm.trt, canopy, sep = "_", remove = FALSE) %>%
  mutate(gm.canopy = factor(gm.canopy, levels = c("weeded_open",
                                                  "ambient_open",
                                                  "weeded_closed",
                                                  "ambient_closed")),
         trt = str_c(canopy, "_", gm.trt))
head(df)

# helper fxn to change "NaN" to "NA"
NaN_to_NA <- function(x) ifelse(is.nan(x), NA, x)

## Read and subset soil dataset
df.soil <- df %>%
  group_by(plot, composite, gm.trt, canopy) %>%
  summarize_at(.vars = vars(phosphate_ppm:inorg_n_ppm),
               .funs = mean, na.rm = TRUE) %>%
  mutate(gm.trt = factor(gm.trt, levels = c("weeded", "ambient")),
         canopy = factor(canopy, levels = c("open", "closed")),
         trt = str_c(canopy, "_", gm.trt)) %>%
  mutate(across(nitrate_ppm:inorg_n_ppm, .fns = NaN_to_NA),
         np.ratio = inorg_n_ppm/phosphate_ppm)

## Read daily soil moisture dataset
df.sm <- read.csv("../data/TT23_tomst_probe_sm_daily.csv") %>%
  mutate(gm.trt = factor(gm.trt, levels = c( "weeded", "ambient")))

## Remove outliers
df.soil$np.ratio[54] <- NA
df$anet[91] <- NA
df$gsw[c(53, 91)] <- NA
df$SPAD[111] <- NA
df$vcmax25[180] <- NA
df$jmax25[c(142, 180)] <- NA
df$jmax.vcmax[c(113, 181, 218, 219)] <- NA

## Create models for soil data
nitrate <- lmer(
  nitrate_ppm ~ gm.trt * canopy + (1 | plot), data = df.soil)

ammonium <- lmer(
  log(ammonium_ppm) ~ gm.trt * canopy + (1 | plot), data = df.soil)

phosphate <- lmer(
  phosphate_ppm ~ gm.trt * canopy + (1 | plot), data = df.soil)

plant_availableN <- lmer(
  log(inorg_n_ppm) ~ gm.trt * canopy + (1 | plot), data = df.soil)

n_to_p_ratio <- lmer(
  log(np.ratio) ~ gm.trt * canopy + (1 | plot), data = df.soil)

sm_model <- lmer(daily_vwc ~ gm.trt * doy + (1 | plot),
                 data = df.sm)

## Create models for photosynthesis data
anet.tri <- lmer(
  log(anet) ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))

anet.mai <- lmer(
  anet ~ gm.trt * canopy  + (1 | plot), data = subset(df, spp == "Mai"))

gsw.tri <- lmer(
  log(gsw) ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))

gsw.mai <- lmer(
  gsw ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))

l.tri <- lmer(
  l ~ gm.trt * canopy  + (1 | plot), data = subset(df, spp == "Tri" & l > 0))

l.mai <- lmer(
  log(l) ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai" & l > 0))

spad.tri <- lmer(
  SPAD ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))

spad.mai <- lmer(
  log(SPAD) ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))

vcmax.tri <- lmer(
  log(vcmax25) ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))

vcmax.mai <- lmer(
  vcmax25 ~ gm.trt * canopy  + (1 | plot), data = subset(df, spp == "Mai"))

jmax.tri <- lmer(
  log(jmax25) ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))

jmax.mai <- lmer(
  log(jmax25) ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))

jmax.vcmax.tri <- lmer(
  jmax.vcmax ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))

jmax.vcmax.mai <- lmer(
  jmax.vcmax ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))

## Add code for facet labels
facet.labs <- c("Trillium spp.", "M. racemosum")
names(facet.labs) <- c("Tri", "Mai")

## Color palettes
gm.colors <- c("#00B2BE", "#F1B700")

##############################################################################
## Soil nitrate availability 
##############################################################################
# Prep
nitrate_results <- cld(emmeans(nitrate, pairwise~canopy*gm.trt),
                      reversed = TRUE, Letters = LETTERS) %>%
  mutate(.group = trimws(.group, "both"))

# Plot
nitrate_plot <- ggplot(data = df.soil,
                       aes(x = canopy, y = nitrate_ppm, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = nitrate_results, 
            aes(y = 40, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open", "closed")) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 10)) +
  labs(x = "Tree canopy status",
       y = expression(bold("Soil NO"["3"]*"-N (ppm)")),
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.grid.minor.y = element_blank())
nitrate_plot

##############################################################################
## Soil ammonium availability 
##############################################################################
# Prep
ammonium_results <- cld(emmeans(ammonium, pairwise~canopy*gm.trt),
                       reversed = TRUE, Letters = LETTERS) %>%
  mutate(.group = trimws(.group, "both"))

# Plot
ammonium_plot <- ggplot(data = df.soil,
                       aes(x = canopy, y = ammonium_ppm, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = ammonium_results, 
            aes(y = 3, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open", "closed")) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  labs(x = "Tree canopy status",
       y = expression(bold("Soil NH"["4"]*"-N (ppm)")),
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.grid.minor.y = element_blank())
ammonium_plot

##############################################################################
## Soil N availability 
##############################################################################
# Prep
inorgN_results <- cld(emmeans(plant_availableN, pairwise~canopy*gm.trt),
                      reversed = TRUE, Letters = LETTERS) %>%
  mutate(.group = trimws(.group, "both"))

# Plot
nitrogen_plot <- ggplot(data = df.soil,
                        aes(x = canopy, y = inorg_n_ppm, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = inorgN_results, 
            aes(y = 40, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open", "closed")) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 10)) +
  labs(x = "Tree canopy status",
       y = "Soil inorg. N (ppm)",
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(face = "bold"),
        panel.grid.minor.y = element_blank())
nitrogen_plot

##############################################################################
## Phosphate figure  
##############################################################################
# Prep
phosphate_results <- cld(emmeans(phosphate, pairwise~canopy*gm.trt), 
                         Letters = LETTERS, reversed = TRUE) %>%
  mutate(.group = trimws(.group, "both"))

# Plot
phosphate_plot <- ggplot(data = df.soil,
                         aes(x = canopy, y = phosphate_ppm, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = phosphate_results, 
            aes(y = 2, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open", "closed")) +
  scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.5)) +
  labs(x = "Tree canopy status",
       y = "Soil phosphate (ppm)",
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(face = "bold"),
        panel.grid.minor.y = element_blank())
phosphate_plot

##############################################################################
## Soil N:P figure  
##############################################################################
Anova(n_to_p_ratio)

# Prep
soil_np_results <- cld(emmeans(n_to_p_ratio, pairwise~canopy*gm.trt), 
                       Letters = LETTERS, reversed = TRUE) %>%
  mutate(.group = trimws(.group, which = "both"))

# Plot
soil_np_plot <- ggplot(data = df.soil,
                       aes(x = canopy, y = np.ratio, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = soil_np_results, 
            aes(y = 40, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open", "closed")) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 10)) +
  labs(x = "Tree canopy status",
       y = "Soil N:P ratio (unitless)",
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(face = "bold"),
        panel.grid.minor.y = element_blank())
soil_np_plot

###########################################################
## Soil moisture
###########################################################
Anova(sm_model)

# Prep
sm_means <- df.sm %>%
  group_by(day, doy, gm.trt) %>%
  summarize(vwc_mean = mean(daily_vwc))

sm_results <- data.frame(
  emmeans(sm_model, ~gm.trt, "doy",
          at = list(doy = seq(116, 181, 1))))

# Plot
sm_plot <- ggplot(data = sm_means, aes(x = doy, y = vwc_mean)) +
  geom_line(aes(color = gm.trt)) +
  geom_point(aes(fill = gm.trt), shape = 21, size = 3) +
  geom_smooth(data = sm_results, 
              aes(x = doy, y = emmean, color = gm.trt),
              se = FALSE, linewidth = 1.5) +
  geom_ribbon(data = sm_results, 
              aes(x = doy, y = emmean, ymin = emmean - SE,
                  ymax = emmean + SE, fill = gm.trt),
              alpha = 0.1) +
  scale_color_manual(values = gm.colors) +
  scale_fill_manual(values = gm.colors) +
  scale_y_continuous(limits = c(0.1, 0.42), 
                     breaks = seq(0.1, 0.4, 0.1)) +
  scale_x_continuous(breaks = seq(116, 181, 13),
                     labels = c("April 26", "May 9", "May 22", 
                                "June 4", "June 17", "June 30")) +
  labs(x = "Date", y = "Volumetric soil moisture (%)",
       color = expression(bolditalic("Alliaria")*bold(" treatment")),
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        panel.grid.minor.y = element_blank())
sm_plot

##############################################################################
## Net photosynthesis - Trillium
##############################################################################
Anova(anet.tri)

anet_tri_results <- cld(emmeans(anet.tri, ~canopy*gm.trt, type = "response"), 
                Letters = LETTERS, reversed = TRUE) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

anet_tri_plot <- ggplot(data = subset(df, spp == "Tri"),
                        aes(x = canopy, y = anet, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = anet_tri_results, 
            aes(y = 18, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open", "closed")) +
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
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank())
anet_tri_plot

##############################################################################
## Net photosynthesis - Maianthemum
##############################################################################
Anova(anet.mai)

anet_mai_results <- cld(emmeans(anet.mai, ~canopy*gm.trt, type = "response"), 
                        Letters = LETTERS, reversed = TRUE) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

anet_mai_plot <- ggplot(data = subset(df, spp == "Mai"),
                        aes(x = canopy, y = anet, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = anet_mai_results, 
            aes(y = 18, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open", "closed")) +
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
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank())
anet_mai_plot

##############################################################################
## Stomatal conductance - Tri
##############################################################################
Anova(gsw.tri)

gsw_tri_results <- cld(emmeans(gsw.tri, ~canopy*gm.trt, type = "response"), 
                        Letters = LETTERS, reversed = TRUE) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

gsw_tri_plot <- ggplot(data = subset(df, spp == "Tri"),
                  aes(x = canopy, y = gsw, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = gsw_tri_results, 
            aes(y = 0.3, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open", "closed")) +
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
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank())
gsw_tri_plot

##############################################################################
## Stomatal conductance - Maianthemum
##############################################################################
Anova(gsw.mai)

gsw_mai_results <- cld(emmeans(gsw.mai, ~canopy*gm.trt, type = "response"), 
                        Letters = LETTERS, reversed = TRUE, alpha = 0.06) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

gsw_mai_plot <- ggplot(data = subset(df, spp == "Mai"),
                       aes(x = canopy, y = gsw, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = gsw_mai_results, 
            aes(y = 0.3, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open", "closed")) +
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
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank())
gsw_mai_plot

##############################################################################
## Stomatal limitation - Tri
##############################################################################
Anova(l.tri)

l_tri_results <- cld(emmeans(l.tri, ~canopy*gm.trt, type = "response"), 
                       Letters = LETTERS, reversed = TRUE) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

l_tri_plot <- ggplot(data = subset(df, spp == "Tri"),
                           aes(x = canopy, y = l, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = l_tri_results, 
            aes(y = 1, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open", "closed")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Tree canopy status",
       y = "Stom. limitation (unitless)",
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank())
l_tri_plot

##############################################################################
## Stomatal limitation - Mai
##############################################################################
Anova(l.mai)

l_mai_results <- cld(emmeans(l.mai, ~canopy*gm.trt, type = "response"), 
                       Letters = LETTERS) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

l_mai_plot <- ggplot(data = subset(df, spp == "Mai"),
                           aes(x = canopy, y = l, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = l_mai_results, 
            aes(y = 1, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open", "canopy")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Tree canopy status",
       y = "Stom. limitation (unitless)",
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank()) +
  guides(fill = "none")
l_mai_plot

##############################################################################
## Vcmax - Tri
##############################################################################
Anova(vcmax.tri)

vcmax_tri_results <- cld(emmeans(vcmax.tri, pairwise~canopy*gm.trt, type = "response"), 
                         Letters = LETTERS, reversed = TRUE) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

vcmax_tri_plot <- ggplot(data = subset(df, spp == "Tri"),
                         aes(x = canopy, y = vcmax25, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = vcmax_tri_results, 
            aes(x = canopy, y = 200, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open", "closed")) +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, 50)) +
  labs(x = "Tree canopy status",
       y = expression(bold(italic("V")["cmax25"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank())
vcmax_tri_plot

##############################################################################
## Vcmax - Mai
##############################################################################
Anova(vcmax.mai)

vcmax_mai_results <- cld(emmeans(vcmax.mai, ~gm.trt*canopy, type = "response"), 
                         Letters = LETTERS, reversed = TRUE) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

vcmax_mai_plot <- ggplot(data = subset(df, spp == "Mai"),
                         aes(x = canopy, y = vcmax25, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = vcmax_mai_results, 
            aes(x = canopy, y = 100, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open", "closed")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(x = "Tree canopy status",
       y = expression(bold(italic("V")["cmax25"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank())
vcmax_mai_plot

##############################################################################
## Jmax - Tri
##############################################################################
Anova(jmax.tri)

jmax_tri_results <- cld(emmeans(jmax.tri, ~gm.trt*canopy, type = "response"), 
                         Letters = LETTERS, reversed = TRUE) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

jmax_tri_plot <- ggplot(data = subset(df, spp == "Tri"),
                         aes(x = canopy, y = jmax25, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = jmax_tri_results, 
            aes(x = canopy, y = 300, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open", "closed")) +
  scale_y_continuous(limits = c(0, 300), breaks = seq(0, 300, 100)) +
  labs(x = "Tree canopy status",
       y = expression(bold(italic("J")["max25"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank())
jmax_tri_plot

##############################################################################
## Jmax - Mai
##############################################################################
Anova(jmax.mai)

jmax_mai_results <- cld(emmeans(jmax.mai, ~gm.trt*canopy, type = "response"), 
                         Letters = LETTERS, reversed = TRUE) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

jmax_mai_plot <- ggplot(data = subset(df, spp == "Mai"),
                         aes(x = canopy, y = jmax25, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = jmax_mai_results, 
            aes(x = canopy, y = 150, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open", "closed")) +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 50)) +
  labs(x = "Tree canopy status",
       y = expression(bold(italic("J")["max25"]*" ("*mu*"mol m"^"-2"*" s"^"-1"*")")),
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank())
jmax_mai_plot

##############################################################################
## Jmax: Vcmax - Tri
##############################################################################
Anova(jmax.vcmax.tri)

jvmax_tri_results <- cld(emmeans(jmax.vcmax.tri, ~gm.trt*canopy, type = "response"), 
                        Letters = LETTERS, reversed = TRUE) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

jvmax_tri_plot <- ggplot(data = subset(df, spp == "Tri"),
                        aes(x = canopy, y = jmax.vcmax, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = jvmax_tri_results, 
            aes(x = canopy, y = 2.3, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open", "closed")) +
  scale_y_continuous(limits = c(1.4, 2.3), breaks = seq(1.4, 2.2, 0.2)) +
  labs(x = "Tree canopy status",
       y = expression(bold(italic("J")["max25"]*":"*italic("V")["cmax25"]*" (unitless)")),
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank())
jvmax_tri_plot

##############################################################################
## Jmax:Vcmax - Mai
##############################################################################
Anova(jmax.vcmax.mai)

jvmax_mai_results <- cld(emmeans(jmax.vcmax.mai, ~gm.trt*canopy, type = "response"), 
                        Letters = LETTERS, reversed = TRUE) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

jvmax_mai_plot <- ggplot(data = subset(df, spp == "Mai"),
                        aes(x = canopy, y = jmax.vcmax, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = jvmax_mai_results, 
            aes(x = canopy, y = 2.3, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open", "closed")) +
  scale_y_continuous(limits = c(1.4, 2.3), breaks = seq(1.4, 2.2, 0.2)) +
  labs(x = "Tree canopy status",
       y = expression(bold(italic("J")["max25"]*":"*italic("V")["cmax25"]*" (unitless)")),
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank())
jvmax_mai_plot

##############################################################################
## SPAD - Tri
##############################################################################
Anova(spad.tri)

## Canopy plot
spad_tri_results <- cld(emmeans(spad.tri, ~canopy*gm.trt, type = "response"),
                               Letters = LETTERS) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

spad_tri_plot <- ggplot(data = subset(df, spp == "Tri"),
                       aes(x = canopy, y = SPAD, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = spad_tri_results, 
            aes(y = 60, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open", "closed")) +
  scale_y_continuous(limits = c(15, 60), breaks = seq(15, 60, 15)) +
  labs(x = "Tree canopy status",
       y = "SPAD (unitless)",
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank())
spad_tri_plot

##############################################################################
## SPAD - Mai
##############################################################################
Anova(spad.mai)

## Canopy plot
spad_mai_results <- cld(emmeans(spad.mai, ~canopy*gm.trt, type = "response"),
                               Letters = LETTERS) %>% 
  data.frame() %>% mutate(.group = trimws(.group, "both"))

spad_mai_plot <- ggplot(data = subset(df, spp == "Mai"),
                               aes(x = canopy, y = SPAD, fill = gm.trt)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(position = position_dodge(0.75),
               width = 0.5, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, 
                                             jitter.width = 0.1),
             alpha = 0.5, size = 2.5, shape = 21) +
  geom_text(data = spad_mai_results, 
            aes(y = 60, group = gm.trt, label = .group),
            position = position_dodge(width = 0.75), 
            fontface = "bold", size = 6) +
  scale_fill_manual(values = gm.colors,
                    labels = c("weeded", "ambient")) +
  scale_x_discrete(labels = c("open", "closed")) +
  scale_y_continuous(limits = c(15, 60), breaks = seq(15, 60, 15)) +
  labs(x = "Tree canopy status",
       y = "SPAD (unitless)",
       fill = expression(bolditalic("Alliaria")*bold(" treatment"))) +
  facet_grid(~spp, labeller = labeller(spp = facet.labs)) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", size = 18),
        panel.grid.minor.y = element_blank())
spad_mai_plot

##############################################################################
## Figure 1: Soil nutrients 
##############################################################################
png("../drafts/figs/TT23_fig1_soilNutrients.png", width = 12, height = 4.5,
    units = "in", res = 600)
ggarrange(nitrogen_plot, phosphate_plot, soil_np_plot, ncol = 3, nrow = 1, 
          hjust = 0, common.legend = TRUE, legend = "bottom",
          align = "hv", labels = c("(a)", "(b)", "(c)", "(d)"), 
          font.label = list(size = 18))
dev.off()

##############################################################################
## Figure 1: Soil moisture
##############################################################################
png("../drafts/figs/TT23_fig2_soilMoisture.png",
    width = 8, height = 4.5, units = "in", res = 600)
sm_plot
dev.off()

##############################################################################
## Figure 2: Gas exchange
##############################################################################
png("../drafts/figs/TT23_fig3_gasExchange.png", width = 8, height = 12,
    units = "in", res = 600)
ggarrange(anet_tri_plot, anet_mai_plot, gsw_tri_plot, gsw_mai_plot,
          l_tri_plot, l_mai_plot, common.legend = TRUE, 
          hjust = 0, legend = "bottom", ncol = 2, nrow = 3, align = "hv",
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), 
          font.label = list(size = 18))
dev.off()

##############################################################################
## Figure 3: Photosynthetic capacity
##############################################################################
png("../drafts/figs/TT23_fig4_photoCapacity.png", 
    width = 8, height = 12, units = "in", res = 600)
ggarrange(vcmax_tri_plot, vcmax_mai_plot, jmax_tri_plot, jmax_mai_plot, 
          jvmax_tri_plot, jvmax_mai_plot, common.legend = TRUE, 
          legend = "bottom", ncol = 2, nrow = 3, align = "hv", hjust = 0,
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), 
          font.label = list(size = 18))
dev.off()

##############################################################################
## Figure S1: Soil nitrogen components
##############################################################################
png("../drafts/figs/TT23_figS1_nitrate_ammonium.png", 
    width = 10, height = 4.5, units = "in", res = 600)
ggarrange(nitrate_plot, ammonium_plot,
          common.legend = TRUE, legend = "right", ncol = 2, nrow = 1, 
          align = "hv", font.label = list(size = 18), hjust = 0,
          labels = c("(a)", "(b)"))
dev.off()

##############################################################################
## Figure S2: Chlorophyll fluorescence
##############################################################################
png("../drafts/figs/TT23_figS2_chlorophyll.png", 
    width = 10, height = 4.5, units = "in", res = 600)
ggarrange(spad_tri_plot, spad_mai_plot,
          common.legend = TRUE, legend = "right", ncol = 2, nrow = 1, 
          align = "hv", font.label = list(size = 18), hjust = 0,
          labels = c("(a)", "(b)"))
dev.off()

##############################################################################
## ESA talk figure: net photosynthesis
##############################################################################
png("../drafts/figs/TT23_ESAtalk_anet.png", 
    width = 10, height = 4.5, units = "in", res = 600)
ggarrange(anet_tri_plot, anet_mai_plot,
          common.legend = TRUE, legend = "right", ncol = 2, nrow = 1, 
          align = "hv", font.label = list(size = 18), hjust = 0,
          labels = c("(a)", "(b)"))
dev.off()

##############################################################################
## ESA talk figure: 
##############################################################################
png("../drafts/figs/TT23_ESAtalk_tri_gasEx.png", 
    width = 10, height = 4.5, units = "in", res = 600)
ggarrange(vcmax_tri_plot, jmax_tri_plot,
          common.legend = TRUE, legend = "right", ncol = 2, nrow = 1, 
          align = "hv", font.label = list(size = 18), hjust = 0,
          labels = c("(a)", "(b)"))
dev.off()

##############################################################################
## ESA talk figure: 
##############################################################################
png("../drafts/figs/TT23_ESAtalk_mai_gasEx.png", 
    width = 10, height = 4.5, units = "in", res = 600)
ggarrange(gsw_mai_plot, l_mai_plot,
          common.legend = TRUE, legend = "right", ncol = 2, nrow = 1, 
          align = "hv", font.label = list(size = 18), hjust = 0,
          labels = c("(a)", "(b)"))
dev.off()

##############################################################################
## ESA talk figure: soil resources
##############################################################################
png("../drafts/figs/TT23_ESAtalk_soil_conclusions.png", 
    width = 6, height = 9, units = "in", res = 600)
ggarrange(nitrogen_plot, sm_plot,
          common.legend = TRUE, legend = "bottom", ncol = 1, nrow = 2, 
          align = "hv", hjust = 0)
dev.off()

