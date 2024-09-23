## Load libraries
library(tidyverse)
library(ggpubr)
library(lme4)
library(car)
library(emmeans)
library(multcomp)
library(MuMIn)

## Read compiled data file
df <- read.csv("../data/TT23_data.csv") %>%
  mutate(gm.trt = factor(gm.trt, levels = c("ambient", "weeded")),
         canopy = factor(canopy, levels = c("open", "closed")))
head(df)

# Turn off digit rounding in emmean args
emm_options(opt.digits = FALSE)

# Helper fxn to change "NaN" to "NA"
NaN_to_NA <- function(x) ifelse(is.nan(x), NA, x)

## Read and subset soil dataset
df.soil <- df %>%
  group_by(plot, composite, gm.trt, canopy) %>%
  summarize_at(.vars = vars(phosphate_ppm:inorg_n_ppm),
               .funs = mean, na.rm = TRUE) %>%
  mutate(gm.trt = factor(gm.trt, levels = c("ambient", "weeded")),
         canopy = factor(canopy, levels = c("open", "closed"))) %>%
  mutate(across(nitrate_ppm:inorg_n_ppm, .fns = NaN_to_NA),
         np.ratio = inorg_n_ppm/phosphate_ppm)

## Read daily soil moisture dataset
df.sm <- read.csv("../data/TT23_tomst_probe_sm_daily.csv")

##############################################################################
## N availability (nitrate + ammonium)
##############################################################################
plant_availableN <- lmer(
  log(inorg_n_ppm) ~ gm.trt * canopy + (1 | plot), data = df.soil)

# Check model assumptions
plot(plant_availableN)
qqnorm(residuals(plant_availableN))
qqline(residuals(plant_availableN))
densityPlot(residuals(plant_availableN))
shapiro.test(residuals(plant_availableN))
outlierTest(plant_availableN)

# Model output
summary(plant_availableN)
Anova(plant_availableN)
r.squaredGLMM(plant_availableN)

# Pairwise comparisons
emmeans(plant_availableN, pairwise~canopy, type = "response")
emmeans(plant_availableN, pairwise~gm.trt, type = "response")

# % change canopy
(4.10 - 16.97) / 16.97 * 100

# % change gm.trt
(9.00 - 7.74) / 7.74 * 100

##############################################################################
## Phosphate
##############################################################################
phosphate <- lmer(
  phosphate_ppm ~ gm.trt * canopy + (1 | plot), data = df.soil)

# Check model assumptions
plot(phosphate)
qqnorm(residuals(phosphate))
qqline(residuals(phosphate))
densityPlot(residuals(phosphate))
shapiro.test(residuals(phosphate))
outlierTest(phosphate)

# Model output
summary(phosphate)
Anova(phosphate)
r.squaredGLMM(phosphate)

# Pairwise comparisons
emmeans(phosphate, pairwise~canopy)
emmeans(phosphate, pairwise~gm.trt)

# % change canopy
(0.813 - 1.095) / 1.095 * 100

# % change gm.trt
(0.88 - 1.03) / 1.03 * 100

##############################################################################
## Nitrate
##############################################################################
nitrate <- lmer(
  nitrate_ppm ~ gm.trt * canopy + (1 | plot), data = df.soil)

# Check model assumptions
plot(nitrate)
qqnorm(residuals(nitrate))
qqline(residuals(nitrate))
densityPlot(residuals(nitrate))
shapiro.test(residuals(nitrate))
outlierTest(nitrate)

# Model output
summary(nitrate)
Anova(nitrate)
r.squaredGLMM(nitrate)

# Pairwise comparisons
emmeans(nitrate, pairwise~canopy)

# % change nitrate with canopy status
(5.47 - 18.57) / 18.57 * 100

##############################################################################
## Ammonium
##############################################################################
ammonium <- lmer(
  sqrt(ammonium_ppm) ~ gm.trt * canopy + (1 | plot), data = df.soil)

# Check model assumptions
plot(ammonium)
qqnorm(residuals(ammonium))
qqline(residuals(ammonium))
densityPlot(residuals(ammonium))
shapiro.test(residuals(ammonium))
outlierTest(ammonium)

# Model output
summary(ammonium)
Anova(ammonium)
r.squaredGLMM(ammonium)

# Pairwise comparisons
cld(emmeans(ammonium, pairwise~canopy*gm.trt, type = "response"))

##############################################################################
## Soil N:P
##############################################################################
df.soil$np.ratio[54] <- NA

n_to_p_ratio <- lmer(
  log(np.ratio) ~ gm.trt * canopy + (1 | plot), data = df.soil)

# Check model assumptions
plot(n_to_p_ratio)
qqnorm(residuals(n_to_p_ratio))
qqline(residuals(n_to_p_ratio))
densityPlot(residuals(n_to_p_ratio))
shapiro.test(residuals(n_to_p_ratio))
outlierTest(n_to_p_ratio)

# Model output
summary(n_to_p_ratio)
Anova(n_to_p_ratio)
r.squaredGLMM(n_to_p_ratio)

# Pairwise comparisons
emmeans(n_to_p_ratio, pairwise~gm.trt, type = "response")
emmeans(n_to_p_ratio, pairwise~canopy, type = "response")

# % change due to gm.trt
(11.76 - 7.8) / 7.8 * 100

# % change canopy
(11.8 - 7.8) / 7.8 * 100

##############################################################################
## Soil moisture (time series) 
##############################################################################
sm_model <- lmer(daily_vwc ~ gm.trt * doy + (1 | plot),
                 data = df.sm)

# Check model assumptions
plot(sm_model)
qqnorm(residuals(sm_model))
qqline(residuals(sm_model))
densityPlot(residuals(sm_model))
shapiro.test(residuals(sm_model))
outlierTest(sm_model)

## Model results
summary(sm_model)
Anova(sm_model)

## Post hoc tests
test(emtrends(sm_model, ~1, "doy"))
emmeans(sm_model, pairwise~trt)

##############################################################################
## Anet - Tri
##############################################################################
anet.tri <- lmer(
  log(anet) ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))

# Check model assumptions
plot(anet.tri)
qqnorm(residuals(anet.tri))
qqline(residuals(anet.tri))
densityPlot(residuals(anet.tri))
shapiro.test(residuals(anet.tri))
outlierTest(anet.tri)

# Model output
summary(anet.tri)
Anova(anet.tri)
r.squaredGLMM(anet.tri)

# Pairwise comparisons
cld(emmeans(anet.tri, pairwise~canopy*gm.trt, type = "response"))
emmeans(anet.tri, pairwise~canopy, type = "response")
emmeans(anet.tri, pairwise~gm.trt, type = "response")

# % change canopy
(4.51 - 12.52) / 12.52 * 100

# % change gm.trt
(7.19 - 7.86) / 7.86 * 100
(4.175 - 4.877) / 4.877  # post canopy 

##############################################################################
## Anet - Mai
##############################################################################
df$anet[91] <- NA

anet.mai <- lmer(
  anet ~ gm.trt * canopy  + (1 | plot), data = subset(df, spp == "Mai"))

# Check model assumptions
plot(anet.mai)
qqnorm(residuals(anet.mai))
qqline(residuals(anet.mai))
densityPlot(residuals(anet.mai))
shapiro.test(residuals(anet.mai))
outlierTest(anet.mai)

# Model output
summary(anet.mai)
Anova(anet.mai)
r.squaredGLMM(anet.mai)

# Pairwise comparisons
emmeans(anet.mai, pairwise~canopy)
emmeans(anet.mai, pairwise~gm.trt)

emmeans(anet.mai, pairwise~canopy * gm.trt)

# % change canopy
(4.02 - 9.85) / 9.85 * 100

# % change gm.trt
(6.25 - 7.62) / 7.62 * 100

# % change gm.trt between open canopy and closed canopy
(9.182 - 10.519) / 10.519 * 100 # open canopy
(3.317 - 4.725) / 4.725 * 100 # closed canopy

##############################################################################
## gs - Tri
##############################################################################
df$gsw[53] <- NA

gsw.tri <- lmer(
  log(gsw) ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))

# Check model assumptions
plot(gsw.tri)
qqnorm(residuals(gsw.tri))
qqline(residuals(gsw.tri))
densityPlot(residuals(gsw.tri))
shapiro.test(residuals(gsw.tri))
outlierTest(gsw.tri)

# Model output
summary(gsw.tri)
Anova(gsw.tri)
r.squaredGLMM(gsw.tri)

# Pairwise comparisons
emmeans(gsw.tri, pairwise~canopy, type = "response")

# Canopy % change
(0.106 - 0.136) / 0.136 * 100

##############################################################################
## gs - Mai
##############################################################################
df$gsw[91] <- NA

gsw.mai <- lmer(
  gsw ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))

# Check model assumptions
plot(gsw.mai)
qqnorm(residuals(gsw.mai))
qqline(residuals(gsw.mai))
densityPlot(residuals(gsw.mai))
shapiro.test(residuals(gsw.mai))
outlierTest(gsw.mai)

# Model output
summary(gsw.mai)
Anova(gsw.mai)
r.squaredGLMM(gsw.mai)

# Pairwise comparisons
emmeans(gsw.mai, pairwise~canopy)
emmeans(gsw.mai, pairwise~gm.trt)

# Canopy % change
(0.0581 - 0.1536) / 0.1536 * 100

# % change gm.trt
(0.0889 - 0.1227) / 0.1227 * 100

##############################################################################
## stomatal limitation - Tri
##############################################################################
l.tri <- lmer(
  l ~ gm.trt * canopy  + (1 | plot), data = subset(df, spp == "Tri" & l > 0))

# Check model assumptions
plot(l.tri)
qqnorm(residuals(l.tri))
qqline(residuals(l.tri))
densityPlot(residuals(l.tri))
shapiro.test(residuals(l.tri))
outlierTest(l.tri)

# Model output
summary(l.tri)
Anova(l.tri)
r.squaredGLMM(l.tri)

# Pairwise comparisons
emmeans(l.tri, pairwise~canopy, type = "response")

# % change canopy
(0.228 - 0.510) / 0.510 * 100

##############################################################################
## stomatal limitation - Mai
##############################################################################
l.mai <- lmer(
  log(l) ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai" & l > 0))

# Check model assumptions
plot(l.mai)
qqnorm(residuals(l.mai))
qqline(residuals(l.mai))
densityPlot(residuals(l.mai))
shapiro.test(residuals(l.mai))
outlierTest(l.mai)

# Model output
summary(l.mai)
Anova(l.mai)
r.squaredGLMM(l.mai)

# Pairwise comparisons
cld(emmeans(l.mai, pairwise~gm.trt*canopy, type = "response"))
emmeans(l.mai, pairwise~canopy, type = "response")
emmeans(l.mai, pairwise~gm.trt, type = "response")

# % change canopy
(0.328 - 0.376) / 0.376 * 100

# % change gm.trt
(0.398 - 0.310) / 0.310 * 100

##############################################################################
## SPAD - Tri
##############################################################################
spad.tri <- lmer(
  SPAD ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))

# Check model assumptions
plot(spad.tri)
qqnorm(residuals(spad.tri))
qqline(residuals(spad.tri))
densityPlot(residuals(spad.tri))
shapiro.test(residuals(spad.tri))
outlierTest(spad.tri)

# Model output
summary(spad.tri)
Anova(spad.tri)
r.squaredGLMM(spad.tri)

# Pairwise comparisons
emmeans(spad.tri, pairwise~canopy)

# % change canopy
(44.587 - 35.52) / 35.52 * 100

##############################################################################
## SPAD - Mai
##############################################################################
df$SPAD[111] <- NA

spad.mai <- lmer(
  log(SPAD) ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))

# Check model assumptions
plot(spad.mai)
qqnorm(residuals(spad.mai))
qqline(residuals(spad.mai))
densityPlot(residuals(spad.mai))
shapiro.test(residuals(spad.mai))
outlierTest(spad.mai)

# Model output
summary(spad.mai)
Anova(spad.mai)
r.squaredGLMM(spad.mai)

# Pairwise comparisons
emmeans(spad.mai, pairwise~canopy, type = "response")

# % change canopy
(39.68 - 26.28) / 26.28 * 100

##############################################################################
## Vcmax - Tri
##############################################################################
df$vcmax25[180] <- NA

vcmax.tri <- lmer(
  log(vcmax25) ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))

# Check model assumptions
plot(vcmax.tri)
qqnorm(residuals(vcmax.tri))
qqline(residuals(vcmax.tri))
densityPlot(residuals(vcmax.tri))
shapiro.test(residuals(vcmax.tri))
outlierTest(vcmax.tri)

# Model output
summary(vcmax.tri)
Anova(vcmax.tri)
r.squaredGLMM(vcmax.tri)

# Pairwise comparisons
emmeans(vcmax.tri, pairwise~canopy, type = "response")

# % change canopy
(23.690 - 99.356) / 99.356 * 100


# What is the mean +/- SD of Tri Vcmax25?
df %>%
  filter(spp == "Tri") %>%
  summarize(vcmax.mean = mean(vcmax25, na.rm = TRUE),
            vcmax.stdev = sd(vcmax25, na.rm = TRUE))

##############################################################################
## Vcmax - Mai
##############################################################################
vcmax.mai <- lmer(
  vcmax25 ~ gm.trt * canopy  + (1 | plot), data = subset(df, spp == "Mai"))

# Check model assumptions
plot(vcmax.mai)
qqnorm(residuals(vcmax.mai))
qqline(residuals(vcmax.mai))
densityPlot(residuals(vcmax.mai))
shapiro.test(residuals(vcmax.mai))
outlierTest(vcmax.mai)

# Model output
summary(vcmax.mai)
Anova(vcmax.mai)
r.squaredGLMM(vcmax.mai)

# Pairwise comparisons
emmeans(vcmax.mai, pairwise~canopy, type = "response")
emmeans(vcmax.mai, pairwise~gm.trt, type = "response")

# % change canopy
(25.081 - 58.061) / 58.061 * 100

# What is the mean +/- SD of Tri Vcmax25?
df %>%
  filter(spp == "Mai") %>%
  summarize(vcmax.mean = mean(vcmax25, na.rm = TRUE),
            vcmax.stdev = sd(vcmax25, na.rm = TRUE))

##############################################################################
## Jmax - Tri
##############################################################################
df$jmax25[c(142, 180)] <- NA

jmax.tri <- lmer(
  log(jmax25) ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))

# Check model assumptions
plot(jmax.tri)
qqnorm(residuals(jmax.tri))
qqline(residuals(jmax.tri))
densityPlot(residuals(jmax.tri))
shapiro.test(residuals(jmax.tri))
outlierTest(jmax.tri)

# Model output
summary(jmax.tri)
Anova(jmax.tri)
r.squaredGLMM(jmax.tri)

# Pairwise comparisons
emmeans(jmax.tri, pairwise~canopy, type = "response")
emmeans(jmax.tri, pairwise~gm.trt, type = "response")
cld(emmeans(jmax.tri, pairwise~gm.trt*canopy, type = "response"))

# % change canopy
(44.570 - 176.151) / 176.151 * 100

# % change GM trt
(85.055 - 92.306) / 92.306 * 100

##############################################################################
## Jmax - Mai
##############################################################################
jmax.mai <- lmer(
  log(jmax25) ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))

# Check model assumptions
plot(jmax.mai)
qqnorm(residuals(jmax.mai))
qqline(residuals(jmax.mai))
densityPlot(residuals(jmax.mai))
shapiro.test(residuals(jmax.mai))
outlierTest(jmax.mai)

# Model output
summary(jmax.mai)
Anova(jmax.mai)
r.squaredGLMM(jmax.mai)

# Pairwise comparisons
emmeans(jmax.mai, pairwise~canopy, type = "response")
emmeans(jmax.mai, pairwise~gm.trt, type = "response")

cld(emmeans(jmax.mai, pairwise~canopy*gm.trt))


# % change canopy
(43.887 - 102.424) / 102.424 * 100

# % change gm.trt
(66.3 - 67.8) / 67.8 * 100

##############################################################################
## Jmax : Vcmax - Tri
##############################################################################
df$jmax.vcmax[181] <- NA

jmax.vcmax.tri <- lmer(
  jmax.vcmax ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Tri"))

# Check model assumptions
plot(jmax.vcmax.tri)
qqnorm(residuals(jmax.vcmax.tri))
qqline(residuals(jmax.vcmax.tri))
densityPlot(residuals(jmax.vcmax.tri))
shapiro.test(residuals(jmax.vcmax.tri))
outlierTest(jmax.vcmax.tri)

# Model output
summary(jmax.vcmax.tri)
Anova(jmax.vcmax.tri)
r.squaredGLMM(jmax.vcmax.tri)

# Pairwise comparisons
emmeans(jmax.vcmax.tri, pairwise~canopy)

# % change canopy
(1.85 - 1.78) / 1.78 * 100

##############################################################################
## Jmax : Vcmax - Mai
##############################################################################
df$jmax.vcmax[c(113, 218, 219)] <- NA

jmax.vcmax.mai <- lmer(
  jmax.vcmax ~ gm.trt * canopy + (1 | plot), data = subset(df, spp == "Mai"))

# Check model assumptions
plot(jmax.vcmax.mai)
qqnorm(residuals(jmax.vcmax.mai))
qqline(residuals(jmax.vcmax.mai))
densityPlot(residuals(jmax.vcmax.mai))
shapiro.test(residuals(jmax.vcmax.mai))
outlierTest(jmax.vcmax.mai)


# Model output
summary(jmax.vcmax.mai)
Anova(jmax.vcmax.mai)
r.squaredGLMM(jmax.vcmax.mai)

# Pairwise comparisons
emmeans(jmax.vcmax.mai, pairwise~canopy)
emmeans(jmax.vcmax.mai, pairwise~gm.trt)
cld(emmeans(jmax.vcmax.mai, pairwise~canopy*gm.trt))

# % change canopy
(1.702 - 1.808) / 1.808 * 100 

# % change gm.trt
(1.73 - 1.78) / 1.78 * 100

##############################################################################
## Write Table 1: Soil nutrients
##############################################################################
soil.nitrogen.table <- data.frame(Anova(plant_availableN)) %>%
  mutate(treatment = row.names(.),
         chisq_soilN = Chisq,
         p_soilN = Pr..Chisq.,
         across(chisq_soilN:p_soilN, \(x) round(x, digits = 3)),
         chisq_soilN = ifelse(chisq_soilN <0.001 & chisq_soilN >= 0, 
                                "<0.001", chisq_soilN),
         p_soilN = ifelse(p_soilN <0.001 & p_soilN >= 0, 
                               "<0.001", p_soilN)) %>%
  dplyr::select(treatment, Df, chisq_soilN, p_soilN)

soil.nitrate.table <- data.frame(Anova(nitrate)) %>%
  mutate(treatment = row.names(.),
         chisq_nitrate = Chisq,
         p_nitrate = Pr..Chisq.,
         across(chisq_nitrate:p_nitrate, \(x) round(x, digits = 3)),
         chisq_nitrate = ifelse(chisq_nitrate < 0.001 & chisq_nitrate >= 0, 
                              "<0.001", chisq_nitrate),
         p_nitrate = ifelse(p_nitrate <0.001 & p_nitrate >= 0, 
                          "<0.001", p_nitrate)) %>%
  dplyr::select(treatment, chisq_nitrate, p_nitrate)

soil.ammonium.table <- data.frame(Anova(ammonium)) %>%
  mutate(treatment = row.names(.),
         chisq_ammonium = Chisq,
         p_ammonium = Pr..Chisq.,
         across(chisq_ammonium:p_ammonium, \(x) round(x, digits = 3)),
         chisq_ammonium = ifelse(chisq_ammonium < 0.001 & chisq_ammonium >= 0, 
                                "<0.001", chisq_ammonium),
         p_ammonium = ifelse(p_ammonium <0.001 & p_ammonium >= 0, 
                            "<0.001", p_ammonium)) %>%
  dplyr::select(treatment, chisq_ammonium, p_ammonium)

soil.phosphate.table <- data.frame(Anova(phosphate)) %>%
  mutate(treatment = row.names(.),
         chisq_phosphate = Chisq,
         p_phosphate = Pr..Chisq.,
         across(chisq_phosphate:p_phosphate, \(x) round(x, digits = 3)),
         chisq_phosphate = ifelse(chisq_phosphate < 0.001 & chisq_phosphate >= 0, 
                                 "<0.001", chisq_phosphate),
         p_phosphate = ifelse(p_phosphate <0.001 & p_phosphate >= 0, 
                             "<0.001", p_phosphate)) %>%
  dplyr::select(treatment, chisq_phosphate, p_phosphate)

table1 <- soil.nitrogen.table %>% full_join(soil.nitrate.table) %>% 
  full_join(soil.ammonium.table) %>% full_join(soil.phosphate.table)

write.csv(table1, "../drafts/tables/TT23_tableS1_soil_nutrients.csv",
          row.names = FALSE)

##############################################################################
## Write Table 2: Gas exchange
##############################################################################
anet.tri.table <- data.frame(Anova(anet.tri)) %>%
  mutate(treatment = row.names(.),
         chisq_anet.tri = Chisq,
         p_anet.tri = Pr..Chisq.,
         across(chisq_anet.tri:p_anet.tri, \(x) round(x, digits = 3)),
         chisq_anet.tri = ifelse(chisq_anet.tri < 0.001 & chisq_anet.tri >= 0, 
                              "<0.001", chisq_anet.tri),
         p_anet.tri = ifelse(p_anet.tri <0.001 & p_anet.tri >= 0, 
                          "<0.001", p_anet.tri)) %>%
  dplyr::select(treatment, Df, chisq_anet.tri, p_anet.tri)

anet.mai.table <- data.frame(Anova(anet.mai)) %>%
  mutate(treatment = row.names(.),
         chisq_anet.mai = Chisq,
         p_anet.mai = Pr..Chisq.,
         across(chisq_anet.mai:p_anet.mai, \(x) round(x, digits = 3)),
         chisq_anet.mai = ifelse(chisq_anet.mai < 0.001 & chisq_anet.mai >= 0, 
                                 "<0.001", chisq_anet.mai),
         p_anet.mai = ifelse(p_anet.mai <0.001 & p_anet.mai >= 0, 
                             "<0.001", p_anet.mai)) %>%
  dplyr::select(treatment, chisq_anet.mai, p_anet.mai)

gsw.tri.table <- data.frame(Anova(gsw.tri)) %>%
  mutate(treatment = row.names(.),
         chisq_gsw.tri = Chisq,
         p_gsw.tri = Pr..Chisq.,
         across(chisq_gsw.tri:p_gsw.tri, \(x) round(x, digits = 3)),
         chisq_gsw.tri = ifelse(chisq_gsw.tri < 0.001 & chisq_gsw.tri >= 0, 
                                 "<0.001", chisq_gsw.tri),
         p_gsw.tri = ifelse(p_gsw.tri <0.001 & p_gsw.tri >= 0, 
                            "<0.001", p_gsw.tri)) %>%
  dplyr::select(treatment, chisq_gsw.tri, p_gsw.tri)

gsw.mai.table <- data.frame(Anova(gsw.mai)) %>%
  mutate(treatment = row.names(.),
         chisq_gsw.mai = Chisq,
         p_gsw.mai = Pr..Chisq.,
         across(chisq_gsw.mai:p_gsw.mai, \(x) round(x, digits = 3)),
         chisq_gsw.mai = ifelse(chisq_gsw.mai < 0.001 & chisq_gsw.mai >= 0, 
                                "<0.001", chisq_gsw.mai),
         p_gsw.mai = ifelse(p_gsw.mai <0.001 & p_gsw.mai >= 0, 
                            "<0.001", p_gsw.mai)) %>%
  dplyr::select(treatment, chisq_gsw.mai, p_gsw.mai)


l.tri.table <- data.frame(Anova(l.tri)) %>%
  mutate(treatment = row.names(.),
         chisq_l.tri = Chisq,
         p_l.tri = Pr..Chisq.,
         across(chisq_l.tri:p_l.tri, \(x) round(x, digits = 3)),
         chisq_l.tri = ifelse(chisq_l.tri < 0.001 & 
                                chisq_l.tri >= 0, 
                              "<0.001", chisq_l.tri),
         p_l.tri = ifelse(p_l.tri <0.001 & p_l.tri >= 0, 
                          "<0.001", p_l.tri)) %>%
  dplyr::select(treatment, chisq_l.tri, p_l.tri)

l.mai.table <- data.frame(Anova(l.mai)) %>%
  mutate(treatment = row.names(.),
         chisq_l.mai = Chisq,
         p_l.mai = Pr..Chisq.,
         across(chisq_l.mai:p_l.mai, \(x) round(x, digits = 3)),
         chisq_l.mai = ifelse(chisq_l.mai < 0.001 & 
                                chisq_l.mai >= 0, 
                              "<0.001", chisq_l.mai),
         p_l.mai = ifelse(p_l.mai <0.001 & p_l.mai >= 0, 
                          "<0.001", p_l.mai)) %>%
  dplyr::select(treatment, chisq_l.mai, p_l.mai)

spad.tri.table <- data.frame(Anova(spad.tri)) %>%
  mutate(treatment = row.names(.),
         chisq_spad.tri = Chisq,
         p_spad.tri = Pr..Chisq.,
         across(chisq_spad.tri:p_spad.tri, \(x) round(x, digits = 3)),
         chisq_spad.tri = ifelse(chisq_spad.tri < 0.001 & chisq_spad.tri >= 0, 
                                 "<0.001", chisq_spad.tri),
         p_spad.tri = ifelse(p_spad.tri <0.001 & p_spad.tri >= 0, 
                             "<0.001", p_spad.tri)) %>%
  dplyr::select(treatment, Df, chisq_spad.tri, p_spad.tri)

spad.mai.table <- data.frame(Anova(spad.mai)) %>%
  mutate(treatment = row.names(.),
         chisq_spad.mai = Chisq,
         p_spad.mai = Pr..Chisq.,
         across(chisq_spad.mai:p_spad.mai, \(x) round(x, digits = 3)),
         chisq_spad.mai = ifelse(chisq_spad.mai < 0.001 & chisq_spad.mai >= 0, 
                                 "<0.001", chisq_spad.mai),
         p_spad.mai = ifelse(p_spad.mai <0.001 & p_spad.mai >= 0, 
                             "<0.001", p_spad.mai)) %>%
  dplyr::select(treatment, chisq_spad.mai, p_spad.mai)

table2 <- anet.tri.table %>% full_join(gsw.tri.table) %>% 
  full_join(l.tri.table) %>% full_join(spad.tri.table) %>% 
  full_join(anet.mai.table) %>% full_join(gsw.mai.table) %>% 
  full_join(l.mai.table) %>% full_join(spad.mai.table)

write.csv(table2, "../drafts/tables/TT23_table2_gas_exchange.csv",
          row.names = FALSE)

##############################################################################
## Write Table 3: Indices of photosynthetic capacity
##############################################################################
vcmax.tri <- data.frame(Anova(vcmax.tri)) %>%
  mutate(treatment = row.names(.),
         chisq_vcmax.tri = Chisq,
         p_vcmax.tri = Pr..Chisq.,
         across(chisq_vcmax.tri:p_vcmax.tri, \(x) round(x, digits = 3)),
         chisq_vcmax.tri = ifelse(chisq_vcmax.tri < 0.001 & chisq_vcmax.tri >= 0, 
                                 "<0.001", chisq_vcmax.tri),
         p_vcmax.tri = ifelse(p_vcmax.tri <0.001 & p_vcmax.tri >= 0, 
                             "<0.001", p_vcmax.tri)) %>%
  dplyr::select(treatment, Df, chisq_vcmax.tri, p_vcmax.tri)

vcmax.mai <- data.frame(Anova(vcmax.mai)) %>%
  mutate(treatment = row.names(.),
         chisq_vcmax.mai = Chisq,
         p_vcmax.mai = Pr..Chisq.,
         across(chisq_vcmax.mai:p_vcmax.mai, \(x) round(x, digits = 3)),
         chisq_vcmax.mai = ifelse(chisq_vcmax.mai < 0.001 & chisq_vcmax.mai >= 0, 
                                 "<0.001", chisq_vcmax.mai),
         p_vcmax.mai = ifelse(p_vcmax.mai <0.001 & p_vcmax.mai >= 0, 
                             "<0.001", p_vcmax.mai)) %>%
  dplyr::select(treatment, chisq_vcmax.mai, p_vcmax.mai)

jmax.tri <- data.frame(Anova(jmax.tri)) %>%
  mutate(treatment = row.names(.),
         chisq_jmax.tri = Chisq,
         p_jmax.tri = Pr..Chisq.,
         across(chisq_jmax.tri:p_jmax.tri, \(x) round(x, digits = 3)),
         chisq_jmax.tri = ifelse(chisq_jmax.tri < 0.001 & chisq_jmax.tri >= 0, 
                                  "<0.001", chisq_jmax.tri),
         p_jmax.tri = ifelse(p_jmax.tri <0.001 & p_jmax.tri >= 0, 
                              "<0.001", p_jmax.tri)) %>%
  dplyr::select(treatment, Df, chisq_jmax.tri, p_jmax.tri)

jmax.mai <- data.frame(Anova(jmax.mai)) %>%
  mutate(treatment = row.names(.),
         chisq_jmax.mai = Chisq,
         p_jmax.mai = Pr..Chisq.,
         across(chisq_jmax.mai:p_jmax.mai, \(x) round(x, digits = 3)),
         chisq_jmax.mai = ifelse(chisq_jmax.mai < 0.001 & chisq_jmax.mai >= 0, 
                                  "<0.001", chisq_jmax.mai),
         p_jmax.mai = ifelse(p_jmax.mai <0.001 & p_jmax.mai >= 0, 
                              "<0.001", p_jmax.mai)) %>%
  dplyr::select(treatment, chisq_jmax.mai, p_jmax.mai)

jvmax.tri <- data.frame(Anova(jmax.vcmax.tri)) %>%
  mutate(treatment = row.names(.),
         chisq_jvmax.tri = Chisq,
         p_jvmax.tri = Pr..Chisq.,
         across(chisq_jvmax.tri:p_jvmax.tri, \(x) round(x, digits = 3)),
         chisq_jvmax.tri = ifelse(chisq_jvmax.tri < 0.001 & chisq_jvmax.tri >= 0, 
                                 "<0.001", chisq_jvmax.tri),
         p_jvmax.tri = ifelse(p_jvmax.tri <0.001 & p_jvmax.tri >= 0, 
                             "<0.001", p_jvmax.tri)) %>%
  dplyr::select(treatment, Df, chisq_jvmax.tri, p_jvmax.tri)

jvmax.mai <- data.frame(Anova(jmax.vcmax.mai)) %>%
  mutate(treatment = row.names(.),
         chisq_jvmax.mai = Chisq,
         p_jvmax.mai = Pr..Chisq.,
         across(chisq_jvmax.mai:p_jvmax.mai, \(x) round(x, digits = 3)),
         chisq_jvmax.mai = ifelse(chisq_jvmax.mai < 0.001 & chisq_jvmax.mai >= 0, 
                                 "<0.001", chisq_jvmax.mai),
         p_jvmax.mai = ifelse(p_jvmax.mai <0.001 & p_jvmax.mai >= 0, 
                             "<0.001", p_jvmax.mai)) %>%
  dplyr::select(treatment, chisq_jvmax.mai, p_jvmax.mai)

table3 <- vcmax.tri %>% full_join(jmax.tri) %>%  full_join(jvmax.tri) %>% 
  full_join(vcmax.mai) %>%  full_join(jmax.mai) %>% full_join(jvmax.mai)

write.csv(table3, "../drafts/tables/TT23_table3_photoCapacity.csv", 
          row.names = FALSE)