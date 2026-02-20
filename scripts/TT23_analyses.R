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
  dplyr::select(id, spp = spp.x, plot:inorg_n_ppm, n_meas)

# Turn off digit rounding in emmean args
emm_options(opt.digits = FALSE)

## Read and subset soil dataset
df.soil <- df %>%
  distinct(plot, composite, canopy, .keep_all = TRUE) %>%
  mutate(gm.trt = factor(gm.trt, levels = c("ambient", "weeded")),
         canopy = factor(canopy, levels = c("open", "closed")),
         np.ratio = inorg_n_ppm/phosphate_ppm)

## Read daily soil moisture dataset
df.sm <- read.csv("../data/TT23_tomst_probe_sm_daily.csv")

## How many Trillium and Maianthemum individuals?
unique(subset(df2, spp == "Tri")$id)
unique(subset(df2, spp == "Mai")$id)

df2 %>% group_by(spp, plot, gm.trt) %>%
  summarize(n_spp = length(id))

## Mean Vcmax25 for each spp (for discussion section)
df2 %>%
  group_by(spp) %>%
  summarize(vcmax25_mean = mean(vcmax25, na.rm = TRUE),
            vcmax25_sd = sd(vcmax25, na.rm = TRUE))

##############################################################################
## N availability (nitrate + ammonium)
##############################################################################
plant_availableN <- lmer(
  inorg_n_ppm ~ gm.trt * canopy + (1 | plot), data = df.soil)

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
emmeans(plant_availableN, pairwise~canopy)
emmeans(plant_availableN, pairwise~gm.trt)

# % change canopy
(5.81 - 18.57) / 18.57 * 100

# % change gm.trt
(13.417 - 10.965) / 10.965 * 100

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
(0.882 - 1.026) / 1.026 * 100

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
emmeans(nitrate, pairwise~gm.trt)

# % change nitrate with canopy status
(5.465 - 18.573) / 18.573 * 100

# % change nitrate with gm.trt

(13.254 - 10.784) / 10.784 * 100


##############################################################################
## Ammonium
##############################################################################
ammonium <- lmer(
  log(ammonium_ppm) ~ gm.trt * canopy + (1 | plot), data = df.soil)

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
n_to_p_ratio <- lmer(log(np.ratio) ~ gm.trt * canopy + (1 | plot), data = df.soil)

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
(11.792 - 8.251) / 8.251 * 100

# % change canopy
(16.746 - 5.810) / 5.810 * 100

##############################################################################
## Soil moisture (time series) 
##############################################################################
sm_model <- lmer(daily_vwc ~ gm.trt * doy + (1 | plot), data = df.sm)

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
emmeans(sm_model, pairwise~gm.trt)

##############################################################################
## Anet - Tri
##############################################################################
anet.tri <- lmer(log(anet) ~ gm.trt * canopy + (1 | plot) + (1 | id), 
  data = subset(df2, spp == "Tri"))

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
(4.444 - 12.363) / 12.363 * 100

# % change gm.trt
(7.19 - 7.86) / 7.86 * 100
(4.128 - 4.784) / 4.784  # post canopy 

##############################################################################
## gsw - Tri
##############################################################################
gsw.tri <- lmer(log(gsw) ~ gm.trt * canopy + (1 | plot) + (1 | id), 
                data = subset(df2, spp == "Tri"))

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

# % change canopy
(0.102 - 0.128) / 0.128 * 100

##############################################################################
## stomatal limitation - Tri
##############################################################################
l.tri <- lmer(l ~ gm.trt * canopy + (1 | plot) + (1 | id), 
              data = subset(df2, spp == "Tri"))

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
emmeans(l.tri, pairwise~canopy)

# % change canopy
(0.239 - 0.531) / 0.531 * 100

##############################################################################
## Vcmax25 - Tri
##############################################################################
vcmax25.tri <- lmer(log(vcmax25) ~ gm.trt * canopy + (1 | plot) + (1 | id), 
                    data = subset(df2, spp == "Tri" & id != "4414" & id != "2988"))

# Check model assumptions
plot(vcmax25.tri)
qqnorm(residuals(vcmax25.tri))
qqline(residuals(vcmax25.tri))
densityPlot(residuals(vcmax25.tri))
shapiro.test(residuals(vcmax25.tri))
outlierTest(vcmax25.tri)

# Model output
summary(vcmax25.tri)
Anova(vcmax25.tri)
r.squaredGLMM(vcmax25.tri)

# Pairwise comparisons
cld(emmeans(vcmax25.tri, pairwise~gm.trt*canopy, type = "response"))
emmeans(vcmax25.tri, pairwise~gm.trt, type = "response")
emmeans(vcmax25.tri, pairwise~canopy, type = "response")

# % change gm.trt
(47.452 - 51.209) / 51.209 * 100

# % change canopy
(24.023 - 101.152) / 101.152 * 100

##############################################################################
## Jmax25 - Tri
##############################################################################
jmax25.tri <- lmer(log(jmax25) ~ gm.trt * canopy + (1 | plot) + (1 | id), 
                   data = subset(df2, spp == "Tri" & id != "4414" & id != "2988"))

# Check model assumptions
plot(jmax25.tri)
qqnorm(residuals(jmax25.tri))
qqline(residuals(jmax25.tri))
densityPlot(residuals(jmax25.tri))
shapiro.test(residuals(jmax25.tri))
outlierTest(jmax25.tri)

# Model output
summary(jmax25.tri)
Anova(jmax25.tri)
r.squaredGLMM(jmax25.tri)

# Pairwise comparisons
cld(emmeans(jmax25.tri, pairwise~gm.trt*canopy, type = "response"))
emmeans(jmax25.tri, pairwise~gm.trt, type = "response")
emmeans(jmax25.tri, pairwise~canopy, type = "response")

# % change gm.trt
(85.409 - 94.635) / 94.635 * 100

# % change canopy
(45.298 - 178.433) / 178.433 * 100

##############################################################################
## Jmax25:Vcmax25 - Tri
##############################################################################
jmax25_vcmax25.tri <- lmer(jmax.vcmax ~ gm.trt * canopy + (1 | plot) + (1 | id), 
                           data = subset(df2, spp == "Tri"  & id != "4431"))

# Check model assumptions
plot(jmax25_vcmax25.tri)
qqnorm(residuals(jmax25_vcmax25.tri))
qqline(residuals(jmax25_vcmax25.tri))
densityPlot(residuals(jmax25_vcmax25.tri))
shapiro.test(residuals(jmax25_vcmax25.tri))
outlierTest(jmax25_vcmax25.tri)

# Model output
summary(jmax25_vcmax25.tri)
Anova(jmax25_vcmax25.tri)
r.squaredGLMM(jmax25_vcmax25.tri)

# Pairwise comparisons
emmeans(jmax25_vcmax25.tri, pairwise~gm.trt)
emmeans(jmax25_vcmax25.tri, pairwise~canopy)

# % change gm.trt
(1.790 - 1.852) / 1.852 * 100

# % change canopy
(1.871 - 1.770) / 1.770 * 100

##############################################################################
## SPAD - Tri
##############################################################################
spad.tri <- lmer(SPAD ~ gm.trt * canopy + (1 | plot) + (1 | id), 
                 data = subset(df2, spp == "Tri"))

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
cld(emmeans(spad.tri, pairwise~gm.trt*canopy))
emmeans(spad.tri, pairwise~gm.trt)
emmeans(spad.tri, pairwise~canopy)

# % change canopy
(45.530 - 34.714) / 34.714 * 100

##############################################################################
## Anet - Mai
##############################################################################
anet.mai <- lmer(anet ~ gm.trt * canopy + (1 | plot) + (1 | id), 
                 data = subset(df2, spp == "Mai" & id != "5069"))

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
emmeans(anet.mai, pairwise~canopy, type = "response")
emmeans(anet.mai, pairwise~gm.trt, type = "response")
cld(emmeans(anet.mai, pairwise~canopy*gm.trt, type = "response"))

# % change canopy
(3.806 - 9.96) / 9.96 * 100

# % change gm.trt
(5.50 - 6.89) / 6.89 * 100

##############################################################################
## gsw - Mai
##############################################################################
gsw.mai <- lmer(gsw ~ gm.trt * canopy + (1 | plot) + (1 | id), 
                data = subset(df2, spp == "Mai" & id != "5069"))

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
emmeans(gsw.mai, pairwise~gm.trt)
emmeans(gsw.mai, pairwise~canopy)
emmeans(gsw.mai, pairwise~canopy*gm.trt)


# % change canopy
(0.058 - 0.156) / 0.156 * 100

# % change gm.trt
(0.094 - 0.120) / 0.120 * 100

##############################################################################
## stomatal limitation - Mai
##############################################################################
l.mai <- lmer(log(l) ~ gm.trt * canopy + (1 | plot) + (1 | id), 
              data = subset(df2, spp == "Mai"))

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
emmeans(l.mai, pairwise~gm.trt, type = "response")
emmeans(l.mai, pairwise~canopy, type = "response")

# % change gm.trt
(0.398 - 0.321) / 0.321 * 100

# % change canopy
(0.383 - 0.335) / 0.335 * 100

##############################################################################
## Vcmax25 - Mai
##############################################################################
vcmax25.mai <- lmer(log(vcmax25) ~ gm.trt * canopy + (1 | plot) + (1 | id), 
                    data = subset(df2, spp == "Mai"))

# Check model assumptions
plot(vcmax25.mai)
qqnorm(residuals(vcmax25.mai))
qqline(residuals(vcmax25.mai))
densityPlot(residuals(vcmax25.mai))
shapiro.test(residuals(vcmax25.mai))
outlierTest(vcmax25.mai)

# Model output
summary(vcmax25.mai)
Anova(vcmax25.mai)
r.squaredGLMM(vcmax25.mai)

# Pairwise comparisons
emmeans(vcmax25.mai, pairwise~canopy, type = "response")

# % change canopy
(26.032 - 59.038) / 59.038 * 100

##############################################################################
## Jmax25 - Tri
##############################################################################
jmax25.mai <- lmer(log(jmax25) ~ gm.trt * canopy + (1 | plot) + (1 | id), 
                   data = subset(df2, spp == "Mai"))

# Check model assumptions
plot(jmax25.mai)
qqnorm(residuals(jmax25.mai))
qqline(residuals(jmax25.mai))
densityPlot(residuals(jmax25.mai))
shapiro.test(residuals(jmax25.mai))
outlierTest(jmax25.mai)

# Model output
summary(jmax25.mai)
Anova(jmax25.mai)
r.squaredGLMM(jmax25.mai)

# Pairwise comparisons
emmeans(jmax25.mai, pairwise~canopy, type = "response")

# % change canopy
(45.161 - 105.718) / 105.718 * 100

##############################################################################
## Jmax25:Vcmax25 - Tri
##############################################################################
jmax25_vcmax25.mai <- lmer(jmax.vcmax ~ gm.trt * canopy + (1 | plot) + (1 | id), 
                           data = subset(df2, spp == "Mai"))

# Check model assumptions
plot(jmax25_vcmax25.mai)
qqnorm(residuals(jmax25_vcmax25.mai))
qqline(residuals(jmax25_vcmax25.mai))
densityPlot(residuals(jmax25_vcmax25.mai))
shapiro.test(residuals(jmax25_vcmax25.mai))
outlierTest(jmax25_vcmax25.mai)

# Model output
summary(jmax25_vcmax25.mai)
Anova(jmax25_vcmax25.mai)
r.squaredGLMM(jmax25_vcmax25.mai)

##############################################################################
## SPAD - Mai
##############################################################################
spad.mai <- lmer(SPAD ~ gm.trt * canopy + (1 | plot) + (1 | id), 
                 data = subset(df2, spp == "Mai"))

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
emmeans(spad.mai, pairwise~canopy)

# % change canopy
(39.769 - 26.955) / 26.955 * 100

##############################################################################
## Write Table 3: Soil nutrients
##############################################################################
# Soil inorganic nitrogen
soil.nitrogen.table <- data.frame(Anova(plant_availableN)) %>%
  dplyr::select(chisq = Chisq, df = Df, p_val = Pr..Chisq.) %>%
  mutate(trait = "soil_inorg_N",
         treatment = row.names(.),
         across(chisq:p_val, \(x) round(x, digits = 3)),
         chisq = ifelse(chisq < 0.001 & chisq >= 0, "<0.001", chisq),
         p_val = ifelse(p_val < 0.001 & p_val >= 0, "<0.001", p_val)) %>%
  pivot_wider(names_from = treatment, values_from = chisq:p_val) %>%
  dplyr::select(trait, df = df_gm.trt, chisq_gm.trt, p_val_gm.trt,
                chisq_canopy, p_val_canopy, chisq_int = `chisq_gm.trt:canopy`, 
                p_val_int = `p_val_gm.trt:canopy`)

# Soil nitrate availability
soil.nitrate.table <- data.frame(Anova(nitrate)) %>%
  dplyr::select(chisq = Chisq, df = Df, p_val = Pr..Chisq.) %>%
  mutate(trait = "soil_nitrate",
         treatment = row.names(.),
         across(chisq:p_val, \(x) round(x, digits = 3)),
         chisq = ifelse(chisq < 0.001 & chisq >= 0, "<0.001", chisq),
         p_val = ifelse(p_val < 0.001 & p_val >= 0, "<0.001", p_val)) %>%
  pivot_wider(names_from = treatment, values_from = chisq:p_val) %>%
  dplyr::select(trait, df = df_gm.trt, chisq_gm.trt, p_val_gm.trt,
                chisq_canopy, p_val_canopy, chisq_int = `chisq_gm.trt:canopy`, 
                p_val_int = `p_val_gm.trt:canopy`)

# Soil ammonium availability
soil.ammonium.table <- data.frame(Anova(ammonium)) %>%
  dplyr::select(chisq = Chisq, df = Df, p_val = Pr..Chisq.) %>%
  mutate(trait = "soil_ammonium",
         treatment = row.names(.),
         across(chisq:p_val, \(x) round(x, digits = 3)),
         chisq = ifelse(chisq < 0.001 & chisq >= 0, "<0.001", chisq),
         p_val = ifelse(p_val < 0.001 & p_val >= 0, "<0.001", p_val)) %>%
  pivot_wider(names_from = treatment, values_from = chisq:p_val) %>%
  dplyr::select(trait, df = df_gm.trt, chisq_gm.trt, p_val_gm.trt,
                chisq_canopy, p_val_canopy, chisq_int = `chisq_gm.trt:canopy`, 
                p_val_int = `p_val_gm.trt:canopy`)

# Soil phosphate availability
soil.phosphate.table <- data.frame(Anova(phosphate)) %>%
  dplyr::select(chisq = Chisq, df = Df, p_val = Pr..Chisq.) %>%
  mutate(trait = "soil_phosphate",
         treatment = row.names(.),
         across(chisq:p_val, \(x) round(x, digits = 3)),
         chisq = ifelse(chisq < 0.001 & chisq >= 0, "<0.001", chisq),
         p_val = ifelse(p_val < 0.001 & p_val >= 0, "<0.001", p_val)) %>%
  pivot_wider(names_from = treatment, values_from = chisq:p_val) %>%
  dplyr::select(trait, df = df_gm.trt, chisq_gm.trt, p_val_gm.trt,
                chisq_canopy, p_val_canopy, chisq_int = `chisq_gm.trt:canopy`, 
                p_val_int = `p_val_gm.trt:canopy`)

# Soil N:P
soil.np.table <- data.frame(Anova(n_to_p_ratio)) %>%
  dplyr::select(chisq = Chisq, df = Df, p_val = Pr..Chisq.) %>%
  mutate(trait = "soil_np",
         treatment = row.names(.),
         across(chisq:p_val, \(x) round(x, digits = 3)),
         chisq = ifelse(chisq < 0.001 & chisq >= 0, "<0.001", chisq),
         p_val = ifelse(p_val < 0.001 & p_val >= 0, "<0.001", p_val)) %>%
  pivot_wider(names_from = treatment, values_from = chisq:p_val) %>%
  dplyr::select(trait, df = df_gm.trt, chisq_gm.trt, p_val_gm.trt,
                chisq_canopy, p_val_canopy, chisq_int = `chisq_gm.trt:canopy`, 
                p_val_int = `p_val_gm.trt:canopy`)

# Soil moisture
soil.moisture.table <- data.frame(Anova(sm_model)) %>%
  dplyr::select(chisq = Chisq, df = Df, p_val = Pr..Chisq.) %>%
  mutate(trait = "soil_moisture",
         treatment = row.names(.),
         across(chisq:p_val, \(x) round(x, digits = 3)),
         chisq = ifelse(chisq < 0.001 & chisq >= 0, "<0.001", chisq),
         p_val = ifelse(p_val < 0.001 & p_val >= 0, "<0.001", p_val)) %>%
  pivot_wider(names_from = treatment, values_from = chisq:p_val) %>%
  dplyr::select(trait, df = df_gm.trt, chisq_gm.trt, p_val_gm.trt,
                chisq_canopy = chisq_doy, p_val_canopy = p_val_doy, chisq_int = `chisq_gm.trt:doy`, 
                p_val_int = `p_val_gm.trt:doy`)

# Write Table 3
table3 <- rbind(soil.nitrogen.table, soil.nitrate.table,
                soil.ammonium.table, soil.phosphate.table,
                soil.np.table, soil.moisture.table)
write.csv(table3, "../tables/TT23_table3.csv", row.names = FALSE)

##############################################################################
## Write Table 4: Gas exchange
##############################################################################

# Net photosynthesis (Trillium)
anet.tri.table <- data.frame(Anova(anet.tri)) %>%
  mutate(spp = "Tri",
         treatment = row.names(.),
         chisq_anet = Chisq,
         p_anet = Pr..Chisq.,
         across(chisq_anet:p_anet, \(x) round(x, digits = 3)),
         chisq_anet = ifelse(chisq_anet < 0.001 & chisq_anet >= 0, 
                              "<0.001", chisq_anet),
         p_anet = ifelse(p_anet < 0.001 & p_anet >= 0, 
                          "<0.001", p_anet)) %>%
  dplyr::select(treatment, spp, Df, chisq_anet, p_anet)

# Net photosynthesis (Maianthemum)
anet.mai.table <- data.frame(Anova(anet.mai)) %>%
  mutate(spp = "Mai",
         treatment = row.names(.),
         chisq_anet = Chisq,
         p_anet = Pr..Chisq.,
         across(chisq_anet:p_anet, \(x) round(x, digits = 3)),
         chisq_anet = ifelse(chisq_anet < 0.001 & chisq_anet >= 0, 
                                 "<0.001", chisq_anet),
         p_anet = ifelse(p_anet < 0.001 & p_anet >= 0, 
                             "<0.001", p_anet)) %>%
  dplyr::select(treatment, spp, Df, chisq_anet, p_anet)

# Stomatal conductance (Trillium)
gsw.tri.table <- data.frame(Anova(gsw.tri)) %>%
  mutate(spp = "Tri",
         treatment = row.names(.),
         chisq_gsw = Chisq,
         p_gsw = Pr..Chisq.,
         across(chisq_gsw:p_gsw, \(x) round(x, digits = 3)),
         chisq_gsw = ifelse(chisq_gsw < 0.001 & chisq_gsw >= 0, 
                             "<0.001", chisq_gsw),
         p_gsw = ifelse(p_gsw < 0.001 & p_gsw >= 0, 
                         "<0.001", p_gsw)) %>%
  dplyr::select(treatment, spp, chisq_gsw, p_gsw)

# Stomatal conductance (Maianthemum)
gsw.mai.table <- data.frame(Anova(gsw.mai)) %>%
  mutate(spp = "Mai",
         treatment = row.names(.),
         chisq_gsw = Chisq,
         p_gsw = Pr..Chisq.,
         across(chisq_gsw:p_gsw, \(x) round(x, digits = 3)),
         chisq_gsw = ifelse(chisq_gsw < 0.001 & chisq_gsw >= 0, 
                            "<0.001", chisq_gsw),
         p_gsw = ifelse(p_gsw < 0.001 & p_gsw >= 0, 
                        "<0.001", p_gsw)) %>%
  dplyr::select(treatment, spp, chisq_gsw, p_gsw)

# Stomatal limitation (Trillium)
l.tri.table <- data.frame(Anova(l.tri)) %>%
  mutate(spp = "Tri",
         treatment = row.names(.),
         chisq_l = Chisq,
         p_l = Pr..Chisq.,
         across(chisq_l:p_l, \(x) round(x, digits = 3)),
         chisq_l = ifelse(chisq_l < 0.001 & 
                                chisq_l >= 0, 
                              "<0.001", chisq_l),
         p_l = ifelse(p_l <0.001 & p_l >= 0, 
                          "<0.001", p_l)) %>%
  dplyr::select(treatment, spp, chisq_l, p_l)

# Stomatal limitation (Maianthemum)
l.mai.table <- data.frame(Anova(l.mai)) %>%
  mutate(spp = "Mai",
         treatment = row.names(.),
         chisq_l = Chisq,
         p_l = Pr..Chisq.,
         across(chisq_l:p_l, \(x) round(x, digits = 3)),
         chisq_l = ifelse(chisq_l < 0.001 & 
                            chisq_l >= 0, 
                          "<0.001", chisq_l),
         p_l = ifelse(p_l <0.001 & p_l >= 0, 
                      "<0.001", p_l)) %>%
  dplyr::select(treatment, spp, chisq_l, p_l)

# SPAD (Trillium)
spad.tri.table <- data.frame(Anova(spad.tri)) %>%
  mutate(spp = "Tri",
         treatment = row.names(.),
         chisq_spad = Chisq,
         p_spad = Pr..Chisq.,
         across(chisq_spad:p_spad, \(x) round(x, digits = 3)),
         chisq_spad = ifelse(chisq_spad < 0.001 & chisq_spad >= 0, 
                                 "<0.001", chisq_spad),
         p_spad = ifelse(p_spad <0.001 & p_spad >= 0, 
                             "<0.001", p_spad)) %>%
  dplyr::select(treatment, spp, chisq_spad, p_spad)

# SPAD (Maianthemum)
spad.mai.table <- data.frame(Anova(spad.mai)) %>%
  mutate(spp = "Mai",
         treatment = row.names(.),
         chisq_spad = Chisq,
         p_spad = Pr..Chisq.,
         across(chisq_spad:p_spad, \(x) round(x, digits = 3)),
         chisq_spad = ifelse(chisq_spad < 0.001 & chisq_spad >= 0, 
                             "<0.001", chisq_spad),
         p_spad = ifelse(p_spad <0.001 & p_spad >= 0, 
                         "<0.001", p_spad)) %>%
  dplyr::select(treatment, spp, chisq_spad, p_spad)

# Create Table 3
anet_rbind <- rbind(anet.tri.table, anet.mai.table)
gsw_rbind <- rbind(gsw.tri.table, gsw.mai.table)
l_rbind <- rbind(l.tri.table, l.mai.table)
spad_rbind <- rbind(spad.tri.table, spad.mai.table)

table4 <- anet_rbind %>% full_join(gsw_rbind) %>%
  full_join(l_rbind) %>% full_join(spad_rbind)
# write.csv(table4, "../tables/TT23_table4.csv", row.names = FALSE)

##############################################################################
## Write Table 45: Indices of photosynthetic capacity
##############################################################################

# Temp. standardized maximum rate of Rubisco carboyxlation (Trillium)
vcmax.tri <- data.frame(Anova(vcmax25.tri)) %>%
  mutate(spp = "Tri",
         treatment = row.names(.),
         chisq_vcmax = Chisq,
         p_vcmax = Pr..Chisq.,
         across(chisq_vcmax:p_vcmax, \(x) round(x, digits = 3)),
         chisq_vcmax = ifelse(chisq_vcmax < 0.001 & chisq_vcmax >= 0, 
                                 "<0.001", chisq_vcmax),
         p_vcmax = ifelse(p_vcmax < 0.001 & p_vcmax >= 0, 
                             "<0.001", p_vcmax)) %>%
  dplyr::select(treatment, spp, Df, chisq_vcmax, p_vcmax)

# Temp. standardized maximum rate of Rubisco carboxylation (Maianthemum)
vcmax.mai <- data.frame(Anova(vcmax25.mai)) %>%
  mutate(spp = "Mai",
         treatment = row.names(.),
         chisq_vcmax = Chisq,
         p_vcmax = Pr..Chisq.,
         across(chisq_vcmax:p_vcmax, \(x) round(x, digits = 3)),
         chisq_vcmax = ifelse(chisq_vcmax < 0.001 & chisq_vcmax >= 0, 
                              "<0.001", chisq_vcmax),
         p_vcmax = ifelse(p_vcmax < 0.001 & p_vcmax >= 0, 
                          "<0.001", p_vcmax)) %>%
  dplyr::select(treatment, spp, Df, chisq_vcmax, p_vcmax)

# Temp. standardized maximum rate of electron transport for RuBP 
# regeneration (Trillium)
jmax.tri <- data.frame(Anova(jmax25.tri)) %>%
  mutate(spp = "Tri",
         treatment = row.names(.),
         chisq_jmax = Chisq,
         p_jmax = Pr..Chisq.,
         across(chisq_jmax:p_jmax, \(x) round(x, digits = 3)),
         chisq_jmax = ifelse(chisq_jmax < 0.001 & chisq_jmax >= 0, 
                              "<0.001", chisq_jmax),
         p_jmax = ifelse(p_jmax < 0.001 & p_jmax >= 0, 
                          "<0.001", p_jmax)) %>%
  dplyr::select(treatment, spp, chisq_jmax, p_jmax)

# Temp. standardized maximum rate of electron transport for RuBP 
# regeneration (Maianthemum)
jmax.mai <- data.frame(Anova(jmax25.mai)) %>%
  mutate(spp = "Mai",
         treatment = row.names(.),
         chisq_jmax = Chisq,
         p_jmax = Pr..Chisq.,
         across(chisq_jmax:p_jmax, \(x) round(x, digits = 3)),
         chisq_jmax = ifelse(chisq_jmax < 0.001 & chisq_jmax >= 0, 
                             "<0.001", chisq_jmax),
         p_jmax = ifelse(p_jmax < 0.001 & p_jmax >= 0, 
                         "<0.001", p_jmax)) %>%
  dplyr::select(treatment, spp, chisq_jmax, p_jmax)

# The ratio of Jmax25 to Vcmax25 (Trillium) 
jvmax.tri <- data.frame(Anova(jmax25_vcmax25.tri)) %>%
  mutate(spp = "Tri",
         treatment = row.names(.),
         chisq_jvmax = Chisq,
         p_jvmax = Pr..Chisq.,
         across(chisq_jvmax:p_jvmax, \(x) round(x, digits = 3)),
         chisq_jvmax = ifelse(chisq_jvmax < 0.001 & chisq_jvmax >= 0, 
                             "<0.001", chisq_jvmax),
         p_jvmax = ifelse(p_jvmax < 0.001 & p_jvmax >= 0, 
                         "<0.001", p_jvmax)) %>%
  dplyr::select(treatment, spp, chisq_jvmax, p_jvmax)

# The ratio of Jmax25 to Vcmax25 (Maianthemum)
jvmax.mai <- data.frame(Anova(jmax25_vcmax25.mai)) %>%
  mutate(spp = "Mai",
         treatment = row.names(.),
         chisq_jvmax = Chisq,
         p_jvmax = Pr..Chisq.,
         across(chisq_jvmax:p_jvmax, \(x) round(x, digits = 3)),
         chisq_jvmax = ifelse(chisq_jvmax < 0.001 & chisq_jvmax >= 0, 
                              "<0.001", chisq_jvmax),
         p_jvmax = ifelse(p_jvmax < 0.001 & p_jvmax >= 0, 
                          "<0.001", p_jvmax)) %>%
  dplyr::select(treatment, spp, chisq_jvmax, p_jvmax)

# Concatenate Table 5
vcmax_rbind <- rbind(vcmax.tri, vcmax.mai)
jmax_rbind <- rbind(jmax.tri, jmax.mai)
jvmax_rbind <- rbind(jvmax.tri, jvmax.mai)

table5 <- vcmax_rbind %>%
  full_join(jmax_rbind) %>% full_join(jvmax_rbind)
# write.csv(table5, "../tables/TT23_table4.csv", row.names = FALSE)
