setwd("/Users/katherinehulting/Documents/MSU/SRS pollination founder plant/review-SRS-2009_founder")

# loading packages
library(tidyverse)
library(glmmTMB)
library(car)
library(emmeans)
library(DHARMa)
library(ggeffects)
library(performance)
library(scales)
library(svglite)
library(ggpubr)
 
# load data - excludes center patch (not used in analyses) 
pollination_2009 <- read.csv("founder_plant_2009.csv") # joined seed count data (structure level) with seed predation/flowering structure count data (plant level)
# 1 row = 1 flowering structure -- viable/nonviable seed counts specific for each structure, other data duplicated for each structure on the same plant (live, reproductive, seed predation, etc.)

pollination_2009 <- pollination_2009 %>%
  filter(patch != "A") # removing patch A (center patch) from analyses

pollination_2009 <- pollination_2009 %>%
  mutate(plant_size = length * width * height) %>%
  mutate(log_plant_size = log(plant_size)) %>%
  filter(!is.na(reproductive)) %>% # removing dead plants and plants with no repro status recorded
  filter(flags != 3) %>% # removing repro plants with no flowering structures recorded
  mutate(no_axil_fl = ifelse(is.na(no_axil_fl), 0, no_axil_fl)) %>% # NAs for AV, AB, SS - no axilary flowering stalks on these species
  mutate(no_basal_fl = ifelse(is.na(no_basal_fl), 0, no_basal_fl)) %>% # NAs = no basal flowering stalks on these species
  mutate(no.structures = no_basal_fl + no_axil_fl)

# subsetting by species - keeping all plants, NA for reproduction = dead plant, excluded from analysis
AB_reproduction <- pollination_2009 %>%
  filter(species == "Aristida") %>%
  distinct(plant_ID, .keep_all = TRUE) # only keeping one observation/plant for reproduction analysis-- there are three observations/plant because of 3 recorded structures/plant

AV_reproduction <- pollination_2009 %>%
  filter(species == "Anthaenantia") %>%
  distinct(plant_ID, .keep_all = TRUE)

CB_reproduction <- pollination_2009 %>%
  filter(species == "Carphephorus") %>%
  distinct(plant_ID, .keep_all = TRUE)

LE_reproduction <- pollination_2009 %>%
  filter(species == "Liatris") %>%
  distinct(plant_ID, .keep_all = TRUE)

SS_reproduction <- pollination_2009 %>%
  filter(species == "Sorghastrum") %>%
  distinct(plant_ID, .keep_all = TRUE)

###### reproductive status models #####
# Aristida beyrichiana
AB_repr_mod <- glmmTMB(no.structures ~ ptype + dist_num + log_plant_size + (1|block/patch/corner),
                       ziformula = ~ ptype + dist_num + log_plant_size + (1|block/patch/corner), 
                       data = AB_reproduction,
                       family = "truncated_nbinom2")
summary(AB_repr_mod)
plot(simulateResiduals(AB_repr_mod))
# Anova
Anova(AB_repr_mod, type = "III", component = "cond")
Anova(AB_repr_mod, type = "III", component = "zi")
# Pairwise comparisons
AB_repr_mod.posthoc <- emmeans(AB_repr_mod, ~ptype, comp = "cond")
pairs(AB_repr_mod.posthoc)
AB_repr_mod.posthoc <- emmeans(AB_repr_mod, ~ptype, comp = "zi")
pairs(AB_repr_mod.posthoc)

# Anthaenantia villosa
AV_repr_mod <- glmmTMB(no.structures ~ ptype + dist_num + log_plant_size + (1|block/patch/corner),
                       ziformula = ~ ptype + dist_num + log_plant_size + (1|block/patch), 
                       data = AV_reproduction,
                       family = "truncated_nbinom2")
summary(AV_repr_mod)
plot(simulateResiduals(AV_repr_mod))
Anova(AV_repr_mod, type = "III", component = "cond")
Anova(AV_repr_mod, type = "III", component = "zi")
# Pairwise comparisons
AV_repr_mod.posthoc <- emmeans(AV_repr_mod, ~ptype, comp = "cond")
pairs(AV_repr_mod.posthoc)
AV_repr_mod.posthoc <- emmeans(AV_repr_mod, ~ptype, comp = "zi")
pairs(AV_repr_mod.posthoc)

# Carphephorus bellidifolius 
CB_repr_mod <- glmmTMB(no.structures ~ ptype + dist_num + log_plant_size + (1|block/patch/corner),
                       ziformula = ~ ptype + dist_num + log_plant_size + (1|block/patch/corner), 
                       data = CB_reproduction,
                       family = "truncated_nbinom2")
summary(CB_repr_mod)
plot(simulateResiduals(CB_repr_mod))
Anova(CB_repr_mod, type = "III", component = "cond")
Anova(CB_repr_mod, type = "III", component = "zi")
# Pairwise comparisons
CB_repr_mod.posthoc <- emmeans(CB_repr_mod, ~ptype, comp = "cond")
pairs(CB_repr_mod.posthoc)
CB_repr_mod.posthoc <- emmeans(CB_repr_mod, ~ptype, comp = "zi")
pairs(CB_repr_mod.posthoc)

# Liatris earlei 
LE_repr_mod <- glmmTMB(no.structures ~ ptype + dist_num + log_plant_size + (1|block/patch/corner),
                       ziformula = ~ ptype + dist_num + log_plant_size + (1|block/patch/corner), 
                       data = LE_reproduction,
                       family = "truncated_nbinom2")
summary(LE_repr_mod)
plot(simulateResiduals(LE_repr_mod))
Anova(LE_repr_mod, type = "III", component = "cond")
Anova(LE_repr_mod, type = "III", component = "zi")
# Pairwise comparisons
LE_repr_mod.posthoc <- emmeans(LE_repr_mod, ~ptype, comp = "cond")
pairs(LE_repr_mod.posthoc)
LE_repr_mod.posthoc <- emmeans(LE_repr_mod, ~ptype, comp = "zi")
pairs(LE_repr_mod.posthoc)

#Sorghastrum secundum
SS_repr_mod <- glmmTMB(no.structures ~ ptype + dist_num + log_plant_size + (1|block/patch/corner),
                       ziformula = ~ ptype + dist_num + log_plant_size + (1|block/patch/corner), 
                       data = SS_reproduction,
                       family = "truncated_nbinom2")
summary(SS_repr_mod)
plot(simulateResiduals(SS_repr_mod))
Anova(SS_repr_mod, type = "III", component = "cond")
Anova(SS_repr_mod, type = "III", component = "zi")
# Pairwise comparisons
SS_repr_mod.posthoc <- emmeans(SS_repr_mod, ~ptype, comp = "cond")
pairs(SS_repr_mod.posthoc)
SS_repr_mod.posthoc <- emmeans(SS_repr_mod, ~ptype, comp = "zi")
pairs(SS_repr_mod.posthoc)


###### reproductive status plotting: Figure 2a #####
# model predictions for plotting - zero-count process
predictAB_repr_mod <- ggpredict(AB_repr_mod, terms=c("dist_num", "ptype"), type = "zi_prob", back.transform = T, allow.new.levels=TRUE)
predictAB_repr_mod <- predictAB_repr_mod %>%
  mutate(species = "Aristida")

predictAV_repr_mod <- ggpredict(AV_repr_mod, terms=c("dist_num", "ptype"), type = "zi_prob", back.transform = T, allow.new.levels=TRUE)
predictAV_repr_mod <- predictAV_repr_mod %>%
  mutate(species = "Anthaenantia")

predictCB_repr_mod <- ggpredict(CB_repr_mod, terms=c("dist_num", "ptype"), type = "zi_prob", back.transform = T, allow.new.levels=TRUE)
predictCB_repr_mod <- predictCB_repr_mod %>%
  mutate(species = "Carphephorus")

predictLE_repr_mod <- ggpredict(LE_repr_mod, terms=c("dist_num", "ptype"), type = "zi_prob", back.transform = T, allow.new.levels=TRUE)
predictLE_repr_mod <- predictLE_repr_mod %>%
  mutate(species = "Liatris")

predictSS_repr_mod <- ggpredict(SS_repr_mod, terms=c("dist_num", "ptype"), type = "zi_prob", back.transform = T, allow.new.levels=TRUE)
predictSS_repr_mod <- predictSS_repr_mod %>%
  mutate(species = "Sorghastrum")

repr_predict_grass <- rbind(predictAB_repr_mod, predictAV_repr_mod, predictSS_repr_mod)
repr_predict_grass$plant_type = "Poaceae (wind pollinated)"
repr_predict_forb <- rbind(predictCB_repr_mod, predictLE_repr_mod)
repr_predict_forb$plant_type = "Asteraceae (insect pollinated)"


# plotting
figure2a_grass <- repr_predict_grass %>%
  ggplot(aes(x = x, y = predicted)) +
  geom_line(aes(x, predicted, linetype = group), linewidth = 0.7) + 
  theme_bw()+
  #geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
  #            alpha = 0.15) +
  labs(title = NULL,
       x = "Distance from Edge",
       y = "Likelihood of Vegetative Status") +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title = element_text(size = 18)) +
  theme(legend.position = "none") +
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1)) +
  facet_grid(rows = vars(species), cols = vars(plant_type), scales = "fixed", drop = T) +
  theme(strip.text.y = element_text(size = 16, colour = "black", angle = 270)) +
  theme(strip.text.x = element_text(size = 16, colour = "black"))
figure2a_grass

figure2a_forb <- repr_predict_forb %>%
  ggplot(aes(x = x, y = predicted)) +
  geom_line(aes(x, predicted, linetype = group), linewidth = 0.7) + 
  theme_bw()+
  #geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
  #            alpha = 0.15) +
  labs(title = NULL,
       x = "Distance from Edge",
       y = NULL) +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title = element_text(size = 18)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text =element_text(size = 14)) +
  theme(plot.margin=unit(c(0.2,1,3,0.5),"cm"))+
  theme(legend.position = c(0.2, -0.3)) +
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1)) +
  guides(linetype = guide_legend(title = "Patch Type")) +
  theme(legend.key.width=unit(2,"cm")) + 
  facet_grid(rows = vars(species), cols = vars(plant_type), scales = "free", drop = T) +
  theme(strip.text.y = element_text(size = 16, colour = "black", angle = 270)) +
  theme(strip.text.x = element_text(size = 16, colour = "black"))
figure2a_forb



left_column_repra <- cowplot::plot_grid(figure2a_forb, NULL, nrow = 2, ncol = 1, rel_heights = c(2.44,0.4))
fig2a <- cowplot::plot_grid(figure2a_grass, left_column_repra, nrow=1, ncol=2, align = "none", axis = "t", rel_widths = c(0.95,1))
fig2a

###### number of flowering structures plotting: Figure 2b #####
# model predictions for plotting - zero-count process
predictAB_repr_mod <- ggpredict(AB_repr_mod, terms=c("dist_num", "ptype"), type = "zi_random", back.transform = T, allow.new.levels=TRUE)
predictAB_repr_mod <- predictAB_repr_mod %>%
  mutate(species = "Aristida")

predictAV_repr_mod <- ggpredict(AV_repr_mod, terms=c("dist_num", "ptype"), type = "zi_random", back.transform = T, allow.new.levels=TRUE)
predictAV_repr_mod <- predictAV_repr_mod %>%
  mutate(species = "Anthaenantia")

predictCB_repr_mod <- ggpredict(CB_repr_mod, terms=c("dist_num", "ptype"), type = "zi_random", back.transform = T, allow.new.levels=TRUE)
predictCB_repr_mod <- predictCB_repr_mod %>%
  mutate(species = "Carphephorus") %>%
  select(!std.error)

predictLE_repr_mod <- ggpredict(LE_repr_mod, terms=c("dist_num", "ptype"), type = "zi_random", back.transform = T, allow.new.levels=TRUE)
predictLE_repr_mod <- predictLE_repr_mod %>%
  mutate(species = "Liatris")

predictSS_repr_mod <- ggpredict(SS_repr_mod, terms=c("dist_num", "ptype"), type = "zi_random", back.transform = T, allow.new.levels=TRUE)
predictSS_repr_mod <- predictSS_repr_mod %>%
  mutate(species = "Sorghastrum")

repr_predict_grass <- rbind(predictAB_repr_mod, predictAV_repr_mod, predictSS_repr_mod)
repr_predict_grass$plant_type = "Poaceae (wind pollinated)"
repr_predict_forb <- rbind(predictCB_repr_mod, predictLE_repr_mod)
repr_predict_forb$plant_type = "Asteraceae (insect pollinated)"


# plotting
figure2b_grass <- repr_predict_grass %>%
  ggplot(aes(x = x, y = predicted)) +
  geom_line(aes(x, predicted, linetype = group), linewidth = 0.7) + 
  theme_bw()+
 # geom_jitter(aes(dist_num, no.structures), data = SS_reproduction, width = 0.3, height = 0.2, alpha = 0.2) +
  #geom_jitter(aes(dist_num, no.structures), data = AB_reproduction, width = 0.3, height = 0.2, alpha = 0.2) +
 # geom_jitter(aes(dist_num, no.structures), data = AV_reproduction, width = 0.3, height = 0.2, alpha = 0.2) +
  labs(title = NULL,
       x = "Distance from Edge",
       y = "Number of Flowering Structures") +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title = element_text(size = 18)) +
  theme(legend.position = "none") +
  #scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1)) +
  facet_grid(rows = vars(species), cols = vars(plant_type), scales = "free", drop = T) +
  theme(strip.text.y = element_text(size = 16, colour = "black", angle = 270)) +
  theme(strip.text.x = element_text(size = 16, colour = "black"))
figure2b_grass

figure2b_forb <- repr_predict_forb %>%
  ggplot(aes(x = x, y = predicted)) +
  geom_line(aes(x, predicted, linetype = group), linewidth = 0.7) + 
  theme_bw()+
  #geom_jitter(aes(dist_num, no.structures), data = CB_reproduction, width = 0.3, height = 0.2, alpha = 0.2) +
  #geom_jitter(aes(dist_num, no.structures), data = LE_reproduction, width = 0.3, height = 0.2, alpha = 0.2) +
  #geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
  #            alpha = 0.15) +
  labs(title = NULL,
       x = "Distance from Edge",
       y = NULL) +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title = element_text(size = 18)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text =element_text(size = 14)) +
  theme(plot.margin=unit(c(0.2,1,3,0.5),"cm"))+
  theme(legend.position = c(0.2, -0.3)) +
  #scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1)) +
  guides(linetype = guide_legend(title = "Patch Type")) +
  theme(legend.key.width=unit(2,"cm")) + 
  facet_grid(rows = vars(species), cols = vars(plant_type), scales = "free", drop = T) +
  theme(strip.text.y = element_text(size = 16, colour = "black", angle = 270)) +
  theme(strip.text.x = element_text(size = 16, colour = "black"))
figure2b_forb



left_column_reprb <- cowplot::plot_grid(figure2b_forb, NULL, nrow = 2, ncol = 1, rel_heights = c(2.44,0.4))
fig2b <- cowplot::plot_grid(figure2b_grass, left_column_reprb, nrow=1, ncol=2, align = "none", axis = "t", rel_widths = c(0.95,1))
fig2b




figure2 <- cowplot::plot_grid(NULL, fig2a, NULL, fig2b, nrow=4, ncol=1, labels = c("","A","", "B"), label_size = 25,
                              label_x = 0,
                              label_y = 1.05,
                              rel_heights = c(0.05, 1, 0.05, 1))
pdf(file = "Figure2.pdf", width = 12, height = 15)
figure2
dev.off()





#### pollination data wrangling: keeping only reproductive plants ####
#count_repro_plants <- pollination_2009 %>% # calculating # of flowering structures per species in each patch
#  mutate(no_axil_fl = ifelse(is.na(no_axil_fl), 0, no_axil_fl)) %>%
#  mutate(flowering_structures = no_basal_fl + no_axil_fl) %>%
#  mutate(patch_ID = paste(block, patch, species, sep = ".")) %>%
#  group_by(patch_ID) %>%
#  mutate(count_repro_plants = sum(flowering_structures, na.rm = TRUE)) %>%
#  select(c("patch_ID", "count_repro_plants")) %>%
#  distinct(patch_ID, .keep_all = TRUE) # keeping one per patch/species

pollination_2009 <- pollination_2009 %>%
  filter(flags == 0) %>% # only using data with no flags for pollination/seed predation data
  filter(reproductive == 1) %>% # removing dead/non-reproductive plants from pollination analysis
  filter(!is.na(no.nonviable_seeds)) %>% # removing empty rows - from plants with < 3 structures collected
  filter(no.nonviable_seeds != "") %>% # removing 2 empty rows - from plants with < 3 structures collected
  filter(!is.na(no.viable_seeds)) %>% # removing 1 row with missing viable seed count
  filter(dispersed_percent < 50 | is.na(dispersed_percent)) %>%  # removing observations with over 50% seed dispersal prior to collection
  filter(!is.na(predisp_seedpred)) # removing observations with no predispersal seed predation data

pollination_2009$predisp_seedpred <- as.numeric(pollination_2009$predisp_seedpred) # converting character to numeric
pollination_2009$no.viable_seeds <- as.numeric(pollination_2009$no.viable_seeds) # converting character to numeric
pollination_2009$no.nonviable_seeds <- as.numeric(pollination_2009$no.nonviable_seeds) # converting character to numeric

#pollination_2009$predisp_seedpred[is.na(pollination_2009$predisp_seedpred)] <- 0 # replacing NAs with 0
#pollination_2009$no.viable_seeds[is.na(pollination_2009$no.viable_seeds)] <- 0 # replacing NAs with 0 - 1 observation with viable seeds = NA, but has data for non-viable seeds

pollination_2009 <- pollination_2009 %>% # calculating pollination rate
  mutate(no.pollinated_seeds = no.viable_seeds + predisp_seedpred * no.nonviable_seeds) %>% # number of pollinated seeds = # of viable seeds corrected for pre dispersal seed predation
  mutate(total_est_fl = no.viable_seeds + no.nonviable_seeds) %>% # total # of flowers = number of viable seeds + number of nonviable seeds
  mutate(pollination = no.pollinated_seeds/total_est_fl) %>% # pollination = # of pollinated flowers / total # of flowers
  filter(!is.nan(pollination)) # removing rows that had 0 viable seeds and 0 nonviable seeds - no seeds formed on structures, no estimate of pollination

pollination_2009$pollination <- as.numeric(pollination_2009$pollination) # converting character to numeric
pollination_2009$total_est_fl <- as.numeric(pollination_2009$total_est_fl) # converting character to numeric

# averaging pollination, # viable seeds, # pollinated seeds per plant for analysis - looking at plant-level, not structure-level
pollination_2009 <- pollination_2009 %>%
  group_by(plant_ID) %>% # grouping by plant
  mutate(pollination_avg = mean(pollination)) %>% # averaging pollination per plant from three structures
  mutate(viable_seed_avg = mean(no.viable_seeds)) %>% # averaging the # of viable seeds per plant from three structures
  mutate(no.pollinated_seeds_avg = mean(no.pollinated_seeds)) %>% # averaging the # of pollinated seeds per plant from three structures
  mutate(total_est_fl_avg = mean(total_est_fl)) %>% # averaging the # of flowers per structure per plant
  ungroup() %>% # ungrouping by plant
  distinct(plant_ID, .keep_all = TRUE) %>% # keeping all columns, one row per plant
  dplyr::select(!c("pollination", "total_est_fl", "no.pollinated_seeds", "no.viable_seeds", "no.nonviable_seeds")) # removing non-summarized data

pollination_2009$no_axil_fl <- as.numeric(pollination_2009$no_axil_fl) # converting character to numeric
pollination_2009$no_basal_fl <- as.numeric(pollination_2009$no_basal_fl) # converting character to numeric

pollination_2009 <- pollination_2009 %>%
  mutate(no_axil_fl = ifelse(is.na(no_axil_fl), 0, no_axil_fl)) %>% # NAs for AV, AB, SS - no axilary flowering stalks on these species
  mutate(total_no_structures = no_basal_fl + no_axil_fl) %>% # calculating total # of structures on a plant = basal flowering stalks + axilary flowering stalks
  mutate(plant_seed_prod = total_no_structures * no.pollinated_seeds_avg) %>% # plant seed production = # of structures on plant * # of pollinated seeds per structure
  mutate(plant_seed_prod = round(plant_seed_prod)) %>% # rounding for analysis
  mutate(total_est_fl_avg = round(total_est_fl_avg))  # rounding for analysis

#pollination_2009 <- pollination_2009 %>%
#  mutate(patch_ID = paste(block, patch, species, sep = ".")) %>%
#  left_join(count_repro_plants, by = "patch_ID")

# separating out species into separate dataframes for analysis
AB_pollination <- pollination_2009 %>%
  filter(species == "Aristida")

AV_pollination <- pollination_2009 %>%
  filter(species == "Anthaenantia")

CB_pollination <- pollination_2009 %>%
  filter(species == "Carphephorus")

LE_pollination <- pollination_2009 %>%
  filter(species == "Liatris")

SS_pollination <- pollination_2009 %>%
  filter(species == "Sorghastrum")


mean(AB_pollination$pollination_avg)
mean(AV_pollination$pollination_avg)
mean(CB_pollination$pollination_avg)
mean(LE_pollination$pollination_avg)
mean(SS_pollination$pollination_avg)


###### pollination models ######
# Aristida beyrichiana
AB_mod <- glmmTMB(pollination_avg ~ ptype + dist_num + log_plant_size + (1|block/patch/corner), 
                  data = AB_pollination,
                  family = betabinomial(),
                  weights = total_est_fl_avg)
summary(AB_mod)
plot(simulateResiduals(AB_mod)) # not the best
Anova(AB_mod, type = "III")
AB_mod.posthoc <- emmeans(AB_mod, "ptype")
pairs(AB_mod.posthoc)


# Anthaenantia villosa
AV_mod <- glmmTMB(pollination_avg ~ ptype + dist_num + log_plant_size + (1|block/patch/corner), 
                  data = AV_pollination,
                  family = betabinomial(),
                  weights = total_est_fl_avg)
summary(AV_mod)
plot(simulateResiduals(AV_mod))
Anova(AV_mod, type = "III")
AB_mod.posthoc <- emmeans(AB_mod, "ptype")
pairs(AB_mod.posthoc)


# Carphephorus bellidifolius 
CB_mod <- glmmTMB(pollination_avg ~ ptype + dist_num + log_plant_size + (1|block/patch/corner), 
                  data = CB_pollination,
                  family = betabinomial(),
                  weights = total_est_fl_avg)
summary(CB_mod)
Anova(CB_mod, type = "III")
plot(simulateResiduals(CB_mod))
CB_mod.posthoc <- emmeans(CB_mod, "ptype")
pairs(CB_mod.posthoc)


# Liatris earlei 
LE_mod <- glmmTMB(pollination_avg ~ ptype + dist_num + log_plant_size + (1|block/patch/corner), 
                  data = LE_pollination,
                  family = betabinomial(),
                  weights = total_est_fl_avg)
summary(LE_mod)
Anova(LE_mod, type = "III")
plot(simulateResiduals(LE_mod))
LE_mod.posthoc <- emmeans(LE_mod, "ptype")
pairs(LE_mod.posthoc)


# Sorghastrum secundum
SS_mod <- glmmTMB(pollination_avg ~ ptype + dist_num + log_plant_size + (1|block/patch/corner), 
                  data = SS_pollination,
                  family = betabinomial(),
                  weights = total_est_fl_avg)
summary(SS_mod)
plot(simulateResiduals(SS_mod))
Anova(SS_mod, type = "III")
SS_mod.posthoc <- emmeans(SS_mod, "ptype")
pairs(SS_mod.posthoc)


#### pollination plotting: Figure S1 #######
predictAB_rate_mod <- ggpredict(AB_mod, terms=c("dist_num", "ptype"), back.transform = T, allow.new.levels=TRUE)
predictAB_rate_mod <- predictAB_rate_mod %>%
  mutate(species = "Aristida")

predictAV_rate_mod <- ggpredict(AV_mod, terms=c("dist_num", "ptype"), back.transform = T, allow.new.levels=TRUE)
predictAV_rate_mod <- predictAV_rate_mod %>%
  mutate(species = "Anthaenantia")

predictCB_rate_mod <- ggpredict(CB_mod, terms=c("dist_num", "ptype"), back.transform = T, allow.new.levels=TRUE)
predictCB_rate_mod <- predictCB_rate_mod %>%
  mutate(species = "Carphephorus")

predictLE_rate_mod <- ggpredict(LE_mod, terms=c("dist_num", "ptype"), back.transform = T, allow.new.levels=TRUE)
predictLE_rate_mod <- predictLE_rate_mod %>%
  mutate(species = "Liatris")

predictSS_rate_mod <- ggpredict(SS_mod, terms=c("dist_num", "ptype"), back.transform = T, allow.new.levels=TRUE)
predictSS_rate_mod <- predictSS_rate_mod %>%
  mutate(species = "Sorghastrum")


rate_predict_grass <- rbind(predictAB_rate_mod, predictAV_rate_mod, predictSS_rate_mod)
rate_predict_grass$plant_type = "Poaceae (wind pollinated)"
rate_predict_forb <- rbind(predictCB_rate_mod, predictLE_rate_mod)
rate_predict_forb$plant_type = "Asteraceae (insect pollinated)"


figureS1_grass <- rate_predict_grass %>%
  ggplot(aes(x = x, y = predicted)) +
  geom_line(aes(x, predicted, linetype = group)) +
  geom_point(aes(x = dist_num, y = pollination_avg), data = AB_pollination, alpha = 0.2) +
  geom_point(aes(x = dist_num, y = pollination_avg), data = SS_pollination, alpha = 0.2) +
  geom_point(aes(x = dist_num, y = pollination_avg), data = AV_pollination, alpha = 0.2) +
  theme_bw()+
  labs(title = NULL,
       x = "Distance from Edge",
       y = "Pollination rate") +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title = element_text(size = 18)) +
  theme(legend.position = "none") +
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1)) +
  facet_grid(rows = vars(species), cols = vars(plant_type), scales = "free", drop = T) +
  theme(strip.text.y = element_text(size = 16, colour = "black", angle = 270)) +
  theme(strip.text.x = element_text(size = 16, colour = "black"))

figureS1_forb <- rate_predict_forb %>%
  ggplot(aes(x = x, y = predicted)) +
  geom_line(aes(x, predicted, linetype = group)) +
  geom_point(aes(x = dist_num, y = pollination_avg), data = CB_pollination, alpha = 0.2) +
  geom_point(aes(x = dist_num, y = pollination_avg), data = LE_pollination, alpha = 0.2) +
  theme_bw()+
  labs(title = NULL,
       x = "Distance from Edge",
       y = NULL) +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title = element_text(size = 18)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text =element_text(size = 14)) +
  theme(plot.margin=unit(c(0.2,1,3,0.5),"cm"))+
  theme(legend.position = c(0.2, -0.3)) +
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1)) +
  guides(linetype = guide_legend(title = "Patch Type")) +
  theme(legend.key.width=unit(2,"cm")) +
  facet_grid(rows = vars(species), cols = vars(plant_type), scales = "free", drop = T) +
  theme(strip.text.y = element_text(size = 16, colour = "black", angle = 270)) +
  theme(strip.text.x = element_text(size = 16, colour = "black"))

left_column_rate <- cowplot::plot_grid(figureS1_forb, NULL, nrow = 2, ncol = 1, rel_heights = c(2.44,0.4))
pdf(file = "FigureS1.pdf", width = 12, height = 8)
cowplot::plot_grid(figureS1_grass, left_column_rate, nrow=1, ncol=2, align = "none", axis = "t", rel_widths = c(1.05,1))
dev.off()



#### seed production ####
# Aristida beyrichiana
#AB_pollination_no.zero <- AB_pollination %>%
#  filter(plant_seed_prod != 0)
#AB_mod1 <- glmmTMB(plant_seed_prod ~ ptype + dist_num + (1|block/patch/corner),
#                   data = AB_pollination_no.zero,
#                   family = "nbinom2")

AB_mod1 <- glmmTMB(plant_seed_prod ~ ptype + dist_num + log_plant_size + (1|block/patch/corner), 
                   ziformula = ~ ptype + dist_num + log_plant_size + (1|block/patch/corner),
                   data = AB_pollination,
                   family = "truncated_nbinom2")
summary(AB_mod1)
plot(simulateResiduals(AB_mod1))
#check_overdispersion(AB_mod1)
#check_zeroinflation(AB_mod1)
Anova(AB_mod1, type = "III", component = "cond")
Anova(AB_mod1, type = "III", component = "zi")
# Post hoc
AB_mod1.posthoc <- emmeans(AB_mod1, ~ ptype, comp = "cond") # posthoc
pairs(AB_mod1.posthoc)
AB_mod1.posthoc.zero <- emmeans(AB_mod1, ~ ptype, comp = "zi") # posthoc zero count process
pairs(AB_mod1.posthoc.zero)



# Anthaenantia villosa
AV_mod1 <- glmmTMB(plant_seed_prod ~ ptype + dist_num + log_plant_size + (1|block/patch/corner), 
                          ziformula = ~ ptype + dist_num + (1|block/patch/corner),
                          data = AV_pollination,
                          family = "truncated_nbinom2")
summary(AV_mod1)
plot(simulateResiduals(AV_mod1))
#check_overdispersion(AV_mod1)
#check_zeroinflation(AV_mod1) 
Anova(AV_mod1, type = "III", component = "cond")
Anova(AV_mod1, type = "III", component = "zi")
(exp(-0.4273)-1) *100

## plotting Anthaenantia patch type differences
AV_mod1.stat.test <- tibble::tribble( # adding p-values above graph
  ~group1, ~group2,   ~p.adj,
  "Connected",     "Rectangle", "0.12",
  "Connected",     "Winged", "0.99",
  "Rectangle",     "Winged", "0.07"
)
predictAV_mod1 <- ggpredict(AV_mod1, terms=c("ptype [all]"), back.transform = T)
plotAV_mod1 <- predictAV_mod1 %>% ggplot(aes(x = x, y = predicted)) +
  geom_jitter(aes(x = ptype, y = plant_seed_prod), data = AV_pollination, alpha = 0.2, width = 0.1, height = 0.15, size = 3)+ 
  geom_point(size = 2)+ 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, linewidth = 1.5) +
  theme_minimal()+
  stat_pvalue_manual(
    AV_mod1.stat.test, 
    size = 4,
    bracket.size = 0.7,
    y.position = 10000, step.increase = 0.1,
    label = "p.adj"
  ) +
  labs(title = NULL,
       x = "Patch Type",
       y = "Anthaenantia Seed Production") +
  theme(axis.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 22))
  #theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1.5))
plotAV_mod1
# exporting
pdf(file = "FigureS2.pdf", width = 10, height = 7)
plotAV_mod1
dev.off()


# Carphephorus bellidifolius 
CB_mod1 <- glmmTMB(plant_seed_prod ~ ptype + dist_num + log_plant_size + (1|block/patch/corner), 
                   ziformula = ~ ptype + dist_num + log_plant_size + (1|block/patch/corner),
                   data = CB_pollination,
                   family = "truncated_nbinom2")
summary(CB_mod1)
plot(simulateResiduals(CB_mod1))
#check_overdispersion(CB_mod1)
#check_zeroinflation(CB_mod1) 
Anova(CB_mod1, type = "III", component = "cond")
Anova(CB_mod1, type = "III", component = "zi")
CB_mod1.posthoc <- emmeans(CB_mod1, "ptype")
pairs(CB_mod1.posthoc)


# Liatris earlei 
LE_mod1 <- glmmTMB(plant_seed_prod ~ ptype + dist_num + log_plant_size + (1|block/patch/corner), 
                   ziformula = ~ ptype + dist_num + (1|block/patch/corner),
                   data = LE_pollination,
                   family = "truncated_nbinom2")
summary(LE_mod1)
plot(simulateResiduals(LE_mod1))
#check_overdispersion(LE_mod1)
#check_zeroinflation(LE_mod1)
Anova(LE_mod1, type = "III", component = "cond")
Anova(LE_mod1, type = "III", component = "zi")
LE_mod1.posthoc <- emmeans(LE_mod1, "ptype")
pairs(LE_mod1.posthoc)


# Sorghastrum secundum
SS_mod1 <- glmmTMB(plant_seed_prod ~ ptype + dist_num + log_plant_size + (1|block/patch/corner), 
                   ziformula = ~ ptype + dist_num + log_plant_size + (1|block/patch/corner),
                   data = SS_pollination,
                   family = "truncated_nbinom2")
summary(SS_mod1)
plot(simulateResiduals(SS_mod1))
#check_overdispersion(SS_mod1)
#check_zeroinflation(SS_mod1)
Anova(SS_mod1, type = "III", component = "cond")
Anova(SS_mod1, type = "III", component = "zi")
SS_mod1.posthoc <- emmeans(SS_mod1, "ptype")
pairs(SS_mod1.posthoc)


###### seed production plotting: Figure 3 ######
predictAB_mod1 <- ggpredict(AB_mod1, terms=c("dist_num"), type = "fe", back.transform = T, allow.new.levels=TRUE)
predictAB_mod1 <- predictAB_mod1 %>%
  mutate(species = "Aristida")

predictAV_mod1 <- ggpredict(AV_mod1, terms=c("dist_num"), type = "fe", back.transform = T, allow.new.levels=TRUE)
predictAV_mod1 <- predictAV_mod1 %>%
  mutate(species = "Anthaenantia")

predictCB_mod1 <- ggpredict(CB_mod1, terms=c("dist_num"), type = "fe", back.transform = T, allow.new.levels=TRUE)
predictCB_mod1 <- predictCB_mod1 %>%
  mutate(species = "Carphephorus")

predictLE_mod1 <- ggpredict(LE_mod1, terms=c("dist_num"), type = "fe", back.transform = T, allow.new.levels=TRUE)
predictLE_mod1 <- predictLE_mod1 %>%
  mutate(species = "Liatris")

predictSS_mod1 <- ggpredict(SS_mod1, terms=c("dist_num"), type = "fe", back.transform = T, allow.new.levels=TRUE)
predictSS_mod1 <- predictSS_mod1 %>%
  mutate(species = "Sorghastrum")

mod1_predict_grass <- rbind(predictAB_mod1, predictAV_mod1, predictSS_mod1)
mod1_predict_grass$plant_type = "Poaceae (wind pollinated)"
mod1_predict_forb <- rbind(predictCB_mod1, predictLE_mod1)
mod1_predict_forb$plant_type = "Asteraceae (insect pollinated)"



SS_non_zero <- SS_pollination %>%
  filter(plant_seed_prod != 0)
AB_non_zero <- AB_pollination %>%
  filter(plant_seed_prod != 0)
AV_non_zero <- AV_pollination %>%
  filter(plant_seed_prod != 0)
CB_non_zero <- CB_pollination %>%
  filter(plant_seed_prod != 0)
LE_non_zero <- LE_pollination %>%
  filter(plant_seed_prod != 0)

figure3_grass <- mod1_predict_grass %>% ggplot(aes(x = x, y = predicted)) +
  geom_jitter(aes(x = dist_num, y = plant_seed_prod), data = SS_non_zero, alpha = 0.3, width = 0.4, height = 0.2) + 
  geom_jitter(aes(x = dist_num, y = plant_seed_prod), data = AB_non_zero, alpha = 0.3, width = 0.4, height = 0.2) + 
  geom_jitter(aes(x = dist_num, y = plant_seed_prod), data = AV_non_zero, alpha = 0.3, width = 0.4, height = 0.2) + 
  geom_line(aes(x, predicted, linetype = group)) +
  geom_line(aes(x, predicted), color = "red", linewidth = 1.5) +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2) +
  theme_bw() +
  labs(title = NULL,
       x = "Distance from Edge",
       y = "Seed Production") +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title = element_text(size = 18)) +
  theme(legend.position = "none") +
  facet_grid(rows = vars(species), cols = vars(plant_type), scales = "free", drop = T) +
  theme(strip.text.y = element_text(size = 16, colour = "black", angle = 270)) +
  theme(strip.text.x = element_text(size = 16, colour = "black"))
figure3_grass

figure3_forb <- mod1_predict_forb %>% ggplot(aes(x = x, y = predicted)) +
  geom_jitter(aes(x = dist_num, y = plant_seed_prod), data = CB_non_zero, alpha = 0.3, width = 0.4, height = 0.2) + 
  geom_jitter(aes(x = dist_num, y = plant_seed_prod), data = LE_non_zero, alpha = 0.3, width = 0.4, height = 0.2) + 
  geom_line(aes(x, predicted, linetype = group)) +
  geom_line(aes(x, predicted), color = "red", linewidth = 1.5) +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2) +
  theme_bw() +
  labs(title = NULL,
       x = "Distance from Edge",
       y = NULL) +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title = element_text(size = 18)) +
  theme(legend.position = "none") +
  facet_grid(rows = vars(species), cols = vars(plant_type), scales = "free", drop = T) +
  theme(strip.text.y = element_text(size = 16, colour = "black", angle = 270)) +
  theme(strip.text.x = element_text(size = 16, colour = "black"))
figure3_forb

left_column <- cowplot::plot_grid(figure3_forb, NULL, nrow = 2, ncol = 1, rel_heights = c(2.44,1))
pdf(file = "Figure3.pdf", width = 12, height = 8)
cowplot::plot_grid(figure3_grass, left_column, nrow=1, ncol=2, align = "none", axis = "t", rel_widths = c(1.05,1))
dev.off()




####### Plant size ######
# testing if plant size of reproductive plants is affected by edge or patch type
AB_size <- glmmTMB(log_plant_size ~ ptype * dist_num + (1|block/patch/corner),
                   data = AB_reproduction,
                   family = "gaussian")
summary(AB_size)
plot(simulateResiduals(AB_size))


AV_size <- glmmTMB(log_plant_size ~ ptype + dist_num + (1|block/patch/corner),
                   data = AV_reproduction,
                   family = "gaussian")
summary(AV_size)
plot(simulateResiduals(AV_size))

CB_size <- glmmTMB(log_plant_size ~ ptype + dist_num + (1|block/patch/corner),
                   data = CB_reproduction,
                   family = "gaussian")
summary(CB_size)
plot(simulateResiduals(CB_size))

LE_size <- glmmTMB(log_plant_size ~ ptype + dist_num + (1|block/patch/corner),
                   data = LE_reproduction,
                   family = "gaussian")
summary(LE_size)
plot(simulateResiduals(LE_size))

SS_size <- glmmTMB(log_plant_size ~ ptype + dist_num + (1|block/patch/corner),
                   data = SS_reproduction,
                   family = "gaussian")
summary(SS_size)
plot(simulateResiduals(SS_size))

