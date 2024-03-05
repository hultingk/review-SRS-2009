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
pollination_2009 <- read.csv("founder_plant_2009.csv") # joined seed count data (structure level) with seed predation/structure count data (plant level)

pollination_2009 <- pollination_2009 %>%
  filter(patch != "A") # removing patch A (center patch) from analyses

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
AB_repr_mod <- glmmTMB(reproductive ~ ptype + dist_num + (1|block/patch/corner),
                       data = AB_reproduction,
                       family = binomial())
summary(AB_repr_mod)
plot(simulateResiduals(AB_repr_mod))
Anova(AB_repr_mod, type = "III")
#AB_repr_mod.posthoc <- emtrends(AB_repr_mod, pairwise ~ ptype, var = "dist_num")
#AB_repr_mod.posthoc
AB_repr_mod.posthoc <- emmeans(AB_repr_mod, "ptype")
pairs(AB_repr_mod.posthoc)


# Anthaenantia villosa
AV_repr_mod <- glmmTMB(reproductive ~ ptype + dist_num + (1|block/patch/corner),
                       data = AV_reproduction,
                       family = binomial())
summary(AV_repr_mod)
plot(simulateResiduals(AV_repr_mod))
Anova(AV_repr_mod, type = "III")
#AV_repr_mod.posthoc <- emtrends(AV_repr_mod, pairwise ~ ptype, var = "dist_num")
#AV_repr_mod.posthoc
AV_repr_mod.posthoc <- emmeans(AV_repr_mod, "ptype")
pairs(AV_repr_mod.posthoc)

# Carphephorus bellidifolius 
CB_repr_mod <- glmmTMB(reproductive ~ ptype + dist_num + (1|block/patch/corner),
                       data = CB_reproduction,
                       family = binomial())
summary(CB_repr_mod)
plot(simulateResiduals(CB_repr_mod))
Anova(CB_repr_mod, type = "III")
#CB_repr_mod.posthoc <- emtrends(CB_repr_mod, pairwise ~ ptype, var = "dist_num")
#CB_repr_mod.posthoc
CB_repr_mod.posthoc <- emmeans(CB_repr_mod, "ptype")
pairs(CB_repr_mod.posthoc)

# Liatris earlei 
LE_repr_mod <- glmmTMB(reproductive ~ ptype + dist_num + (1|block/patch/corner),
                       data = LE_reproduction,
                       family = binomial())
summary(LE_repr_mod)
plot(simulateResiduals(LE_repr_mod))
Anova(LE_repr_mod, type = "III")
#LE_repr_mod.posthoc <- emtrends(LE_repr_mod, pairwise ~ ptype, var = "dist_num")
#LE_repr_mod.posthoc
LE_repr_mod.posthoc <- emmeans(LE_repr_mod, "ptype")
pairs(LE_repr_mod.posthoc)

#Sorghastrum secundum
SS_repr_mod <- glmmTMB(reproductive ~ ptype + dist_num + (1|block/patch/corner),
                       data = SS_reproduction,
                       family = binomial())
summary(SS_repr_mod)
plot(simulateResiduals(SS_repr_mod))
Anova(SS_repr_mod, type = "III")
#SS_repr_mod.posthoc <- emtrends(SS_repr_mod, pairwise ~ ptype, var = "dist_num")
#SS_repr_mod.posthoc
SS_repr_mod.posthoc <- emmeans(SS_repr_mod, "ptype")
pairs(SS_repr_mod.posthoc)


###### reproductive status plotting: Figure 2 #####
# model predictions for plotting
predictAB_repr_mod <- ggpredict(AB_repr_mod, terms=c("dist_num", "ptype"), back.transform = T, allow.new.levels=TRUE)
predictAB_repr_mod <- predictAB_repr_mod %>%
  mutate(species = "Aristida")

predictAV_repr_mod <- ggpredict(AV_repr_mod, terms=c("dist_num", "ptype"), back.transform = T, allow.new.levels=TRUE)
predictAV_repr_mod <- predictAV_repr_mod %>%
  mutate(species = "Anthaenantia")

predictCB_repr_mod <- ggpredict(CB_repr_mod, terms=c("dist_num", "ptype"), back.transform = T, allow.new.levels=TRUE)
predictCB_repr_mod <- predictCB_repr_mod %>%
  mutate(species = "Carphephorus")

predictLE_repr_mod <- ggpredict(LE_repr_mod, terms=c("dist_num", "ptype"), back.transform = T, allow.new.levels=TRUE)
predictLE_repr_mod <- predictLE_repr_mod %>%
  mutate(species = "Liatris")

predictSS_repr_mod <- ggpredict(SS_repr_mod, terms=c("dist_num", "ptype"), back.transform = T, allow.new.levels=TRUE)
predictSS_repr_mod <- predictSS_repr_mod %>%
  mutate(species = "Sorghastrum")

total_repr_predict <- rbind(predictAB_repr_mod, predictAV_repr_mod, predictCB_repr_mod,
                            predictLE_repr_mod, predictSS_repr_mod)

# plotting
figure2 <- total_repr_predict %>%
  ggplot(aes(x = x, y = predicted)) +
  geom_line(aes(x, predicted, linetype = group), linewidth = 0.7) + 
  theme_bw()+
  #geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
  #            alpha = 0.15) +
  labs(title = NULL,
       x = "Distance from Edge",
       y = "Likelihood of flowering") +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title = element_text(size = 18)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text =element_text(size = 14)) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  guides(linetype = guide_legend(title = "Patch Type")) +
  theme(legend.key.width=unit(2,"cm")) #+
#scale_colour_manual(values = c("#648FFF", "#DC267F", "#FFB000")) +
#scale_fill_manual(values = c("#648FFF", "#DC267F", "#FFB000"))

# faceting
figure2 <- figure2 + facet_grid(rows = vars(species), scales = "free") +
  theme(strip.text.y = element_text(size = 16, colour = "black", angle = 270))
figure2

# exporting
svglite(file = "Figure2.svg", width = 9, height = 11)
plot(figure2)
dev.off()


#### pollination data wrangling: keeping only reproductive plants ####
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



###### pollination models ######
# Aristida beyrichiana
AB_mod <- glmmTMB(pollination_avg ~ ptype + dist_num + (1|block/patch/corner), 
                  data = AB_pollination,
                  family = betabinomial(),
                  weights = total_est_fl_avg)
summary(AB_mod)
plot(simulateResiduals(AB_mod)) # not the best
Anova(AB_mod, type = "III")
AB_mod.posthoc <- emmeans(AB_mod, "ptype")
pairs(AB_mod.posthoc)


# Anthaenantia villosa
AV_mod <- glmmTMB(pollination_avg ~ ptype + dist_num + (1|block/patch/corner), 
                  data = AV_pollination,
                  family = betabinomial(),
                  weights = total_est_fl_avg)
summary(AV_mod)
plot(simulateResiduals(AV_mod))
Anova(AV_mod, type = "III")
AB_mod.posthoc <- emmeans(AB_mod, "ptype")
pairs(AB_mod.posthoc)


# Carphephorus bellidifolius 
CB_mod <- glmmTMB(pollination_avg ~ ptype + dist_num + (1|block/patch/corner), 
                  data = CB_pollination,
                  family = betabinomial(),
                  weights = total_est_fl_avg)
summary(CB_mod)
Anova(CB_mod, type = "III")
plot(simulateResiduals(CB_mod))
CB_mod.posthoc <- emmeans(CB_mod, "ptype")
pairs(CB_mod.posthoc)


# Liatris earlei 
LE_mod <- glmmTMB(pollination_avg ~ ptype + dist_num + (1|block/patch/corner), 
                  data = LE_pollination,
                  family = betabinomial(),
                  weights = total_est_fl_avg)
summary(LE_mod)
Anova(LE_mod, type = "III")
plot(simulateResiduals(LE_mod))
LE_mod.posthoc <- emmeans(LE_mod, "ptype")
pairs(LE_mod.posthoc)


# Sorghastrum secundum
SS_mod <- glmmTMB(pollination_avg ~ ptype + dist_num + (1|block/patch/corner), 
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

total_rate_predict <- rbind(predictAB_rate_mod, predictAV_rate_mod, predictCB_rate_mod, predictLE_rate_mod,predictSS_rate_mod)

figureS1 <- total_rate_predict %>%
  ggplot(aes(x = x, y = predicted)) +
  geom_line(aes(x, predicted, linetype = group)) +
  #geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
  #            alpha = 0.15) +
  geom_point(aes(x = dist_num, y = pollination_avg), data = AB_pollination, alpha = 0.2) +
  geom_point(aes(x = dist_num, y = pollination_avg), data = CB_pollination, alpha = 0.2) +
  geom_point(aes(x = dist_num, y = pollination_avg), data = LE_pollination, alpha = 0.2) +
  geom_point(aes(x = dist_num, y = pollination_avg), data = SS_pollination, alpha = 0.2) +
  geom_point(aes(x = dist_num, y = pollination_avg), data = AV_pollination, alpha = 0.2) +
  theme_bw()+
  labs(title = NULL,
       x = "Distance from Edge",
       y = "Pollination rate") +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title = element_text(size = 18)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text =element_text(size = 14)) +
  guides(linetype = guide_legend(title = "Patch Type")) +
  theme(legend.key.width=unit(2,"cm"))

figureS1 <- figureS1 + facet_grid(rows = vars(species), scales = "free") + 
  theme(strip.text.y = element_text(size = 16, colour = "black", angle = 270)) #+
# theme(aspect.ratio=9/20)
figureS1

# exporting
svglite(file = "FigureS1.svg", width = 9, height = 11)
plot(figureS1)
dev.off()


#### seed production ####
# Aristida beyrichiana
AB_mod1 <- glmmTMB(plant_seed_prod ~ ptype + dist_num + (1|block/patch/corner), 
                   data = AB_pollination,
                   family = "nbinom2",
                   ziformula = ~1)
summary(AB_mod1)
plot(simulateResiduals(AB_mod1))
#check_overdispersion(AB_mod1)
#check_zeroinflation(AB_mod1)
Anova(AB_mod1, type = "III")
AB_mod1.posthoc <- emmeans(AB_mod1, "ptype")
pairs(AB_mod1.posthoc)



# Anthaenantia villosa
AV_mod1 <- glmmTMB(plant_seed_prod ~ ptype + dist_num + (1|block/patch/corner), 
                   data = AV_pollination,
                   family = "nbinom2",
                   ziformula = ~1) # bad residuals
summary(AV_mod1)
plot(simulateResiduals(AV_mod1))
#check_overdispersion(AV_mod1)
#check_zeroinflation(AV_mod1) 
AV_mod1.posthoc <- emmeans(AV_mod1, "ptype")
pairs(AV_mod1.posthoc)
Anova(AV_mod1, type = "III")
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
  theme(axis.title = element_text(size = 22)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1.5))
plotAV_mod1
# exporting
pdf(file = "FigureS2.pdf", width = 10, height = 7)
plotAV_mod1
dev.off()


# Carphephorus bellidifolius 
CB_mod1 <- glmmTMB(plant_seed_prod ~ ptype + dist_num + (1|block/patch/corner), 
                   data = CB_pollination,
                   family = "nbinom2",
                   ziformula = ~1)
summary(CB_mod1)
plot(simulateResiduals(CB_mod1))
#check_overdispersion(CB_mod1)
#check_zeroinflation(CB_mod1) 
Anova(CB_mod1, type = "III")
CB_mod1.posthoc <- emmeans(CB_mod1, "ptype")
pairs(CB_mod1.posthoc)


# Liatris earlei 
LE_mod1 <- glmmTMB(plant_seed_prod ~ ptype + dist_num + (1|block/patch/corner), 
                   data = LE_pollination,
                   family = "nbinom2")
summary(LE_mod1)
plot(simulateResiduals(LE_mod1))
#check_overdispersion(LE_mod1)
#check_zeroinflation(LE_mod1)
Anova(LE_mod1, type = "III")
LE_mod1.posthoc <- emmeans(LE_mod1, "ptype")
pairs(LE_mod1.posthoc)


# Sorghastrum secundum
SS_mod1 <- glmmTMB(plant_seed_prod ~ ptype + dist_num + (1|block/patch/corner), 
                   data = SS_pollination,
                   family = "nbinom2",
                   ziformula = ~1)
summary(SS_mod1)
plot(simulateResiduals(SS_mod1))
#check_overdispersion(SS_mod1)
#check_zeroinflation(SS_mod1)
Anova(SS_mod1, type = "III")
SS_mod1.posthoc <- emmeans(SS_mod1, "ptype")
pairs(SS_mod1.posthoc)


###### seed production plotting: Figure 3 ######
predictAB_mod1 <- ggpredict(AB_mod1, terms=c("dist_num"), back.transform = T, allow.new.levels=TRUE)
predictAB_mod1 <- predictAB_mod1 %>%
  mutate(species = "Aristida")

predictAV_mod1 <- ggpredict(AV_mod1, terms=c("dist_num"), back.transform = T, allow.new.levels=TRUE)
predictAV_mod1 <- predictAV_mod1 %>%
  mutate(species = "Anthaenantia")

predictCB_mod1 <- ggpredict(CB_mod1, terms=c("dist_num"), back.transform = T, allow.new.levels=TRUE)
predictCB_mod1 <- predictCB_mod1 %>%
  mutate(species = "Carphephorus")

predictLE_mod1 <- ggpredict(LE_mod1, terms=c("dist_num"), back.transform = T, allow.new.levels=TRUE)
predictLE_mod1 <- predictLE_mod1 %>%
  mutate(species = "Liatris")

predictSS_mod1 <- ggpredict(SS_mod1, terms=c("dist_num"), back.transform = T, allow.new.levels=TRUE)
predictSS_mod1 <- predictSS_mod1 %>%
  mutate(species = "Sorghastrum")

total_mod1_predict <- rbind(predictAB_mod1, predictAV_mod1, predictCB_mod1,
                            predictLE_mod1, predictSS_mod1)

figure3 <- total_mod1_predict %>% ggplot(aes(x = x, y = predicted)) +
  geom_point(aes(x = dist_num, y = plant_seed_prod), data = SS_pollination, alpha = 0.3) + 
  geom_point(aes(x = dist_num, y = plant_seed_prod), data = AB_pollination, alpha = 0.3) + 
  geom_point(aes(x = dist_num, y = plant_seed_prod), data = AV_pollination, alpha = 0.3) + 
  geom_point(aes(x = dist_num, y = plant_seed_prod), data = CB_pollination, alpha = 0.3) + 
  geom_point(aes(x = dist_num, y = plant_seed_prod), data = LE_pollination, alpha = 0.3) + 
  geom_line(aes(x, predicted, linetype = group)) +
  geom_line(aes(x, predicted), color = "red", size = 1.5) +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2) +
  theme_bw() +
  labs(title = NULL,
       x = "Distance from Edge",
       y = "Seed Production") +
  theme(axis.text = element_text(size = 16)) +
  theme(axis.title = element_text(size = 18)) +
  theme(legend.position = "none")
figure3 <- figure3 + facet_grid(rows = vars(species), scales = "free") +
  theme(strip.text.y = element_text(size = 16, colour = "black", angle = 270))
figure3

svglite(file = "Figure3.svg", width = 7, height = 12)
plot(figure3)
dev.off()





####### number of flowering structures ###########
# if there is a significant or marginally significant response of seed production to edge proximity, 
### testing response of number of flowering structures

# Aristida beyrichiana
AB_mod_fl <- glmmTMB(total_no_structures ~ ptype + dist_num + (1|block/patch/corner), 
                     data = AB_pollination,
                     family = "nbinom2")
summary(AB_mod_fl)
plot(simulateResiduals(AB_mod_fl))
Anova(AB_mod_fl, type = "III")
AB_mod_fl.posthoc <- emmeans(AB_mod_fl, "ptype")
pairs(AB_mod_fl.posthoc)

# Anthaenantia villosa
AV_mod_fl <- glmmTMB(total_no_structures ~ ptype + dist_num + (1|block/patch/corner), 
                     data = AV_pollination,
                     family = "nbinom2")
summary(AV_mod_fl)
Anova(AV_mod_fl, type = "III")
plot(simulateResiduals(AV_mod_fl))


# Carphephorus bellidifolius 
CB_mod_fl <- glmmTMB(total_no_structures ~ ptype + dist_num + (1|block/patch/corner), 
                     data = CB_pollination,
                     family = "nbinom2")
summary(CB_mod_fl)
plot(simulateResiduals(CB_mod_fl))
Anova(CB_mod_fl, type = "III")


# Sorghastrum secundum
SS_mod_fl <- glmmTMB(total_no_structures ~ ptype + dist_num + (1|block/patch/corner), 
                     data = SS_pollination,
                     family = "nbinom2")
summary(SS_mod_fl)
plot(simulateResiduals(SS_mod_fl))
Anova(SS_mod_fl, type = "III")












