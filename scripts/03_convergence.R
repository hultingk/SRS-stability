########
## SCRIPT NAME: 03_convergence.R
## AUTHOR: Katherine Hulting
## PURPOSE: Calculate convergence/divergence between patch type communities across time, repeat within dispersal mode groups
## PRODUCTS: Figure2.pdf, TableS1.html
#########

librarian::shelf(tidyverse, vegan, ape, BiodiversityR, glmmTMB, AICcmodavg, 
                 DHARMa, emmeans, car, ggeffects, performance, cowplot, here, kableExtra)

source(here::here("02_pcoa_permanova.R"))
source(here::here("00_functions.R"))

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% 
  filter(transplant != TRUE) %>% # removing experimentally planted species  
  filter(patch_type != "Center") # removing center patch from analysis


#########################
#### ALL SPECIES ####
#########################

##### Convergence/divergence between patch types based on distance #####
# iterate over blocks, for each patch pair within a block, compute jaccard dissimilarity for each year
# splitting into blocks, applying function, putting back together
convergence_jaccard <- srs_data %>%
  count(block, patch, patch_type, unique_id, year, time, sppcode) %>%
  group_by(block) %>%
  group_split() %>%
  lapply(compute_convergence_jaccard) %>%
  bind_rows() # putting together into a dataframe

# removing same patch type comparisons and time 0 (only for 52 and 57)
convergence_jaccard <- convergence_jaccard %>%
  filter(!patch_pair %in% c("Rectangular-Rectangular", "Winged-Winged")) %>%
  filter(time != 0) %>%
  mutate(dispersal_mode = "All Species")
convergence_jaccard$s.time <- as.numeric(scale(convergence_jaccard$time)) # scaling time


##### Models #####
# linear model
m.converge <- glmmTMB(jaccard ~ patch_pair * s.time + (1|block),
              data = convergence_jaccard)

# quadratic model
m.converge_quad <- glmmTMB(jaccard ~ patch_pair * s.time + patch_pair * I(s.time^2) + (1|block),
                         data = convergence_jaccard)
# null model
m.converge_null <- glmmTMB(jaccard ~ 1 + (1|block), # null model
                         data = convergence_jaccard)
# AIC comparison
a <- list(m.converge, m.converge_quad, m.converge_null)
aictab(a) # quadratic much better fit

## model checking
summary(m.converge_quad)

# percent change in dissimilarity from year 1-21 (20 years)
# time 1 = -1.5320869
# time 21 = 1.4326678
# intercept + time_estimate (time) + time^2_estimate (time)
# at time 1
0.386812 + 0.014180*(-1.5320869) + 0.004602*((-1.5320869)^2) # 0.3758892
# at time 21
0.386812 + 0.014180*(1.4326678) + 0.004602*((1.4326678)^2) # 0.416573

# percent change = time 21 - time 1 / time 1 * 100
(0.416573 - 0.3758892)/0.3758892 * 100 #10.82335 % increase

# 95% CI, percent change
confint(m.converge_quad)
# lower 95% CI
# at time 1
0.3632975047 + 0.0088328776*(-1.5320869) + -0.0008100583*((-1.5320869)^2) # 0.3478633
# at time 21
0.3632975047 + 0.0088328776*(1.4326678) + -0.0008100583*((1.4326678)^2) # 0.3742894
# percent change = time 21 - time 1 / time 1 * 100
(0.3742894 - 0.3478633)/0.3478633 * 100 #7.596691 % increase

# upper 95% CI
# at time 1
0.410326911 + 0.019526838*(-1.5320869) + 0.0100148923*((-1.5320869)^2) # 0.403918
# at time 21
0.410326911 + 0.019526838*(1.4326678) + 0.010014892*((1.4326678)^2) # 0.4588583
# percent change = time 21 - time 1 / time 1 * 100
(0.4588583 - 0.403918)/0.403918 * 100 #13.60184 % increase


# emmeans(m.converge_quad, ~s.time+I(s.time^2), at = list(s.time = c(-1.5320869, 1.4326678)),
#         type = "response")
# (0.382 - 0.398) / 0.398 * 100 # 8.28877 % increase over 20 years

# Anova
anova.converge <- Anova(m.converge_quad, type = "III")
# model checking
plot(simulateResiduals(m.converge_quad))
check_model(m.converge_quad)
performance::r2(m.converge_quad)

## posthoc comparisons
m.converge_posthoc <- emmeans(m.converge_quad, ~ patch_pair*s.time+ patch_pair * I(s.time^2))
m.converge_pairs <- pairs(m.converge_posthoc, simple = "patch_pair")
m.converge_pairs


# % increase in dissimilarity from (connected-winged) to (connected-rectangular)
(0.387-0.361)/0.361 * 100 # connected patches are %7.202216 more similar to winged patches than rectangular patches across time
#95% CI
(0.00535 / 0.361) * 100 #1.481994
7.2 +1.96 *1.481994 # upper = 10.10471
7.2 -1.96 *1.481994 # lower = 4.295292

(0.387-0.375)/0.375 * 100 # rectangular patches are %3.2 more similar to winged patches than connected patches across time
(0.375-0.361)/0.361 * 100 # winged patches are %3.878116 more similar to connected patches than rectangular patches across time





###################################
#### ANIMAL DISPERSED SPECIES ####
###################################

##### Convergence/divergence between patch types based on distance #####
# iterate over blocks, for each patch pair within a block, compute jaccard dissimilarity for each year
# splitting into blocks, applying function, putting back together
animal_convergence_jaccard <- srs_data %>%
  filter(dispersal_mode == "Animal") %>%
  count(block, patch, patch_type, unique_id, year, time, sppcode) %>%
  group_by(block) %>%
  group_split() %>%
  lapply(compute_convergence_jaccard) %>%
  bind_rows() # putting together into a dataframe

# removing same patch type comparisons and time 0 (only for 52 and 57)
animal_convergence_jaccard <- animal_convergence_jaccard %>%
  filter(!patch_pair %in% c("Rectangular-Rectangular", "Winged-Winged")) %>%
  filter(time != 0) %>%
  mutate(dispersal_mode = "Animal")

# scaling time
animal_convergence_jaccard$s.time <- as.numeric(scale(animal_convergence_jaccard$time))

##### Models #####
# animal dispersed linear
m.converge_animal_linear <- glmmTMB(jaccard ~ patch_pair*s.time + (1|block),
                                  data = animal_convergence_jaccard)

# animal dispersed quadratic
m.converge_animal_quad <- glmmTMB(jaccard ~ patch_pair*s.time + patch_pair*I(s.time^2) + (1|block),
                             data = animal_convergence_jaccard)

# animal dispersed null
m.converge_animal_null <- glmmTMB(jaccard ~ 1 + (1|block),
                                    data = animal_convergence_jaccard)

# AIC comparison
a <- list(m.converge_animal_linear, m.converge_animal_quad, m.converge_animal_null)
aictab(a) # linear better fit


# model checking
summary(m.converge_animal_linear)
plot(simulateResiduals(m.converge_animal_linear))
check_model(m.converge_animal_linear)
anova.animal.converge <- Anova(m.converge_animal_linear, type = "III")

# posthoc
m.converge_animal_posthoc <- emmeans(m.converge_animal_linear, ~ patch_pair*s.time)
m.converge_animal_pairs <- pairs(m.converge_animal_posthoc, simple = "patch_pair")
m.converge_animal_pairs





###################################
#### GRAVITY DISPERSED SPECIES ####
###################################

##### Convergence/divergence between patch types based on distance #####
# iterate over blocks, for each patch pair within a block, compute jaccard dissimilarity for each year
# splitting into blocks, applying function, putting back together
gravity_convergence_jaccard <- srs_data %>%
  filter(dispersal_mode == "Gravity") %>%
  count(block, patch, patch_type, unique_id, year, time, sppcode) %>%
  group_by(block) %>%
  group_split() %>%
  lapply(compute_convergence_jaccard) %>%
  bind_rows() # putting together into a dataframe

# removing same patch type comparisons and time 0 (only for 52 and 57)
gravity_convergence_jaccard <- gravity_convergence_jaccard %>%
  filter(!patch_pair %in% c("Rectangular-Rectangular", "Winged-Winged")) %>%
  filter(time != 0) %>%
  mutate(dispersal_mode = "Gravity")


# scaling time
gravity_convergence_jaccard$s.time <- as.numeric(scale(gravity_convergence_jaccard$time))

##### Models #####
# gravity dispersed linear
m.converge_gravity_linear <- glmmTMB(jaccard ~ patch_pair*s.time + (1|block),
                                    data = gravity_convergence_jaccard)

# gravity dispersed quadratic
m.converge_gravity_quad <- glmmTMB(jaccard ~ patch_pair*s.time + patch_pair*I(s.time^2) + (1|block),
                                  data = gravity_convergence_jaccard)

# gravity dispersed null
m.converge_gravity_null <- glmmTMB(jaccard ~ 1 + (1|block),
                                  data = gravity_convergence_jaccard)

# AIC comparison
a <- list(m.converge_gravity_linear, m.converge_gravity_quad, m.converge_gravity_null)
aictab(a) # quadratic better fit


# model checking
summary(m.converge_gravity_quad)
plot(simulateResiduals(m.converge_gravity_quad))
check_model(m.converge_gravity_quad)
anova.gravity.converge <- Anova(m.converge_gravity_quad, type = "III")

# posthoc
m.converge_gravity_posthoc <- emmeans(m.converge_gravity_quad, ~ patch_pair*s.time + patch_pair * I(s.time^2))
m.converge_gravity_pairs <- pairs(m.converge_gravity_posthoc, simple = "patch_pair")
m.converge_gravity_pairs






###################################
#### WIND DISPERSED SPECIES ####
###################################

##### Convergence/divergence between patch types based on distance #####
# iterate over blocks, for each patch pair within a block, compute jaccard dissimilarity for each year
# splitting into blocks, applying function, putting back together
wind_convergence_jaccard <- srs_data %>%
  filter(dispersal_mode == "Wind") %>%
  count(block, patch, patch_type, unique_id, year, time, sppcode) %>%
  group_by(block) %>%
  group_split() %>%
  lapply(compute_convergence_jaccard) %>%
  bind_rows() # putting together into a dataframe

# removing same patch type comparisons and time 0 (only for 52 and 57)
wind_convergence_jaccard <- wind_convergence_jaccard %>%
  filter(!patch_pair %in% c("Rectangular-Rectangular", "Winged-Winged")) %>%
  filter(time != 0) %>%
  mutate(dispersal_mode = "Wind")

# scaling time
wind_convergence_jaccard$s.time <- as.numeric(scale(wind_convergence_jaccard$time))

##### Models #####
# wind dispersed linear
m.converge_wind_linear <- glmmTMB(jaccard ~ patch_pair*s.time + (1|block),
                                     data = wind_convergence_jaccard)

# wind dispersed quadratic
m.converge_wind_quad <- glmmTMB(jaccard ~ patch_pair*s.time + patch_pair*I(s.time^2) + (1|block),
                                   data = wind_convergence_jaccard)

# wind dispersed null
m.converge_wind_null <- glmmTMB(jaccard ~ 1 + (1|block),
                                   data = wind_convergence_jaccard)

# AIC comparison
a <- list(m.converge_wind_linear, m.converge_wind_quad, m.converge_wind_null)
aictab(a) # quadratic better fit


# model checking
summary(m.converge_wind_quad)
plot(simulateResiduals(m.converge_wind_quad))
check_model(m.converge_wind_quad)
anova.wind.converge <- Anova(m.converge_wind_quad, type = "III")

# posthoc
m.converge_wind_posthoc <- emmeans(m.converge_wind_quad, ~ patch_pair*s.time + patch_pair * I(s.time^2))
m.converge_wind_pairs <- pairs(m.converge_wind_posthoc, simple = "patch_pair")
m.converge_wind_pairs




#########################
#### TABLES ####
#########################
# Anova table 
anova.converge_df <- as.data.frame(anova.converge)
anova.converge_df <- anova.converge_df %>%
  rownames_to_column("model_term") %>%
  mutate(`Dispersal mode` = "All Species") %>%
  mutate(`Top Model` = "Quadratic")

anova.animal.converge_df <- as.data.frame(anova.animal.converge)
anova.animal.converge_df <- anova.animal.converge_df %>%
  rownames_to_column("model_term") %>%
  mutate(`Dispersal mode` = "Animal-Dispersed") %>%
  mutate(`Top Model` = "Linear")

anova.gravity.converge_df <- as.data.frame(anova.gravity.converge)
anova.gravity.converge_df <- anova.gravity.converge_df %>%
  rownames_to_column("model_term") %>%
  mutate(`Dispersal mode` = "Gravity-Dispersed") %>%
  mutate(`Top Model` = "Quadratic")

anova.wind.converge_df <- as.data.frame(anova.wind.converge)
anova.wind.converge_df <- anova.wind.converge_df %>%
  rownames_to_column("model_term") %>%
  mutate(`Dispersal mode` = "Wind-Dispersed") %>%
  mutate(`Top Model` = "Quadratic")

m.converge_anova_all <- rbind(
  anova.converge_df, anova.animal.converge_df, anova.gravity.converge_df, anova.wind.converge_df
)

rename_variable_anova <- tibble(model_term = c("(Intercept)", "patch_pair", "s.time", "I(s.time^2)", "patch_pair:s.time", "patch_pair:I(s.time^2)"),
                                Variable = c("Intercept", "Patch Pair", "Time", "Time^2", "Patch Pair:Time", "Patch Pair:Time^2"))


m.converge_anova_all <- m.converge_anova_all %>%
  #filter(model_term != "(Intercept)") %>%
  left_join(rename_variable_anova, by = "model_term") %>%
  dplyr::select(`Dispersal mode`, `Top Model`, Variable, Chisq, Df, `Pr(>Chisq)`) %>%
  rename(p.value = `Pr(>Chisq)`, df = Df)

tableS1 <- m.converge_anova_all %>%
  kbl(digits = 3) %>%
  kable_classic(full_width = T) %>%
  kable_styling(html_font = "Times New Roman",
                font_size = 16) %>%
  collapse_rows(columns = c(1, 2), target = 1) %>%
  row_spec(0, extra_css = "border-bottom: 5px double;") %>%
  row_spec(1:nrow(m.converge_anova_all), extra_css = "border-bottom: 1px solid;") %>%
  row_spec(0:nrow(m.converge_anova_all), extra_css = "padding-bottom: 5px;")
tableS1

# exporting
#save_kable(tableS1, file = file.path("plots", "tableS1.html"))


# emmeans posthoc tables
# creating dataframes of results
# all species
m.converge_pairs_df <- as.data.frame(m.converge_pairs)
m.converge_pairs_df <- m.converge_pairs_df %>%
  mutate(`Dispersal mode` = "All Species") 

# animal dispersed
m.converge_animal_pairs_df <- as.data.frame(m.converge_animal_pairs)
m.converge_animal_pairs_df <- m.converge_animal_pairs_df %>%
  mutate(`Dispersal mode` = "Animal-Dispersed") 

# gravity dispersed
m.converge_gravity_pairs_df <- as.data.frame(m.converge_gravity_pairs)
m.converge_gravity_pairs_df <- m.converge_gravity_pairs_df %>%
  mutate(`Dispersal mode` = "Gravity-Dispersed")

# wind dispersed
m.converge_wind_pairs_df <- as.data.frame(m.converge_wind_pairs)
m.converge_wind_pairs_df <- m.converge_wind_pairs_df %>%
  mutate(`Dispersal mode` = "Wind-Dispersed") 

# putting all together
m.converge_emmeans_all <- rbind(
  m.converge_pairs_df, m.converge_animal_pairs_df, m.converge_gravity_pairs_df, m.converge_wind_pairs_df
)

# table
tableS2 <- m.converge_emmeans_all %>% 
  dplyr::select(`Dispersal mode`, contrast, estimate, SE, df, z.ratio, p.value) %>%
  kbl(digits = 3) %>%
  kable_classic(full_width = T) %>%
  kable_styling(html_font = "Times New Roman",
                font_size = 16) %>%
  collapse_rows(columns = 1) %>%
  row_spec(0, extra_css = "border-bottom: 5px double;") %>%
  row_spec(1:nrow(m.converge_emmeans_all), extra_css = "border-bottom: 1px solid;") %>%
  row_spec(0:nrow(m.converge_emmeans_all), extra_css = "padding-bottom: 5px;")
tableS2

# exporting
#save_kable(tableS2, file = file.path("plots", "tableS2.html"))




#########################
#### PLOTS ####
#########################
## All species dispersed plot
# model predictions
m.converge.predict <- ggpredict(m.converge_quad, terms=c("s.time [all]", "patch_pair [all]"), back_transform = T)
m.converge.predict <- as.data.frame(m.converge.predict)
m.converge.predict$dispersal_mode <- "All Species"

## animal dispersed plot
# model predictions
m.converge_animal.predict <- ggpredict(m.converge_animal_linear, terms=c("s.time [all]", "patch_pair [all]"), back_transform = T)
m.converge_animal.predict <- as.data.frame(m.converge_animal.predict)
m.converge_animal.predict$dispersal_mode <- "Animal"

## gravity dispersed plot
# model predictions
m.converge_gravity.predict <- ggpredict(m.converge_gravity_quad, terms=c("s.time [all]", "patch_pair [all]"), back_transform = T)
m.converge_gravity.predict <- as.data.frame(m.converge_gravity.predict)
m.converge_gravity.predict$dispersal_mode <- "Gravity"

## wind dispersed plot
# model predictions
m.converge_wind.predict <- ggpredict(m.converge_wind_quad, terms=c("s.time [all]", "patch_pair [all]"), back_transform = T)
m.converge_wind.predict <- as.data.frame(m.converge_wind.predict)
m.converge_wind.predict$dispersal_mode <- "Wind"


# FACET BY ROWS - total and animal together and gravity and wind together
# joining together predictions
predict_converge_1 <- rbind(
  m.converge.predict, m.converge_animal.predict
)
predict_converge_2 <- rbind(
  m.converge_gravity.predict, m.converge_wind.predict
)


# creating key of scaled times to join to predictions for easy visualization
scaled_time_key <- convergence_jaccard %>%
  count(time, s.time) %>%
  dplyr::select(-n) %>%
  mutate(s.time = round(s.time, 2))

# joining with time
predict_converge_1 <- predict_converge_1 %>%
  left_join(scaled_time_key, by = c("x" = "s.time"))
predict_converge_2 <- predict_converge_2 %>%
  left_join(scaled_time_key, by = c("x" = "s.time"))

predict_converge_1$dispersal_mode <- factor(predict_converge_1$dispersal_mode, levels = c("All Species", "Animal"))
predict_converge_2$dispersal_mode <- factor(predict_converge_2$dispersal_mode, levels = c("Gravity", "Wind"))

# joining together data points
dispersal_mode_convergence_1 <- rbind(
  convergence_jaccard, animal_convergence_jaccard
)
dispersal_mode_convergence_1$dispersal_mode <- factor(dispersal_mode_convergence_1$dispersal_mode, levels = c("All Species", "Animal"))

dispersal_mode_convergence_2 <- rbind(
  gravity_convergence_jaccard, wind_convergence_jaccard
)
dispersal_mode_convergence_2$dispersal_mode <- factor(dispersal_mode_convergence_2$dispersal_mode, levels = c("Gravity", "Wind"))

# first set of plots
converge_plot_1 <- predict_converge_1 %>%
  ggplot() +
  geom_point(aes(time, jaccard, color = patch_pair), size = 3, alpha = 0.05, data = dispersal_mode_convergence_1) +
  geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4) +
  geom_line(aes(time, predicted, color = group), linewidth = 1.4) +
  facet_wrap(~dispersal_mode, scales = "free", labeller = as_labeller(c("All Species" = "(A) All species", "Animal" = "(B) Animal-dispersed"))) +
  theme_minimal(base_size = 20) +
  theme(panel.border = element_rect(color = "darkgrey", fill=NA, linewidth=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "darkgrey", linewidth = 0.5),
        strip.text.x = element_text(hjust = -0.05)) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
  xlab(NULL) +
  ylab(expression(paste("Spatial ", beta, " diversity (Jaccard)"))) +
  guides(fill=guide_legend(ncol=1)) +
  guides(color=guide_legend(ncol=1)) +
  theme(axis.text = element_text(size = 14)) +
  theme(legend.position = "none") 
converge_plot_1

# second set of plots
converge_plot_2 <- predict_converge_2 %>%
  ggplot() +
  geom_point(aes(time, jaccard, color = patch_pair), size = 3, alpha = 0.05, data = dispersal_mode_convergence_2) +
  geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4) +
  geom_line(aes(time, predicted, color = group), linewidth = 1.4) +
  facet_wrap(~dispersal_mode, scales = "free", labeller = as_labeller(c("Gravity" = "(C) Gravity-dispersed", "Wind" = "(D) Wind-dispersed"))) +
  theme_minimal(base_size = 20) +
  theme(panel.border = element_rect(colour = "darkgrey", fill=NA, linewidth=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "darkgrey", linewidth = 0.5),
        strip.text.x = element_text(hjust = -0.05)) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
  xlab("Years since site creation") +
  ylab(expression(paste("Spatial ", beta, " diversity (Jaccard)"))) +
  guides(fill=guide_legend(ncol=1)) +
  guides(color=guide_legend(ncol=1)) +
  theme(axis.text = element_text(size = 14)) +
  theme(legend.position = "none") 
converge_plot_2

# get legend
pL <- m.converge_wind.predict %>%
  left_join(scaled_time_key, by = c("x" = "s.time")) %>%
  ggplot() +
  geom_point(aes(time, jaccard, color = patch_pair), size = 4, alpha = 0.05, data = wind_convergence_jaccard) +
  geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.5) +
  geom_line(aes(time, predicted, color = group), linewidth = 1.5) +
  theme_minimal(base_size = 20) +
  theme(panel.border = element_rect(colour = "darkgrey", fill=NA, linewidth=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison")
l <- get_legend(pL)

# put together
figure2 <- cowplot::plot_grid(converge_plot_1, l, converge_plot_2, 
                                      ncol = 2, nrow = 2, rel_widths = c(1, 0.4), rel_heights = c(1, 1.1),
                                      label_size = 20, label_x = 0.2, label_y = 0.95)
figure2

# exporting
# pdf(file = file.path("plots", "figure2.pdf"), width = 11.5, height = 8.7)
# figure2
# dev.off()




### individual total plot

### plotting model predictions
# creating key of scaled times to join to predictions for easy visualization
# scaled_time_key <- convergence_jaccard %>%
#   count(time, s.time) %>%
#   dplyr::select(-n) %>%
#   mutate(s.time = round(s.time, 2))
# 
# # model predictions
# m.converge.predict <- ggpredict(m.converge_quad, terms=c("s.time [all]", "patch_pair [all]"), back_transform = T)
# m.converge.predict <- as.data.frame(m.converge.predict)
# m.converge.predict$dispersal_mode <- "Total"
# # plotting
# convergence_plot <- m.converge.predict %>%
#   left_join(scaled_time_key, by = c("x" = "s.time")) %>%
#   ggplot() +
#   geom_point(aes(time, jaccard, color = patch_pair), size = 4, alpha = 0.05, data = convergence_jaccard) +
#   geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.5) +
#   geom_line(aes(time, predicted, color = group), linewidth = 1.5) +
#   theme_minimal(base_size = 22) +
#   scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
#   scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
#   #xlab("Time since site creation (years)") +
#   ylab(expression(paste("Spatial ", beta, " diversity (Jaccard)"))) +
#   xlab(NULL) +
#   theme(legend.position = "none")
#  # annotate("text", x = 18, y=0.59, label = expression(paste('R'^2*' = 0.321')), size=7)
# convergence_plot

# pdf(file = file.path("plots", "convergence_plot.pdf"), width = 12, height = 8)
# convergence_plot
# dev.off()

