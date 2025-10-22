librarian::shelf(tidyverse, vegan, ape, BiodiversityR, glmmTMB, AICcmodavg, 
                 DHARMa, emmeans, car, ggeffects, performance, cowplot)

source(here::here("02_pcoa_permanova.R"))
source(here::here("00_functions.R"))

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE) %>%
  #filter(!block %in% c("75W", "75E")) %>%
  filter(patch_type != "center")


# # pivot to wider format
# srs_data_wider <- srs_data %>%
#  # filter(!time %in% c("0", "4")) %>%
#   dplyr::count(unique_id, time, year, sppcode, soil_moisture, year_since_fire) %>%
#   pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) # wide format
# 
# # make factor
# srs_data_wider$time <- as.factor(srs_data_wider$time)
# srs_data_wider$unique_id <- as.factor(srs_data_wider$unique_id)
# srs_data_wider$year <- as.factor(srs_data_wider$year)
# srs_data_wider$year_since_fire <- as.numeric(srs_data_wider$year_since_fire)
# 
# # patch data
# patch_info <- srs_data_wider %>% 
#   arrange(unique_id, time) %>%
#   select(unique_id, time, year, soil_moisture, year_since_fire)
# 
# patch_info <- patch_info %>%
#   separate(unique_id, into = c("block", "patch_rep", "patch_type"), sep = "-", remove = F) %>%
#   mutate(patch_time = paste(patch_type, time, sep = "-"))
# 
# 
# # species matrix
# sp_info <- srs_data_wider %>%
#   arrange(unique_id, time) %>%
#   mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
#   column_to_rownames("unique_id_year") %>%
#   select(!c("unique_id", "time", "year", "soil_moisture", "year_since_fire"))
# 
# 
# 
# 
# 
# 
# # calculate distance to centroid within a patch type
# bd_patch_list <- betadisper(vegdist(sp_info, method = "jaccard"), patch_info$patch_time, type = "centroid")   
# 
# # make a dataframe
# bd_patch <- as.data.frame(bd_patch_list$group.distances, row.names = rownames(bd_patch_list$group))
# bd_patch <- bd_patch %>%
#   rownames_to_column("patch_time") %>%
#   separate(patch_time, into = c("patch_type", "time"), sep = "-")
# 
# # renaming column
# bd_patch$distance_centroids <- bd_patch$`bd_patch_list$group.distances`
# bd_patch$`bd_patch_list$group.distances` <- NULL
# bd_patch$time <- as.numeric(bd_patch$time)
# 
# ##### trying with individual patch distances to group centroid
# # make a dataframe
# dist_patch <- as.data.frame(bd_patch_list$distances)
# dist_patch <- dist_patch %>%
#   rownames_to_column("patch_time") %>%
#   separate(patch_time, into = c("block", "patch_rep", "patch_type", "time", "year"), sep = "-")
# 
# # renaming column
# dist_patch$distance_centroids <- dist_patch$`bd_patch_list$distances`
# dist_patch$`bd_patch_list$distances` <- NULL
# dist_patch$time <- as.numeric(dist_patch$time)
# 
# 
# 
# 
# #############################################################
# ######## convergence/divergence WITHIN a patch type #########
# #### plot distance to patch type centroid (convergence within a patch type)
# bd_patch %>%
#   ggplot(aes(time, distance_centroids, color = patch_type, fill = patch_type)) +
#   geom_point() +
#   geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.2) +
#   theme_minimal(base_size = 16) +
#   xlab("Time since site creation (years)")+
#   ylab("Average distance to patch type centroid") +
#   scale_fill_brewer(palette = "Set2") +
#   scale_color_brewer(palette ="Set2")
# 
# # alternative - broken up by block, not an average
# conv_within_patch_plot <- dist_patch %>%
#   mutate(patch_type = dplyr::case_when(
#     patch_type %in% c("connected") ~ "Connected",
#     patch_type %in% c("rectangle") ~ "Rectangular",
#     patch_type %in% c("wing") ~ "Winged"
#   )) %>%
#   ggplot(aes(time, distance_centroids, color = patch_type, fill = patch_type)) +
#   geom_point(size = 3, alpha = 0.3) +
#   geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.5, linewidth = 2) +
#   theme_minimal(base_size = 28) +
#   xlab("Time since site creation (years)")+
#   ylab("Distance to patch type centroid") +
#   scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
#   scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type")
# conv_within_patch_plot
# 
# 
# # 
# # pdf(file = file.path("plots", "NO75_conv_within_patch_plot.pdf"), width = 13, height =8)
# # conv_within_patch_plot
# # dev.off()
# 
# 
# # models
# m_centroid <- glmmTMB(distance_centroids ~ patch_type + time + (1|block/patch_rep),
#                       data = dist_patch)
# summary(m_centroid)
# 
# m_centroid_quad <- glmmTMB(distance_centroids ~ patch_type + time + I(time^2) + (1|block/patch_rep),
#                       data = dist_patch)
# summary(m_centroid_quad)
# 
# m_centroid_null <- glmmTMB(distance_centroids ~ 1,
#                            data = dist_patch)
# summary(m_centroid_null)
# 
# # AIC comparison
# a <- list(m_centroid, m_centroid_quad, m_centroid_null)
# aictab(a) # quadratic much better fit
# 
# 
# 
# 
# #############################################################
# ####### calculating distance BETWEEN patch type centroids #######
# # pulling out centroids
# # centroids <- as.data.frame(bd_patch_list[["centroids"]])
# # centroids <- centroids[,1:2] # only keeping first 2 axes
# # centroids <- centroids %>%
# #   rownames_to_column("patch_time") %>%
# #   separate(patch_time, into = c("patch_type", "time"), sep = "-", remove = F)
# # 
# # # separating by patch type
# # centroids_c <- centroids %>% filter(patch_type == "connected")
# # centroids_w <- centroids %>% filter(patch_type == "wing")
# # centroids_r <- centroids %>% filter(patch_type == "rectangle")
# # 
# # # putting all together
# # centroids_all <- merge(centroids_c, centroids_w, by = "time", all.x = TRUE)
# # centroids_all <- merge(centroids_all, centroids_r, by = "time", all.x = TRUE)
# # centroids_all <- centroids_all %>% # renaming columns for clarity
# #   rename(c_PCoA1 = PCoA1.x, c_PCoA2 = PCoA2.x,
# #          w_PCoA1 = PCoA1.y, w_PCoA2 = PCoA2.y,
# #          r_PCoA1 = PCoA1, r_PCoA2 = PCoA2)
# # 
# # # calculate distance between centroids in pairs of patch types
# # centroids_all <- centroids_all %>%
# #   mutate(dist_c_w = sqrt((w_PCoA1 - c_PCoA1)^2 + (w_PCoA2 - c_PCoA2)^2)) %>%
# #   mutate(dist_c_r = sqrt((r_PCoA1 - c_PCoA1)^2 + (r_PCoA2 - c_PCoA2)^2)) %>%
# #   mutate(dist_r_w = sqrt((w_PCoA1 - r_PCoA1)^2 + (w_PCoA2 - r_PCoA2)^2))
# # 
# # # putting together
# # bw_group_dist <- centroids_all %>%
# #   select(time, dist_c_w, dist_c_r, dist_r_w)
# # bw_group_dist <- bw_group_dist %>%
# #   pivot_longer(cols = c("dist_c_w", "dist_c_r", "dist_r_w"), names_to = "patch_pair", values_to = "distance")
# # bw_group_dist$time <- as.numeric(bw_group_dist$time)
# # 
# # # plotting
# # bw_group_dist %>%
# #   ggplot(aes(time, distance, color = patch_pair, fill = patch_pair)) +
# #   geom_point() +
# #   geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.2) +
# #   theme_minimal(base_size = 16) +
# #   scale_fill_brewer(palette = "Set2") +
# #   scale_color_brewer(palette ="Set2") +
# #   xlab("Time since site creation (years)") +
# #   ylab("Average distance between patch type centroids")
# # 
# # # models
# # m_centroid_pt <- glmmTMB(distance ~ time * patch_pair,
# #                          data = bw_group_dist)
# # summary(m_centroid_pt)
# # 
# # 
# #### calculate distance between patch types within a block for every year
# # use PCoA axes
# pcoa_dist_bw_patch <- pcoa_axes %>%
#   select(!c("soil_moisture", "year_since_fire")) %>%
#   separate(unique_id, into = c("block", "patch_rep", "patch_type"), remove = F) %>%
#   #filter(!block %in% c("75W", "75E")) %>%
#   mutate(block_time = paste(block, time, sep = "-"))
# 
# 
# # separating by patch replicate
# dist_bw_b <- pcoa_dist_bw_patch %>% filter(patch_rep == "B") %>%
#   rename(PCoA1.B = Axis.1, PCoA2.B = Axis.2)
# dist_bw_c <- pcoa_dist_bw_patch %>% filter(patch_rep == "C") %>%
#   rename(PCoA1.C = Axis.1, PCoA2.C = Axis.2)
# dist_bw_d <- pcoa_dist_bw_patch %>% filter(patch_rep == "D") %>%
#   rename(PCoA1.D = Axis.1, PCoA2.D = Axis.2)
# dist_bw_e <- pcoa_dist_bw_patch %>% filter(patch_rep == "E") %>%
#   rename(PCoA1.E = Axis.1, PCoA2.E = Axis.2)
# 
# # joining all together
# dist_bw_all <- dist_bw_b %>%
#   left_join(dist_bw_c, by = c("block_time", "block", "time", "year")) %>%
#   left_join(dist_bw_d, by = c("block_time", "block", "time", "year")) %>%
#   left_join(dist_bw_e, by = c("block_time", "block", "time", "year"))
# 
# # calculate distance using pythagorean theorem
# dist_bw_all <- dist_bw_all %>%
#   mutate(dist_b_c = sqrt((PCoA1.C - PCoA1.B)^2 + (PCoA2.C - PCoA2.B)^2)) %>%
#   mutate(dist_b_d = sqrt((PCoA1.D - PCoA1.B)^2 + (PCoA2.D - PCoA2.B)^2)) %>%
#   mutate(dist_b_e = sqrt((PCoA1.E - PCoA1.B)^2 + (PCoA2.E - PCoA2.B)^2)) %>%
#   mutate(dist_c_d = sqrt((PCoA1.D - PCoA1.C)^2 + (PCoA2.D - PCoA2.C)^2)) %>%
#   mutate(dist_c_e = sqrt((PCoA1.E - PCoA1.C)^2 + (PCoA2.E - PCoA2.C)^2)) %>%
#   mutate(dist_d_e = sqrt((PCoA1.E - PCoA1.D)^2 + (PCoA2.E - PCoA2.D)^2))
# 
# 
# # putting together
# dist_bw_all <- dist_bw_all %>%
#   select(block, time, dist_b_c, dist_b_d, dist_b_e, dist_c_d, dist_c_e, dist_d_e)
# dist_bw_all <- dist_bw_all %>%
#   pivot_longer(cols = c("dist_b_c", "dist_b_d", "dist_b_e", "dist_c_d", "dist_c_e", "dist_d_e"), 
#                names_to = "patch_pair", values_to = "distance")
# dist_bw_all$time <- as.numeric(dist_bw_all$time)
# 
# # reading in key to pairs of patches
# patch_pair_ID <- read.csv(file = file.path("data", "L2_summarized", "patch_pair_ID.csv"))
# # joining to data 
# dist_bw_all <- dist_bw_all %>%
#   left_join(patch_pair_ID, by = c("block", "patch_pair")) %>%
#   filter(!patch_pair_ID %in% c("Winged-Winged", "Rectangular-Rectangular"))
#   
# 
# dist_bw_all %>%
#   count(block)
# conv_bw_patch_plot <- dist_bw_all %>%
#   ggplot(aes(time, distance, color = patch_pair_ID, fill = patch_pair_ID)) +
#   geom_point(size = 3, alpha = 0.3) +
#   geom_smooth(method = "lm", alpha = 0.5, linewidth = 2) +
#   theme_minimal(base_size = 24) +
#   scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
#   scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
#   xlab("Time since site creation (years)") +
#   ylab("Distance between patch type communities") +
#   annotate("text", x = 4, y=0.35, label = expression(paste('R'^2*' = 0.264')), size=7)
# conv_bw_patch_plot
# 
# pdf(file = file.path("plots", "conv_bw_patch_plot.pdf"), width = 13, height =8)
# conv_bw_patch_plot
# dev.off()
# 
# # models
# m_centroid_pt <- glmmTMB(distance ~ patch_pair_ID * time + (1|block),
#                          data = dist_bw_all)
# summary(m_centroid_pt)
# plot(simulateResiduals(m_centroid_pt))
# Anova(m_centroid_pt, type = "III")
# 0.0005185 * 10 # 0.005185
# (0.005185 / 0.1095479) * 100 # 9.646922%
# performance::r2(m_centroid_pt)
# ## CI 
# 0.0010568 * 10 # 0.010568
# (0.010568 / 0.1095479) * 100 # 4.733089%
# 9.646922-4.733089 # 4.913833
# 9.646922+4.733089 # 14.38001
# 
# m_centroid_pt_quad <- glmmTMB(distance ~ patch_pair_ID * time + patch_pair_ID * I(time^2) + (1|block),
#                          data = dist_bw_all)
# summary(m_centroid_pt_quad)
# 
# m_centroid_pt_null <- glmmTMB(distance ~ 1 + (1|block),
#                               data = dist_bw_all)
# summary(m_centroid_pt_null)
# 
# # AIC comparison
# a <- list(m_centroid_pt, m_centroid_pt_quad, m_centroid_pt_null)
# aictab(a) # quadratic much better fit
# 
# 
# # posthoc
# m_centroid_pt_posthoc <- emmeans(m_centroid_pt, ~ patch_pair_ID*time)
# pairs(m_centroid_pt_posthoc, simple = "patch_pair_ID")
# m_centroid_pt_posthoc
# # how much more similar are rectangular patches to winged patches compared to connected patches
# 0.122 - 0.108 # 0.014
# (0.014/0.108)*100 # 12.96296 % more similar to winged patches
# 
# # CI
# 0.0105 - 0.0103 # 2e-04
# (2e-04/0.0103)*100 # 1.941748
# 12.96296 - 1.941748
# 12.96296 + 1.941748
# 
# pairs(emtrends(m_centroid_pt, "patch_pair_ID", var = "time"))
# 
# 
# 
# 
# 
# 







#### convergence/divergence between patch types based on distance ####
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
  filter(!patch_pair %in% c("rectangle-rectangle", "wing-wing", "Rectangular-Rectangular", "Winged-Winged")) %>%
  filter(time != 0) %>%
  mutate(dispersal_mode = "Total")
convergence_jaccard$s.time <- as.numeric(scale(convergence_jaccard$time)) # scaling time


## models
# linear model
m.converge <- glmmTMB(jaccard ~ patch_pair * s.time + (1|block),
              data = convergence_jaccard)
summary(m.converge)
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

Anova(m.converge_quad)
plot(simulateResiduals(m.converge_quad))
check_model(m.converge_quad)
performance::r2(m.converge_quad)
## posthoc comparisons
m.converge_posthoc <- emmeans(m.converge_quad, ~ patch_pair*s.time+ patch_pair * I(s.time^2))
pairs(m.converge_posthoc, simple = "patch_pair")
m.converge_posthoc
m.converge_posthoc <- emtrends(m.converge_quad, "patch_pair", var = "s.time")
pairs(m.converge_posthoc)

# % increase in dissimilarity from (connected-winged) to (connected-rectangular)
(0.387-0.361)/0.361 * 100 # connected patches are %7.202216 more similar to winged patches than rectangular patches across time
#95% CI
(0.00535 / 0.361) * 100 #1.481994
7.2 +1.96 *1.481994 # upper = 10.10471
7.2 -1.96 *1.481994 # lower = 4.295292

(0.387-0.375)/0.375 * 100 # rectangular patches are %3.2 more similar to winged patches than connected patches across time
(0.375-0.361)/0.361 * 100 # winged patches are %3.878116 more similar to connected patches than rectangular patches across time


### plotting model predictions
# creating key of scaled times to join to predictions for easy visualization
scaled_time_key <- convergence_jaccard %>%
  count(time, s.time) %>%
  dplyr::select(-n) %>%
  mutate(s.time = round(s.time, 2))

# model predictions
m.converge.predict <- ggpredict(m.converge_quad, terms=c("s.time [all]", "patch_pair [all]"), back_transform = T)
m.converge.predict <- as.data.frame(m.converge.predict)
m.converge.predict$dispersal_mode <- "Total"
# plotting
convergence_plot <- m.converge.predict %>%
  left_join(scaled_time_key, by = c("x" = "s.time")) %>%
  ggplot() +
  geom_point(aes(time, jaccard, color = patch_pair), size = 4, alpha = 0.05, data = convergence_jaccard) +
  geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.5) +
  geom_line(aes(time, predicted, color = group), linewidth = 1.5) +
  theme_minimal(base_size = 22) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
  #xlab("Time since site creation (years)") +
  ylab(expression(paste("Spatial ", beta, " diversity (Jaccard)"))) +
  xlab(NULL) +
  theme(legend.position = "none")
 # annotate("text", x = 18, y=0.59, label = expression(paste('R'^2*' = 0.321')), size=7)
convergence_plot

pdf(file = file.path("plots", "convergence_plot.pdf"), width = 12, height = 8)
convergence_plot
dev.off()





#### Dispersal mode convergence/divergence between patch types ####
##### animal dispersed convergence/divergence #####
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

# ##### gravity dispersed convergence/divergence #####
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


# ##### wind dispersed convergence/divergence #####
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
animal_convergence_jaccard$s.time <- as.numeric(scale(animal_convergence_jaccard$time))
gravity_convergence_jaccard$s.time <- as.numeric(scale(gravity_convergence_jaccard$time))
wind_convergence_jaccard$s.time <- as.numeric(scale(wind_convergence_jaccard$time))
dispersal_mode_convergence <- rbind(
  convergence_jaccard, animal_convergence_jaccard, gravity_convergence_jaccard, wind_convergence_jaccard
)


### individual dispersal mode models
# animal dispersed
m.converge_animal <- glmmTMB(jaccard ~ patch_pair*s.time + patch_pair*I(s.time^2) + (1|block),
                             data = animal_convergence_jaccard)
summary(m.converge_animal)
m.converge_animal_posthoc <- emmeans(m.converge_animal, ~ patch_pair*s.time + patch_pair * I(s.time^2))
pairs(m.converge_animal_posthoc, simple = "patch_pair")

# gravity dispersed
m.converge_gravity <- glmmTMB(jaccard ~ patch_pair*s.time + patch_pair*I(s.time^2) + (1|block),
                              data = gravity_convergence_jaccard)
summary(m.converge_gravity)
m.converge_gravity_posthoc <- emmeans(m.converge_gravity, ~ patch_pair*s.time + patch_pair * I(s.time^2))
pairs(m.converge_gravity_posthoc, simple = "patch_pair")
m.converge_gravity_posthoc <- emtrends(m.converge_gravity, "patch_pair", var = "s.time")
pairs(m.converge_gravity_posthoc)

# wind
m.converge_wind <- glmmTMB(jaccard ~ patch_pair*s.time + patch_pair*I(s.time^2) + (1|block),
                           data = wind_convergence_jaccard)
summary(m.converge_wind)
m.converge_wind_posthoc <- emmeans(m.converge_wind, ~ patch_pair*s.time + patch_pair * I(s.time^2), at = list(s.time = 1.88))
pairs(m.converge_wind_posthoc, simple = "patch_pair")
m.converge_wind_posthoc <- emtrends(m.converge_wind, "patch_pair", var = "I(s.time^2)")
pairs(m.converge_wind_posthoc)

## animal dispersed plot
# model predictions
m.converge_animal.predict <- ggpredict(m.converge_animal, terms=c("s.time [all]", "patch_pair [all]"), back_transform = T)
m.converge_animal.predict <- as.data.frame(m.converge_animal.predict)
m.converge_animal.predict$dispersal_mode <- "Animal"

## gravity dispersed plot
# model predictions
m.converge_gravity.predict <- ggpredict(m.converge_gravity, terms=c("s.time [all]", "patch_pair [all]"), back_transform = T)
m.converge_gravity.predict <- as.data.frame(m.converge_gravity.predict)
m.converge_gravity.predict$dispersal_mode <- "Gravity"

## wind dispersed plot
# model predictions
m.converge_wind.predict <- ggpredict(m.converge_wind, terms=c("s.time [all]", "patch_pair [all]"), back_transform = T)
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
# joining with time
predict_converge_1 <- predict_converge_1 %>%
  left_join(scaled_time_key, by = c("x" = "s.time"))
predict_converge_2 <- predict_converge_2 %>%
  left_join(scaled_time_key, by = c("x" = "s.time"))

predict_converge_1$dispersal_mode <- factor(predict_converge_1$dispersal_mode, levels = c("Total", "Animal"))
predict_converge_2$dispersal_mode <- factor(predict_converge_2$dispersal_mode, levels = c("Gravity", "Wind"))

# joining together data points
dispersal_mode_convergence_1 <- rbind(
  convergence_jaccard, animal_convergence_jaccard
)
dispersal_mode_convergence_1$dispersal_mode <- factor(dispersal_mode_convergence_1$dispersal_mode, levels = c("Total", "Animal"))

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
  facet_wrap(~dispersal_mode, scales = "free") +
  theme_minimal(base_size = 20) +
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
  facet_wrap(~dispersal_mode, scales = "free") +
  theme_minimal(base_size = 20) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
  xlab("Time since site creation (years)") +
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
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison")
l <- get_legend(pL)

# put together
convergence_all <- cowplot::plot_grid(converge_plot_1, l, converge_plot_2, 
                                      ncol = 2, nrow = 2, rel_widths = c(1, 0.4), rel_heights = c(1, 1.1),
                                      label_size = 20, label_x = 0.2, label_y = 0.95)
convergence_all
# exporting
pdf(file = file.path("plots", "convergence_plot_all.pdf"), width = 11.5, height = 9)
convergence_all
dev.off()



### individual plots
# animal_converge_plot <- m.converge_animal.predict %>%
#   left_join(scaled_time_key, by = c("x" = "s.time")) %>%
#   ggplot() +
#   geom_point(aes(time, jaccard, color = patch_pair), size = 4, alpha = 0.05, data = animal_convergence_jaccard) +
#   geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.5) +
#   geom_line(aes(time, predicted, color = group), linewidth = 1.5) +
#   theme_minimal(base_size = 22) +
#   scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
#   scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
#   xlab(NULL) +
#   ylab(NULL) +
#   theme(legend.position = "none")
#  # theme(legend.position = "right") +
#   #theme(legend.justification.bottom = "left")
# animal_converge_plot
# 
# gravity_converge_plot <- m.converge_gravity.predict %>%
#   left_join(scaled_time_key, by = c("x" = "s.time")) %>%
#   ggplot() +
#   geom_point(aes(time, jaccard, color = patch_pair), size = 4, alpha = 0.05, data = gravity_convergence_jaccard) +
#   geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.5) +
#   geom_line(aes(time, predicted, color = group), linewidth = 1.5) +
#   theme_minimal(base_size = 22) +
#   scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
#   scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
#   xlab("Time since site creation (years)") +
#   ylab(expression(paste("Spatial ", beta, " diversity (Jaccard)"))) +
#   theme(legend.position = "none") 
# gravity_converge_plot
# 
# wind_converge_plot <- m.converge_wind.predict %>%
#   left_join(scaled_time_key, by = c("x" = "s.time")) %>%
#   ggplot() +
#   geom_point(aes(time, jaccard, color = patch_pair), size = 4, alpha = 0.05, data = wind_convergence_jaccard) +
#   geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.5) +
#   geom_line(aes(time, predicted, color = group), linewidth = 1.5) +
#   theme_minimal(base_size = 22) +
#   scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
#   scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
#   xlab("Time since site creation (years)") +
#   ylab(NULL) +
#   theme(legend.position = "none") 
# wind_converge_plot
# 
# pL <- m.converge_wind.predict %>%
#   left_join(scaled_time_key, by = c("x" = "s.time")) %>%
#   ggplot() +
#   geom_point(aes(time, jaccard, color = patch_pair), size = 4, alpha = 0.05, data = wind_convergence_jaccard) +
#   geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.5) +
#   geom_line(aes(time, predicted, color = group), linewidth = 1.5) +
#   theme_minimal(base_size = 20) +
#   scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
#   scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison")
# l <- get_legend(pL)
# 
# 
# convergence_all <- cowplot::plot_grid(convergence_plot, animal_converge_plot, l, gravity_converge_plot, wind_converge_plot,
#                                       ncol = 3, nrow = 2, rel_widths = c(1.1, 1, 0.8),
#                                       labels = c("(a)", "(b)", NA, "(c)", "(d)"),
#                                       label_size = 20, label_x = 0.2, label_y = 0.95)
# convergence_all
# # exporting
# pdf(file = file.path("plots", "convergence_plot_all.pdf"), width = 11.5, height = 9)
# convergence_all
# dev.off()


# # joining together
# predict_converge_disp_all <- rbind(
#   m.converge.predict, m.converge_animal.predict, m.converge_gravity.predict, m.converge_wind.predict
# )
# # joining with time
# predict_converge_disp_all <- predict_converge_disp_all %>%
#   left_join(scaled_time_key, by = c("x" = "s.time"))
# 
# predict_converge_disp_all$dispersal_mode <- factor(predict_converge_disp_all$dispersal_mode, levels = c("Total", "Animal", "Gravity", "Wind"))
# dispersal_mode_convergence$dispersal_mode <- factor(dispersal_mode_convergence$dispersal_mode, levels = c("Total", "Animal", "Gravity", "Wind"))
# 
# convergence_plot_all <- predict_converge_disp_all %>%
#   ggplot() +
#   geom_point(aes(time, jaccard, color = patch_pair), size = 3, alpha = 0.05, data = dispersal_mode_convergence) +
#   geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4) +
#   geom_line(aes(time, predicted, color = group), linewidth = 1.4) +
#   facet_wrap(~dispersal_mode) +
#   theme_minimal(base_size = 22) +
#   scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
#   scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
#   xlab("Time since site creation (years)") +
#   ylab(expression(paste("Spatial ", beta, " diversity (Jaccard)"))) +
#   theme(panel.spacing.x = unit(1.4, "lines")) +
#   guides(fill=guide_legend(ncol=1)) +
#   guides(color=guide_legend(ncol=1)) +
#   theme(legend.title = element_text(size = 13), 
#         legend.text = element_text(size = 12)) +
#   theme(legend.position = "right") +
#   theme(legend.justification.bottom = "left")
# convergence_plot_all
# 
# # exporting
# pdf(file = file.path("plots", "convergence_plot_all.pdf"), width = 10, height = 8)
# convergence_plot_all
# dev.off()
# 
# 
# 
