### community trajectory analysis - segement lengths
librarian::shelf(tidyverse, vegan, ecotraj, glmmTMB, DHARMa, emmeans, ggeffects, AICcmodavg) # Install missing packages and load needed libraries

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE) %>%
  #filter(!block %in% c("75W", "75E")) %>%
  filter(patch_type != "center")


# pivot to wider format
srs_data_wider <- srs_data %>%
  dplyr::count(unique_id, time, year, sppcode) %>%
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) # wide format


# make factor
srs_data_wider$time <- as.numeric(srs_data_wider$time)
srs_data_wider$unique_id <- as.factor(srs_data_wider$unique_id)
srs_data_wider$year <- as.factor(srs_data_wider$year)

# patch data
patch_info <- srs_data_wider %>% 
  arrange(unique_id, time) %>%
  select(unique_id, time, year)

# species matrix
sp_info <- srs_data_wider %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  select(!c("unique_id", "time", "year"))


# Jaccard distance matrix
jaccard_dist <- vegdist(sp_info, method = "jaccard")

# defining trajectories
srs_trajectory <- defineTrajectories(jaccard_dist, sites = patch_info$unique_id, surveys = patch_info$time)


# segment lengths of trajectories between consectutive years
segment_lengths <- trajectoryLengths(srs_trajectory)
segment_lengths <- segment_lengths %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  pivot_longer(cols = S1:S22, names_to = "time", values_to = "distance") %>%
  mutate(time = as.numeric(sub("S", "", time))) %>%
  filter(!is.na(distance)) %>%
  dplyr::select(!time)

# creating time info to join with segment lengths - some surveys were not consecutive years
time_surveys <- patch_info %>%
  filter(year != 2001) %>% # removing first year for sites created in 2000
  filter(time!= 0) # removing first survey for sites created in 2007

# joining with segment lengths
segment_lengths <- cbind(segment_lengths, time_surveys)

# segment lengths model 
m_length <- glmmTMB(distance ~ time + patch_type + (1|block/patch),
                    data = segment_lengths)
m_length_quad <- glmmTMB(distance ~ patch_type + time + I(time^2) + (1|block/patch),
                    data = segment_lengths)
#m_length <- nls(distance ~ SSasymp(time, Asym, R0, lrc), data=segment_lengths) 
m_length_null <- glmmTMB(distance ~ 1, # null model
                         data = segment_lengths)
# AIC comparison
a <- list(m_length, m_length_quad, m_length_null)
aictab(a) # quadratic much better fit


summary(m_length_quad)
#plot(simulateResiduals(m_length_quad))

# posthoc
m_length_posthoc <- emmeans(m_length_quad, ~ patch_type)
pairs(m_length_posthoc)

# prediction plot
# m_length_predict <- ggpredict(m_length, terms = c("time [all]", "patch_type"))
# m_length_predict %>%
#   ggplot() +
#   geom_point(aes(time, distance, color = patch_type), data = segment_lengths) +
#   geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
#   geom_line(aes(x, predicted, color = group), linewidth = 2) +
#   scale_color_brewer(palette = "Set2", name = "Patch Type") +
#   scale_fill_brewer(palette = "Set2", name = "Patch Type") +
#   theme_bw()



# plotting segment lengths 
segment_lengths_plot <- segment_lengths %>%
  #filter(!block %in% c("75W", "75E")) %>%
  mutate(patch_type = dplyr::case_when(
    patch_type %in% c("connected") ~ "Connected",
    patch_type %in% c("rectangle") ~ "Rectangular",
    patch_type %in% c("wing") ~ "Winged"
  )) %>%
  ggplot(aes(time, distance, color = patch_type, fill = patch_type)) +
  geom_point(size = 4.5, alpha = 0.3) +
  theme_minimal(base_size = 28) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.5, linewidth = 2) +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  
  #scale_color_brewer(palette = "Set2", name = "Patch Type") +
 # scale_fill_brewer(palette = "Set2", name = "Patch Type") +
  ylab(expression(atop("Trajectory distance", paste("between consecutive surveys")))) +
  #ylab("Trajectory distance between consecutive surveys") +
  xlab("Time since site creation (years)")
segment_lengths_plot

# pdf(file = file.path("plots", "segment_lengths.pdf"), width = 12, height = 8)
# segment_lengths_plot
# dev.off()







