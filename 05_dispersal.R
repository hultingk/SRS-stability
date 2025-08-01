### looking at dispersal mode info
librarian::shelf(tidyverse, vegan, ecotraj, glmmTMB, DHARMa, emmeans, ggeffects, AICcmodavg) # Install missing packages and load needed libraries

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE) %>%
  #filter(rare == 1) %>%
  #filter(!block %in% c("75W", "75E")) %>%
  filter(patch_type != "center")

srs_dispersal_prop <- srs_data %>%
  count(block, patch, patch_type, unique_id, year, time, dispersal_mode) %>%
  pivot_wider(names_from = dispersal_mode, values_from = n) %>%
  dplyr::select(!c("NA")) %>%
  mutate(sp_richness = Animal + Gravity + Wind) %>%
  mutate(prop_animal = Animal/sp_richness,
         prop_gravity = Gravity/sp_richness,
         prop_wind = Wind/sp_richness)

srs_dispersal_prop <- srs_dispersal_prop %>%
  pivot_longer(cols = 11:13, names_to = "dispersal_mode", values_to = "proportion") 

m1 <- glmmTMB(proportion ~ patch_type*time*dispersal_mode + (1|block),
              data = srs_dispersal_prop,
              family = "binomial",
              weights = sp_richness)
summary(m1)
plot(simulateResiduals(m1))


srs_dispersal_prop %>%
  ggplot(aes(time, proportion, color = patch_type, fill = patch_type)) +
  facet_wrap(~dispersal_mode) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  theme_minimal(base_size = 12) +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") 

srs_dispersal_prop %>%
  ggplot(aes(time, sp_richness, color = patch_type, fill = patch_type)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  theme_minimal(base_size = 12) +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") 

