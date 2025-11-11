## supplemental exploratory plots

librarian::shelf(tidyverse, vegan, ecotraj, glmmTMB, DHARMa, emmeans, ggeffects, 
                 AICcmodavg, performance, cowplot) # Install missing packages and load needed libraries

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE) %>%
  #filter(rare == 1) %>%
  #filter(!block %in% c("75W", "75E")) %>%
  filter(patch_type != "center")



### species richness over time
srs_richness <- srs_data %>%
  count(unique_id, time) 

# sp rich by dispersal mode
srs_richness_dispersal <- srs_data %>%
  count(unique_id, time, dispersal_mode) %>%
  pivot_wider(names_from = dispersal_mode, values_from = n) %>%
  select(-`NA`) %>%
  left_join(srs_richness, by = c("unique_id", "time")) %>% # joining with total richness
  rename(Total = n) %>%
  pivot_longer(cols = c("Animal", "Gravity", "Wind", "Total"), names_to = "dispersal_mode", 
               values_to = "richness") %>%
  separate(unique_id, into = c("block", "patch_rep", "patch_type"), sep = "-")

# ordering factors
srs_richness_dispersal$dispersal_mode <- factor(srs_richness_dispersal$dispersal_mode, levels = c("Total", "Animal", "Gravity", "Wind"))

## sp richness plot
richness_plot <- srs_richness_dispersal %>%
  ggplot(aes(time, richness, color = patch_type, fill = patch_type)) +
  geom_point(alpha = 0.15, size = 3) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.5, linewidth = 1.5) +
  facet_wrap(~dispersal_mode, scales = "free") +
  theme_minimal(base_size = 20) +
  xlab("Time since site creation (years)") +
  ylab("Species richness") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type", labels = c("Connected", "Rectangular", "Winged")) +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type", labels = c("Connected", "Rectangular", "Winged"))
richness_plot

pdf(file = file.path("plots", "richness_plot.pdf"), width = 10.5, height = 9.5)
richness_plot
dev.off()


## proportion of dispersal mode over time
srs_dispersal_prop <- srs_richness_dispersal %>%
  pivot_wider(names_from = dispersal_mode, values_from = richness) %>%
  mutate(prop_animal = Animal/Total,
         prop_gravity = Gravity/Total,
         prop_wind = Wind/Total) %>%
  select(-Total, -Animal, -Gravity, -Wind) %>%
  pivot_longer(cols = c("prop_animal", "prop_gravity", "prop_wind"), names_to = "dispersal_mode", values_to = "proportion") %>%
  mutate(patch_type = recode(as.factor(patch_type), "connected" = "Connected", 
                                           "rectangle" = "Rectangular", 
                                           "wing" = "Winged"))

# proportion plot over time
dispersal_prop_plot <- srs_dispersal_prop %>%
  ggplot(aes(time, proportion, color = dispersal_mode, fill = dispersal_mode)) +
  geom_point(alpha = 0.15, size = 3) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.5, linewidth = 1.5) +
  facet_wrap(~patch_type) +
  theme_minimal(base_size = 20) +
  xlab("Time since site creation (years)") +
  ylab("Proportion") +
  scale_fill_manual(values = c("#E1BE6A", "#8c510a", "#40B0A6"), name = "Dispersal mode", labels = c("Animal", "Gravity", "Wind")) +
  scale_color_manual(values = c("#E1BE6A", "#8c510a", "#40B0A6"), name = "Dispersal mode", labels = c("Animal", "Gravity", "Wind"))
dispersal_prop_plot

# pdf(file = file.path("plots", "dispersal_prop.pdf"), width = 10.5, height = 5.5)
# dispersal_prop_plot
# dev.off()





