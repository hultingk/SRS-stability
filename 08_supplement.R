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

# pdf(file = file.path("plots", "richness_plot.pdf"), width = 10.5, height = 9.5)
# richness_plot
# dev.off()


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



#### change in richness between consecutive years
richness_change <- srs_richness %>%
  arrange(unique_id, time) %>%                              
  group_by(unique_id) %>%                                  
  mutate(
    richness_change = n - lag(n),
    year_diff = time - lag(time)                   
  ) %>%
  ungroup() %>%
  separate(unique_id, into = c("block", "patch_rep", "patch_type"), sep = "-")


richness_change %>%
  ggplot(aes(time, richness_change, color = patch_type, fill = patch_type)) +
  geom_point(alpha = 0.15, size = 3) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.5, linewidth = 1.5) +
  theme_minimal(base_size = 20) +
  xlab("Time since site creation (years)") +
  ylab("Change in richness between consecutive years") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type", labels = c("Connected", "Rectangular", "Winged")) +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type", labels = c("Connected", "Rectangular", "Winged"))




# raw composition change between consecutive years


compute_composition_change <- function(df) {
  # convert data to correct format
  df_wide <- df %>% 
    pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>% # wide format
    arrange(unique_id, time) %>%  # Ensure years are sorted properly
    mutate(unique_id = paste(unique_id, time, sep = "-")) %>%
    dplyr::select(-block, -patch_type, -patch, -year) %>%
    column_to_rownames("unique_id")
  
  # calculate losses and gains between consecutive years
  turnover <- df_wide %>%
    arrange(time) %>%
    mutate(
      gains = rowSums(across(-time, ~ (.x == 1 & lag(.x) == 0)), na.rm = TRUE),     # 0 → 1
      losses = rowSums(across(-time, ~ (.x == 0 & lag(.x) == 1)), na.rm = TRUE),    # 1 → 0
      stayed_present = rowSums(across(-time, ~ (.x == 1 & lag(.x) == 1)), na.rm = TRUE),  # 1 → 1
      changed_total = gains + losses
    ) %>%
    select(time, gains, losses, stayed_present, changed_total)
  

  return(turnover)
}


species_changes <- srs_data %>%
    count(block, patch, patch_type, unique_id, year, time, sppcode) %>%
    group_by(block, patch) %>%
    group_split() %>%
    lapply(compute_composition_change) %>%
    bind_rows() %>% # putting together into a dataframe
    rownames_to_column("unique_id") %>%
    separate(unique_id, into = c("block", "patch_rep", "patch_type"), sep = "-")

species_changes <- species_changes %>%
  filter(!time %in% c("0", "1")) %>%
  pivot_longer(5:8, names_to = "type", values_to = "change")

gains <- species_changes %>%
  filter(type %in% c("gains"))
m1 <- glmmTMB(change ~ patch_type * time + patch_type * I(time^2) + (1|block),
              data = gains,
              family = "nbinom2")
summary(m1)
losses <- species_changes %>%
  filter(type %in% c("losses"))

# all together
species_changes %>%
  ggplot(aes(time, change, color = patch_type, fill = patch_type)) +
  geom_point(alpha = 0.15, size = 3) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.5, linewidth = 1.5) +
  theme_minimal(base_size = 20) +
  facet_wrap(~type, scales = "free") +
  ylab(expression(atop("Number of gains and losses", paste("between consecutive surveys")))) +
  xlab("Time since site creation (years)") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type", labels = c("Connected", "Rectangular", "Winged")) +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type", labels = c("Connected", "Rectangular", "Winged")) +
  theme(legend.position = "none")


# amount of species turnoever
changed_total_plot <- species_changes %>%
  filter(type %in% c("changed_total")) %>%
  ggplot(aes(time, change, color = patch_type, fill = patch_type)) +
  geom_point(alpha = 0.15, size = 3) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.5, linewidth = 1.5) +
  theme_minimal(base_size = 20) +
  ylim(0, 150) +
  ylab(expression(atop("Number of gains and losses", paste("between consecutive surveys")))) +
  xlab("Time since site creation (years)") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type", labels = c("Connected", "Rectangular", "Winged")) +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type", labels = c("Connected", "Rectangular", "Winged")) +
  theme(legend.position = "none")
changed_total_plot

# number of consistent species 
stayed_present_plot <- species_changes %>%
  filter(type %in% c("stayed_present")) %>%
  ggplot(aes(time, change, color = patch_type, fill = patch_type)) +
  geom_point(alpha = 0.15, size = 3) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.5, linewidth = 1.5) +
  theme_minimal(base_size = 20) +
  ylim(0, 150) +
  ylab(expression(atop("Number of species consistent", paste("between consecutive surveys")))) +
  xlab("Time since site creation (years)") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type", labels = c("Connected", "Rectangular", "Winged")) +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type", labels = c("Connected", "Rectangular", "Winged"))
stayed_present_plot

# all together
comp_change_plot <- cowplot::plot_grid(changed_total_plot, stayed_present_plot, rel_widths = c(1, 1.44),
                                       labels = c("(A)", "(B)"), label_size = 16)
comp_change_plot

# pdf(file = file.path("plots", "comp_change_plot.pdf"), width = 12, height = 5)
# comp_change_plot
# dev.off()











