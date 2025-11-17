########
## SCRIPT NAME: 05_supplement.R
## AUTHOR: Katherine Hulting
## PURPOSE: Supplemental exploratory plots: species richness over time, proportion of dispersal mode, turnover
## PRODUCTS: figureS2.pdf, figureS3.pdf, figureS6.pdf
#########

# loading libraries
librarian::shelf(tidyverse, vegan, ecotraj, glmmTMB, DHARMa, emmeans, ggeffects, 
                 AICcmodavg, performance, cowplot) # Install missing packages and load needed libraries

source(here::here("00_functions.R")) # loading functions

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE) %>%
  #filter(rare == 1) %>%
  #filter(!block %in% c("75W", "75E")) %>%
  filter(patch_type != "Center")



### species richness over time
srs_richness <- srs_data %>%
  count(unique_id, time) 

# species richness by dispersal mode
srs_richness_dispersal <- srs_data %>%
  count(unique_id, time, dispersal_mode) %>%
  pivot_wider(names_from = dispersal_mode, values_from = n) %>%
  select(-`NA`) %>% # removing the few observations that dont have dispersal mode data
  left_join(srs_richness, by = c("unique_id", "time")) %>% # joining with total richness
  rename(`All Species` = n) %>% # renaming the count of total richness as "all species"
  pivot_longer(cols = c("Animal", "Gravity", "Wind", "All Species"), names_to = "dispersal_mode", # pivoting longer to making column of dispersal mode
               values_to = "richness") %>%
  separate(unique_id, into = c("block", "patch_rep", "patch_type"), sep = "-") # seperating unique ID

# ordering factors
srs_richness_dispersal$dispersal_mode <- factor(srs_richness_dispersal$dispersal_mode, levels = c("All Species", "Animal", "Gravity", "Wind"))

## species richness plot
figureS2 <- srs_richness_dispersal %>%
  ggplot(aes(time, richness, color = patch_type, fill = patch_type)) +
  geom_point(alpha = 0.15, size = 3) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.5, linewidth = 1.5) + # allowing quadratic line
  facet_wrap(~dispersal_mode, scales = "free") +
  theme_minimal(base_size = 20) +
  xlab("Years since site creation") +
  ylab("Species richness") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type")
figureS2

# exporting
# pdf(file = file.path("plots", "figureS2.pdf"), width = 10.5, height = 9.5)
# figureS2
# dev.off()


## proportion of dispersal mode over time
srs_dispersal_prop <- srs_richness_dispersal %>%
  pivot_wider(names_from = dispersal_mode, values_from = richness) %>%
  mutate(prop_animal = Animal/`All Species`, # calculating proportion of each dispersal mode by dividing the sp richnes of that dispersal mode by the total species richness at that time
         prop_gravity = Gravity/`All Species`,
         prop_wind = Wind/`All Species`) %>%
  select(-`All Species`, -Animal, -Gravity, -Wind) %>%
  pivot_longer(cols = c("prop_animal", "prop_gravity", "prop_wind"), names_to = "dispersal_mode", values_to = "proportion")

# proportion plot over time
figureS3 <- srs_dispersal_prop %>%
  ggplot(aes(time, proportion, color = dispersal_mode, fill = dispersal_mode)) +
  geom_point(alpha = 0.15, size = 3) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.5, linewidth = 1.5) +
  facet_wrap(~patch_type) +
  theme_minimal(base_size = 20) +
  xlab("Years since site creation") +
  ylab("Proportion") +
  scale_fill_manual(values = c("#E1BE6A", "#8c510a", "#40B0A6"), name = "Dispersal mode", labels = c("Animal", "Gravity", "Wind")) +
  scale_color_manual(values = c("#E1BE6A", "#8c510a", "#40B0A6"), name = "Dispersal mode", labels = c("Animal", "Gravity", "Wind"))
figureS3

# exporting
# pdf(file = file.path("plots", "figureS3.pdf"), width = 10.5, height = 5.5)
# figureS3
# dev.off()




##### compositional changes #####
# calculating raw changes in composition: # of gains, # of losses, # of species staying the same from one year to the next
species_changes <- srs_data %>%
    count(block, patch, patch_type, unique_id, year, time, sppcode) %>%
    group_by(block, patch) %>%
    group_split() %>%
    lapply(compute_composition_change) %>%
    bind_rows() %>% # putting together into a dataframe
    rownames_to_column("unique_id") %>%
    separate(unique_id, into = c("block", "patch_rep", "patch_type", "time2"), sep = "-") %>%
    dplyr::select(-time2) # removing duplicate time column

# removing the first year 
species_changes <- species_changes %>%
  filter(!time %in% c("0")) %>% # removing time 0 - only recorded for 52 and 57, first year of sampling so no turnover from previous year
  filter(!(time == 1 & block %in% c("08", "10", "53N", "53S", "54N", "54S", "75E", "75W"))) %>% # first year of sampling for these blocks so no turnover from previous year
  pivot_longer(5:8, names_to = "type", values_to = "change")


# all together
# species_changes %>%
#   ggplot(aes(time, change, color = patch_type, fill = patch_type)) +
#   geom_point(alpha = 0.15, size = 3) +
#   geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.5, linewidth = 1.5) +
#   theme_minimal(base_size = 20) +
#   facet_wrap(~type, scales = "free") +
#   ylab(expression(atop("Number of gains and losses", paste("between consecutive surveys")))) +
#   xlab("Years since site creation") +
#   scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
#   scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
#   theme(legend.position = "none")


# amount of species turnover
changed_total_plot <- species_changes %>%
  filter(type %in% c("changed_total")) %>%
  ggplot(aes(time, change, color = patch_type, fill = patch_type)) +
  geom_point(alpha = 0.15, size = 3) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.5, linewidth = 1.5) +
  theme_minimal(base_size = 20) +
  ylim(0, 150) +
  ylab(expression(atop("Number of gains and losses", paste("between consecutive surveys")))) +
  xlab("Years since site creation") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
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
  xlab("Years since site creation") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type", labels = c("Connected", "Rectangular", "Winged")) +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type", labels = c("Connected", "Rectangular", "Winged"))
stayed_present_plot

# all together
figureS6 <- cowplot::plot_grid(changed_total_plot, stayed_present_plot, rel_widths = c(1, 1.44),
                                       labels = c("(A)", "(B)"), label_size = 16)
figureS6


# exporting
# pdf(file = file.path("plots", "figureS6.pdf"), width = 12, height = 5)
# figureS6
# dev.off()











