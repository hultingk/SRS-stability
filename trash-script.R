# loading libraries
# Load `librarian` package
library(librarian)
# Install missing packages and load needed libraries
shelf(tidyverse, vegan, codyn, glmmTMB, ggeffects, performance, DHARMa, emmeans, ecotraj)

# loading data from EDI 
srs_data_url <- "https://pasta.lternet.edu/package/data/eml/edi/414/1/2429e7fc1b33cefb59bab8451aaa8327"
srs_data_raw <- read.csv(file = srs_data_url)

# loading species info from EDI
srs_dispersal_url <- "https://portal.edirepository.org/nis/dataviewer?packageid=edi.414.1&entityid=8de4a490a6ac6b05d2406c975d25b649"
srs_dispersal <- read.csv(file = srs_dispersal_url)

# loading site info about replicate patches
site_info <- read_csv("site_info.csv")


srs_data_raw <- srs_data_raw %>%
  rowwise() %>% # grouping by rows
  mutate(occurance = sum(across(starts_with("X")), na.rm = T)) # summing across rows - if there is a 1 in any of the plots the species occurred in that patch

srs_data <- srs_data_raw %>%
  select(!X1:X8) %>% # removing plot info - redundant
  filter(occurance == 1) %>% # only including observations where the species was present
  left_join(site_info, by = c("site", "EU", "PatchType", "Year")) %>% # joining patch rep info
  mutate(unique_ID = paste(EU, patch_rep, PatchType, sep = "_")) %>% # creating unique ID for each rep
  count(unique_ID, Time, SppCode) %>% # summarizing what species occured in each patch x year
  mutate(Year = Time) %>% # renaming
  select(!c("Time"))

srs_data_split <- srs_data %>% # splitting dataframe into each patch as a list for function
  group_by(unique_ID) %>%
  group_split()

##### Jaccard ####
# writing function to calculate jaccard's dissimilarity iteratively between consectutive years
compute_jaccard <- function(df) {
  # convert data to correct format
  df_wide <- df %>% 
    pivot_wider(names_from = SppCode, values_from = n, values_fill = 0) %>% # wide format
    dplyr::select(!c("unique_ID")) %>% # remove unique ID column
    arrange(Year) %>% # Ensure years are sorted properly
    column_to_rownames("Year") #convert years to rownames
  
  # compute jaccard dissimilarity iteratively between consecutive years
  jaccard_values <- sapply(1:(nrow(df_wide) - 1), function(i) {
    vegdist(df_wide[i:(i+1), ], method = "jaccard")
  })
  
  # store results with year pairs
  result <- data.frame(
    unique_ID = unique(df$unique_ID),
    year_pair = paste(sort(unique(df$Year))[-length(unique(df$Year))], 
                      sort(unique(df$Year))[-1], sep = " - "),
    jaccard_dissimilarity = jaccard_values
  )
  
  return(result)
}


jaccard_results <- srs_data_split %>% # applying function to data
  lapply(compute_jaccard) %>%
  bind_rows() # putting together into a dataframe



year_factor <- jaccard_results %>% # converting year pairs into a usable format for modeling/plotting
  count(year_pair) %>%
  mutate(time = as.factor(year_pair)) %>% # converting to factors
  mutate(time = recode_factor(time, `1 - 2` = "2", `2 - 3` = "3", `3 - 4` = "4", `4 - 5` = "5",
                `4 - 6` = "6", `4 - 8` = "8", `5 - 6` = "6", `6 - 7` = "7",
                `6 - 8` = "8", `7 - 8` = "8", `8 - 9` = "9", `9 - 10` = "10",
                `10 - 11` = "11", `11 - 12` = "12", `12 - 13` = "13", `13 - 14` = "14",
                `13 - 15` = "15", `14 - 15` = "15", `15 - 16` = "16", `16 - 17` = "17",
                `17 - 18` = "18")) # reassigning each year pair as the later year

jaccard_results <- jaccard_results %>%
  separate(unique_ID, into = c("EU", "patch_rep", "patch"), sep = "_", remove = F) %>% # seperating unique ID into columns
  mutate(year_pair = as.factor(year_pair)) %>% # making a factor to prepare for joining
  left_join(year_factor, by = c("year_pair" = "year_pair")) %>% # joining to year info
  mutate(time = as.numeric(time)) # making time numeric


#site_info <- srs_data_raw %>%
  #filter(EU == "EU53S") %>%
 # count(site, EU, PatchType, Year)
#write_csv(site_info, file = "site_info.csv")

jaccard_plot <- jaccard_results %>% # plotting jaccard dissimilarity across year 
  ggplot(aes(time, jaccard_dissimilarity, color = patch)) + 
  facet_wrap(~patch) + 
  geom_point(size = 2, alpha = 0.7) + 
  #geom_line(aes(time, jaccard_dissimilarity, group = unique_ID), linewidth = 1, alpha = 0.3) +
  theme_bw() + 
  scale_color_brewer(palette = "Set2") +
  stat_smooth(method = "lm", se = F, linewidth = 3)
jaccard_plot

pdf(file = "Figure1.pdf", width = 8, height = 5)
jaccard_plot
dev.off()

hist(jaccard_results$jaccard_dissimilarity) # histogram of dissimilarity values

# how does dissimilarity change across time for each patch type?
m1 <- glmmTMB(jaccard_dissimilarity ~ patch * time + (1|EU/patch_rep), 
              data = jaccard_results)
summary(m1) # model summary
plot(simulateResiduals(m1)) # looks not the best but not the worst
check_model(m1) # okay for now

m1_posthoc <- emtrends(m1, specs = pairwise ~ patch, var = "time") # posthoc test for differences between slopes
m1_posthoc # connected different than rectangular, but not winged



#### nestedness and turnover ####
compute_jaccard2 <- function(df) {
  # convert data to correct format
  df_wide <- df %>% 
    pivot_wider(names_from = SppCode, values_from = n, values_fill = 0) %>% # wide format
    dplyr::select(!c("unique_ID")) %>% # remove unique ID column
    arrange(Year) %>% # Ensure years are sorted properly
    column_to_rownames("Year") #convert years to rownames
  
  # compute jaccard dissimilarity, nestedness, and turnover iteratively between consecutive years
  jaccard_values2 <- sapply(1:(nrow(df_wide) - 1), function(i) {
    nestedbetajac(df_wide[i:(i+1), ])
  })
  
  # store results with year pairs
  result <- data.frame(
    unique_ID = unique(df$unique_ID),
    year_pair = paste(sort(unique(df$Year))[-length(unique(df$Year))], 
                      sort(unique(df$Year))[-1], sep = " - "),
    jaccard_dissimilarity = jaccard_values2[row.names(jaccard_values2) %in% c("jaccard"),],
    turnover_values = jaccard_values2[row.names(jaccard_values2) %in% c("turnover"),],
    nestedness_values = jaccard_values2[row.names(jaccard_values2) %in% c("nestedness"),]
  )
  
  return(result)
}


jaccard_results2 <- srs_data_split %>% # applying function to data
  lapply(compute_jaccard2) %>%
  bind_rows() # putting together into a dataframe


jaccard_results2 <- jaccard_results2 %>%
  separate(unique_ID, into = c("EU", "patch_rep", "patch"), sep = "_", remove = T) %>% # seperating unique ID into columns
  mutate(year_pair = as.factor(year_pair)) %>% # making a factor to prepare for joining
  left_join(year_factor, by = c("year_pair" = "year_pair")) %>% # joining to year info
  mutate(time = as.numeric(time)) # making time numeric


colors <- c("Jaccard dissimilarity" = "#8DA0CB", "Turnover" = "#FC8D62", "Nestedness" = "#66C2A5")
jaccard_plot2 <- jaccard_results2 %>% # plotting jaccard dissimilarity across year 
  ggplot() + 
  facet_wrap(~patch) + 
  geom_point(aes(time, nestedness_values, color = "Nestedness"), size = 2, alpha = 0.7) + 
  geom_point(aes(time, turnover_values, color = "Turnover"), size = 2, alpha = 0.7) + 
  geom_point(aes(time, jaccard_dissimilarity, color = "Jaccard dissimilarity"), size = 2, alpha = 0.7) + 
  theme_bw() + 
  stat_smooth(aes(time, nestedness_values, color = "Nestedness"), method = "lm", se = F, linewidth = 3) +
  stat_smooth(aes(time, turnover_values, color = "Turnover"), method = "lm", se = F, linewidth = 3)+ 
  stat_smooth(aes(time, jaccard_dissimilarity, color = "Jaccard dissimilarity"), method = "lm", se = F, linewidth = 3) +
  labs(x = "Time",
       y = "",
       color = "Legend") +
  scale_color_manual(values = colors)
jaccard_plot2

#pdf(file = "Figure2.pdf", width = 10, height = 5)
#jaccard_plot2
#dev.off()

# how does dissimilarity change across time for each patch type?
m2 <- glmmTMB(turnover_values ~ patch * time + (1|EU/patch_rep), 
              data = jaccard_results2)
summary(m2) # model summary
plot(simulateResiduals(m2)) # looks not the best but not the worst
check_model(m2) # okay for now

m2_posthoc <- emtrends(m2, pairwise ~ patch, var = "time") # posthoc test for differences between slopes
m2_posthoc # connected different than rectangular, but not winged



#### wind dispersed ####
# looking at only wind dispersed species
wind_jaccard <- srs_data %>%
  left_join(srs_dispersal, by = c("SppCode")) %>%
  filter(DispersalMode == "Wind") %>%
  dplyr::select(!c("Genus", "Species", "DispersalMode", "LongleafPineSpecies")) %>%
  group_by(unique_ID) %>%
  group_split() %>%
  lapply(compute_jaccard) %>%
  bind_rows() # putting together into a dataframe

wind_jaccard <- wind_jaccard %>%
  separate(unique_ID, into = c("EU", "patch_rep", "patch"), sep = "_") %>% # seperating unique ID into columns
  mutate(year_pair = as.factor(year_pair)) %>% # making a factor to prepare for joining
  left_join(year_factor, by = c("year_pair" = "year_pair")) %>% # joining to year info
  mutate(time = as.numeric(time)) # making time numeric

wind_plot <- wind_jaccard %>% # plotting jaccard dissimilarity across year 
  ggplot(aes(time, jaccard_dissimilarity, color = patch)) + 
  facet_wrap(~patch) + 
  geom_point(size = 2, alpha = 0.7) + 
  theme_bw() + 
  scale_color_brewer(palette = "Set2") +
  ylab("jaccard_dissimilarity, wind dispered") +
  stat_smooth(method = "lm", se = F, linewidth = 3)
wind_plot

# how does dissimilarity change across time for each patch type?
m.wind <- glmmTMB(jaccard_dissimilarity ~ patch * time + (1|EU/patch_rep), 
              data = wind_jaccard)
summary(m.wind) # model summary
plot(simulateResiduals(m.wind)) # looks not the best but not the worst
check_model(m.wind) # okay for now

m.wind_posthoc <- emtrends(m.wind, pairwise ~ patch, var = "time") # posthoc test for differences between slopes
m.wind_posthoc # connected different than rectangular, but not winged

#### gravity dispersed ####
# looking at only gravity dispersed species
gravity_jaccard <- srs_data %>%
  left_join(srs_dispersal, by = c("SppCode")) %>%
  filter(DispersalMode == "Gravity") %>%
  dplyr::select(!c("Genus", "Species", "DispersalMode", "LongleafPineSpecies")) %>%
  group_by(unique_ID) %>%
  group_split() %>%
  lapply(compute_jaccard) %>%
  bind_rows() # putting together into a dataframe

gravity_jaccard <- gravity_jaccard %>%
  separate(unique_ID, into = c("EU", "patch_rep", "patch"), sep = "_") %>% # seperating unique ID into columns
  mutate(year_pair = as.factor(year_pair)) %>% # making a factor to prepare for joining
  left_join(year_factor, by = c("year_pair" = "year_pair")) %>% # joining to year info
  mutate(time = as.numeric(time)) # making time numeric

gravity_plot <- gravity_jaccard %>% # plotting jaccard dissimilarity across year 
  ggplot(aes(time, jaccard_dissimilarity, color = patch)) + 
  facet_wrap(~patch) + 
  geom_point(size = 2, alpha = 0.7) + 
  theme_bw() + 
  scale_color_brewer(palette = "Set2") +
  ylab("jaccard_dissimilarity, gravity dispered") +
  stat_smooth(method = "lm", se = F, linewidth = 3)
gravity_plot

# how does dissimilarity change across time for each patch type?
m.gravity <- glmmTMB(jaccard_dissimilarity ~ patch * time + (1|EU/patch_rep), 
                  data = gravity_jaccard)
summary(m.gravity) # model summary
plot(simulateResiduals(m.gravity)) # looks not the best but not the worst
check_model(m.gravity) # okay for now

m.gravity_posthoc <- emtrends(m.gravity, pairwise ~ patch, var = "time") # posthoc test for differences between slopes
m.gravity_posthoc # connected different than rectangular, but not winged


#### animal dispersed ####
# looking at only animal dispersed species
animal_jaccard <- srs_data %>%
  left_join(srs_dispersal, by = c("SppCode")) %>%
  filter(DispersalMode == "Animal") %>%
  dplyr::select(!c("Genus", "Species", "DispersalMode", "LongleafPineSpecies")) %>%
  group_by(unique_ID) %>%
  group_split() %>%
  lapply(compute_jaccard) %>%
  bind_rows() # putting together into a dataframe

animal_jaccard <- animal_jaccard %>%
  separate(unique_ID, into = c("EU", "patch_rep", "patch"), sep = "_") %>% # seperating unique ID into columns
  mutate(year_pair = as.factor(year_pair)) %>% # making a factor to prepare for joining
  left_join(year_factor, by = c("year_pair" = "year_pair")) %>% # joining to year info
  mutate(time = as.numeric(time)) # making time numeric

animal_plot <- animal_jaccard %>% # plotting jaccard dissimilarity across year 
  ggplot(aes(time, jaccard_dissimilarity, color = patch)) + 
  facet_wrap(~patch) + 
  geom_point(size = 2, alpha = 0.7) + 
  theme_bw() + 
  scale_color_brewer(palette = "Set2") +
  ylab("jaccard_dissimilarity, animal dispered") +
  stat_smooth(method = "lm", se = F, linewidth = 3)
animal_plot

# how does dissimilarity change across time for each patch type?
m.animal <- glmmTMB(jaccard_dissimilarity ~ patch * time + (1|EU/patch_rep), 
                     data = animal_jaccard)
summary(m.animal) # model summary
plot(simulateResiduals(m.animal)) # looks not the best but not the worst
check_model(m.animal) # okay for now

m.animal_posthoc <- emtrends(m.animal, pairwise ~ patch, var = "time") # posthoc test for differences between slopes
m.animal_posthoc # connected different than rectangular, but not winged


# putting plots together
cowplot::plot_grid(jaccard_plot, wind_plot, gravity_plot, animal_plot, 
                   label_size =30, nrow=2, ncol=2, label_x = 0.11, label_y = 0.92, align = "hv")

pdf(file = "Figure3.pdf", width = 12, height = 8)
cowplot::plot_grid(jaccard_plot, wind_plot, gravity_plot, animal_plot, 
                   label_size =30, nrow=2, ncol=2, label_x = 0.11, label_y = 0.92, align = "hv")
dev.off()


#### longleaf species ####
# looking at only longleaf pine species
longleaf_jaccard <- srs_data %>%
  left_join(srs_dispersal, by = c("SppCode")) %>%
  filter(LongleafPineSpecies == 1) %>%
  dplyr::select(!c("Genus", "Species", "DispersalMode", "LongleafPineSpecies")) %>%
  group_by(unique_ID) %>%
  group_split() %>%
  lapply(compute_jaccard) %>%
  bind_rows() # putting together into a dataframe

longleaf_jaccard <- longleaf_jaccard %>%
  separate(unique_ID, into = c("EU", "patch_rep", "patch"), sep = "_") %>% # seperating unique ID into columns
  mutate(year_pair = as.factor(year_pair)) %>% # making a factor to prepare for joining
  left_join(year_factor, by = c("year_pair" = "year_pair")) %>% # joining to year info
  mutate(time = as.numeric(time)) # making time numeric

longleaf_plot <- longleaf_jaccard %>% # plotting jaccard dissimilarity across year 
  ggplot(aes(time, jaccard_dissimilarity, color = patch)) + 
  facet_wrap(~patch) + 
  geom_point() + 
  theme_bw() + 
  scale_color_brewer(palette = "Set2") +
  ylab("jaccard_dissimilarity, longleaf species") +
  stat_smooth(method = "lm", se = F, linewidth = 2)
longleaf_plot

# how does dissimilarity change across time for each patch type?
m.longleaf <- glmmTMB(jaccard_dissimilarity ~ patch * time + (1|EU/patch_rep), 
                    data = longleaf_jaccard)
summary(m.longleaf) # model summary
plot(simulateResiduals(m.longleaf)) # looks not the best but not the worst
check_model(m.longleaf) # okay for now

m.longleaf_posthoc <- emtrends(m.longleaf, pairwise ~ patch, var = "time") # posthoc test for differences between slopes
m.longleaf_posthoc # connected different than rectangular, but not winged




#### common species ####
# repeating but excluding rare species (species with less than 10 detections over time)
common_spp <- srs_data %>%
  count(SppCode) %>%
  arrange(n) %>%
  filter(n > 10) %>%
  mutate(rare = 0) %>%
  select(!c("n"))

common_spp_jaccard <- srs_data %>%
  left_join(common_spp, by = c("SppCode")) %>%
  filter(rare == 0) %>%
  dplyr::select(!c("rare")) %>%
  group_by(unique_ID) %>%
  group_split() %>%
  lapply(compute_jaccard) %>%
  bind_rows() # putting together into a dataframe

common_spp_jaccard <- common_spp_jaccard %>%
  separate(unique_ID, into = c("EU", "patch_rep", "patch"), sep = "_") %>% # seperating unique ID into columns
  mutate(year_pair = as.factor(year_pair)) %>% # making a factor to prepare for joining
  left_join(year_factor, by = c("year_pair" = "year_pair")) %>% # joining to year info
  mutate(time = as.numeric(time)) # making time numeric

common_spp_plot <- common_spp_jaccard %>% # plotting jaccard dissimilarity across year 
  ggplot(aes(time, jaccard_dissimilarity, color = patch)) + 
  facet_wrap(~patch) + 
  geom_point() + 
  theme_bw() + 
  scale_color_brewer(palette = "Set2") +
  ylab("jaccard_dissimilarity, common species") +
  stat_smooth(method = "lm", se = F, linewidth = 2)
common_spp_plot

# how does dissimilarity change across time for each patch type?
m.common_spp <- glmmTMB(jaccard_dissimilarity ~ patch * time + (1|EU/patch_rep), 
                      data = common_spp_jaccard,
                      family = "beta_family")
summary(m.common_spp) # model summary
plot(simulateResiduals(m.common_spp)) # looks not the best but not the worst
check_model(m.common_spp) # okay for now

m.common_spp_posthoc <- emtrends(m.common_spp, pairwise ~ patch, var = "time") # posthoc test for differences between slopes
m.common_spp_posthoc # connected different than rectangular, but not winged






#### ordination ####

cta_trial <- srs_data %>%
  filter(unique_ID %in% c("EU53S_B_C")) 
unique_sites <- unique(cta_trial$unique_ID)

cta_trial_sites <- cta_trial %>%
  count(unique_ID, Year)

cta_trial_wider <- cta_trial %>%
  pivot_wider(names_from = "SppCode", values_from = n, values_fill = 0) %>%
  dplyr::select(!c("unique_ID")) %>% # remove unique ID column
  arrange(Year) %>% # Ensure years are sorted properly
  column_to_rownames("Year") #convert years to rownames


d1 <- vegan::vegdist(cta_trial_wider, "jaccard")

#segment_lengths <- trajectoryLengths(d1, sites = cta_trial_sites$unique_ID, surveys = cta_trial_sites$Year)
#segment_length_dataframe <- data.frame(
#  segment = seq_along(segment_lengths),
#  length = segment_lengths
#)
#color <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728")
#pdf(file = "53S_E_W.pdf", width = 6, height = 5)
trajectoryPCoA(d1, sites = cta_trial_sites$unique_ID, surveys = cta_trial_sites$Year,
               survey.labels = T, length = 0.1, lwd = 2)
#dev.off()


#####
 

calculate_lengths <- function(df) {
  df_wide <- df %>% 
    pivot_wider(names_from = "SppCode", values_from = n, values_fill = 0) %>%
    dplyr::select(!c("unique_ID")) %>% # remove unique ID column
    arrange(Year) %>% # Ensure years are sorted properly
    column_to_rownames("Year") #convert years to rownames
  
  site_names <- df %>%
    count(unique_ID, Year)
  
  distances <- vegan::vegdist(df_wide, "jaccard")
  
  segment_lengths <- trajectoryLengths(distances, 
                                       sites = site_names$unique_ID, 
                                       surveys = site_names$Year)
  
  segment_lengths <- segment_lengths %>%
    pivot_longer(cols = S1:Trajectory,
                 names_to = "time_step",
                 values_to = "length") %>%
    mutate(unique_ID = unique(df$unique_ID))

}

cta_lengths <- srs_data_split %>%
  lapply(calculate_lengths) %>%
  bind_rows()

cta_trajectory_length <- cta_lengths %>%
  filter(time_step %in% c("Trajectory")) %>%
  separate(unique_ID, into = c("EU", "patch_rep", "patch"), sep = "_")

cta_tajectory_timestep <- cta_lengths %>%
  filter(!time_step %in% c("Trajectory")) %>%
  separate(unique_ID, into = c("EU", "patch_rep", "patch"), sep = "_") %>%
  mutate(time = as.numeric(str_sub(time_step, 2)))


m.trajectory <- glmmTMB(length ~ patch + (1|EU/patch_rep),
                        data = cta_trajectory_length, 
                        family = "gaussian")
summary(m.trajectory)
plot(simulateResiduals(m.trajectory))

m.trajectory.posthoc <- emmeans(m.trajectory, specs = "patch")
pairs(m.trajectory.posthoc)

m.trajectory.predict <- ggpredict(m.trajectory, terms=c("patch"))
m.trajectory.predict %>%
  ggplot() +
  geom_point(aes(x = x, y = predicted), size = 4.5, data = m.trajectory.predict,  position = position_dodge(0.5)) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),data = m.trajectory.predict, width = 0.4, linewidth = 2,  position = position_dodge(0.5)) +
  theme_classic() +
  geom_jitter(aes(x = patch, y = length), data = cta_trajectory_length, alpha = 0.2, size = 3.8, width = 0.2, height = 0) 

cta_trajectory_length %>%
  ggplot() +
  theme_bw() +
  geom_boxplot(aes(patch, length)) +
  geom_jitter(aes(patch, length))

m.trajectory.time <- glmmTMB(length ~ patch*time + (1|EU/patch_rep),
                             data = cta_tajectory_timestep, 
                             family = "beta_family")
summary(m.trajectory.time)
plot(simulateResiduals(m.trajectory.time))

cta_tajectory_timestep %>% 
  ggplot(aes(time, length, color = patch)) + 
  facet_wrap(~patch) + 
  geom_point() + 
  theme_bw() + 
  scale_color_brewer(palette = "Set2") +
  ylab("length") +
  stat_smooth(method = "lm", se = F, linewidth = 2)






#site_names <- unique(srs_data$unique_ID)

#cta_data <- srs_data_split %>%
 # lapply(prepare_matrix)
#names(cta_data) <- site_names
#cta_data_dist <- lapply(vegdist(cta_data, "jaccard"))

#d1 <- vegan::vegdist(cta_trial_wider, "jaccard")
#names(cta_data["EU08_B_C"])

#segment_lengths <- trajectoryLengths(d1, sites = cta_trial_sites$unique_ID, surveys = cta_trial_sites$Year)
#segment_length_dataframe <- data.frame(
#  segment = seq_along(segment_lengths),
#  length = segment_lengths
#)




#### stability of sp richness ####
srs_data %>%
  filter(unique_ID == "EU08_B_C") %>%
  pivot_wider(names_from = "SppCode", values_from = "n", values_fill = 0)

srs_data_wider <- srs_data %>%
  pivot_wider(names_from = "SppCode", values_from = "n", values_fill = 0)

library(plyr)
srs_richness <- ddply(srs_data_wider, c("unique_ID", "Year"), function(x) {
     data.frame(RICHNESS=sum(x[-2]>0))
   })

srs_richness_plot <- srs_richness %>%
  separate(unique_ID, into = c("EU", "patch_rep", "patch"), sep = "_") # seperating unique ID into columns

srs_richness_plot %>% 
  ggplot(aes(Year, RICHNESS, color = patch)) +
  geom_point() +
  stat_smooth(method = "lm", se = F, linewidth = 2)
summary(glmmTMB(RICHNESS ~ patch * Year + (1|EU/patch_rep), 
        data = srs_richness_plot,
        family = "poisson"))


# coefficient of variation
CV<-function(x){
  return(sd(x,na.rm=T)/mean(x,na.rm=T))
}

# stability
stability <- function(x){
  1/CV(x)
}

srs_cv <- srs_richness %>%
  dplyr::group_by(unique_ID) %>%
  dplyr::summarise(cv = CV(RICHNESS),
                   stability = stability(RICHNESS))

srs_cv <- srs_cv %>%
  separate(unique_ID, into = c("EU", "patch_rep", "patch"), sep = "_") # seperating unique ID into columns

m_stability <- glmmTMB(stability ~ patch + (1|EU/patch_rep),
                       data = srs_cv, 
                       family = "gaussian")
summary(m_stability)
check_model(m_stability)


##########
# playing around with some other metrics
srs_turnover <- turnover(df=srs_data,  
                      time.var = "Year", 
                      species.var = "SppCode",
                      abundance.var = "n", 
                      replicate.var="unique_ID")
srs_appearances <- turnover(df=srs_data,  
                         time.var = "Year", 
                         species.var = "SppCode",
                         abundance.var = "n", 
                         replicate.var="unique_ID",
                         metric="appearance")
srs_disappearances <- turnover(df=srs_data,  
                            time.var = "Year", 
                            species.var = "SppCode",
                            abundance.var = "n", 
                            replicate.var="unique_ID",
                            metric="disappearance")
srs_turnover <- srs_turnover %>%
  left_join(srs_appearances, by = c("Year", "unique_ID")) %>%
  left_join(srs_disappearances, by = c("Year", "unique_ID")) %>%
  separate(unique_ID, into = c("EU", "patch_rep", "patch"), sep = "_")
  
srs_turnover %>%
  ggplot() +
  geom_point(aes(Year, total, color = patch))
srs_turnover %>%
  ggplot() +
  geom_point(aes(Year, disappearance, color = patch))
srs_turnover %>%
 # separate(unique_ID, into = c("EU", "patch_rep", "patch"), sep = "_") %>%
  ggplot() +
  geom_point(aes(Year, appearance, color = patch))

m1 <- glmmTMB(total ~ patch * Year + (1|EU/patch_rep),
              data = srs_turnover,
              family = "beta_family")
summary(m1)
plot(simulateResiduals(m1)) # looks good
check_model(m1)

srs_turnover %>%
  ggplot(aes(Year, appearance, color = patch)) + facet_wrap(~patch) + 
  geom_point() + theme_bw() + stat_smooth(method = "lm", se = F, size = 2)

rate.res <- rate_change(srs_data,  
                        time.var = "Year", 
                        species.var = "SppCode",
                        abundance.var = "n", 
                        replicate.var="unique_ID")
comm.res <- rate_change_interval(srs_data,  
                        time.var = "Year", 
                        species.var = "SppCode",
                        abundance.var = "n", 
                        replicate.var="unique_ID")


comm.res %>%
  separate(unique_ID, into = c("EU", "patch_rep", "patch"), sep = "_") %>%
  ggplot(aes(interval, distance, color = patch)) + facet_wrap(~patch) + 
  geom_point() + theme_bw() + stat_smooth(method = "lm", se = F, size = 2)








#### TRIAL ####
setwd("~/Documents/MSU/SRS community stability/SRS-stability")
# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE)




#### raup crick function ####
format_raup <- function(df) {
  # convert data to correct format
  df_wide <- df %>% 
    pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>% # wide format
    dplyr::select(!c("unique_id")) %>% # remove unique ID column
    arrange(time) %>% # Ensure years are sorted properly
    column_to_rownames("time") #convert years to rownames
  
  # compute raup
  rc_dist <- RC.pc(df_wide, weighted = F, taxo.metric = "jaccard")
  
  # format raup
  raup_crick_matrix <- as.matrix(rc_dist[["index"]])
  consecutive_dissimilarity <- sapply(1:(nrow(raup_crick_matrix) - 1), function(i) {
    raup_crick_matrix[i, i + 1]
  })
  
  
  result <- data.frame(
    unique_id = unique(df$unique_id),
    year_pair = paste(sort(unique(df$time))[-length(unique(df$time))], 
                      sort(unique(df$time))[-1], sep = " - "),
    raup_dissimilarity = consecutive_dissimilarity
  )
  
  return(result)
}


### apply to data
srs_data_split <- srs_data %>% 
  dplyr::count(unique_id, time, sppcode) %>%
  group_by(unique_id) %>%
  group_split() 

raup_results <- srs_data_split %>%
  lapply(format_raup)

# format
raup_results2 <- raup_results %>%
  bind_rows() %>%
  separate(unique_id, into = c("block", "patch", "patch_type"), sep = "-", remove = F) %>% # separating unique ID into columns
  separate(year_pair, into = c("time1", "time2"), sep = " - ", remove = FALSE) %>% # converting year pairs into a usable format for modeling/plotting
  mutate(year_pair = as.factor(year_pair)) %>% # converting to factor
  mutate(time1 = as.numeric(time1)) %>% # converting to numeric
  mutate(time2 = as.numeric(time2)) # converting to numeric

# write results
#write.csv(raup_results2, file = file.path("data", "L2_summarized", "raup_values.csv"), row.names = F)

# plotting
raup_results2 %>%
  #mutate(raup_dissimilarity = (raup_dissimilarity-.5)*2) %>%
  filter(!patch_type %in% c("center")) %>%
  ggplot(aes(time2, raup_dissimilarity, color = patch_type)) + 
  facet_wrap(~patch_type) + 
  geom_jitter(size = 2, alpha = 0.7) + 
  #geom_line(aes(time2, raup_dissimilarity, group = patch_type), linewidth = 1, alpha = 0.3) +
  theme_bw() + 
  scale_color_brewer(palette = "Set2") +
  stat_smooth(method = "lm", se = F, linewidth = 3)

raup_results2 %>%
  ggplot(aes(patch_type, raup_dissimilarity, color = patch_type)) + 
  geom_boxplot()

raup_results3 <- raup_results2 %>%
  filter(raup_dissimilarity != -1)

hist(raup_results3$raup_dissimilarity)
m1 <- glmmTMB(raup_dissimilarity ~ patch_type * time2 + (1|block/patch),
              data = raup_results2, 
              family = gaussian)
summary(m1)
plot(simulateResiduals(m1))





###### LOCAL STABILITY #######
# loading libraries
library(librarian) # Load `librarian` package
shelf(tidyverse, vegan) # Install missing packages and load needed libraries

source("00_functions.R")
# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE)



#### ALL SPECIES ####
# splitting dataframe into each patch as a list for function
srs_data_split <- srs_data %>% 
  dplyr::count(unique_id, time, sppcode) %>%
  group_by(unique_id) %>%
  group_split()

# applying function to data
jaccard_results <- srs_data_split %>% 
  lapply(compute_jaccard) %>%
  bind_rows() # putting together into a dataframe

## sorensen values
sorensen_results <- srs_data_split %>% 
  lapply(compute_sorensen) %>%
  bind_rows() # putting together into a dataframe


#### COMMON SPECIES ####
# repeating excluding rare species (species with less than 10 occurrences across all surveys)

# figuring out which species are common, assigning a value to those to keep later 
common_spp <- srs_data %>% 
  count(sppcode) %>% # counting # of occurrences per species across all surveys
  filter(n > 10) %>% # removing species that had less than 10 occurrences
  mutate(rare = 0) %>% #  assigning a 0 to common species, rare will be NA later when joining
  dplyr::select(!c("n")) # removing counts of occurrences

commom_spp_data <- srs_data %>%
  left_join(common_spp, by = c("sppcode")) %>% # joining data on which species are common
  filter(rare == 0) %>% # only keeping common species
  dplyr::select(!c("rare")) # removing column

common_spp_jaccard <- commom_spp_data %>%
  dplyr::count(unique_id, time, sppcode) %>% # preparing data for function
  group_by(unique_id) %>%
  group_split() %>% # splitting dataframe into each patch as a list for function
  lapply(compute_jaccard) %>% # applying function to common species
  bind_rows() %>% # putting together into a dataframe
  rename(COMMON_jaccard = jaccard_dissimilarity, # renaming columns to be different than calculations with all species
         COMMON_turnover = turnover_values, 
         COMMON_nestedness = nestedness_values,
         COMMON_alpha = mean_alpha_diversity)

#### LONGLEAF SPECIES (excluding rare) ####
# IGNORE FOR NOW 

#### DISPERSAL MODE (excluding rare) ####
## gravity dispersed
gravity_jaccard <- commom_spp_data %>%
  filter(dispersal_mode == "Gravity") %>%
  dplyr::count(unique_id, time, sppcode) %>%
  group_by(unique_id) %>%
  group_split() %>%
  lapply(compute_jaccard) %>%
  bind_rows() %>% # putting together into a dataframe
  rename(GRAVITY_jaccard = jaccard_dissimilarity, # renaming columns to be different than calculations with all species
         GRAVITY_turnover = turnover_values, 
         GRAVITY_nestedness = nestedness_values,
         GRAVITY_alpha = mean_alpha_diversity)

## wind dispersed 
wind_jaccard <- commom_spp_data %>%
  filter(dispersal_mode == "Wind") %>%
  dplyr::count(unique_id, time, sppcode) %>%
  group_by(unique_id) %>%
  group_split() %>%
  lapply(compute_jaccard) %>%
  bind_rows() %>% # putting together into a dataframe
  rename(WIND_jaccard = jaccard_dissimilarity, # renaming columns to be different than calculations with all species
         WIND_turnover = turnover_values, 
         WIND_nestedness = nestedness_values,
         WIND_alpha = mean_alpha_diversity)

## animal dispersed 
animal_jaccard <- commom_spp_data %>%
  filter(dispersal_mode == "Animal") %>%
  dplyr::count(unique_id, time, sppcode) %>%
  group_by(unique_id) %>%
  group_split() %>%
  lapply(compute_jaccard) %>%
  bind_rows() %>% # putting together into a dataframe
  rename(ANIMAL_jaccard = jaccard_dissimilarity, # renaming columns to be different than calculations with all species
         ANIMAL_turnover = turnover_values, 
         ANIMAL_nestedness = nestedness_values,
         ANIMAL_alpha = mean_alpha_diversity)




#### Temporal gamma diversity (spp pool size of each patch?)

# gamma diversity across time
gamma_diversity <- commom_spp_data %>%
  count(unique_id, sppcode) %>%
  mutate(n = 1) %>% # changing values to 1 for occurence
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>%
  rowwise() %>% # grouping by rows
  mutate(gamma_diversity = sum(across(-1), na.rm = F)) %>% # summing across rows - if there is a 1 in any of the plots the species occurred in that patch
  dplyr::select(unique_id, gamma_diversity)

# gamma of gravity dispersed
GRAVITY_gamma <- commom_spp_data %>%
  filter(dispersal_mode == "Gravity") %>%
  count(unique_id, sppcode) %>%
  mutate(n = 1) %>% # changing values to 1 for occurence
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>%
  rowwise() %>% # grouping by rows
  mutate(GRAVITY_gamma = sum(across(-1), na.rm = F)) %>% # summing across rows - if there is a 1 in any of the plots the species occurred in that patch
  dplyr::select(unique_id, GRAVITY_gamma)

# gamma of wind dispersed
WIND_gamma <- commom_spp_data %>%
  filter(dispersal_mode == "Wind") %>%
  count(unique_id, sppcode) %>%
  mutate(n = 1) %>% # changing values to 1 for occurence
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>%
  rowwise() %>% # grouping by rows
  mutate(WIND_gamma = sum(across(-1), na.rm = F)) %>% # summing across rows - if there is a 1 in any of the plots the species occurred in that patch
  dplyr::select(unique_id, WIND_gamma)

# gamma of animal dispersed
ANIMAL_gamma <- commom_spp_data %>%
  filter(dispersal_mode == "Animal") %>%
  count(unique_id, sppcode) %>%
  mutate(n = 1) %>% # changing values to 1 for occurence
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>%
  rowwise() %>% # grouping by rows
  mutate(ANIMAL_gamma = sum(across(-1), na.rm = F)) %>% # summing across rows - if there is a 1 in any of the plots the species occurred in that patch
  dplyr::select(unique_id, ANIMAL_gamma)


#####
#jaccard_results <- jaccard_results %>%
#  mutate(dispersal_mode = "all")

#wind_jaccard <- wind_jaccard %>%
#  mutate(dispersal_mode = "wind")

#gravity_jaccard <- gravity_jaccard %>%
#  mutate(dispersal_mode = "gravity")

#animal_jaccard <- animal_jaccard %>%
#  mutate(dispersal_mode = "animal")

#jaccard_agg <- rbind(
#  jaccard_results,
# wind_jaccard,
#  gravity_jaccard,
#  animal_jaccard
#)

#jaccard_agg <- jaccard_agg %>%
#  separate(unique_id, into = c("block", "patch", "patch_type"), sep = "-", remove = F) %>% # separating unique ID into columns
#  separate(year_pair, into = c("time1", "time2"), sep = " - ", remove = FALSE) %>% # converting year pairs into a usable format for modeling/plotting
#  mutate(year_pair = as.factor(year_pair)) %>% # converting to factor
#  mutate(time1 = as.numeric(time1)) %>% # converting to numeric
# mutate(time2 = as.numeric(time2)) %>% # converting to numeric
# mutate(years_bw_surveys = time2 - time1) 



#### joining all together ####
jaccard_results <- jaccard_results %>%
  left_join(common_spp_jaccard, by = c("unique_id", "year_pair")) %>% # joining with other data
  left_join(gravity_jaccard, by = c("unique_id", "year_pair")) %>% # joining with other data
  left_join(wind_jaccard, by = c("unique_id", "year_pair")) %>% # joining with other data
  left_join(animal_jaccard, by = c("unique_id", "year_pair")) %>% # joining with other data
  left_join(sorensen_results, by = c("unique_id", "year_pair")) %>% # joining with other data
  separate(unique_id, into = c("block", "patch", "patch_type"), sep = "-", remove = F) %>% # separating unique ID into columns
  separate(year_pair, into = c("time1", "time2"), sep = " - ", remove = FALSE) %>% # converting year pairs into a usable format for modeling/plotting
  mutate(year_pair = as.factor(year_pair)) %>% # converting to factor
  mutate(time1 = as.numeric(time1)) %>% # converting to numeric
  mutate(time2 = as.numeric(time2)) %>% # converting to numeric
  mutate(years_bw_surveys = time2 - time1) %>% # calculating # of years between surveys (most are annual)
  left_join(gamma_diversity, by = c("unique_id")) %>%
  left_join(GRAVITY_gamma, by = c("unique_id")) %>%
  left_join(WIND_gamma, by = c("unique_id")) %>%
  left_join(ANIMAL_gamma, by = c("unique_id"))


# adding soil moisture
soil_moisture <- srs_data %>%
  count(block, patch, soil_moisture) %>%
  dplyr::select(!n)

jaccard_results <- jaccard_results %>%
  left_join(soil_moisture, by = c("block", "patch"))

# writing summarized file
write_csv(jaccard_results, file = file.path("data", "L2_summarized", "patch_jaccard.csv"))


##### METACOMMUNITY STABILITY #######
# loading libraries
library(librarian) # Load `librarian` package
shelf(tidyverse, vegan) # Install missing packages and load needed libraries

source("00_functions.R")
# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE)


srs_data_sub <- srs_data %>%
  filter(patch %in% c("A", "B", "C", "D")) %>%
  filter(!year %in% c("2004", "2013", "2024")) # years that center patch was not sampled in any block


srs_metacommunity <- srs_data_sub %>%
  mutate(metacommunity = if_else(str_detect(unique_id, "A-center|B-connected"), 
                                 "connected-metacommunity", "unconnected-metacommunity")) %>%
  mutate(unique_id = paste(block, metacommunity, sep = "-"))


srs_metacommunity_split <- srs_metacommunity %>%
  dplyr::count(unique_id, time, sppcode) %>%
  mutate(n = 1) %>%
  group_by(unique_id) %>%
  group_split()

metacommunity_jaccard <- srs_metacommunity_split %>%
  lapply(compute_jaccard) %>%
  bind_rows()

metacommunity_jaccard <- metacommunity_jaccard %>%
  separate(unique_id, into = c("block", "patch_type", "metacommunity"), sep = "-", remove = F) %>% # separating unique ID into columns
  separate(year_pair, into = c("time1", "time2"), sep = " - ", remove = FALSE) %>% # converting year pairs into a usable format for modeling/plotting
  mutate(year_pair = as.factor(year_pair)) %>% # converting to factor
  mutate(time1 = as.numeric(time1)) %>% # converting to numeric
  mutate(time2 = as.numeric(time2)) %>% # converting to numeric
  mutate(years_bw_surveys = time2 - time1)


metacommunity_jaccard %>% # plotting jaccard dissimilarity across year 
  ggplot(aes(time2, jaccard_dissimilarity, color = patch_type)) + 
  facet_wrap(~patch_type) + 
  geom_point(size = 2, alpha = 0.7) + 
  #geom_line(aes(time, jaccard_dissimilarity, group = unique_id), linewidth = 1, alpha = 0.3) +
  theme_bw() + 
  scale_color_brewer(palette = "Set2") +
  stat_smooth(method = "lm", se = F, linewidth = 3)

#### regressing against alpha diversity ####
m_resid <- glmmTMB(jaccard_dissimilarity ~ mean_alpha_diversity,
                   data = metacommunity_jaccard)
summary(m_resid)
resid <- residuals(m_resid, type = "response")
metacommunity_jaccard$resid_beta_corrected <- resid


m_metacommunity_beta <- glmmTMB(resid_beta_corrected ~ patch_type*time2 + (1|block),
                                data = metacommunity_jaccard)
summary(m_metacommunity_beta)
plot(simulateResiduals(m_metacommunity_beta))




###### TIME LAG ########
librarian::shelf(tidyverse, vegan, glmmTMB) # Install missing packages and load needed libraries


#### computing time lag jaccard values ####
compute_time_lag <- function(df) {
  # convert data to correct format
  df_wide <- df %>% 
    pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>% # wide format
    dplyr::select(!c("unique_id")) %>% # remove unique ID column
    arrange(time) %>% # Ensure years are sorted properly
    column_to_rownames("time") #convert years to rownames
  
  # Calculate pairwise dissimilarities
  beta_dist <- vegdist(df_wide, method = "jaccard") 
  
  # Create time lag matrix
  n <- nrow(df_wide)
  time_lags <- abs(outer(1:n, 1:n, "-"))
  time_lag_values <- time_lags[lower.tri(time_lags)]
  beta_values <- beta_dist[lower.tri(as.matrix(beta_dist))]
  
  
  # store results with year pairs
  result <- data.frame(
    unique_id = unique(df$unique_id), # unique id
    time_lags = time_lag_values,
    beta_values = beta_values
  )
  
  return(result)
}


### applying function
# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE)



#### ALL SPECIES ####
# splitting dataframe into each patch as a list for function
srs_data_split <- srs_data %>% 
  dplyr::count(unique_id, time, sppcode) %>%
  group_by(unique_id) %>%
  group_split()

# applying function to data
time_lag_results <- srs_data_split %>% 
  lapply(compute_time_lag) %>%
  bind_rows() # putting together into a dataframe

time_lag_results <- time_lag_results %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  filter(patch_type != "center")

m1 <- glmmTMB(beta_values ~ patch_type*time_lags + (1|block/patch),
              data = time_lag_results)
summary(m1)

time_lag_results %>%
  ggplot() +
  geom_point(aes(time_lags, beta_values, color = patch_type)) +
  geom_smooth(aes(time_lags, beta_values, color = patch_type), method = "lm") +
  theme_classic(base_size = 16) + 
  facet_wrap(~patch_type) +
  scale_color_brewer(palette = "Set2")




###### JACCARD ANALYSIS #######
# loading libraries
library(librarian) # Load `librarian` package
shelf(tidyverse, glmmTMB, performance, DHARMa, emmeans, ggeffects, car) # Install missing packages and load needed libraries


# loading data
jaccard_patch <- read_csv(file = file.path("data", "L2_summarized", "patch_jaccard.csv"))

jaccard_patch <- jaccard_patch %>%
  filter(patch_type != "center") %>% # removing center patch from analysis
  mutate(patch_type = as.factor(patch_type))

#### Q1: How is connectivity related to community stability across succession? ####
m1 <- glmmTMB(COMMON_jaccard ~ patch_type * time2 + (1|block/patch),
              data = jaccard_patch)
summary(m1)
plot(simulateResiduals(m1))
Anova(m1)

# posthoc for interaction
m1_posthoc <- emtrends(m1, specs = pairwise ~ patch_type, var = "time2")
m1_posthoc

# posthoc for patch type
m1_posthoc <- emmeans(m1, specs = pairwise ~ patch_type * time2)
pairs(m1_posthoc, simple = "patch_type")

# plotting
m1_predict <- ggpredict(m1, terms = c("time2", "patch_type"))
m1_plot <- m1_predict %>%
  rename(patch_type = group) %>%
  ggplot() +
  geom_line(aes(x, predicted, color = patch_type), linewidth = 2) +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high, fill = patch_type), alpha = 0.3) + 
  geom_point(aes(time2, jaccard_dissimilarity, color = patch_type), data = jaccard_patch, size = 2.5, alpha = 0.3) +
  theme_classic(base_size = 16) + 
  facet_wrap(~patch_type) +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  xlab("Year of succession") +
  ylab("Jaccard dissimilarity") +
  labs(color = "Patch type", fill = "Patch type") +
  theme(panel.spacing = unit(1.2, "lines"))

pdf(file = "jaccard.pdf", width = 10, height = 6)
m1_plot
dev.off()
#jaccard_plot <- jaccard_patch %>% # plotting jaccard dissimilarity across year 
#  ggplot(aes(time2, jaccard_dissimilarity, color = patch_type)) + 
#  facet_wrap(~patch_type) + 
#  geom_point(size = 2, alpha = 0.7) + 
#geom_line(aes(time, jaccard_dissimilarity, group = unique_ID), linewidth = 1, alpha = 0.3) +
#  theme_bw() + 
#  scale_color_brewer(palette = "Set2") +
# stat_smooth(method = "lm", se = F, linewidth = 3)
#jaccard_plot


# How does temporal gamma diversity play into this? 
# Would expect that a higher temporal beta diversity would increase gamma diversity



# initial turnover did not differ between patch types 
initial_jaccard <- jaccard_patch %>%
  filter(time2 == 2)

m2 <- glmmTMB(jaccard_dissimilarity ~ patch_type + (1|block/patch),
              data = initial_jaccard)
summary(m2)



#### dispersal mode ####
# gravity dispersed 
m.gravity <- glmmTMB(GRAVITY_jaccard ~ patch_type * time2 + (1|block/patch),
                     data = jaccard_patch)
summary(m.gravity)
# posthoc for interaction
m.gravity_posthoc <- emtrends(m.gravity, specs = pairwise ~ patch_type, var = "time2")
m.gravity_posthoc

# posthoc for patch type
m.gravity_posthoc <- emmeans(m.gravity, specs = pairwise ~ patch_type * time2)
pairs(m.gravity_posthoc, simple = "patch_type")


# wind dispersed 
m.wind <- glmmTMB(WIND_jaccard ~ patch_type * time2 + (1|block/patch),
                  data = jaccard_patch)
summary(m.wind)
# posthoc for interaction
m.wind_posthoc <- emtrends(m.wind, specs = pairwise ~ patch_type, var = "time2")
m.wind_posthoc

# posthoc for patch type
m.wind_posthoc <- emmeans(m.wind, specs = pairwise ~ patch_type * time2)
pairs(m.wind_posthoc, simple = "patch_type")




# animal dispersed 
m.animal <- glmmTMB(ANIMAL_jaccard ~ patch_type * time2 + (1|block/patch),
                    data = jaccard_patch)
summary(m.animal)
# posthoc for interaction
m.animal_posthoc <- emtrends(m.animal, specs = pairwise ~ patch_type, var = "time2")
m.animal_posthoc

# posthoc for patch type
m.animal_posthoc <- emmeans(m.animal, specs = pairwise ~ patch_type * time2)
pairs(m.animal_posthoc, simple = "patch_type")





#### regressing beta diversity against mean alpha diversity to control for alpha diversity ####
m_beta_corrected <- glmmTMB(COMMON_jaccard ~ COMMON_alpha,
                            data = jaccard_patch)
summary(m_beta_corrected)
resid <- residuals(m_beta_corrected, type = "response")
jaccard_patch$resid_beta_corrected <- resid


m_beta <- glmmTMB(resid_beta_corrected ~ patch_type + time2 + (1|block/patch),
                  data = jaccard_patch)
summary(m_beta)
Anova(m_beta)

# posthoc for interaction
m_beta_posthoc <- emtrends(m_beta, specs = pairwise ~ patch_type, var = "time2")
m_beta_posthoc

# posthoc for patch type
m_beta_posthoc <- emmeans(m_beta, specs = pairwise ~ patch_type)
pairs(m_beta_posthoc, simple = "patch_type")

# plotting
# plotting
m_resid_predict <- ggpredict(m_beta, terms = c("time2", "patch_type"))
m_resid_plot <- m_resid_predict %>%
  rename(patch_type = group) %>%
  ggplot() +
  geom_line(aes(x, predicted, color = patch_type), linewidth = 2) +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high, fill = patch_type), alpha = 0.3) + 
  geom_point(aes(time2, resid_beta_corrected, color = patch_type), data = jaccard_patch, size = 2.5, alpha = 0.3) +
  theme_classic(base_size = 16) + 
  facet_wrap(~patch_type) +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  xlab("Year of succession") +
  ylab("Jaccard dissimilarity (residuals)") +
  labs(color = "Patch type", fill = "Patch type") +
  theme(panel.spacing = unit(1.2, "lines"))
m_resid_plot

pdf(file = "jaccard_resid.pdf", width = 10, height = 6)
m_resid_plot
dev.off()

jaccard_results %>% # plotting jaccard dissimilarity across year 
  ggplot(aes(time2, resid_beta_corrected, color = patch_type)) + 
  facet_wrap(~patch_type) + 
  geom_point(size = 2, alpha = 0.7) + 
  #geom_line(aes(time, jaccard_dissimilarity, group = unique_id), linewidth = 1, alpha = 0.3) +
  theme_bw() + 
  scale_color_brewer(palette = "Set2") +
  stat_smooth(method = "lm", se = F, linewidth = 3)




#### directionality over time ####
# 1. for every three consecutive surveys (e.g., 1-3, 2-4, 3-5)
# ---- calculate jaccard distance
# ---- define trajectories
# ---- calculate directionality
# ---- use last year of three year period in analysis as time
srs_data_wider %>%
  count(time) %>%
  View()

calculate_directionality <- function(df) {
  patch_info <- df %>% 
    arrange(unique_id, time) %>%
    dplyr::select(unique_id, time, year)
  
  # species matrix
  sp_info <- df %>%
    arrange(unique_id, time) %>%
    mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
    column_to_rownames("unique_id_year") %>%
    dplyr::select(!c("unique_id", "time", "year"))
  
  jaccard <- vegdist(sp_info, method = "jaccard")
  
  trajectory <- defineTrajectories(jaccard, sites = patch_info$unique_id, surveys = patch_info$time)
  
  segment_direction <- trajectoryDirectionality(trajectory)
  segment_direction <- data.frame(segment_direction)
  segment_direction <- segment_direction %>%
    rownames_to_column("unique_id") %>%
    separate(unique_id, into = c("block", "patch", "patch_type")) %>%
    #mutate(time = "Year 13-24") %>%
    rename(directionality = segment_direction)
  
  return(segment_direction)
}

patch_info <- srs_data_wider %>%
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
sp_info <- srs_data_wider %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))








calculate_trajectory_directionality <- function(comm_matrix, time_vector) {
  # comm_matrix: site × species matrix, with rows as observations across time
  # time_vector: vector of time points corresponding to each row in comm_matrix
  
  directionality_results <- data.frame(
    end_year = numeric(),
    directionality = numeric(),
    stringsAsFactors = FALSE
  )
  
  patch_info <- comm_matrix %>%
    arrange(unique_id, time) %>%
    dplyr::select(unique_id, time, year)
  
  sp_info <- comm_matrix %>%
    arrange(unique_id, time) %>%
    mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
    column_to_rownames("unique_id_year") %>%
    dplyr::select(!c("unique_id", "time", "year"))
  
  n_years <- max(patch_info$time)
  
  for (i in 1:(n_years - 2)) {
    # Subset for 3 consecutive years
    sub_matrix <- sp_info[i:(i+2), , drop = FALSE]
    sub_patch <- patch_info[i:(i+2), ]
    time_value <- sub_patch$time
    
    # Skip if any NA in years or community data
    if (any(is.na(time_value)) || any(is.na(sub_matrix))) next
    
    if (length(unique(sub_patch)) == 3 && all(diff(sub_patch$time) > 0)) {
      # Step 1: Calculate Jaccard distances
      jaccard_dist <- vegdist(sub_matrix, method = "jaccard", binary = TRUE)
      
      # Step 2: Define trajectory
      traj <- defineTrajectories(jaccard_dist, sites = sub_patch$unique_id, surveys = sub_patch$time)
      
      # Step 3: Calculate directionality
      dir_val <- trajectoryDirectionality(traj)
      
      # Step 4: Store the result with the last year in the 3-year window
      directionality_results <- rbind(directionality_results, data.frame(
        end_year = max(time_value),
        directionality = dir_val
      ))
    }
  }
  
  return(directionality_results)
}

srs_direction_split <- srs_data_wider %>%
  #separate(unique_id, into = c("block", "patch_rep", "patch_type"), sep = "-", remove = F) %>%
  #filter(!block %in% c("75E", "75W")) %>%
  #select(!c("block", "patch_rep", "patch_type")) %>%
  group_by(unique_id) %>%
  group_split()

directionality <- srs_direction_split %>%
  lapply(calculate_trajectory_directionality)


directionality <- directionality %>%
  bind_rows() # putting together into a dataframe

directionality <- directionality %>%
  rownames_to_column("unique_ID") %>%
  separate(unique_ID, into = c("block", "patch_rep", "patch_type_time"), sep = "-") %>%
  dplyr::select(!patch_type_time)


patch_type <- srs_data %>%
  count(block, patch, patch_type) %>%
  dplyr::select(!n)

directionality <- directionality %>%
  left_join(patch_type, by = c("block" = "block",
                               "patch_rep" = "patch"))


directionality %>%
  ggplot(aes(end_year, directionality, color = patch_type, fill = patch_type)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
  theme_minimal() +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type")

m1 <- glmmTMB(directionality ~ patch_type * end_year + (1|block),
              data = directionality)
summary(m1)
plot(simulateResiduals(m1))
m1.posthoc <- emtrends(m1, "patch_type", var = "end_year")
pairs(m1.posthoc)

hist(directionality$directionality)







