# loading libraries
# Load `librarian` package
library(librarian)
# Install missing packages and load needed libraries
shelf(tidyverse, vegan, codyn, glmmTMB, performance, DHARMa, emmeans)

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
  separate(unique_ID, into = c("EU", "patch_rep", "patch"), sep = "_") %>% # seperating unique ID into columns
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
  geom_point() + 
  theme_bw() + 
  scale_color_brewer(palette = "Set2") +
  stat_smooth(method = "lm", se = F, linewidth = 2)
jaccard_plot

hist(jaccard_results$jaccard_dissimilarity) # histogram of dissimilarity values

# how does dissimilarity change across time for each patch type?
m1 <- glmmTMB(jaccard_dissimilarity ~ patch * time + (1|EU/patch_rep), 
              data = jaccard_results)
summary(m1) # model summary
plot(simulateResiduals(m1)) # looks not the best but not the worst
check_model(m1) # okay for now

m1_posthoc <- emtrends(m1, pairwise ~ patch, var = "time") # posthoc test for differences between slopes
m1_posthoc # connected different than rectangular, but not winged



########


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
  separate(unique_ID, into = c("EU", "patch_rep", "patch"), sep = "_") %>% # seperating unique ID into columns
  mutate(year_pair = as.factor(year_pair)) %>% # making a factor to prepare for joining
  left_join(year_factor, by = c("year_pair" = "year_pair")) %>% # joining to year info
  mutate(time = as.numeric(time)) # making time numeric


colors <- c("Jaccard dissimilarity" = "#8DA0CB", "Turnover" = "#FC8D62", "Nestedness" = "#66C2A5")
jaccard_plot2 <- jaccard_results2 %>% # plotting jaccard dissimilarity across year 
  ggplot() + 
  facet_wrap(~patch) + 
  geom_point(aes(time, nestedness_values, color = "Nestedness")) + 
  geom_point(aes(time, turnover_values, color = "Turnover")) + 
  geom_point(aes(time, jaccard_dissimilarity, color = "Jaccard dissimilarity")) + 
  theme_bw() + 
  stat_smooth(aes(time, nestedness_values, color = "Nestedness"), method = "lm", se = F, linewidth = 2) +
  stat_smooth(aes(time, turnover_values, color = "Turnover"), method = "lm", se = F, linewidth = 2)+ 
  stat_smooth(aes(time, jaccard_dissimilarity, color = "Jaccard dissimilarity"), method = "lm", se = F, linewidth = 2) +
  labs(x = "Time",
       y = "",
       color = "Legend") +
  scale_color_manual(values = colors)
  
jaccard_plot2

# how does dissimilarity change across time for each patch type?
m2 <- glmmTMB(turnover_values ~ patch * time + (1|EU/patch_rep), 
              data = jaccard_results2)
summary(m2) # model summary
plot(simulateResiduals(m2)) # looks not the best but not the worst
check_model(m2) # okay for now

m2_posthoc <- emtrends(m2, pairwise ~ patch, var = "time") # posthoc test for differences between slopes
m2_posthoc # connected different than rectangular, but not winged




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
  geom_point() + 
  theme_bw() + 
  scale_color_brewer(palette = "Set2") +
  ylab("jaccard_dissimilarity, wind dispered") +
  stat_smooth(method = "lm", se = F, linewidth = 2)
wind_plot

# how does dissimilarity change across time for each patch type?
m.wind <- glmmTMB(jaccard_dissimilarity ~ patch * time + (1|EU/patch_rep), 
              data = wind_jaccard)
summary(m.wind) # model summary
plot(simulateResiduals(m.wind)) # looks not the best but not the worst
check_model(m.wind) # okay for now

m.wind_posthoc <- emtrends(m.wind, pairwise ~ patch, var = "time") # posthoc test for differences between slopes
m.wind_posthoc # connected different than rectangular, but not winged


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
  geom_point() + 
  theme_bw() + 
  scale_color_brewer(palette = "Set2") +
  ylab("jaccard_dissimilarity, gravity dispered") +
  stat_smooth(method = "lm", se = F, linewidth = 2)
gravity_plot

# how does dissimilarity change across time for each patch type?
m.gravity <- glmmTMB(jaccard_dissimilarity ~ patch * time + (1|EU/patch_rep), 
                  data = gravity_jaccard)
summary(m.gravity) # model summary
plot(simulateResiduals(m.gravity)) # looks not the best but not the worst
check_model(m.gravity) # okay for now

m.gravity_posthoc <- emtrends(m.gravity, pairwise ~ patch, var = "time") # posthoc test for differences between slopes
m.gravity_posthoc # connected different than rectangular, but not winged



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
  geom_point() + 
  theme_bw() + 
  scale_color_brewer(palette = "Set2") +
  ylab("jaccard_dissimilarity, animal dispered") +
  stat_smooth(method = "lm", se = F, linewidth = 2)
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








