# loading libraries
# Load `librarian` package
library(librarian)
# Install missing packages and load needed libraries
shelf(tidyverse, vegan, codyn, glmmTMB, performance, DHARMa, emmeans)

# loading data from EDI
srs_data_url <- "https://pasta.lternet.edu/package/data/eml/edi/414/1/2429e7fc1b33cefb59bab8451aaa8327"
srs_data_raw <- read.csv(file = srs_data_url)

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

jaccard_results %>% # plotting jaccard dissimilarity across year 
  ggplot(aes(time, jaccard_dissimilarity, color = patch)) + 
  facet_wrap(~patch) + 
  geom_point() + 
  theme_bw() + 
  stat_smooth(method = "lm", se = F, size = 2)

hist(jaccard_results$jaccard_dissimilarity) # histogram of dissimilarity values

# how does dissimilarity change across time for each patch type?
m1 <- glmmTMB(jaccard_dissimilarity ~ patch * time + (1|EU/patch_rep), 
              data = jaccard_results)
summary(m1) # model summary
plot(simulateResiduals(m1)) # looks not the best but not the worst
check_model(m1) # okay for now

m1_posthoc <- emtrends(m1, pairwise ~ patch, var = "time") # posthoc test for differences between slopes
m1_posthoc # connected different than rectangular, but not winged







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








