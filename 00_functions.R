# loading libraries
librarian::shelf(tidyverse, vegan, iCAMP) # Install missing packages and load needed libraries

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











