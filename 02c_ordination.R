# loading libraries
librarian::shelf(tidyverse, vegan, ecotraj) # Install missing packages and load needed libraries

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE) %>%
  filter(!year %in% c("2004", "2013", "2024")) # years that center patch was not sampled in any block


srs_data_wider <- srs_data %>%
  dplyr::count(unique_id, time, sppcode) %>%
  mutate(unique_id = paste(unique_id, time, sep = "-")) %>%
  dplyr::select(!time) %>%
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>% # wide format
  column_to_rownames("unique_id") #convert years to rownames
  
srs_metadata <- srs_data_wider %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type", "time")) %>%
  mutate(time = as.numeric(time))


dissim <- vegdist(srs_data_wider, method = "jaccard")  
adonis2(dissim ~ patch_type * time, data = srs_metadata, permutations = 999,
        strata = srs_metadata$block)



set.seed(300)
# Choosing number of dimensions
# Specifying the raw data models for all dimensions k
mods.r <- list() # an empty list to hold model objects
stress.r <- numeric() # an empty numeric vector to hold final stress values 
conv.r <- character() # an empty character vector to hold convergence messages 
for(i in 1:4){
  mod <- metaMDS(srs_data_wider, distance = 'jaccard', k = i, try = 20, trymax = 100)
  mods.r[[as.character(i)]] <- mod 
  stress.r[i] <- mod$stress 
  conv.r[i] <- mod$converged
}
# scree plot
palette(c('red', 'black'))
plot(x = 1:4, y = stress.r, main = 'Raw Data', xlab = 'Dimensionality',
     ylab = 'Stress',
     pch = 19, col = factor(conv.r))



# Distance matrix
srs_dist <- vegdist(srs_data_wider, method = "jaccard")
# Compute NMDS 
srs_nmds <- metaMDS(srs_dist, k = 3, trymax = 100)
srs_nmds$stress

adonis_result <- adonis2(srs_dist ~ patch_type * time, data = srs_metadata, permutations = 999,
                         strata = srs_metadata$block)
print(adonis_result)

disp <- betadisper(srs_dist, srs_metadata$patch_type)
anova(disp)  # Tests if group variances are equal

plot(disp)   # Visualize group spread
boxplot(disp)


# Extract scores and merge with metadata
scores <- scores(srs_nmds)
site_scores <- as.data.frame(scores[[1]])
site_scores <- site_scores %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type", "time")) %>%
  mutate(time = as.numeric(time))

hist(site_scores$NMDS2)

model1 <- glmmTMB(NMDS1 ~ patch_type * time + (1|block/patch), 
                 data = site_scores)
summary(model1)
plot(simulateResiduals(model1))



# Plot with ggplot2
site_scores %>%
  ggplot(aes(x = NMDS1, y = NMDS2, color = block)) +
  facet_wrap(~patch_type) + 
  geom_path(aes(group = time)) +   # trajectories
  geom_point() + # points
  theme_minimal()





