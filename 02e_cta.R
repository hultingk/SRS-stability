### community trajectory analysis
librarian::shelf(tidyverse, vegan, ecotraj) # Install missing packages and load needed libraries

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE) %>%
  filter(patch_type != "center")




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




## trash analysis 
d1 <- vegan::vegdist(srs_data_wider, "jaccard")
cta_trial_sites <- srs_data_wider %>%
  rownames_to_column("unique_id") %>%
  select(unique_id) %>%
  separate(unique_id, into = c("block", "patch", "patch_type", "time")) %>%
  mutate(time = as.numeric(time)) %>%
  mutate(unique_id = paste(block, patch, patch_type, sep = "-")) %>%
  mutate(color = dplyr::case_when(
    patch_type %in% c("connected") ~ "yellow", 
    patch_type %in% c("center") ~ "blue", 
    patch_type %in% c("rectangle") ~ "purple", 
    patch_type %in% c("wing") ~ "red"
  ))


trajectoryPCoA(d1, sites = cta_trial_sites$unique_id, surveys = cta_trial_sites$time,
               traj.colors = cta_trial_sites$color, survey.labels = T, length = 0.1, lwd = 2)

segment_lengths <- trajectoryLengths(d1, 
                                     sites = cta_trial_sites$unique_id, 
                                     surveys = cta_trial_sites$time)
segment_lengths <- segment_lengths %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type"))

segment_lengths <- segment_lengths %>%
  pivot_longer(cols = S1:S22, names_to = "time", values_to = "distance") %>%
  mutate(time = as.numeric(sub("S", "", time)))

segment_lengths <- segment_lengths %>%
  filter(patch_type != "center")
m_total <- glmmTMB(distance ~ patch_type*time + (1|block/patch),
                   data = segment_lengths,
                   family = gaussian)
summary(m_total)
plot(simulateResiduals(m_total))
# posthoc for patch type
m_total_posthoc <- emmeans(m_total, specs = pairwise ~ patch_type)
pairs(m_total_posthoc)







x <- defineTrajectories(d1, sites = cta_trial_sites$unique_id, surveys = cta_trial_sites$time)
srs_convergence <- trajectoryConvergence(x)
srs_convergance_longer <- as.data.frame(srs_convergence[["tau"]])
srs_convergance_longer <- srs_convergance_longer %>%
  rownames_to_column(var = "unique_id") %>%
  pivot_longer(cols = 2:51, names_to = "unique_id2", values_to = "convergence")

srs_convergance <- srs_convergance_longer %>%
  separate(unique_id, into = c("block", "patch", "patch_type"), sep = "-") %>%
  separate(unique_id2, into = c("block2", "patch2", "patch_type2"), sep = "-") %>%
  filter(block == block2)

srs_convergance <- srs_convergance %>%
  mutate(patch_pair = paste(patch, patch2, sep = "-")) %>%
  mutate(patch_type_pair = paste(patch_type, patch_type2, sep = "-")) %>%
  mutate(connected_unconnected = dplyr::case_when(
    patch_type_pair %in% c("center-connected", "connected-center") ~ "connected pair",
    patch_type_pair %in% c("rectangle-wing", "wing-rectangle", "wing-wing", "rectangle-rectangle") ~ "unconnected pair",
    .default = NA
  ))


srs_convergance %>%
  ggplot() +
  geom_boxplot(aes(connected_unconnected, convergence)) 

m1 <- glmmTMB(convergence ~ connected_unconnected + (1|block/patch_pair),
              data = srs_convergance)
summary(m1) # most unconnected pairs within a block are diverging in composition, while the connected and center patch are convering more




