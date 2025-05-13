
librarian::shelf(tidyverse, vegan, ape, BiodiversityR)

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE) %>%
  filter(patch_type != "center")


# pivot to wider format
srs_data_wider <- srs_data %>%
  dplyr::count(unique_id, time, year, sppcode) %>%
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) # wide format


# make factor
srs_data_wider$time <- as.factor(srs_data_wider$time)
srs_data_wider$unique_id <- as.factor(srs_data_wider$unique_id)
srs_data_wider$year <- as.factor(srs_data_wider$year)

# patch data
patch_info <- srs_data_wider %>% 
  arrange(unique_id, time) %>%
  select(unique_id, time, year)

# species matrix
sp_info <- srs_data_wider %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  select(!c("unique_id", "time", "year"))


# Jaccard distance matrix
jaccard_dist <- vegdist(sp_info, method = "jaccard")
jaccard_dist_all <- as.matrix(jaccard_dist)
jaccard_dist_all_df <- as.data.frame(jaccard_dist_all)

jaccard_dist_all_df <- cbind(patch_info, jaccard_dist_all_df) # merge with patch info

# pcoa
pcoa_all <- pcoa(jaccard_dist)
pcoa_all_cmd <- cmdscale(jaccard_dist, eig=TRUE, add=FALSE)  ## Another way of doing the pcoa that gives same results#


# add axes to patch and time
pcoa_axes <- pcoa_all$vectors[,c(1,2)]
pcoa_axes <- cbind(patch_info, pcoa_axes)


#### species scores ####
pcoa_all_cmd_scores <- add.spec.scores(pcoa_all_cmd, sp_info, method = 'cor.scores')

spscores<-as.data.frame(pcoa_all_cmd_scores$cproj)
names<-rownames(spscores)
rownames(spscores) <- NULL

spscores_df <- cbind(names,spscores)


#### permanova ####
srs_wider_permanova <- srs_data_wider %>%
  separate(unique_id, into = c("block", "patch_rep", "patch_type"), sep = "-", remove = F)

permanova_01_03 <- srs_wider_permanova %>%
  filter(year %in% c("2001", "2002", "2003"))
permanova_05_06 <- srs_wider_permanova %>%
  filter(year %in% c("2005", "2006"))
permanova_07 <- srs_wider_permanova %>%
  filter(year %in% c("2007"))
permanova_08_12 <- srs_wider_permanova %>%
  filter(year %in% c("2008", "2009", "2010", "2011", "2012"))
permanova_13 <- srs_wider_permanova %>%
  filter(year %in% c("2013"))
permanova_14_15 <- srs_wider_permanova %>%
  filter(year %in% c("2014", "2015"))
permanova_16 <- srs_wider_permanova %>%
  filter(year %in% c("2016"))
permanova_17 <- srs_wider_permanova %>%
  filter(year %in% c("2017"))
permanova_18_24 <- srs_wider_permanova %>%
  filter(year %in% c("2018", "2019", "2020", "2021", "2022", "2023", "2024"))

loop_permanova <- function(df_wide) {
  results <- adonis2(vegdist(df_wide[,7:331], method = "jaccard") ~ patch_type + block, data = df_wide, method = "jaccard")
  return(data.frame(results))
}


permanova_results_01_03 <- permanova_01_03 %>%
  group_by(year) %>%
  group_split() %>%
  lapply(loop_permanova) %>%
  bind_rows()

names<-c("Patch Type", "Block", "residual", "total")
sourcevariation<-rep(names,length(unique(permanova_results_01_03$.id)))



#### plotting ####
pcoa_axes_plot <- pcoa_axes %>%
  separate(unique_id, into = c("block", "patch_rep", "patch_type"), sep = "-", remove = F)


pcoa_axes_plot %>%
  filter(patch_type != "center") %>%
  ggplot(aes(Axis.1, Axis.2, color = year, shape = patch_type)) +
  geom_point() +
  geom_path(aes(Axis.1, Axis.2, group = patch_type, color = year)) +
  facet_grid(block~patch_rep) +
  scale_color_viridis_d() +
  theme_minimal(base_size = 14)







