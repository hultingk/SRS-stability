
librarian::shelf(tidyverse, vegan, ape, BiodiversityR)

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

srs_data <- srs_data %>% # removing experimentally planted species 
  filter(transplant != TRUE) %>%
  filter(patch_type != "center")


# pivot to wider format
srs_data_wider <- srs_data %>%
  dplyr::count(unique_id, time, year, sppcode, soil_moisture) %>%
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) # wide format


# make factor
srs_data_wider$time <- as.factor(srs_data_wider$time)
srs_data_wider$unique_id <- as.factor(srs_data_wider$unique_id)
srs_data_wider$year <- as.factor(srs_data_wider$year)

# patch data
patch_info <- srs_data_wider %>% 
  arrange(unique_id, time) %>%
  select(unique_id, time, year, soil_moisture)

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
  separate(unique_id, into = c("block", "patch_rep", "patch_type"), sep = "-", remove = F) %>%
  filter(patch_rep %in% c("B", "C", "D", "E")) %>%
  filter(!block %in% c("52", "57"))


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
  results <- adonis2(vegdist(df_wide[,8:332], method = "jaccard") ~ patch_type + soil_moisture + block, data = df_wide, by = "margin", method = "jaccard")
  return(data.frame(results))
}


permanova_results_01_03 <- permanova_01_03 %>%
  group_by(year) %>%
  group_split() %>%
  lapply(loop_permanova) %>%
  bind_rows() %>%
  rownames_to_column("variable") %>%
  mutate(year = rep(2001:2003, each = 5)) %>%
  mutate(variable = dplyr::case_when(
    str_detect(variable, "patch_type") ~ "Patch Type",
    str_detect(variable, "block") ~ "Block",
    str_detect(variable, "soil_moisture") ~ "Soil Moisture",
    str_detect(variable, "Residual") ~ "Residual",
    str_detect(variable, "Total") ~ "Total"
  )) %>%
  select(year, variable, R2)

permanova_results_05_06 <- permanova_05_06 %>%
  group_by(year) %>%
  group_split() %>%
  lapply(loop_permanova) %>%
  bind_rows() %>%
  rownames_to_column("variable") %>%
  mutate(year = rep(2005:2006, each = 5)) %>%
  mutate(variable = dplyr::case_when(
    str_detect(variable, "patch_type") ~ "Patch Type",
    str_detect(variable, "block") ~ "Block",
    str_detect(variable, "soil_moisture") ~ "Soil Moisture",
    str_detect(variable, "Residual") ~ "Residual",
    str_detect(variable, "Total") ~ "Total"
  )) %>%
  select(year, variable, R2)

permanova_results_07 <- permanova_07 %>%
  group_by(year) %>%
  group_split() %>%
  lapply(loop_permanova) %>%
  bind_rows() %>%
  rownames_to_column("variable") %>%
  mutate(year = rep(2007, each = 5)) %>%
  mutate(variable = dplyr::case_when(
    str_detect(variable, "patch_type") ~ "Patch Type",
    str_detect(variable, "block") ~ "Block",
    str_detect(variable, "soil_moisture") ~ "Soil Moisture",
    str_detect(variable, "Residual") ~ "Residual",
    str_detect(variable, "Total") ~ "Total"
  )) %>%
  select(year, variable, R2)

permanova_results_08_12 <- permanova_08_12 %>%
  group_by(year) %>%
  group_split() %>%
  lapply(loop_permanova) %>%
  bind_rows() %>%
  rownames_to_column("variable") %>%
  mutate(year = rep(2008:2012, each = 5)) %>%
  mutate(variable = dplyr::case_when(
    str_detect(variable, "patch_type") ~ "Patch Type",
    str_detect(variable, "block") ~ "Block",
    str_detect(variable, "soil_moisture") ~ "Soil Moisture",
    str_detect(variable, "Residual") ~ "Residual",
    str_detect(variable, "Total") ~ "Total"
  )) %>%
  select(year, variable, R2)

permanova_results_13 <- permanova_13 %>%
  group_by(year) %>%
  group_split() %>%
  lapply(loop_permanova) %>%
  bind_rows() %>%
  rownames_to_column("variable") %>%
  mutate(year = rep(2013, each = 5)) %>%
  mutate(variable = dplyr::case_when(
    str_detect(variable, "patch_type") ~ "Patch Type",
    str_detect(variable, "block") ~ "Block",
    str_detect(variable, "soil_moisture") ~ "Soil Moisture",
    str_detect(variable, "Residual") ~ "Residual",
    str_detect(variable, "Total") ~ "Total"
  )) %>%
  select(year, variable, R2)

permanova_results_14_15 <- permanova_14_15 %>%
  group_by(year) %>%
  group_split() %>%
  lapply(loop_permanova) %>%
  bind_rows() %>%
  rownames_to_column("variable") %>%
  mutate(year = rep(2014:2015, each = 5)) %>%
  mutate(variable = dplyr::case_when(
    str_detect(variable, "patch_type") ~ "Patch Type",
    str_detect(variable, "block") ~ "Block",
    str_detect(variable, "soil_moisture") ~ "Soil Moisture",
    str_detect(variable, "Residual") ~ "Residual",
    str_detect(variable, "Total") ~ "Total"
  )) %>%
  select(year, variable, R2)

permanova_results_16 <- permanova_16 %>%
  group_by(year) %>%
  group_split() %>%
  lapply(loop_permanova) %>%
  bind_rows() %>%
  rownames_to_column("variable") %>%
  mutate(year = rep(2016, each = 5)) %>%
  mutate(variable = dplyr::case_when(
    str_detect(variable, "patch_type") ~ "Patch Type",
    str_detect(variable, "block") ~ "Block",
    str_detect(variable, "soil_moisture") ~ "Soil Moisture",
    str_detect(variable, "Residual") ~ "Residual",
    str_detect(variable, "Total") ~ "Total"
  )) %>%
  select(year, variable, R2)

permanova_results_17 <- permanova_17 %>%
  group_by(year) %>%
  group_split() %>%
  lapply(loop_permanova) %>%
  bind_rows() %>%
  rownames_to_column("variable") %>%
  mutate(year = rep(2017, each = 5)) %>%
  mutate(variable = dplyr::case_when(
    str_detect(variable, "patch_type") ~ "Patch Type",
    str_detect(variable, "block") ~ "Block",
    str_detect(variable, "soil_moisture") ~ "Soil Moisture",
    str_detect(variable, "Residual") ~ "Residual",
    str_detect(variable, "Total") ~ "Total"
  )) %>%
  select(year, variable, R2)

permanova_results_18_24 <- permanova_18_24 %>%
  group_by(year) %>%
  group_split() %>%
  lapply(loop_permanova) %>%
  bind_rows() %>%
  rownames_to_column("variable") %>%
  mutate(year = rep(2018:2024, each = 5)) %>%
  mutate(variable = dplyr::case_when(
    str_detect(variable, "patch_type") ~ "Patch Type",
    str_detect(variable, "block") ~ "Block",
    str_detect(variable, "soil_moisture") ~ "Soil Moisture",
    str_detect(variable, "Residual") ~ "Residual",
    str_detect(variable, "Total") ~ "Total"
  )) %>%
  select(year, variable, R2)

permanova_results <- rbind(
  permanova_results_01_03, permanova_results_05_06, permanova_results_07, permanova_results_08_12, permanova_results_13, permanova_results_14_15, permanova_results_16, permanova_results_17, permanova_results_18_24
)

r2_permanova_plot <- permanova_results %>%
  filter(!variable %in% c("Residual", "Total")) %>%
  ggplot(aes(year, R2)) +
  geom_point(size = 3) + 
  geom_smooth(color = "black") +
  facet_grid(cols = vars(variable)) +
  theme_minimal(base_size = 22) +
  xlab("Year") +
  theme(panel.spacing = unit(1.5, "lines"),
        plot.margin = margin(12,24,12,12)) +
  theme(axis.text.x = element_text(angle = 60,  hjust=1))
r2_permanova_plot

pdf(file = "r2_permanova_plot.pdf", width = 11, height = 7)
r2_permanova_plot
dev.off()




#### plotting ####
pcoa_axes_plot <- pcoa_axes %>%
  separate(unique_id, into = c("block", "patch_rep", "patch_type"), sep = "-", remove = F)


plot_pcoa <- pcoa_axes_plot %>%
  filter(block %in% c("08", "10", "53N", "53S", "54N", "54S")) %>%
  filter(patch_type != "center") %>%
  filter(patch_rep %in% c("B", "C", "D")) %>%
  ggplot(aes(Axis.1, Axis.2, color = time)) +
  geom_point(size = 2) +
  geom_path(aes(Axis.1, Axis.2, group = patch_type, color = time)) +
  facet_grid(block~patch_type) +
  scale_color_viridis_d(option = "plasma") +
  theme_minimal(base_size = 20)
plot_pcoa

pdf(file = "plot_pcoa.pdf", width = 11, height = 7)
plot_pcoa
dev.off()






