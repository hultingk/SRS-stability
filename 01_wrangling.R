# loading libraries
# Load `librarian` package
library(librarian)
# Install missing packages and load needed libraries
shelf(tidyverse, summarytools)


# loading data
srs_all <- read.csv(file = file.path("data", "L0_original", "combined_for_KH_20250415.csv"))
patch_info <- read.csv(file = file.path("data", "L0_original", "corridor_patch_info.csv"))
year_created <- read.csv(file = file.path("data", "L0_original", "year_site_creation.csv"))
dispersal_mode <- read.csv(file = file.path("data", "L0_original", "dispersal_mode_20250418.csv"))

# making column names lower case for joining and preference
colnames(srs_all) <- stringr::str_to_lower(colnames(srs_all))

# changing "8" to "08" to be consistent 
year_created <- year_created %>%
  mutate(block = if_else(block == "8", "08", block))
patch_info <- patch_info %>%
  mutate(block = if_else(block == "8", "08", block))

# getting rid of columns of dispersal mode that I don't need
dispersal_mode <- dispersal_mode %>%
  dplyr::select(c("Species.Code", "USDACode", "DispMode1", "DispMode2"))

# adding patch info, year, and dispersal mode to plant data
srs_all <- srs_all %>%
  left_join(year_created, by = c("block")) %>%
  left_join(patch_info, by = c("block", "patch")) %>%
  left_join(dispersal_mode, by = c("sppcode" = "Species.Code"))


# adding survey time since patch creation
srs_all <- srs_all %>%
  mutate(time = year - year.created)

# checking for missing values and duplicates
summarytools::view(summarytools::dfSummary(srs_all), footnote = NA)


# renaming, removing unneeded columns, and rearranging
srs_all <- srs_all %>%
  mutate(year_created = year.created) %>% # renaming
  mutate(dipsersal_mode = DispMode1) %>% # renaming
  mutate(dipsersal_mode2 = DispMode2) %>% # renaming
  mutate(unique_id = paste(block, patch, patch_type, sep = "-")) %>% # creating unique ID for each patch
  dplyr::select(block, patch, cell, patch_type, unique_id, core, year, year_created, time,
                sppcode, dipsersal_mode, dipsersal_mode2)


# writing cleaned file
write_csv(srs_all, file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))






