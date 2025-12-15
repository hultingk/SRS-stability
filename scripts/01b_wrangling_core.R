# loading packages
librarian::shelf(tidyverse, summarytools)

# sourcing soil moisture and other wrangling from other script
source(here::here("01_wrangling.R"))

# loading data
srs_core <- read.csv(file = file.path("data", "L0_original", "core_for_KH_20250415.csv"))


# making column names lower case for joining and preference
colnames(srs_core) <- stringr::str_to_lower(colnames(srs_core))

# adding patch info, year, and dispersal mode to plant data
srs_core <- srs_core %>%
  left_join(year_created, by = c("block")) %>%
  left_join(patch_info, by = c("block", "patch")) %>%
  left_join(dispersal_mode, by = c("sppcode" = "Species.Code"))


#### spp codes, dispersal modes ####
# resolving spp codes
srs_core <- srs_core %>%
  mutate(sppcode = dplyr::case_when(
    sppcode %in% c("DESPER") ~ "DESGLA", # changing DESPER to DESGLA at Christopher's suggestion -- email 05/24/2025
    sppcode %in% c("SOLPTY", "SOLNIG", "SOLAME") ~ "SOLNIG", # grouping these at Christopher's suggestion -- email 05/24/2025
    sppcode %in% c("DIGSSPP") ~ "DIGSPP", # fixing typo -- 7 letter code DIGSSPP to 6 letter code DIGSPP
    .default = sppcode
  ))

# adding missing dispersal modes 
srs_core <- srs_core %>%
  mutate(DispMode1 = dplyr::case_when(
    sppcode %in% c("ALLCUT", "GALMOL", "MANVIR", "PASSPP", "STRUMB", 
                   "BAPALS", "CHEALB", "DICSPP", "DIGSPP", "RHYSPP") ~ "Gravity",
    sppcode %in% c("ANDSPP", "CATBIG", "CIRHOR", "CUSSPP", 
                   "LIASCG", "OSMREG", "SABDIF", "SENSMA", "TRAGOP",
                   "HYPSPP") ~ "Wind",
    sppcode %in% c("CARSPP", "GALSPP", "LANTAN", "PLUSPP", "SOLNIG", 
                   "DESLRG", "DESGLA", "ILEDEC") ~ "Animal",
    .default = DispMode1
  ))

# grouping dispersal modes into gravity, wind, or animal dispersed
srs_core <- srs_core %>%
  mutate(dispersal_mode = dplyr::case_when(
    DispMode1 %in% c("Adhesive", "Animal", "Ant", "Bird", "Mammal", "MammalBird") ~ "Animal",
    DispMode1 %in% c("Tumbling", "Wind") ~ "Wind",
    DispMode1 %in% c("Ballistic", "Gravity", "Unassisted") ~ "Gravity",
    .default = DispMode1
  ))


# fixing transplant column
srs_core <- srs_core %>%
  mutate(transplant = dplyr::case_when(
    sppcode %in% c("ARIBEY", "SORSEC", "ANTVIL", "CARBEL", "LIAEAR", "PHYAME",
                   "SOLAME", "NOLGEO", "GAYDUM", "RUDHIR", "LANCAM", "ILEVER") ~ TRUE,
    .default = FALSE
  )) %>%
  dplyr::select(!c("transplant.")) # removing old transplant column


# adding survey time since patch creation
srs_core <- srs_core %>%
  mutate(time = year - year.created)

# joining soil moisture 
srs_core <- srs_core %>%
  left_join(soil_moisture, by = c("block", "patch"))

# joining year since fire 
srs_core <- srs_core %>%
  left_join(year_fire, by = c("block", "year"))


#### final dataset ####
# renaming, removing unneeded columns, and rearranging
srs_core <- srs_core %>%
  mutate(year_created = year.created) %>% # renaming
  mutate(unique_id = paste(block, patch, patch_type, sep = "-")) %>% # creating unique ID for each patch
  dplyr::select(block, patch, cell, patch_type, unique_id, soil_moisture, core, year, year_created, time,
                sppcode, transplant, dispersal_mode, year_since_fire)

# checking for missing values and duplicates
summarytools::view(summarytools::dfSummary(srs_core), footnote = NA)

# writing cleaned file
write_csv(srs_core, file = file.path("data", "L1_wrangled", "srs_plant_core.csv"))
