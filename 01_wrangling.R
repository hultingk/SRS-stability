# loading libraries
# Install missing packages and load needed libraries
librarian::shelf(tidyverse, summarytools)


# loading data
srs_all <- read.csv(file = file.path("data", "L0_original", "combined_for_KH_20250415.csv"))
patch_info <- read.csv(file = file.path("data", "L0_original", "corridor_patch_info.csv"))
year_created <- read.csv(file = file.path("data", "L0_original", "year_site_creation.csv"))
dispersal_mode <- read.csv(file = file.path("data", "L0_original", "dispersal_mode_20250418.csv"))
soil_moisture_2003 <- read.csv(file = file.path("data", "L0_original", "soil_moisture_2003.csv"))
soil_moisture_2007 <- read.csv(file = file.path("data", "L0_original", "soil_moisture_demo_plots.csv"))
year_fire <- read.csv(file = file.path("data", "L0_original", "burn_years.csv"))

# making column names lower case for joining and preference
colnames(srs_all) <- stringr::str_to_lower(colnames(srs_all))

# changing "8" to "08" to be consistent 
year_created <- year_created %>% # doing this for creation year
  mutate(block = if_else(block == "8", "08", block))
patch_info <- patch_info %>% # for patch info
  mutate(block = if_else(block == "8", "08", block))
soil_moisture_2007 <- soil_moisture_2007 %>% # and for soil moisture data
  mutate(block = if_else(block == "8", "08", block))
soil_moisture_2003 <- soil_moisture_2003 %>% # and for soil moisture data
  mutate(EU = if_else(EU == "8", "08", EU))
year_fire <- year_fire %>% # and for year since fire
  mutate(block = if_else(block == "8", "08", block))

# getting rid of columns of dispersal mode that I don't need
dispersal_mode <- dispersal_mode %>%
  dplyr::select(c("Species.Code", "USDACode", "DispMode1", "DispMode2")) # keeping only these for now


# adding patch info, year, and dispersal mode to plant data
srs_all <- srs_all %>%
  left_join(year_created, by = c("block")) %>%
  left_join(patch_info, by = c("block", "patch")) %>%
  left_join(dispersal_mode, by = c("sppcode" = "Species.Code"))


srs_all %>%
  count(sppcode)
# matching dispersal mode to EDI dispersal mode to look at missing values #
# loading species info from EDI
#edi_dispersal_url <- "https://pasta.lternet.edu/package/data/eml/edi/414/1/8de4a490a6ac6b05d2406c975d25b649"
#edi_dispersal <- read.csv(file = srs_dispersal_url)

#srs_all %>%
#  filter(is.na(DispMode1)) %>% # looking to see which species don't have dispersal mode info
#  count(sppcode) %>%
#  left_join(edi_dispersal, by = c("sppcode" = "SppCode")) # looking to see if EDI data solves any of these

#### spp codes, dispersal modes ####
# resolving spp codes
srs_all <- srs_all %>%
  mutate(sppcode = dplyr::case_when(
    sppcode %in% c("DESPER") ~ "DESGLA", # changing DESPER to DESGLA at Christopher's suggestion -- email 05/24/2025
    sppcode %in% c("SOLPTY", "SOLNIG", "SOLAME") ~ "SOLNIG", # grouping these at Christopher's suggestion -- email 05/24/2025
    sppcode %in% c("DIGSSPP") ~ "DIGSPP", # fixing typo -- 7 letter code DIGSSPP to 6 letter code DIGSPP
    .default = sppcode
  ))


# adding missing dispersal modes -- doing this manually, going to delete EDI portion later
srs_all <- srs_all %>%
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
srs_all <- srs_all %>%
  mutate(dispersal_mode = dplyr::case_when(
    DispMode1 %in% c("Adhesive", "Animal", "Ant", "Bird", "Mammal", "MammalBird") ~ "Animal",
    DispMode1 %in% c("Tumbling", "Wind") ~ "Wind",
    DispMode1 %in% c("Ballistic", "Gravity", "Unassisted") ~ "Gravity",
    .default = DispMode1
  ))

# fixing transplant column
srs_all <- srs_all %>%
  mutate(transplant = dplyr::case_when(
    sppcode %in% c("ARIBEY", "SORSEC", "ANTVIL", "CARBEL", "LIAEAR", "PHYAME",
                   "SOLAME", "NOLGEO", "GAYDUM", "RUDHIR", "LANCAM", "ILEVER") ~ TRUE,
    .default = FALSE
  )) %>%
  dplyr::select(!c("transplant.")) # removing old transplant column



# adding survey time since patch creation
srs_all <- srs_all %>%
  mutate(time = year - year.created)

#### soil moisture ####
# soil moisture data -- need to decide how to combine -- use 2007 for all blocks except the 75s??
# data from 2003
soil_moisture_2003 <- soil_moisture_2003 %>%
  rename(block = EU, patch = Patch, plot = Plot) %>%
  group_by(block, patch) %>%
  summarise(soil_moisture_2003 = mean(Pct_moisture_by_wt))

# data from 2007
soil_moisture_2007 <- soil_moisture_2007 %>%
  mutate(wet = wet.weight - tin.plus.filter.paper) %>%
  mutate(dry = dry.weight - tin.plus.filter.paper) %>%
  mutate(soil_moisture = (wet-dry)/dry*100) %>%
  filter(wet.weight < 22500.00) %>% # removing outlier
  filter(wet.weight > 116) %>% # removing outlier, dry weight more than wet weight
  group_by(block, patch) %>%
  summarise(soil_moisture_2007 = mean(soil_moisture))

# joining 2003 and 2007 data together
soil_moisture <- soil_moisture_2007 %>%
  full_join(soil_moisture_2003, by = c("block", "patch"))

# looking ar correlations between soil moisture values between 2003 and 2007
soil_moisture2 <- soil_moisture %>%
  filter(!is.na(soil_moisture_2003)) %>%
  filter(!is.na(soil_moisture_2007))
cor(soil_moisture2$soil_moisture_2007, soil_moisture2$soil_moisture_2003) # not that great, ~0.6 correlation

# deciding to use 2007 data, except for the 75s which will use 2003 data
soil_moisture <- soil_moisture %>%
  mutate(soil_moisture = if_else(is.na(soil_moisture_2007), soil_moisture_2003, soil_moisture_2007))

# joining to full dataset
srs_all <- srs_all %>%
  left_join(soil_moisture, by = c("block", "patch"))


#### year since fire ####
year_fire <- year_fire %>%
  dplyr::select(!notes)

srs_all <- srs_all %>%
  left_join(year_fire, by = c("block", "year"))


rare_species <- srs_all %>%
 # filter(!block %in% c("75E", "75W")) %>%
  count(sppcode) %>%
  arrange(n) %>%
  mutate(rare = if_else(n < 10, 0, 1))


#### final dataset ####
# renaming, removing unneeded columns, and rearranging
srs_all <- srs_all %>%
  #filter(!block %in% c("75E", "75W")) %>%
  left_join(rare_species, by = c("sppcode")) %>%
  mutate(year_created = year.created) %>% # renaming
  mutate(unique_id = paste(block, patch, patch_type, sep = "-")) %>% # creating unique ID for each patch
  dplyr::select(block, patch, cell, patch_type, unique_id, soil_moisture, core, year, year_created, time,
                sppcode, transplant, rare, dispersal_mode, year_since_fire)

# checking for missing values and duplicates
summarytools::view(summarytools::dfSummary(srs_all), footnote = NA)

# writing cleaned file
write_csv(srs_all, file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))



