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
soil_moisture_2003 <- read.csv(file = file.path("data", "L0_original", "soil_moisture_2003.csv"))
soil_moisture_2007 <- read.csv(file = file.path("data", "L0_original", "soil_moisture_demo_plots.csv"))

# making column names lower case for joining and preference
colnames(srs_all) <- stringr::str_to_lower(colnames(srs_all))

# changing "8" to "08" to be consistent 
year_created <- year_created %>%
  mutate(block = if_else(block == "8", "08", block))
patch_info <- patch_info %>%
  mutate(block = if_else(block == "8", "08", block))
soil_moisture_2007 <- soil_moisture_2007 %>%
  mutate(block = if_else(block == "8", "08", block))

# getting rid of columns of dispersal mode that I don't need
dispersal_mode <- dispersal_mode %>%
  dplyr::select(c("Species.Code", "USDACode", "DispMode1", "DispMode2"))


# adding patch info, year, and dispersal mode to plant data
srs_all <- srs_all %>%
  left_join(year_created, by = c("block")) %>%
  left_join(patch_info, by = c("block", "patch")) %>%
  left_join(dispersal_mode, by = c("sppcode" = "Species.Code"))




# matching dispersal mode to EDI dispersal mode to look at missing values #
# loading species info from EDI
#edi_dispersal_url <- "https://pasta.lternet.edu/package/data/eml/edi/414/1/8de4a490a6ac6b05d2406c975d25b649"
#edi_dispersal <- read.csv(file = srs_dispersal_url)

#srs_all %>%
#  filter(is.na(DispMode1)) %>% # looking to see which species don't have dispersal mode info
#  count(sppcode) %>%
#  left_join(edi_dispersal, by = c("sppcode" = "SppCode")) # looking to see if EDI data solves any of these


# adding missing dispersal modes
srs_all <- srs_all %>%
  mutate(DispMode1 = dplyr::case_when(
    sppcode %in% c("ALLCUT", "GALMOL", "MANVIR", "PASSPP", "STRUMB") ~ "Gravity",
    sppcode %in% c("ANDSPP", "CATBIG", "CIRHOR", "CUSSPP", 
                   "LIASCG", "OSMREG", "SABDIF", "SENSMA", "TRAGOP") ~ "Wind",
    sppcode %in% c("CARSPP", "GALSPP", "LANTAN", "PLUSPP", "SOLNIG") ~ "Animal",
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


# adding survey time since patch creation
srs_all <- srs_all %>%
  mutate(time = year - year.created)

# soil moisture data -- need to decide how to combine -- use 2007 for all blocks except the 75s??
### NOTE FOR LATER: check to see how different values are between 2003 and 2007
soil_moisture_2003 %>%
  count(EU, Patch)
summarytools::view(summarytools::dfSummary(soil_moisture_2003), footnote = NA)

soil_moisture_2007 %>%
  count(block, patch)
summarytools::view(summarytools::dfSummary(soil_moisture_2007), footnote = NA)

# renaming, removing unneeded columns, and rearranging
srs_all <- srs_all %>%
  mutate(year_created = year.created) %>% # renaming
  mutate(unique_id = paste(block, patch, patch_type, sep = "-")) %>% # creating unique ID for each patch
  dplyr::select(block, patch, cell, patch_type, unique_id, core, year, year_created, time,
                sppcode, dispersal_mode)

# checking for missing values and duplicates
summarytools::view(summarytools::dfSummary(srs_all), footnote = NA)

# writing cleaned file
write_csv(srs_all, file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))






