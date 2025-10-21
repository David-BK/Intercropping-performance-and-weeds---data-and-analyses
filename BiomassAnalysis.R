## Author: David Kottelenberg
## Institute: Wageningen University & Research
## Last edited: 21-10-2025
## Summary: Analyses of 2022-2024 field experiments, biomass

rm(list = ls())

library(tidyverse)
library(readxl)
library(RColorBrewer)
library(lme4)
library(car)
library(scales)
library(caret)
library(gridExtra)
library(segmented)
library(lmtest)
library(stringr)
library(ggbeeswarm)
library(ggpubr)
library(patchwork)
library(emmeans)
library(broom.mixed)
library(purrr)
library(lmerTest)
library(grid)
library(xtable)

# Set directory
# setwd("")

source("WCFBiomassFunctions.R")
source("WCFYieldFunctions.R")

sowing_date2022 <- as.Date("2022-04-19")
sowing_date2023 <- as.Date("2023-03-02")
sowing_date2024 <- as.Date("2023-04-12")

colour_vector <- c(brewer.pal(n = 8, "Set1"), "#FDC086")

model_colors <- c("Cereal" = colour_vector[1],
                  "Legume" = colour_vector[2],
                  "Weed" = colour_vector[4],
                  "T" = colour_vector[1],
                  "Triticale" = colour_vector[1],
                  "F" = colour_vector[2],
                  "Faba" = colour_vector[2],
                  "TF" = colour_vector[3],
                  "C" = colour_vector[1],
                  "L" = colour_vector[2],
                  "CL" = colour_vector[3],
                  "TA" = colour_vector[5],
                  "TM" = colour_vector[9],
                  "FA" = colour_vector[7],
                  "FM" =  colour_vector[8],
                  "Intercrop" = colour_vector[3],
                  "Sole crop" = colour_vector[2],
                  "Sole cereal" = colour_vector[1],
                  "Sole legume" = colour_vector[2],
                  "Sole triticale" = colour_vector[1],
                  "Sole faba" = colour_vector[2],
                  "Sole weeds" = colour_vector[4],
                  "Weed infested" = colour_vector[2],
                  "Herbicide treated" = colour_vector[1],
                  "W" = colour_vector[4],
                  "Y" = colour_vector[4],
                  "N" = colour_vector[5])

#####################
##                 ##
## Experiment 2022 ##
##                 ##
#####################

##
## First harvest
##

# Read data
data_2022_biomass_1 <- read_xlsx("WCF_2022_data.xlsx", sheet = "Biomass1", range = "A1:H77", col_names = TRUE) %>% 
  dplyr::select(Plot, Treatment, "Block", "BiomassLegume", "BiomassCereal") %>% 
  rename(BiomassC = BiomassCereal, BiomassL = BiomassLegume) %>% 
  arrange(Plot) %>% 
  filter(!grepl("Lupine", Treatment)) %>% # Lupine was abandoned in the experiment due to a reaction to the herbicide, and no space left in the no-herbicide parts of the plots
  dplyr::select(Plot, Treatment, Block, BiomassC, BiomassL) %>% 
  mutate(Treatment = case_when(
    Treatment == "Barley" ~ "B",
    Treatment == "Rye" ~ "R",
    Treatment == "Triticale" ~ "T",
    Treatment == "Wheat" ~ "W",
    Treatment == "Pea" ~ "P",
    Treatment == "Faba" ~ "F",
    Treatment == "Barley_Faba" ~ "1B:1F",
    Treatment == "Rye_Faba" ~ "1R:1F",
    Treatment == "Triticale_Faba" ~ "1T:1F",
    Treatment == "Wheat_Faba" ~ "1W:1F",
    Treatment == "Barley_Pea" ~ "1B:1P",
    Treatment == "Rye_Pea" ~ "1R:1P",
    Treatment == "Triticale_Pea" ~ "1T:1P",
    Treatment == "Wheat_Pea" ~ "1W:1P",
    TRUE ~ Treatment
  )) %>% 
  mutate(BiomassC = ifelse(BiomassC == 0, NA, BiomassC),
    BiomassL = ifelse(BiomassL == 0, NA, BiomassL))

##
## Second harvest
##

# Read data
data_2022_biomass_2 <- read_xlsx("WCF_2022_data.xlsx", sheet = "Biomass2", range = "A1:J77", col_names = TRUE) %>% 
  dplyr::select(Plot, Treatment, "Block", "BiomassLegumeNW", "BiomassCerealNW", "BiomassLegumeWW", "BiomassCerealWW")

data_2022_biomass_2NW <- data_2022_biomass_2 %>% 
  dplyr::select(-BiomassLegumeWW, -BiomassCerealWW) %>% 
  rename(BiomassC = BiomassCerealNW, BiomassL = BiomassLegumeNW) %>% 
  arrange(Plot) %>% 
  filter(!grepl("Lupine", Treatment)) %>% # Lupine was abandoned in the experiment due to a reaction to the herbicide, and no space left in the no-herbicide parts of the plots
  dplyr::select(Plot, Treatment, Block, BiomassC, BiomassL) %>% 
  mutate(Treatment = case_when(
    Treatment == "Barley" ~ "B",
    Treatment == "Rye" ~ "R",
    Treatment == "Triticale" ~ "T",
    Treatment == "Wheat" ~ "W",
    Treatment == "Pea" ~ "P",
    Treatment == "Faba" ~ "F",
    Treatment == "Barley_Faba" ~ "1B:1F",
    Treatment == "Rye_Faba" ~ "1R:1F",
    Treatment == "Triticale_Faba" ~ "1T:1F",
    Treatment == "Wheat_Faba" ~ "1W:1F",
    Treatment == "Barley_Pea" ~ "1B:1P",
    Treatment == "Rye_Pea" ~ "1R:1P",
    Treatment == "Triticale_Pea" ~ "1T:1P",
    Treatment == "Wheat_Pea" ~ "1W:1P",
    TRUE ~ Treatment
  )) %>% 
  mutate(BiomassC = ifelse(BiomassC == 0, NA, BiomassC),
    BiomassL = ifelse(BiomassL == 0, NA, BiomassL)) %>% 
  mutate(BiomassC = as.numeric(BiomassC),
         BiomassL = as.numeric(BiomassL))

data_2022_biomass_2WW <- data_2022_biomass_2 %>% 
    dplyr::select(-BiomassLegumeNW, -BiomassCerealNW) %>% 
  rename(BiomassC = BiomassCerealWW, BiomassL = BiomassLegumeWW) %>% 
  arrange(Plot) %>% 
  filter(!grepl("Lupine", Treatment)) %>% # Lupine was abandoned in the experiment due to a reaction to the herbicide, and no space left in the no-herbicide parts of the plots
  dplyr::select(Plot, Treatment, Block, BiomassC, BiomassL) %>% 
  mutate(Treatment = case_when(
    Treatment == "Barley" ~ "B",
    Treatment == "Rye" ~ "R",
    Treatment == "Triticale" ~ "T",
    Treatment == "Wheat" ~ "W",
    Treatment == "Pea" ~ "P",
    Treatment == "Faba" ~ "F",
    Treatment == "Barley_Faba" ~ "1B:1F",
    Treatment == "Rye_Faba" ~ "1R:1F",
    Treatment == "Triticale_Faba" ~ "1T:1F",
    Treatment == "Wheat_Faba" ~ "1W:1F",
    Treatment == "Barley_Pea" ~ "1B:1P",
    Treatment == "Rye_Pea" ~ "1R:1P",
    Treatment == "Triticale_Pea" ~ "1T:1P",
    Treatment == "Wheat_Pea" ~ "1W:1P",
    TRUE ~ Treatment
  )) %>% 
  mutate(BiomassC = ifelse(BiomassC == 0, NA, BiomassC),
    BiomassL = ifelse(BiomassL == 0, NA, BiomassL)) %>% 
  mutate(BiomassC = as.numeric(BiomassC),
         BiomassL = as.numeric(BiomassL))


##
## Final harvest
##

# Read data
data_2022_biomass_3 <- read_xlsx("WCF_2022_data.xlsx", sheet = "FinalHarvest", range = "A1:P77", col_names = TRUE) %>% 
  dplyr::select(Plot, Treatment, "Block", "BiomassLegume", "BiomassCereal") %>% 
  rename(BiomassC = BiomassCereal, BiomassL = BiomassLegume) %>% 
  arrange(Plot) %>% 
  filter(!grepl("Lupine", Treatment)) %>% # Lupine was abandoned in the experiment due to a reaction to the herbicide, and no space left in the no-herbicide parts of the plots
  dplyr::select(Plot, Treatment, Block, BiomassC, BiomassL) %>% 
  mutate(Treatment = case_when(
    Treatment == "Barley" ~ "B",
    Treatment == "Rye" ~ "R",
    Treatment == "Triticale" ~ "T",
    Treatment == "Wheat" ~ "W",
    Treatment == "Pea" ~ "P",
    Treatment == "Faba" ~ "F",
    Treatment == "Barley_Faba" ~ "1B:1F",
    Treatment == "Rye_Faba" ~ "1R:1F",
    Treatment == "Triticale_Faba" ~ "1T:1F",
    Treatment == "Wheat_Faba" ~ "1W:1F",
    Treatment == "Barley_Pea" ~ "1B:1P",
    Treatment == "Rye_Pea" ~ "1R:1P",
    Treatment == "Triticale_Pea" ~ "1T:1P",
    Treatment == "Wheat_Pea" ~ "1W:1P",
    TRUE ~ Treatment
  )) %>% 
  mutate(BiomassC = ifelse(BiomassC == 0, NA, BiomassC),
    BiomassL = ifelse(BiomassL == 0, NA, BiomassL)) %>% 
  mutate(BiomassC = as.numeric(BiomassC),
         BiomassL = as.numeric(BiomassL))


#####################
##                 ##
## Experiment 2023 ##
##                 ##
#####################

##
## First harvest
##

# Read data
data2023_A1 <- read_xlsx("WCF_2023_data.xlsx", sheet = "BiomassA1", range = "A1:M41", col_names = TRUE) %>% 
  rename(Date = HarvestDate) %>% 
  dplyr::select(Plot, Treatment, Block, Date, TriticaleBiomassA, TriticaleBiomassM, TriticaleBiomass, FabaBiomassA, FabaBiomassM, FabaBiomass, WeedBiomass) %>% 
  rename(TA = TriticaleBiomassA, TM = TriticaleBiomassM, T = TriticaleBiomass, FA = FabaBiomassA, FM = FabaBiomassM, F = FabaBiomass, W = WeedBiomass) %>% 
  mutate(T = ifelse(is.na(T) & grepl("T", Treatment), TA + TM, T),
         F = ifelse(is.na(F) & grepl("F", Treatment), FA + FM, F),
         W = as.numeric(W)) %>% 
  arrange(Plot) %>% 
  rename(BiomassC = T, BiomassL = F) %>% 
  dplyr::select(Plot, Treatment, Block, BiomassC, BiomassL) %>% 
  mutate(BiomassC = as.numeric(BiomassC),
         BiomassL = as.numeric(BiomassL))

##
## Second harvest
##

# Read data
data2023_A2 <- read_xlsx("WCF_2023_data.xlsx", sheet = "BiomassA2", range = "A1:M81", col_names = TRUE) %>% 
  rename(Date = HarvestDate) %>% 
  dplyr::select(Plot, Treatment, Block, Weeds, Date, TriticaleBiomassA, TriticaleBiomassM, TriticaleBiomass, FabaBiomassA, FabaBiomassM, FabaBiomass, WeedBiomass) %>% 
  rename(TA = TriticaleBiomassA, TM = TriticaleBiomassM, T = TriticaleBiomass, FA = FabaBiomassA, FM = FabaBiomassM, F = FabaBiomass, W = WeedBiomass) %>% 
  mutate(T = ifelse(is.na(T) & grepl("T", Treatment), TA + TM, T),
         F = ifelse(is.na(F) & grepl("F", Treatment), FA + FM, F),
         W = as.numeric(W)) %>% 
  arrange(Plot) %>% 
  rename(BiomassC = T, BiomassL = F) %>% 
  dplyr::select(Plot, Treatment, Block, Weeds, BiomassC, BiomassL)

data2023_A2WW <- data2023_A2 %>% 
  filter(Weeds == "Y") %>% 
  mutate(BiomassC = as.numeric(BiomassC),
         BiomassL = as.numeric(BiomassL))

data2023_A2NW <- data2023_A2 %>% 
  filter(Weeds == "N") %>% 
  mutate(BiomassC = as.numeric(BiomassC),
         BiomassL = as.numeric(BiomassL))


##
## Final harvest
##

# Read data
data2023_A3 <- read_xlsx("WCF_2023_data.xlsx", sheet = "FinalHarvestA", range = "A1:AJ81", col_names = TRUE) %>%
  dplyr::select(-TreatmentN) %>% 
  mutate(BiomassFabaA = as.numeric(BiomassFabaA),
         BiomassFabaM = as.numeric(BiomassFabaM)) %>% 
  mutate(BiomassTriticale = ifelse(is.na(BiomassTriticale) & grepl("T", Treatment), BiomassTriticaleA + BiomassTriticaleM, BiomassTriticale),
         BiomassFaba = ifelse(is.na(BiomassFaba) & grepl("F", Treatment), BiomassFabaA + BiomassFabaM, BiomassFaba)) %>% 
  arrange(Weeds, Plot) %>% 
  dplyr::select(Plot, Treatment, Block, Weeds, BiomassTriticaleA, BiomassTriticaleM, BiomassTriticale, BiomassFabaA, BiomassFabaM, BiomassFaba, BiomassWeed) %>% 
  rename(TA = BiomassTriticaleA, TM = BiomassTriticaleM, T = BiomassTriticale, FA = BiomassFabaA, FM = BiomassFabaM, F = BiomassFaba, W = BiomassWeed) %>% 
  rename(BiomassC = T, BiomassL = F) %>% 
  dplyr::select(Plot, Treatment, Block, Weeds, BiomassC, BiomassL)
  
data2023_A3WW <- data2023_A3 %>% 
  filter(Weeds == "Y") %>% 
  mutate(BiomassC = as.numeric(BiomassC),
         BiomassL = as.numeric(BiomassL))

data2023_A3NW <- data2023_A3 %>% 
  filter(Weeds == "N") %>% 
  mutate(BiomassC = as.numeric(BiomassC),
         BiomassL = as.numeric(BiomassL))

######################
##                  ##
## Experiment 2023B ##
##                  ##
######################

##
## First harvest
##

# Read data
data2023_B1 <- read_xlsx("WCF_2023_data.xlsx", sheet = "BiomassB1", range = "A1:I45", col_names = TRUE) %>% 
  dplyr::select(Plot, Treatment, Block, BiomassTriticale, BiomassFaba, BiomassWeed) %>% 
  rename(T = BiomassTriticale, F = BiomassFaba, W = BiomassWeed) %>%  
  arrange(Plot) %>% 
  filter(!Treatment %in% c("T-25", "1T:1F-25", "TF-M-25", "F-25")) %>% 
  rename(BiomassC = T, BiomassL = F) %>% 
  dplyr::select(Plot, Treatment, Block, BiomassC, BiomassL) %>% 
  mutate(BiomassC = as.numeric(BiomassC),
         BiomassL = as.numeric(BiomassL))

##
## Second harvest
##

data2023_B2 <- read_xlsx("WCF_2023_data.xlsx", sheet = "BiomassB2", range = "A1:I89", col_names = TRUE) %>% 
  dplyr::select(Plot, Treatment, Block, Weeds, BiomassTriticale, BiomassFaba, BiomassWeed) %>% 
  rename(T = BiomassTriticale, F = BiomassFaba, W = BiomassWeed) %>%  
  arrange(Plot) %>% 
  filter(!Treatment %in% c("T-25", "1T:1F-25", "TF-M-25", "F-25")) %>% 
  rename(BiomassC = T, BiomassL = F) %>% 
  dplyr::select(Plot, Treatment, Block, Weeds, BiomassC, BiomassL)

data2023_B2WW <- data2023_B2 %>% 
  filter(Weeds == "Y") %>% 
  mutate(BiomassC = as.numeric(BiomassC),
         BiomassL = as.numeric(BiomassL))

data2023_B2NW <- data2023_B2 %>% 
  filter(Weeds == "N") %>% 
  mutate(BiomassC = as.numeric(BiomassC),
         BiomassL = as.numeric(BiomassL))

##
## Final harvest
##

# Read data
data2023_B3 <- read_xlsx("WCF_2023_data.xlsx", sheet = "FinalHarvestB", range = "A1:P89", col_names = TRUE) %>%
  dplyr::select(-TreatmentN) %>% 
  mutate(BiomassFaba = as.numeric(BiomassFaba),
         BiomassFaba = as.numeric(BiomassFaba)) %>% 
  arrange(Weeds, Plot) %>% 
  dplyr::select(Plot, Treatment, Block, Weeds, BiomassTriticale, BiomassFaba, BiomassWeed) %>% 
  rename(T = BiomassTriticale, F = BiomassFaba, W = BiomassWeed) %>% 
  filter(!Treatment %in% c("T-25", "1T:1F-25", "TF-M-25", "F-25")) %>% 
  rename(BiomassC = T, BiomassL = F) %>% 
  dplyr::select(Plot, Treatment, Block, Weeds, BiomassC, BiomassL)

data2023_B3WW <- data2023_B3 %>% 
  filter(Weeds == "Y") %>% 
  mutate(BiomassC = as.numeric(BiomassC),
         BiomassL = as.numeric(BiomassL))

data2023_B3NW <- data2023_B3 %>% 
  filter(Weeds == "N") %>% 
  mutate(BiomassC = as.numeric(BiomassC),
         BiomassL = as.numeric(BiomassL))

#####################
##                 ##
## Experiment 2024 ##
##                 ##
#####################

##
## First harvest
##

# Read data
data2024_1 <- read_xlsx("WCF_2024_data.xlsx", sheet = "Biomass1", range = "A1:F41", col_names = TRUE) %>% 
  dplyr::select(Plot, Treatment, Block, Crop, Biomass) %>%
  pivot_wider(names_from = "Crop", values_from = Biomass) %>% 
  arrange(Plot) %>% 
  rename(BiomassC = T, BiomassL = F) %>% 
  dplyr::select(Plot, Treatment, Block, BiomassC, BiomassL) %>% 
  mutate(BiomassC = as.numeric(BiomassC),
         BiomassL = as.numeric(BiomassL))

##
## Second harvest
##

# Read data
data2024_2 <- read_xlsx("WCF_2024_data.xlsx", sheet = "Biomass2", range = "A1:G61", col_names = TRUE) %>% 
  dplyr::select(Plot, Treatment, Block, Weeds, Crop, Biomass) %>%
  pivot_wider(names_from = "Crop", values_from = Biomass) %>% 
  arrange(Plot) %>% 
  rename(BiomassC = T, BiomassL = F) %>% 
  dplyr::select(Plot, Treatment, Block, Weeds, BiomassC, BiomassL)

data2024_2WW <- data2024_2 %>% 
  filter(Weeds == "Y") %>% 
  mutate(BiomassC = as.numeric(BiomassC),
         BiomassL = as.numeric(BiomassL))

data2024_2NW <- data2024_2 %>% 
  filter(Weeds == "N") %>% 
  mutate(BiomassC = as.numeric(BiomassC),
         BiomassL = as.numeric(BiomassL))

##
## Final harvest
##

# Read data
data2024_3 <- read_xlsx("WCF_2024_data.xlsx", sheet = "FinalHarvest", range = "A1:O41", col_names = TRUE) %>%
  arrange(Weeds, Plot) %>% 
  dplyr::select(Plot, Treatment, Block, Weeds, BiomassTriticale, BiomassFaba, BiomassWeed) %>% 
  rename(T = BiomassTriticale, F = BiomassFaba, W = BiomassWeed) %>% 
  rename(BiomassC = T, BiomassL = F) %>% 
  dplyr::select(Plot, Treatment, Block, Weeds, BiomassC, BiomassL)

data2024_3WW <- data2024_3 %>% 
  filter(Weeds == "Y") %>% 
  mutate(BiomassC = as.numeric(BiomassC),
         BiomassL = as.numeric(BiomassL))

data2024_3NW <- data2024_3 %>% 
  filter(Weeds == "N") %>% 
  mutate(BiomassC = as.numeric(BiomassC),
         BiomassL = as.numeric(BiomassL))

#####################
## Create plots     ##
#####################

## 2022, first harvest

data_2022_biomass_groups_1 <- data_2022_biomass_1 %>%
  mutate(
    CerealGroup = case_when(
      str_detect(Treatment, "(^|1)R") ~ "R",
      str_detect(Treatment, "(^|1)B") ~ "B",
      str_detect(Treatment, "(^|1)T") ~ "T",
      str_detect(Treatment, "(^|1)W") ~ "W",
      TRUE ~ NA_character_
    ),
    LegumeGroup = case_when(
      str_detect(Treatment, "^P$") ~ "P",
      str_detect(Treatment, "^F$") ~ "F",
      str_detect(Treatment, ":1P") ~ "P",
      str_detect(Treatment, ":1F") ~ "F",
      TRUE ~ NA_character_
    )
  )

data_2022_biomass_cereal_1 <- data_2022_biomass_groups_1 %>%
  filter(!is.na(CerealGroup)) %>%
  group_by(CerealGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassC, na.rm = TRUE),
    sd_biomass = sd(BiomassC, na.rm = TRUE),
    n = sum(!is.na(BiomassC)),
    .groups = "drop"
  )

data_2022_biomass_legume_1 <- data_2022_biomass_groups_1 %>%
  filter(!is.na(LegumeGroup)) %>%
  group_by(LegumeGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassL, na.rm = TRUE),
    sd_biomass = sd(BiomassL, na.rm = TRUE),
    n = sum(!is.na(BiomassL)),
    .groups = "drop"
  )

cereal_levels <- c(
  "R", "1R:1P", "1R:1F",
  "B", "1B:1P", "1B:1F",
  "T", "1T:1P", "1T:1F",
  "W", "1W:1P", "1W:1F"
)

data_2022_biomass_cereal2_1 <- data_2022_biomass_cereal_1 %>%
  mutate(Treatment = factor(Treatment, levels = cereal_levels))

legume_levels <- c(
  "P", "1R:1P", "1B:1P", "1T:1P", "1W:1P",
  "F", "1R:1F", "1B:1F", "1T:1F", "1W:1F"
)

data_2022_biomass_legume2_1 <- data_2022_biomass_legume_1 %>%
  mutate(Treatment = factor(Treatment, levels = legume_levels))

cereal_cld_2022 <- lapply(unique(data_2022_biomass_cereal2_1$CerealGroup), function(g) {
  treatments <- cereal_levels[grepl(g, cereal_levels)]
  df_sub <- data_2022_biomass_1 %>%
    filter(Treatment %in% treatments & !is.na(BiomassC))
    cld_df <- get_cld(df_sub, "BiomassC", "CerealGroup", treatments)
    cld_df$CerealGroup <- g
    return(cld_df)
}) %>% bind_rows()

legume_cld_2022 <- lapply(unique(data_2022_biomass_legume2_1$LegumeGroup), function(g) {
  if(g == "P") {
    treatments <- legume_levels[legume_levels %in% c("P", "1R:1P", "1B:1P", "1T:1P", "1W:1P")]
  } else if(g == "F") {
    treatments <- legume_levels[legume_levels %in% c("F", "1R:1F", "1B:1F", "1T:1F", "1W:1F")]
  } else {
    treatments <- character(0)
  }
  df_sub <- data_2022_biomass_1 %>%
    filter(Treatment %in% treatments & !is.na(BiomassL))
  if(nrow(df_sub) > 0) {
    cld_df <- get_cld(df_sub, "BiomassL", "LegumeGroup", treatments)
    cld_df$LegumeGroup <- g
    return(cld_df)
  }
}) %>% bind_rows()

data_2022_biomass_cereal_cld_1 <- data_2022_biomass_cereal2_1 %>%
  left_join(cereal_cld_2022, by = c("Treatment", "CerealGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2022_biomass_legume_cld_1 <- data_2022_biomass_legume2_1 %>%
  left_join(legume_cld_2022, by = c("Treatment", "LegumeGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

spacer_cereals <- tibble(
  CerealGroup = c("R", "B", "T"),
  Treatment = factor(c("spacer1", "spacer2", "spacer3"), levels = c(levels(data_2022_biomass_cereal_cld_1$Treatment), "spacer1", "spacer2", "spacer3")),
  mean_biomass = NA_real_,
  sd_biomass = NA_real_,
  n = NA_integer_,
  .group = NA_character_
)

data_2022_biomass_cereal_cld2_1 <- data_2022_biomass_cereal_cld_1 %>%
  bind_rows(spacer_cereals) %>%
  mutate(Treatment = factor(Treatment, levels = c("R", "1R:1F", "1R:1P", "spacer1",
                                                  "B", "1B:1F", "1B:1P", "spacer2",
                                                  "T", "1T:1F", "1T:1P", "spacer3",
                                                  "W", "1W:1F", "1W:1P"))) %>%
  mutate(CerealGroup = factor(CerealGroup, levels = c("R", "B", "T", "W")))

spacer_legumes <- tibble(
  LegumeGroup = c("P"),
  Treatment = factor("spacer1", levels = c(levels(data_2022_biomass_legume_cld_1$Treatment), "spacer1")),
  mean_biomass = NA_real_,
  sd_biomass = NA_real_,
  n = NA_integer_,
  .group = NA_character_
)

data_2022_biomass_legume_cld2_1 <- data_2022_biomass_legume_cld_1 %>%
  bind_rows(spacer_legumes) %>%
  mutate(Treatment = factor(Treatment, levels = c("F", "1R:1F", "1B:1F", "1T:1F", "1W:1F", 
                                                  "spacer1",
                                                  "P", "1R:1P", "1B:1P", "1T:1P", "1W:1P"
                                                  ))) %>%
  mutate(LegumeGroup = factor(LegumeGroup, levels = c("F", "P")))

plot_cereal_2022_1 <- ggplot(data_2022_biomass_cereal_cld2_1, aes(x = Treatment, y = mean_biomass, fill = CerealGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 3), na.rm = TRUE, size = 6) +
  scale_x_discrete(labels = function(x) ifelse(x %in% c("spacer1", "spacer2", "spacer3"), "", x)) +
  scale_fill_manual(
    name = "Cereal species",
    values = c(
      "R" = colour_vector[1],
      "1R:1P" = colour_vector[1],
      "1R:1F" = colour_vector[1],
      "B" = colour_vector[2],
      "1B:1P" = colour_vector[2],
      "1B:1F" = colour_vector[2],
      "T" = colour_vector[3],
      "1T:1P" = colour_vector[3],
      "1T:1F" = colour_vector[3],
      "W" = colour_vector[4],
      "1W:1P" = colour_vector[4],
      "1W:1F" = colour_vector[4]
    ),
    labels = c(
        "R" = "Rye",
        "B" = "Barley",
        "T" = "Triticale",
        "W" = "Wheat"
    )
  ) +
  labs(
    subtitle = "Cereal biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    plot.subtitle = element_text(hjust = 0.5))

plot_legume_2022_1 <- ggplot(data_2022_biomass_legume_cld2_1, aes(x = Treatment, y = mean_biomass, fill = LegumeGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 2), na.rm = TRUE, size = 6) +
  scale_x_discrete(labels = function(x) ifelse(x %in% c("spacer1"), "", x)) +
  scale_fill_manual(
    name = "Legume species",
    values = c(
      "P" = colour_vector[5],
      "1R:1P" = colour_vector[5],
      "1B:1P" = colour_vector[5],
      "1T:1P" = colour_vector[5],
      "1W:1P" = colour_vector[5],
      "F" = colour_vector[7],
      "1R:1F" = colour_vector[7],
      "1B:1F" = colour_vector[7],
      "1T:1F" = colour_vector[7],
      "1W:1F" = colour_vector[7]
    ),
    labels = c(
        "P" = "Pea",
        "F" = "Faba bean"
    )
  ) +
  labs(
    subtitle = "Legume biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    plot.subtitle = element_text(hjust = 0.5))

combined_plot_2022_1 <- plot_cereal_2022_1 / plot_legume_2022_1 +
    plot_layout(ncol = 1
    ) +
    plot_annotation(tag_levels = "a")
combined_plot_2022_1

## 2022, second harvest, no weeds

data_2022_biomass_groups_2NW <- data_2022_biomass_2NW %>%
  mutate(
    CerealGroup = case_when(
      str_detect(Treatment, "(^|1)R") ~ "R",
      str_detect(Treatment, "(^|1)B") ~ "B",
      str_detect(Treatment, "(^|1)T") ~ "T",
      str_detect(Treatment, "(^|1)W") ~ "W",
      TRUE ~ NA_character_
    ),
    LegumeGroup = case_when(
      str_detect(Treatment, "^P$") ~ "P",
      str_detect(Treatment, "^F$") ~ "F",
      str_detect(Treatment, ":1P") ~ "P",
      str_detect(Treatment, ":1F") ~ "F",
      TRUE ~ NA_character_
    )
  )

data_2022_biomass_cereal_2NW <- data_2022_biomass_groups_2NW %>%
  filter(!is.na(CerealGroup)) %>%
  group_by(CerealGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassC, na.rm = TRUE),
    sd_biomass = sd(BiomassC, na.rm = TRUE),
    n = sum(!is.na(BiomassC)),
    .groups = "drop"
  )

data_2022_biomass_legume_2NW <- data_2022_biomass_groups_2NW %>%
  filter(!is.na(LegumeGroup)) %>%
  group_by(LegumeGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassL, na.rm = TRUE),
    sd_biomass = sd(BiomassL, na.rm = TRUE),
    n = sum(!is.na(BiomassL)),
    .groups = "drop"
  )

cereal_levels <- c(
  "R", "1R:1P", "1R:1F",
  "B", "1B:1P", "1B:1F",
  "T", "1T:1P", "1T:1F",
  "W", "1W:1P", "1W:1F"
)

data_2022_biomass_cereal2_2NW <- data_2022_biomass_cereal_2NW %>%
  mutate(Treatment = factor(Treatment, levels = cereal_levels))

legume_levels <- c(
  "P", "1R:1P", "1B:1P", "1T:1P", "1W:1P",
  "F", "1R:1F", "1B:1F", "1T:1F", "1W:1F"
)

data_2022_biomass_legume2_2NW <- data_2022_biomass_legume_2NW %>%
  mutate(Treatment = factor(Treatment, levels = legume_levels))

cereal_cld_2022 <- lapply(unique(data_2022_biomass_cereal2_2NW$CerealGroup), function(g) {
  treatments <- cereal_levels[grepl(g, cereal_levels)]
  df_sub <- data_2022_biomass_2NW %>%
    filter(Treatment %in% treatments & !is.na(BiomassC))
    cld_df <- get_cld(df_sub, "BiomassC", "CerealGroup", treatments)
    cld_df$CerealGroup <- g
    return(cld_df)
}) %>% bind_rows()

legume_cld_2022 <- lapply(unique(data_2022_biomass_legume2_2NW$LegumeGroup), function(g) {
  if(g == "P") {
    treatments <- legume_levels[legume_levels %in% c("P", "1R:1P", "1B:1P", "1T:1P", "1W:1P")]
  } else if(g == "F") {
    treatments <- legume_levels[legume_levels %in% c("F", "1R:1F", "1B:1F", "1T:1F", "1W:1F")]
  } else {
    treatments <- character(0)
  }
  df_sub <- data_2022_biomass_2NW %>%
    filter(Treatment %in% treatments & !is.na(BiomassL))
  if(nrow(df_sub) > 0) {
    cld_df <- get_cld(df_sub, "BiomassL", "LegumeGroup", treatments)
    cld_df$LegumeGroup <- g
    return(cld_df)
  }
}) %>% bind_rows()

data_2022_biomass_cereal_cld_2NW <- data_2022_biomass_cereal2_2NW %>%
  left_join(cereal_cld_2022, by = c("Treatment", "CerealGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2022_biomass_legume_cld_2NW <- data_2022_biomass_legume2_2NW %>%
  left_join(legume_cld_2022, by = c("Treatment", "LegumeGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

spacer_cereals <- tibble(
  CerealGroup = c("R", "B", "T"),
  Treatment = factor(c("spacer1", "spacer2", "spacer3"), levels = c(levels(data_2022_biomass_cereal_cld_2NW$Treatment), "spacer1", "spacer2", "spacer3")),
  mean_biomass = NA_real_,
  sd_biomass = NA_real_,
  n = NA_integer_,
  .group = NA_character_
)

data_2022_biomass_cereal_cld2_2NW <- data_2022_biomass_cereal_cld_2NW %>%
  bind_rows(spacer_cereals) %>%
  mutate(Treatment = factor(Treatment, levels = c("R", "1R:1F", "1R:1P", "spacer1",
                                                  "B", "1B:1F", "1B:1P", "spacer2",
                                                  "T", "1T:1F", "1T:1P", "spacer3",
                                                  "W", "1W:1F", "1W:1P"))) %>%
  mutate(CerealGroup = factor(CerealGroup, levels = c("R", "B", "T", "W")))

spacer_legumes <- tibble(
  LegumeGroup = c("P"),
  Treatment = factor("spacer1", levels = c(levels(data_2022_biomass_legume_cld_2NW$Treatment), "spacer1")),
  mean_biomass = NA_real_,
  sd_biomass = NA_real_,
  n = NA_integer_,
  .group = NA_character_
)

data_2022_biomass_legume_cld2_2NW <- data_2022_biomass_legume_cld_2NW %>%
  bind_rows(spacer_legumes) %>%
  mutate(Treatment = factor(Treatment, levels = c("F", "1R:1F", "1B:1F", "1T:1F", "1W:1F", 
                                                  "spacer1",
                                                  "P", "1R:1P", "1B:1P", "1T:1P", "1W:1P"
                                                  ))) %>%
  mutate(LegumeGroup = factor(LegumeGroup, levels = c("F", "P")))

plot_cereal_2022_2NW <- ggplot(data_2022_biomass_cereal_cld2_2NW, aes(x = Treatment, y = mean_biomass, fill = CerealGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 60), na.rm = TRUE, size = 6) +
  scale_x_discrete(labels = function(x) ifelse(x %in% c("spacer1", "spacer2", "spacer3"), "", x)) +
  scale_fill_manual(
    name = "Cereal species",
    values = c(
      "R" = colour_vector[1],
      "1R:1P" = colour_vector[1],
      "1R:1F" = colour_vector[1],
      "B" = colour_vector[2],
      "1B:1P" = colour_vector[2],
      "1B:1F" = colour_vector[2],
      "T" = colour_vector[3],
      "1T:1P" = colour_vector[3],
      "1T:1F" = colour_vector[3],
      "W" = colour_vector[4],
      "1W:1P" = colour_vector[4],
      "1W:1F" = colour_vector[4]
    ),
    labels = c(
        "R" = "Rye",
        "B" = "Barley",
        "T" = "Triticale",
        "W" = "Wheat"
    )
  ) +
  labs(
    subtitle = "Cereal biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    plot.subtitle = element_text(hjust = 0.5))

plot_legume_2022_2NW <- ggplot(data_2022_biomass_legume_cld2_2NW, aes(x = Treatment, y = mean_biomass, fill = LegumeGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 50), na.rm = TRUE, size = 6) +
  scale_x_discrete(labels = function(x) ifelse(x %in% c("spacer1"), "", x)) +
  scale_fill_manual(
    name = "Legume species",
    values = c(
      "P" = colour_vector[5],
      "1R:1P" = colour_vector[5],
      "1B:1P" = colour_vector[5],
      "1T:1P" = colour_vector[5],
      "1W:1P" = colour_vector[5],
      "F" = colour_vector[7],
      "1R:1F" = colour_vector[7],
      "1B:1F" = colour_vector[7],
      "1T:1F" = colour_vector[7],
      "1W:1F" = colour_vector[7]
    ),
    labels = c(
        "P" = "Pea",
        "F" = "Faba bean"
    )
  ) +
  labs(
    subtitle = "Legume biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    plot.subtitle = element_text(hjust = 0.5))

combined_plot_2022_2NW <- plot_cereal_2022_2NW / plot_legume_2022_2NW +
    plot_layout(ncol = 1
    ) +
    plot_annotation(tag_levels = "a")
combined_plot_2022_2NW

## 2022, second harvest, with weeds
data_2022_biomass_groups_2WW <- data_2022_biomass_2WW %>%
  mutate(
    CerealGroup = case_when(
      str_detect(Treatment, "(^|1)R") ~ "R",
      str_detect(Treatment, "(^|1)B") ~ "B",
      str_detect(Treatment, "(^|1)T") ~ "T",
      str_detect(Treatment, "(^|1)W") ~ "W",
      TRUE ~ NA_character_
    ),
    LegumeGroup = case_when(
      str_detect(Treatment, "^P$") ~ "P",
      str_detect(Treatment, "^F$") ~ "F",
      str_detect(Treatment, ":1P") ~ "P",
      str_detect(Treatment, ":1F") ~ "F",
      TRUE ~ NA_character_
    )
  )

data_2022_biomass_cereal_2WW <- data_2022_biomass_groups_2WW %>%
  filter(!is.na(CerealGroup)) %>%
  group_by(CerealGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassC, na.rm = TRUE),
    sd_biomass = sd(BiomassC, na.rm = TRUE),
    n = sum(!is.na(BiomassC)),
    .groups = "drop"
  )

data_2022_biomass_legume_2WW <- data_2022_biomass_groups_2WW %>%
  filter(!is.na(LegumeGroup)) %>%
  group_by(LegumeGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassL, na.rm = TRUE),
    sd_biomass = sd(BiomassL, na.rm = TRUE),
    n = sum(!is.na(BiomassL)),
    .groups = "drop"
  )

cereal_levels <- c(
  "R", "1R:1P", "1R:1F",
  "B", "1B:1P", "1B:1F",
  "T", "1T:1P", "1T:1F",
  "W", "1W:1P", "1W:1F"
)

data_2022_biomass_cereal2_2WW <- data_2022_biomass_cereal_2WW %>%
  mutate(Treatment = factor(Treatment, levels = cereal_levels))

legume_levels <- c(
  "P", "1R:1P", "1B:1P", "1T:1P", "1W:1P",
  "F", "1R:1F", "1B:1F", "1T:1F", "1W:1F"
)

data_2022_biomass_legume2_2WW <- data_2022_biomass_legume_2WW %>%
  mutate(Treatment = factor(Treatment, levels = legume_levels))

cereal_cld_2022 <- lapply(unique(data_2022_biomass_cereal2_2WW$CerealGroup), function(g) {
  treatments <- cereal_levels[grepl(g, cereal_levels)]
  df_sub <- data_2022_biomass_2WW %>%
    filter(Treatment %in% treatments & !is.na(BiomassC))
    cld_df <- get_cld(df_sub, "BiomassC", "CerealGroup", treatments)
    cld_df$CerealGroup <- g
    return(cld_df)
}) %>% bind_rows()

legume_cld_2022 <- lapply(unique(data_2022_biomass_legume2_2WW$LegumeGroup), function(g) {
  if(g == "P") {
    treatments <- legume_levels[legume_levels %in% c("P", "1R:1P", "1B:1P", "1T:1P", "1W:1P")]
  } else if(g == "F") {
    treatments <- legume_levels[legume_levels %in% c("F", "1R:1F", "1B:1F", "1T:1F", "1W:1F")]
  } else {
    treatments <- character(0)
  }
  df_sub <- data_2022_biomass_2WW %>%
    filter(Treatment %in% treatments & !is.na(BiomassL))
  if(nrow(df_sub) > 0) {
    cld_df <- get_cld(df_sub, "BiomassL", "LegumeGroup", treatments)
    cld_df$LegumeGroup <- g
    return(cld_df)
  }
}) %>% bind_rows()

data_2022_biomass_cereal_cld_2WW <- data_2022_biomass_cereal2_2WW %>%
  left_join(cereal_cld_2022, by = c("Treatment", "CerealGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2022_biomass_legume_cld_2WW <- data_2022_biomass_legume2_2WW %>%
  left_join(legume_cld_2022, by = c("Treatment", "LegumeGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

spacer_cereals <- tibble(
  CerealGroup = c("R", "B", "T"),
  Treatment = factor(c("spacer1", "spacer2", "spacer3"), levels = c(levels(data_2022_biomass_cereal_cld_2WW$Treatment), "spacer1", "spacer2", "spacer3")),
  mean_biomass = NA_real_,
  sd_biomass = NA_real_,
  n = NA_integer_,
  .group = NA_character_
)

data_2022_biomass_cereal_cld2_2WW <- data_2022_biomass_cereal_cld_2WW %>%
  bind_rows(spacer_cereals) %>%
  mutate(Treatment = factor(Treatment, levels = c("R", "1R:1F", "1R:1P", "spacer1",
                                                  "B", "1B:1F", "1B:1P", "spacer2",
                                                  "T", "1T:1F", "1T:1P", "spacer3",
                                                  "W", "1W:1F", "1W:1P"))) %>%
  mutate(CerealGroup = factor(CerealGroup, levels = c("R", "B", "T", "W")))

spacer_legumes <- tibble(
  LegumeGroup = c("P"),
  Treatment = factor("spacer1", levels = c(levels(data_2022_biomass_legume_cld_2WW$Treatment), "spacer1")),
  mean_biomass = NA_real_,
  sd_biomass = NA_real_,
  n = NA_integer_,
  .group = NA_character_
)

data_2022_biomass_legume_cld2_2WW <- data_2022_biomass_legume_cld_2WW %>%
  bind_rows(spacer_legumes) %>%
  mutate(Treatment = factor(Treatment, levels = c("F", "1R:1F", "1B:1F", "1T:1F", "1W:1F", 
                                                  "spacer1",
                                                  "P", "1R:1P", "1B:1P", "1T:1P", "1W:1P"
                                                  ))) %>%
  mutate(LegumeGroup = factor(LegumeGroup, levels = c("F", "P")))

plot_cereal_2022_2WW <- ggplot(data_2022_biomass_cereal_cld2_2WW, aes(x = Treatment, y = mean_biomass, fill = CerealGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 60), na.rm = TRUE, size = 6) +
  scale_x_discrete(labels = function(x) ifelse(x %in% c("spacer1", "spacer2", "spacer3"), "", x)) +
  scale_fill_manual(
    name = "Cereal species",
    values = c(
      "R" = colour_vector[1],
      "1R:1P" = colour_vector[1],
      "1R:1F" = colour_vector[1],
      "B" = colour_vector[2],
      "1B:1P" = colour_vector[2],
      "1B:1F" = colour_vector[2],
      "T" = colour_vector[3],
      "1T:1P" = colour_vector[3],
      "1T:1F" = colour_vector[3],
      "W" = colour_vector[4],
      "1W:1P" = colour_vector[4],
      "1W:1F" = colour_vector[4]
    ),
    labels = c(
        "R" = "Rye",
        "B" = "Barley",
        "T" = "Triticale",
        "W" = "Wheat"
    )
  ) +
  labs(
    subtitle = "Cereal biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    plot.subtitle = element_text(hjust = 0.5))

plot_legume_2022_2WW <- ggplot(data_2022_biomass_legume_cld2_2WW, aes(x = Treatment, y = mean_biomass, fill = LegumeGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 50), na.rm = TRUE, size = 6) +
  scale_x_discrete(labels = function(x) ifelse(x %in% c("spacer1"), "", x)) +
  scale_fill_manual(
    name = "Legume species",
    values = c(
      "P" = colour_vector[5],
      "1R:1P" = colour_vector[5],
      "1B:1P" = colour_vector[5],
      "1T:1P" = colour_vector[5],
      "1W:1P" = colour_vector[5],
      "F" = colour_vector[7],
      "1R:1F" = colour_vector[7],
      "1B:1F" = colour_vector[7],
      "1T:1F" = colour_vector[7],
      "1W:1F" = colour_vector[7]
    ),
    labels = c(
        "P" = "Pea",
        "F" = "Faba bean"
    )
  ) +
  labs(
    subtitle = "Legume biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    plot.subtitle = element_text(hjust = 0.5))

combined_plot_2022_2WW <- plot_cereal_2022_2WW / plot_legume_2022_2WW +
    plot_layout(ncol = 1
    ) +
    plot_annotation(tag_levels = "a")
combined_plot_2022_2WW

## 2022, final harvest
data_2022_biomass_groups_3 <- data_2022_biomass_3 %>%
  mutate(
    CerealGroup = case_when(
      str_detect(Treatment, "(^|1)R") ~ "R",
      str_detect(Treatment, "(^|1)B") ~ "B",
      str_detect(Treatment, "(^|1)T") ~ "T",
      str_detect(Treatment, "(^|1)W") ~ "W",
      TRUE ~ NA_character_
    ),
    LegumeGroup = case_when(
      str_detect(Treatment, "^P$") ~ "P",
      str_detect(Treatment, "^F$") ~ "F",
      str_detect(Treatment, ":1P") ~ "P",
      str_detect(Treatment, ":1F") ~ "F",
      TRUE ~ NA_character_
    )
  )

data_2022_biomass_cereal_3 <- data_2022_biomass_groups_3 %>%
  filter(!is.na(CerealGroup)) %>%
  group_by(CerealGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassC, na.rm = TRUE),
    sd_biomass = sd(BiomassC, na.rm = TRUE),
    n = sum(!is.na(BiomassC)),
    .groups = "drop"
  )

data_2022_biomass_legume_3 <- data_2022_biomass_groups_3 %>%
  filter(!is.na(LegumeGroup)) %>%
  group_by(LegumeGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassL, na.rm = TRUE),
    sd_biomass = sd(BiomassL, na.rm = TRUE),
    n = sum(!is.na(BiomassL)),
    .groups = "drop"
  )

cereal_levels <- c(
  "R", "1R:1P", "1R:1F",
  "B", "1B:1P", "1B:1F",
  "T", "1T:1P", "1T:1F",
  "W", "1W:1P", "1W:1F"
)

data_2022_biomass_cereal2_3 <- data_2022_biomass_cereal_3 %>%
  mutate(Treatment = factor(Treatment, levels = cereal_levels))

legume_levels <- c(
  "P", "1R:1P", "1B:1P", "1T:1P", "1W:1P",
  "F", "1R:1F", "1B:1F", "1T:1F", "1W:1F"
)

data_2022_biomass_legume2_3 <- data_2022_biomass_legume_3 %>%
  mutate(Treatment = factor(Treatment, levels = legume_levels))

cereal_cld_2022 <- lapply(unique(data_2022_biomass_cereal2_3$CerealGroup), function(g) {
  treatments <- cereal_levels[grepl(g, cereal_levels)]
  df_sub <- data_2022_biomass_3 %>%
    filter(Treatment %in% treatments & !is.na(BiomassC))
    cld_df <- get_cld(df_sub, "BiomassC", "CerealGroup", treatments)
    cld_df$CerealGroup <- g
    return(cld_df)
}) %>% bind_rows()

legume_cld_2022 <- lapply(unique(data_2022_biomass_legume2_3$LegumeGroup), function(g) {
  if(g == "P") {
    treatments <- legume_levels[legume_levels %in% c("P", "1R:1P", "1B:1P", "1T:1P", "1W:1P")]
  } else if(g == "F") {
    treatments <- legume_levels[legume_levels %in% c("F", "1R:1F", "1B:1F", "1T:1F", "1W:1F")]
  } else {
    treatments <- character(0)
  }
  df_sub <- data_2022_biomass_3 %>%
    filter(Treatment %in% treatments & !is.na(BiomassL))
  if(nrow(df_sub) > 0) {
    cld_df <- get_cld(df_sub, "BiomassL", "LegumeGroup", treatments)
    cld_df$LegumeGroup <- g
    return(cld_df)
  }
}) %>% bind_rows()

data_2022_biomass_cereal_cld_3 <- data_2022_biomass_cereal2_3 %>%
  left_join(cereal_cld_2022, by = c("Treatment", "CerealGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2022_biomass_legume_cld_3 <- data_2022_biomass_legume2_3 %>%
  left_join(legume_cld_2022, by = c("Treatment", "LegumeGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

spacer_cereals <- tibble(
  CerealGroup = c("R", "B", "T"),
  Treatment = factor(c("spacer1", "spacer2", "spacer3"), levels = c(levels(data_2022_biomass_cereal_cld_3$Treatment), "spacer1", "spacer2", "spacer3")),
  mean_biomass = NA_real_,
  sd_biomass = NA_real_,
  n = NA_integer_,
  .group = NA_character_
)

data_2022_biomass_cereal_cld2_3 <- data_2022_biomass_cereal_cld_3 %>%
  bind_rows(spacer_cereals) %>%
  mutate(Treatment = factor(Treatment, levels = c("R", "1R:1F", "1R:1P", "spacer1",
                                                  "B", "1B:1F", "1B:1P", "spacer2",
                                                  "T", "1T:1F", "1T:1P", "spacer3",
                                                  "W", "1W:1F", "1W:1P"))) %>%
  mutate(CerealGroup = factor(CerealGroup, levels = c("R", "B", "T", "W")))

spacer_legumes <- tibble(
  LegumeGroup = c("P"),
  Treatment = factor("spacer1", levels = c(levels(data_2022_biomass_legume_cld_3$Treatment), "spacer1")),
  mean_biomass = NA_real_,
  sd_biomass = NA_real_,
  n = NA_integer_,
  .group = NA_character_
)

data_2022_biomass_legume_cld2_3 <- data_2022_biomass_legume_cld_3 %>%
  bind_rows(spacer_legumes) %>%
  mutate(Treatment = factor(Treatment, levels = c("F", "1R:1F", "1B:1F", "1T:1F", "1W:1F", 
                                                  "spacer1",
                                                  "P", "1R:1P", "1B:1P", "1T:1P", "1W:1P"
                                                  ))) %>%
  mutate(LegumeGroup = factor(LegumeGroup, levels = c("F", "P")))

plot_cereal_2022_3 <- ggplot(data_2022_biomass_cereal_cld2_3, aes(x = Treatment, y = mean_biomass, fill = CerealGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 100), na.rm = TRUE, size = 6) +
  scale_x_discrete(labels = function(x) ifelse(x %in% c("spacer1", "spacer2", "spacer3"), "", x)) +
  scale_fill_manual(
    name = "Cereal species",
    values = c(
      "R" = colour_vector[1],
      "1R:1P" = colour_vector[1],
      "1R:1F" = colour_vector[1],
      "B" = colour_vector[2],
      "1B:1P" = colour_vector[2],
      "1B:1F" = colour_vector[2],
      "T" = colour_vector[3],
      "1T:1P" = colour_vector[3],
      "1T:1F" = colour_vector[3],
      "W" = colour_vector[4],
      "1W:1P" = colour_vector[4],
      "1W:1F" = colour_vector[4]
    ),
    labels = c(
        "R" = "Rye",
        "B" = "Barley",
        "T" = "Triticale",
        "W" = "Wheat"
    )
  ) +
  labs(
    subtitle = "Cereal biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    plot.subtitle = element_text(hjust = 0.5))

plot_legume_2022_3 <- ggplot(data_2022_biomass_legume_cld2_3, aes(x = Treatment, y = mean_biomass, fill = LegumeGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 75), na.rm = TRUE, size = 6) +
  scale_x_discrete(labels = function(x) ifelse(x %in% c("spacer1"), "", x)) +
  scale_fill_manual(
    name = "Legume species",
    values = c(
      "P" = colour_vector[5],
      "1R:1P" = colour_vector[5],
      "1B:1P" = colour_vector[5],
      "1T:1P" = colour_vector[5],
      "1W:1P" = colour_vector[5],
      "F" = colour_vector[7],
      "1R:1F" = colour_vector[7],
      "1B:1F" = colour_vector[7],
      "1T:1F" = colour_vector[7],
      "1W:1F" = colour_vector[7]
    ),
    labels = c(
        "P" = "Pea",
        "F" = "Faba bean"
    )
  ) +
  labs(
    subtitle = "Legume biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    plot.subtitle = element_text(hjust = 0.5))

combined_plot_2022_3 <- plot_cereal_2022_3 / plot_legume_2022_3 +
    plot_layout(ncol = 1
    ) +
    plot_annotation(tag_levels = "a")
combined_plot_2022_3

## 2023A, first harvest

data_2023A_biomass_NW_1 <- data2023_A1 %>%
  mutate(Treatment = gsub("TF", "1T:1F", Treatment))

data_2023A_biomass_groups_1 <- data_2023A_biomass_NW_1 %>%
  mutate(
    CerealGroup = ifelse(grepl("T", Treatment), "T", NA_character_),
    LegumeGroup = ifelse(grepl("F", Treatment), "F", NA_character_)
  )

data_2023A_biomass_cereal_1 <- data_2023A_biomass_groups_1 %>%
  filter(!is.na(CerealGroup)) %>%
  group_by(CerealGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassC, na.rm = TRUE),
    sd_biomass = sd(BiomassC, na.rm = TRUE),
    n = sum(!is.na(BiomassC)),
    .groups = "drop"
  )

data_2023A_biomass_legume_1 <- data_2023A_biomass_groups_1 %>%
  filter(!is.na(LegumeGroup)) %>%
  group_by(LegumeGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassL, na.rm = TRUE),
    sd_biomass = sd(BiomassL, na.rm = TRUE),
    n = sum(!is.na(BiomassL)),
    .groups = "drop"
  )

cereal_levels <- c(
  "T+", "T", "3T:1F", "1T:1F", "1T:1F-M", "1T:3F"
)

data_2023A_biomass_cereal2_1 <- data_2023A_biomass_cereal_1 %>%
  mutate(Treatment = factor(Treatment, levels = cereal_levels))

legume_levels <- c(
  "F+", "F", "1T:3F", "1T:1F", "1T:1F-M", "3T:1F"
)

data_2023A_biomass_legume2_1 <- data_2023A_biomass_legume_1 %>%
  mutate(Treatment = factor(Treatment, levels = legume_levels))

cereal_cld_2023A_1 <- lapply(unique(data_2023A_biomass_cereal2_1$CerealGroup), function(g) {
  treatments <- cereal_levels[grepl(g, cereal_levels)]
  df_sub <- data_2023A_biomass_NW_1 %>%
    filter(Treatment %in% treatments & !is.na(BiomassC))
    cld_df <- get_cld(df_sub, "BiomassC", "CerealGroup", treatments)
    cld_df$CerealGroup <- g
    return(cld_df)
}) %>% bind_rows()

legume_cld_2023A_1 <- lapply(unique(data_2023A_biomass_legume2_1$LegumeGroup), function(g) {
  treatments <- legume_levels[grepl(g, legume_levels)]
  df_sub <- data_2023A_biomass_NW_1 %>%
    filter(Treatment %in% treatments & !is.na(BiomassL))
  if(nrow(df_sub) > 0) {
    cld_df <- get_cld2(df_sub, "BiomassL", "LegumeGroup", treatments)
    cld_df$LegumeGroup <- g
    return(cld_df)
  }
}) %>% bind_rows()

data_2023A_biomass_cereal_cld_1 <- data_2023A_biomass_cereal2_1 %>%
  left_join(cereal_cld_2023A_1, by = c("Treatment", "CerealGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2023A_biomass_legume_cld_1 <- data_2023A_biomass_legume2_1 %>%
  left_join(legume_cld_2023A_1, by = c("Treatment", "LegumeGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2023A_biomass_cereal_cld2_1 <- data_2023A_biomass_cereal_cld_1 %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "T+", "T", "3T:1F", "1T:1F", "1T:1F-M", "1T:3F"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("T+", "T", "3T:1F", "1T:1F", "1T:1F-M", "1T:3F"))) %>% 
  mutate(CerealGroup = factor(CerealGroup, levels = c("T")))

data_2023A_biomass_legume_cld2_1 <- data_2023A_biomass_legume_cld_1 %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "F+", "F", "1T:3F", "1T:1F", "1T:1F-M", "3T:1F"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("F+", "F", "1T:3F", "1T:1F", "1T:1F-M", "3T:1F"
                                                  ))) %>% 
  mutate(LegumeGroup = factor(LegumeGroup, levels = c("F")))

plot_cereal_2023A_1 <- ggplot(data_2023A_biomass_cereal_cld2_1, aes(x = Treatment, y = mean_biomass, fill = CerealGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 8), na.rm = TRUE, size = 6) +
  scale_x_discrete(labels = function(x) ifelse(x %in% c("spacer1", "spacer2", "spacer3"), "", x)) +
  scale_fill_manual(
    name = "Cereal species",
    values = c(
      "T+" = colour_vector[1],
      "T" = colour_vector[1],
      "3T:1F" = colour_vector[1],
      "1T:1F" = colour_vector[1],
      "1T:1F-M" = colour_vector[1],
      "1T:3F" = colour_vector[1]
    ),
    labels = c(
        "T" = "Triticale"
    )
  ) +
  labs(
    subtitle = "Triticale biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

plot_legume_2023A_1 <- ggplot(data_2023A_biomass_legume_cld2_1, aes(x = Treatment, y = mean_biomass, fill = LegumeGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 3), na.rm = TRUE, size = 6) +
  scale_x_discrete(labels = function(x) ifelse(x %in% c("spacer1"), "", x)) +
  scale_fill_manual(
    name = "Legume species",
    values = c(
      "F+" = colour_vector[2],
      "F" = colour_vector[2],
      "1T:3F" = colour_vector[2],
      "1T:1F" = colour_vector[2],
      "1T:1F-M" = colour_vector[2],
      "3T:1F" = colour_vector[2]
    ),
    labels = c(
        "F" = "Faba bean"
    )
  ) +
  labs(
    subtitle = "Faba bean biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

combined_plot_2023A_1 <- plot_cereal_2023A_1 / plot_legume_2023A_1 +
    plot_layout(ncol = 1
    ) +
    plot_annotation(tag_levels = "a")
combined_plot_2023A_1

## 2023A, second harvest, no weeds

data_2023A_biomass_2NW <- data2023_A2NW %>%
  mutate(Treatment = gsub("TF", "1T:1F", Treatment))

data_2023A_biomass_groups_2NW <- data_2023A_biomass_2NW %>%
  mutate(
    CerealGroup = ifelse(grepl("T", Treatment), "T", NA_character_),
    LegumeGroup = ifelse(grepl("F", Treatment), "F", NA_character_)
  )

data_2023A_biomass_cereal_2NW <- data_2023A_biomass_groups_2NW %>%
  filter(!is.na(CerealGroup)) %>%
  group_by(CerealGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassC, na.rm = TRUE),
    sd_biomass = sd(BiomassC, na.rm = TRUE),
    n = sum(!is.na(BiomassC)),
    .groups = "drop"
  )

data_2023A_biomass_legume_2NW <- data_2023A_biomass_groups_2NW %>%
  filter(!is.na(LegumeGroup)) %>%
  group_by(LegumeGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassL, na.rm = TRUE),
    sd_biomass = sd(BiomassL, na.rm = TRUE),
    n = sum(!is.na(BiomassL)),
    .groups = "drop"
  )

cereal_levels <- c(
  "T+", "T", "3T:1F", "1T:1F", "1T:1F-M", "1T:3F"
)

data_2023A_biomass_cereal2_2NW <- data_2023A_biomass_cereal_2NW %>%
  mutate(Treatment = factor(Treatment, levels = cereal_levels))

legume_levels <- c(
  "F+", "F", "1T:3F", "1T:1F", "1T:1F-M", "3T:1F"
)

data_2023A_biomass_legume2_2NW <- data_2023A_biomass_legume_2NW %>%
  mutate(Treatment = factor(Treatment, levels = legume_levels))

cereal_cld_2023A_2NW <- lapply(unique(data_2023A_biomass_cereal2_2NW$CerealGroup), function(g) {
  treatments <- cereal_levels[grepl(g, cereal_levels)]
  df_sub <- data_2023A_biomass_2NW %>%
    filter(Treatment %in% treatments & !is.na(BiomassC))
    cld_df <- get_cld(df_sub, "BiomassC", "CerealGroup", treatments)
    cld_df$CerealGroup <- g
    return(cld_df)
}) %>% bind_rows()

legume_cld_2023A_2NW <- lapply(unique(data_2023A_biomass_legume2_2NW$LegumeGroup), function(g) {
  treatments <- legume_levels[grepl(g, legume_levels)]
  df_sub <- data_2023A_biomass_2NW %>%
    filter(Treatment %in% treatments & !is.na(BiomassL))
  if(nrow(df_sub) > 0) {
    cld_df <- get_cld2(df_sub, "BiomassL", "LegumeGroup", treatments)
    cld_df$LegumeGroup <- g
    return(cld_df)
  }
}) %>% bind_rows()

data_2023A_biomass_cereal_cld_2NW <- data_2023A_biomass_cereal2_2NW %>%
  left_join(cereal_cld_2023A_2NW, by = c("Treatment", "CerealGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2023A_biomass_legume_cld_2NW <- data_2023A_biomass_legume2_2NW %>%
  left_join(legume_cld_2023A_2NW, by = c("Treatment", "LegumeGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2023A_biomass_cereal_cld2_2NW <- data_2023A_biomass_cereal_cld_2NW %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "T+", "T", "3T:1F", "1T:1F", "1T:1F-M", "1T:3F"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("T+", "T", "3T:1F", "1T:1F", "1T:1F-M", "1T:3F"))) %>% 
  mutate(CerealGroup = factor(CerealGroup, levels = c("T")))

data_2023A_biomass_legume_cld2_2NW <- data_2023A_biomass_legume_cld_2NW %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "F+", "F", "1T:3F", "1T:1F", "1T:1F-M", "3T:1F"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("F+", "F", "1T:3F", "1T:1F", "1T:1F-M", "3T:1F"
                                                  ))) %>% 
  mutate(LegumeGroup = factor(LegumeGroup, levels = c("F")))

plot_cereal_2023A_2NW <- ggplot(data_2023A_biomass_cereal_cld2_2NW, aes(x = Treatment, y = mean_biomass, fill = CerealGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 90), na.rm = TRUE, size = 6) +
  scale_x_discrete(labels = function(x) ifelse(x %in% c("spacer1", "spacer2", "spacer3"), "", x)) +
  scale_fill_manual(
    name = "Cereal species",
    values = c(
      "T+" = colour_vector[1],
      "T" = colour_vector[1],
      "3T:1F" = colour_vector[1],
      "1T:1F" = colour_vector[1],
      "1T:1F-M" = colour_vector[1],
      "1T:3F" = colour_vector[1]
    ),
    labels = c(
        "T" = "Triticale"
    )
  ) +
  labs(
    subtitle = "Triticale biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

plot_legume_2023A_2NW <- ggplot(data_2023A_biomass_legume_cld2_2NW, aes(x = Treatment, y = mean_biomass, fill = LegumeGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 60), na.rm = TRUE, size = 6) +
  scale_x_discrete(labels = function(x) ifelse(x %in% c("spacer1"), "", x)) +
  scale_fill_manual(
    name = "Legume species",
    values = c(
      "F+" = colour_vector[2],
      "F" = colour_vector[2],
      "1T:3F" = colour_vector[2],
      "1T:1F" = colour_vector[2],
      "1T:1F-M" = colour_vector[2],
      "3T:1F" = colour_vector[2]
    ),
    labels = c(
        "F" = "Faba bean"
    )
  ) +
  labs(
    subtitle = "Faba bean biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

combined_plot_2023A_2NW <- plot_cereal_2023A_2NW / plot_legume_2023A_2NW +
    plot_layout(ncol = 1
    ) +
    plot_annotation(tag_levels = "a")
combined_plot_2023A_2NW

## 2023A, second harvest, with weeds

data_2023A_biomass_2WW <- data2023_A2WW %>%
  mutate(Treatment = gsub("TF", "1T:1F", Treatment))

data_2023A_biomass_groups_2WW <- data_2023A_biomass_2WW %>%
  mutate(
    CerealGroup = ifelse(grepl("T", Treatment), "T", NA_character_),
    LegumeGroup = ifelse(grepl("F", Treatment), "F", NA_character_)
  )

data_2023A_biomass_cereal_2WW <- data_2023A_biomass_groups_2WW %>%
  filter(!is.na(CerealGroup)) %>%
  group_by(CerealGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassC, na.rm = TRUE),
    sd_biomass = sd(BiomassC, na.rm = TRUE),
    n = sum(!is.na(BiomassC)),
    .groups = "drop"
  )

data_2023A_biomass_legume_2WW <- data_2023A_biomass_groups_2WW %>%
  filter(!is.na(LegumeGroup)) %>%
  group_by(LegumeGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassL, na.rm = TRUE),
    sd_biomass = sd(BiomassL, na.rm = TRUE),
    n = sum(!is.na(BiomassL)),
    .groups = "drop"
  )

cereal_levels <- c(
  "T+", "T", "3T:1F", "1T:1F", "1T:1F-M", "1T:3F"
)

data_2023A_biomass_cereal2_2WW <- data_2023A_biomass_cereal_2WW %>%
  mutate(Treatment = factor(Treatment, levels = cereal_levels))

legume_levels <- c(
  "F+", "F", "1T:3F", "1T:1F", "1T:1F-M", "3T:1F"
)

data_2023A_biomass_legume2_2WW <- data_2023A_biomass_legume_2WW %>%
  mutate(Treatment = factor(Treatment, levels = legume_levels))

cereal_cld_2023A_2WW <- lapply(unique(data_2023A_biomass_cereal2_2WW$CerealGroup), function(g) {
  treatments <- cereal_levels[grepl(g, cereal_levels)]
  df_sub <- data_2023A_biomass_2WW %>%
    filter(Treatment %in% treatments & !is.na(BiomassC))
    cld_df <- get_cld(df_sub, "BiomassC", "CerealGroup", treatments)
    cld_df$CerealGroup <- g
    return(cld_df)
}) %>% bind_rows()

legume_cld_2023A_2WW <- lapply(unique(data_2023A_biomass_legume2_2WW$LegumeGroup), function(g) {
  treatments <- legume_levels[grepl(g, legume_levels)]
  df_sub <- data_2023A_biomass_2WW %>%
    filter(Treatment %in% treatments & !is.na(BiomassL))
  if(nrow(df_sub) > 0) {
    cld_df <- get_cld2(df_sub, "BiomassL", "LegumeGroup", treatments)
    cld_df$LegumeGroup <- g
    return(cld_df)
  }
}) %>% bind_rows()

data_2023A_biomass_cereal_cld_2WW <- data_2023A_biomass_cereal2_2WW %>%
  left_join(cereal_cld_2023A_2WW, by = c("Treatment", "CerealGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2023A_biomass_legume_cld_2WW <- data_2023A_biomass_legume2_2WW %>%
  left_join(legume_cld_2023A_2WW, by = c("Treatment", "LegumeGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2023A_biomass_cereal_cld2_2WW <- data_2023A_biomass_cereal_cld_2WW %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "T+", "T", "3T:1F", "1T:1F", "1T:1F-M", "1T:3F"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("T+", "T", "3T:1F", "1T:1F", "1T:1F-M", "1T:3F"))) %>% 
  mutate(CerealGroup = factor(CerealGroup, levels = c("T")))

data_2023A_biomass_legume_cld2_2WW <- data_2023A_biomass_legume_cld_2WW %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "F+", "F", "1T:3F", "1T:1F", "1T:1F-M", "3T:1F"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("F+", "F", "1T:3F", "1T:1F", "1T:1F-M", "3T:1F"
                                                  ))) %>% 
  mutate(LegumeGroup = factor(LegumeGroup, levels = c("F")))

plot_cereal_2023A_2WW <- ggplot(data_2023A_biomass_cereal_cld2_2WW, aes(x = Treatment, y = mean_biomass, fill = CerealGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 90), na.rm = TRUE, size = 6) +
  scale_x_discrete(labels = function(x) ifelse(x %in% c("spacer1", "spacer2", "spacer3"), "", x)) +
  scale_fill_manual(
    name = "Cereal species",
    values = c(
      "T+" = colour_vector[1],
      "T" = colour_vector[1],
      "3T:1F" = colour_vector[1],
      "1T:1F" = colour_vector[1],
      "1T:1F-M" = colour_vector[1],
      "1T:3F" = colour_vector[1]
    ),
    labels = c(
        "T" = "Triticale"
    )
  ) +
  labs(
    subtitle = "Triticale biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

plot_legume_2023A_2WW <- ggplot(data_2023A_biomass_legume_cld2_2WW, aes(x = Treatment, y = mean_biomass, fill = LegumeGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 60), na.rm = TRUE, size = 6) +
  scale_x_discrete(labels = function(x) ifelse(x %in% c("spacer1"), "", x)) +
  scale_fill_manual(
    name = "Legume species",
    values = c(
      "F+" = colour_vector[2],
      "F" = colour_vector[2],
      "1T:3F" = colour_vector[2],
      "1T:1F" = colour_vector[2],
      "1T:1F-M" = colour_vector[2],
      "3T:1F" = colour_vector[2]
    ),
    labels = c(
        "F" = "Faba bean"
    )
  ) +
  labs(
    subtitle = "Faba bean biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

combined_plot_2023A_2WW <- plot_cereal_2023A_2WW / plot_legume_2023A_2WW +
    plot_layout(ncol = 1
    ) +
    plot_annotation(tag_levels = "a")
combined_plot_2023A_2WW


## 2023A, final harvest, no weeds
data_2023A_biomass_3NW <- data2023_A3NW %>% 
  mutate(Treatment = gsub("TF", "1T:1F", Treatment))


data_2023A_biomass_groups_3NW <- data_2023A_biomass_3NW %>%
  mutate(
    CerealGroup = ifelse(grepl("T", Treatment), "T", NA_character_),
    LegumeGroup = ifelse(grepl("F", Treatment), "F", NA_character_)
  )

data_2023A_biomass_cereal_3NW <- data_2023A_biomass_groups_3NW %>%
  filter(!is.na(CerealGroup)) %>%
  group_by(CerealGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassC, na.rm = TRUE),
    sd_biomass = sd(BiomassC, na.rm = TRUE),
    n = sum(!is.na(BiomassC)),
    .groups = "drop"
  )

data_2023A_biomass_legume_3NW <- data_2023A_biomass_groups_3NW %>%
  filter(!is.na(LegumeGroup)) %>%
  group_by(LegumeGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassL, na.rm = TRUE),
    sd_biomass = sd(BiomassL, na.rm = TRUE),
    n = sum(!is.na(BiomassL)),
    .groups = "drop"
  )

cereal_levels <- c(
  "T+", "T", "3T:1F", "1T:1F", "1T:1F-M", "1T:3F"
)

data_2023A_biomass_cereal2_3NW <- data_2023A_biomass_cereal_3NW %>%
  mutate(Treatment = factor(Treatment, levels = cereal_levels))

legume_levels <- c(
  "F+", "F", "1T:3F", "1T:1F", "1T:1F-M", "3T:1F"
)

data_2023A_biomass_legume2_3NW <- data_2023A_biomass_legume_3NW %>%
  mutate(Treatment = factor(Treatment, levels = legume_levels))

cereal_cld_2023A_3NW <- lapply(unique(data_2023A_biomass_cereal2_3NW$CerealGroup), function(g) {
  treatments <- cereal_levels[grepl(g, cereal_levels)]
  df_sub <- data_2023A_biomass_3NW %>%
    filter(Treatment %in% treatments & !is.na(BiomassC))
    cld_df <- get_cld(df_sub, "BiomassC", "CerealGroup", treatments)
    cld_df$CerealGroup <- g
    return(cld_df)
}) %>% bind_rows()

legume_cld_2023A_3NW <- lapply(unique(data_2023A_biomass_legume2_3NW$LegumeGroup), function(g) {
  treatments <- legume_levels[grepl(g, legume_levels)]
  df_sub <- data_2023A_biomass_3NW %>%
    filter(Treatment %in% treatments & !is.na(BiomassL))
  if(nrow(df_sub) > 0) {
    cld_df <- get_cld2(df_sub, "BiomassL", "LegumeGroup", treatments)
    cld_df$LegumeGroup <- g
    return(cld_df)
  }
}) %>% bind_rows()

data_2023A_biomass_cereal_cld_3NW <- data_2023A_biomass_cereal2_3NW %>%
  left_join(cereal_cld_2023A_3NW, by = c("Treatment", "CerealGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2023A_biomass_legume_cld_3NW <- data_2023A_biomass_legume2_3NW %>%
  left_join(legume_cld_2023A_3NW, by = c("Treatment", "LegumeGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2023A_biomass_cereal_cld2_3NW <- data_2023A_biomass_cereal_cld_3NW %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "T+", "T", "3T:1F", "1T:1F", "1T:1F-M", "1T:3F"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("T+", "T", "3T:1F", "1T:1F", "1T:1F-M", "1T:3F"))) %>% 
  mutate(CerealGroup = factor(CerealGroup, levels = c("T")))

data_2023A_biomass_legume_cld2_3NW <- data_2023A_biomass_legume_cld_3NW %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "F+", "F", "1T:3F", "1T:1F", "1T:1F-M", "3T:1F"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("F+", "F", "1T:3F", "1T:1F", "1T:1F-M", "3T:1F"
                                                  ))) %>% 
  mutate(LegumeGroup = factor(LegumeGroup, levels = c("F")))

plot_cereal_2023A_3NW <- ggplot(data_2023A_biomass_cereal_cld2_3NW, aes(x = Treatment, y = mean_biomass, fill = CerealGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 90), na.rm = TRUE, size = 6) +
  scale_x_discrete(labels = function(x) ifelse(x %in% c("spacer1", "spacer2", "spacer3"), "", x)) +
  scale_fill_manual(
    name = "Cereal species",
    values = c(
      "T+" = colour_vector[1],
      "T" = colour_vector[1],
      "3T:1F" = colour_vector[1],
      "1T:1F" = colour_vector[1],
      "1T:1F-M" = colour_vector[1],
      "1T:3F" = colour_vector[1]
    ),
    labels = c(
        "T" = "Triticale"
    )
  ) +
  labs(
    subtitle = "Triticale biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

plot_legume_2023A_3NW <- ggplot(data_2023A_biomass_legume_cld2_3NW, aes(x = Treatment, y = mean_biomass, fill = LegumeGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 60), na.rm = TRUE, size = 6) +
  scale_x_discrete(labels = function(x) ifelse(x %in% c("spacer1"), "", x)) +
  scale_fill_manual(
    name = "Legume species",
    values = c(
      "F+" = colour_vector[2],
      "F" = colour_vector[2],
      "1T:3F" = colour_vector[2],
      "1T:1F" = colour_vector[2],
      "1T:1F-M" = colour_vector[2],
      "3T:1F" = colour_vector[2]
    ),
    labels = c(
        "F" = "Faba bean"
    )
  ) +
  labs(
    subtitle = "Faba bean biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

combined_plot_2023A_3NW <- plot_cereal_2023A_3NW / plot_legume_2023A_3NW +
    plot_layout(ncol = 1
    ) +
    plot_annotation(tag_levels = "a")
combined_plot_2023A_3NW

## 2023A, final harvest, with weeds
data_2023A_biomass_3WW <- data2023_A3WW %>% 
  mutate(Treatment = gsub("TF", "1T:1F", Treatment))


data_2023A_biomass_groups_3WW <- data_2023A_biomass_3WW %>%
  mutate(
    CerealGroup = ifelse(grepl("T", Treatment), "T", NA_character_),
    LegumeGroup = ifelse(grepl("F", Treatment), "F", NA_character_)
  )

data_2023A_biomass_cereal_3WW <- data_2023A_biomass_groups_3WW %>%
  filter(!is.na(CerealGroup)) %>%
  group_by(CerealGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassC, na.rm = TRUE),
    sd_biomass = sd(BiomassC, na.rm = TRUE),
    n = sum(!is.na(BiomassC)),
    .groups = "drop"
  )

data_2023A_biomass_legume_3WW <- data_2023A_biomass_groups_3WW %>%
  filter(!is.na(LegumeGroup)) %>%
  group_by(LegumeGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassL, na.rm = TRUE),
    sd_biomass = sd(BiomassL, na.rm = TRUE),
    n = sum(!is.na(BiomassL)),
    .groups = "drop"
  )

cereal_levels <- c(
  "T+", "T", "3T:1F", "1T:1F", "1T:1F-M", "1T:3F"
)

data_2023A_biomass_cereal2_3WW <- data_2023A_biomass_cereal_3WW %>%
  mutate(Treatment = factor(Treatment, levels = cereal_levels))

legume_levels <- c(
  "F+", "F", "1T:3F", "1T:1F", "1T:1F-M", "3T:1F"
)

data_2023A_biomass_legume2_3WW <- data_2023A_biomass_legume_3WW %>%
  mutate(Treatment = factor(Treatment, levels = legume_levels))

cereal_cld_2023A_3WW <- lapply(unique(data_2023A_biomass_cereal2_3WW$CerealGroup), function(g) {
  treatments <- cereal_levels[grepl(g, cereal_levels)]
  df_sub <- data_2023A_biomass_3WW %>%
    filter(Treatment %in% treatments & !is.na(BiomassC))
    cld_df <- get_cld(df_sub, "BiomassC", "CerealGroup", treatments)
    cld_df$CerealGroup <- g
    return(cld_df)
}) %>% bind_rows()

legume_cld_2023A_3WW <- lapply(unique(data_2023A_biomass_legume2_3WW$LegumeGroup), function(g) {
  treatments <- legume_levels[grepl(g, legume_levels)]
  df_sub <- data_2023A_biomass_3WW %>%
    filter(Treatment %in% treatments & !is.na(BiomassL))
  if(nrow(df_sub) > 0) {
    cld_df <- get_cld2(df_sub, "BiomassL", "LegumeGroup", treatments)
    cld_df$LegumeGroup <- g
    return(cld_df)
  }
}) %>% bind_rows()

data_2023A_biomass_cereal_cld_3WW <- data_2023A_biomass_cereal2_3WW %>%
  left_join(cereal_cld_2023A_3WW, by = c("Treatment", "CerealGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2023A_biomass_legume_cld_3WW <- data_2023A_biomass_legume2_3WW %>%
  left_join(legume_cld_2023A_3WW, by = c("Treatment", "LegumeGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2023A_biomass_cereal_cld2_3WW <- data_2023A_biomass_cereal_cld_3WW %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "T+", "T", "3T:1F", "1T:1F", "1T:1F-M", "1T:3F"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("T+", "T", "3T:1F", "1T:1F", "1T:1F-M", "1T:3F"))) %>% 
  mutate(CerealGroup = factor(CerealGroup, levels = c("T")))

data_2023A_biomass_legume_cld2_3WW <- data_2023A_biomass_legume_cld_3WW %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "F+", "F", "1T:3F", "1T:1F", "1T:1F-M", "3T:1F"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("F+", "F", "1T:3F", "1T:1F", "1T:1F-M", "3T:1F"
                                                  ))) %>% 
  mutate(LegumeGroup = factor(LegumeGroup, levels = c("F")))

plot_cereal_2023A_3WW <- ggplot(data_2023A_biomass_cereal_cld2_3WW, aes(x = Treatment, y = mean_biomass, fill = CerealGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 90), na.rm = TRUE, size = 6) +
  scale_x_discrete(labels = function(x) ifelse(x %in% c("spacer1", "spacer2", "spacer3"), "", x)) +
  scale_fill_manual(
    name = "Cereal species",
    values = c(
      "T+" = colour_vector[1],
      "T" = colour_vector[1],
      "3T:1F" = colour_vector[1],
      "1T:1F" = colour_vector[1],
      "1T:1F-M" = colour_vector[1],
      "1T:3F" = colour_vector[1]
    ),
    labels = c(
        "T" = "Triticale"
    )
  ) +
  labs(
    subtitle = "Triticale biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

plot_legume_2023A_3WW <- ggplot(data_2023A_biomass_legume_cld2_3WW, aes(x = Treatment, y = mean_biomass, fill = LegumeGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 60), na.rm = TRUE, size = 6) +
  scale_x_discrete(labels = function(x) ifelse(x %in% c("spacer1"), "", x)) +
  scale_fill_manual(
    name = "Legume species",
    values = c(
      "F+" = colour_vector[2],
      "F" = colour_vector[2],
      "1T:3F" = colour_vector[2],
      "1T:1F" = colour_vector[2],
      "1T:1F-M" = colour_vector[2],
      "3T:1F" = colour_vector[2]
    ),
    labels = c(
        "F" = "Faba bean"
    )
  ) +
  labs(
    subtitle = "Faba bean biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

combined_plot_2023A_3WW <- plot_cereal_2023A_3WW / plot_legume_2023A_3WW +
    plot_layout(ncol = 1
    ) +
    plot_annotation(tag_levels = "a")
combined_plot_2023A_3WW

## 2023B, first harvest
data_2023B_biomass_1 <- data2023_B1 %>% 
  mutate(Treatment = gsub("TF", "1T:1F", Treatment))


data_2023B_biomass_groups_1 <- data_2023B_biomass_1 %>%
  mutate(
    CerealGroup = ifelse(grepl("T", Treatment), "T", NA_character_),
    LegumeGroup = ifelse(grepl("F", Treatment), "F", NA_character_)
  )

data_2023B_biomass_cereal_1 <- data_2023B_biomass_groups_1 %>%
  filter(!is.na(CerealGroup)) %>%
  group_by(CerealGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassC, na.rm = TRUE),
    sd_biomass = sd(BiomassC, na.rm = TRUE),
    n = sum(!is.na(BiomassC)),
    .groups = "drop"
  )

data_2023B_biomass_legume_1 <- data_2023B_biomass_groups_1 %>%
  filter(!is.na(LegumeGroup)) %>%
  group_by(LegumeGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassL, na.rm = TRUE),
    sd_biomass = sd(BiomassL, na.rm = TRUE),
    n = sum(!is.na(BiomassL)),
    .groups = "drop"
  )

cereal_levels <- c(
  "T", "1T:1F", "1T:1F-M", "T-375", "1T:1F-375", "1T:1F-M-375"
)

data_2023B_biomass_cereal2_1 <- data_2023B_biomass_cereal_1 %>%
  mutate(Treatment = factor(Treatment, levels = cereal_levels))

legume_levels <- c(
  "F", "1T:1F", "1T:1F-M", "F-375", "1T:1F-375", "1T:1F-M-375"
)

data_2023B_biomass_legume2_1 <- data_2023B_biomass_legume_1 %>%
  mutate(Treatment = factor(Treatment, levels = legume_levels))

cereal_cld_2023B_1 <- lapply(unique(data_2023B_biomass_cereal2_1$CerealGroup), function(g) {
  treatments <- cereal_levels[grepl(g, cereal_levels)]
  df_sub <- data_2023B_biomass_1 %>%
    filter(Treatment %in% treatments & !is.na(BiomassC))
    cld_df <- get_cld(df_sub, "BiomassC", "CerealGroup", treatments)
    cld_df$CerealGroup <- g
    return(cld_df)
}) %>% bind_rows()

legume_cld_2023B_1 <- lapply(unique(data_2023B_biomass_legume2_1$LegumeGroup), function(g) {
  treatments <- legume_levels[grepl(g, legume_levels)]
  df_sub <- data_2023B_biomass_1 %>%
    filter(Treatment %in% treatments & !is.na(BiomassL))
  if(nrow(df_sub) > 0) {
    cld_df <- get_cld2(df_sub, "BiomassL", "LegumeGroup", treatments)
    cld_df$LegumeGroup <- g
    return(cld_df)
  }
}) %>% bind_rows()

data_2023B_biomass_cereal_cld_1 <- data_2023B_biomass_cereal2_1 %>%
  left_join(cereal_cld_2023B_1, by = c("Treatment", "CerealGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2023B_biomass_legume_cld_1 <- data_2023B_biomass_legume2_1 %>%
  left_join(legume_cld_2023B_1, by = c("Treatment", "LegumeGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2023B_biomass_cereal_cld2_1 <- data_2023B_biomass_cereal_cld_1 %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "T", "1T:1F", "1T:1F-M", "T-375", "1T:1F-375", "1T:1F-M-375"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("T", "1T:1F", "1T:1F-M", "T-375", "1T:1F-375", "1T:1F-M-375"))) %>% 
  mutate(CerealGroup = factor(CerealGroup, levels = c("T")))

data_2023B_biomass_legume_cld2_1 <- data_2023B_biomass_legume_cld_1 %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "F", "F-375", "1T:1F", "1T:1F-M", "1T:1F-375", "1T:1F-M-375"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("F", "1T:1F", "1T:1F-M", "F-375", "1T:1F-375", "1T:1F-M-375"
                                                  ))) %>% 
  mutate(LegumeGroup = factor(LegumeGroup, levels = c("F")))

plot_cereal_2023B_1 <- ggplot(data_2023B_biomass_cereal_cld2_1, aes(x = Treatment, y = mean_biomass, fill = CerealGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 6), na.rm = TRUE, size = 6) +
  scale_fill_manual(
    name = "Cereal species",
    values = c(
      "T+" = colour_vector[1],
      "T" = colour_vector[1],
      "3T:1F" = colour_vector[1],
      "1T:1F" = colour_vector[1],
      "1T:1F-M" = colour_vector[1],
      "1T:3F" = colour_vector[1]
    ),
    labels = c(
        "T" = "Triticale"
    )
  ) +
  labs(
    subtitle = "Triticale biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

plot_legume_2023B_1 <- ggplot(data_2023B_biomass_legume_cld2_1, aes(x = Treatment, y = mean_biomass, fill = LegumeGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 3), na.rm = TRUE, size = 6) +
  scale_fill_manual(
    name = "Legume species",
    values = c(
      "F+" = colour_vector[2],
      "F" = colour_vector[2],
      "1T:3F" = colour_vector[2],
      "1T:1F" = colour_vector[2],
      "1T:1F-M" = colour_vector[2],
      "3T:1F" = colour_vector[2]
    ),
    labels = c(
        "F" = "Faba bean"
    )
  ) +
  labs(
    subtitle = "Faba bean biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

combined_plot_2023B_1 <- plot_cereal_2023B_1 / plot_legume_2023B_1 +
    plot_layout(ncol = 1
    ) +
    plot_annotation(tag_levels = "a")
combined_plot_2023B_1

## 2023B, second harvest, no weeds
data_2023B_biomass_2NW <- data2023_B2NW %>% 
  mutate(Treatment = gsub("TF", "1T:1F", Treatment))


data_2023B_biomass_groups_2NW <- data_2023B_biomass_2NW %>%
  mutate(
    CerealGroup = ifelse(grepl("T", Treatment), "T", NA_character_),
    LegumeGroup = ifelse(grepl("F", Treatment), "F", NA_character_)
  )

data_2023B_biomass_cereal_2NW <- data_2023B_biomass_groups_2NW %>%
  filter(!is.na(CerealGroup)) %>%
  group_by(CerealGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassC, na.rm = TRUE),
    sd_biomass = sd(BiomassC, na.rm = TRUE),
    n = sum(!is.na(BiomassC)),
    .groups = "drop"
  )

data_2023B_biomass_legume_2NW <- data_2023B_biomass_groups_2NW %>%
  filter(!is.na(LegumeGroup)) %>%
  group_by(LegumeGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassL, na.rm = TRUE),
    sd_biomass = sd(BiomassL, na.rm = TRUE),
    n = sum(!is.na(BiomassL)),
    .groups = "drop"
  )

cereal_levels <- c(
  "T", "1T:1F", "1T:1F-M", "T-375", "1T:1F-375", "1T:1F-M-375"
)

data_2023B_biomass_cereal2_2NW <- data_2023B_biomass_cereal_2NW %>%
  mutate(Treatment = factor(Treatment, levels = cereal_levels))

legume_levels <- c(
  "F", "1T:1F", "1T:1F-M", "F-375", "1T:1F-375", "1T:1F-M-375"
)

data_2023B_biomass_legume2_2NW <- data_2023B_biomass_legume_2NW %>%
  mutate(Treatment = factor(Treatment, levels = legume_levels))

cereal_cld_2023B_2NW <- lapply(unique(data_2023B_biomass_cereal2_2NW$CerealGroup), function(g) {
  treatments <- cereal_levels[grepl(g, cereal_levels)]
  df_sub <- data_2023B_biomass_2NW %>%
    filter(Treatment %in% treatments & !is.na(BiomassC))
    cld_df <- get_cld(df_sub, "BiomassC", "CerealGroup", treatments)
    cld_df$CerealGroup <- g
    return(cld_df)
}) %>% bind_rows()

legume_cld_2023B_2NW <- lapply(unique(data_2023B_biomass_legume2_2NW$LegumeGroup), function(g) {
  treatments <- legume_levels[grepl(g, legume_levels)]
  df_sub <- data_2023B_biomass_2NW %>%
    filter(Treatment %in% treatments & !is.na(BiomassL))
  if(nrow(df_sub) > 0) {
    cld_df <- get_cld2(df_sub, "BiomassL", "LegumeGroup", treatments)
    cld_df$LegumeGroup <- g
    return(cld_df)
  }
}) %>% bind_rows()

data_2023B_biomass_cereal_cld_2NW <- data_2023B_biomass_cereal2_2NW %>%
  left_join(cereal_cld_2023B_2NW, by = c("Treatment", "CerealGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2023B_biomass_legume_cld_2NW <- data_2023B_biomass_legume2_2NW %>%
  left_join(legume_cld_2023B_2NW, by = c("Treatment", "LegumeGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2023B_biomass_cereal_cld2_2NW <- data_2023B_biomass_cereal_cld_2NW %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "T", "1T:1F", "1T:1F-M", "T-375", "1T:1F-375", "1T:1F-M-375"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("T", "1T:1F", "1T:1F-M", "T-375", "1T:1F-375", "1T:1F-M-375"))) %>% 
  mutate(CerealGroup = factor(CerealGroup, levels = c("T")))

data_2023B_biomass_legume_cld2_2NW <- data_2023B_biomass_legume_cld_2NW %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "F", "F-375", "1T:1F", "1T:1F-M", "1T:1F-375", "1T:1F-M-375"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("F", "1T:1F", "1T:1F-M", "F-375", "1T:1F-375", "1T:1F-M-375"
                                                  ))) %>% 
  mutate(LegumeGroup = factor(LegumeGroup, levels = c("F")))

plot_cereal_2023B_2NW <- ggplot(data_2023B_biomass_cereal_cld2_2NW, aes(x = Treatment, y = mean_biomass, fill = CerealGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 100), na.rm = TRUE, size = 6) +
  scale_fill_manual(
    name = "Cereal species",
    values = c(
      "T+" = colour_vector[1],
      "T" = colour_vector[1],
      "3T:1F" = colour_vector[1],
      "1T:1F" = colour_vector[1],
      "1T:1F-M" = colour_vector[1],
      "1T:3F" = colour_vector[1]
    ),
    labels = c(
        "T" = "Triticale"
    )
  ) +
  labs(
    subtitle = "Triticale biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

plot_legume_2023B_2NW <- ggplot(data_2023B_biomass_legume_cld2_2NW, aes(x = Treatment, y = mean_biomass, fill = LegumeGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 60), na.rm = TRUE, size = 6) +
  scale_fill_manual(
    name = "Legume species",
    values = c(
      "F+" = colour_vector[2],
      "F" = colour_vector[2],
      "1T:3F" = colour_vector[2],
      "1T:1F" = colour_vector[2],
      "1T:1F-M" = colour_vector[2],
      "3T:1F" = colour_vector[2]
    ),
    labels = c(
        "F" = "Faba bean"
    )
  ) +
  labs(
    subtitle = "Faba bean biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

combined_plot_2023B_2NW <- plot_cereal_2023B_2NW / plot_legume_2023B_2NW +
    plot_layout(ncol = 1
    ) +
    plot_annotation(tag_levels = "a")
combined_plot_2023B_2NW

## 2023B, second harvest, with weeds
data_2023B_biomass_2WW <- data2023_B2WW %>% 
  mutate(Treatment = gsub("TF", "1T:1F", Treatment))


data_2023B_biomass_groups_2WW <- data_2023B_biomass_2WW %>%
  mutate(
    CerealGroup = ifelse(grepl("T", Treatment), "T", NA_character_),
    LegumeGroup = ifelse(grepl("F", Treatment), "F", NA_character_)
  )

data_2023B_biomass_cereal_2WW <- data_2023B_biomass_groups_2WW %>%
  filter(!is.na(CerealGroup)) %>%
  group_by(CerealGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassC, na.rm = TRUE),
    sd_biomass = sd(BiomassC, na.rm = TRUE),
    n = sum(!is.na(BiomassC)),
    .groups = "drop"
  )

data_2023B_biomass_legume_2WW <- data_2023B_biomass_groups_2WW %>%
  filter(!is.na(LegumeGroup)) %>%
  group_by(LegumeGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassL, na.rm = TRUE),
    sd_biomass = sd(BiomassL, na.rm = TRUE),
    n = sum(!is.na(BiomassL)),
    .groups = "drop"
  )

cereal_levels <- c(
  "T", "1T:1F", "1T:1F-M", "T-375", "1T:1F-375", "1T:1F-M-375"
)

data_2023B_biomass_cereal2_2WW <- data_2023B_biomass_cereal_2WW %>%
  mutate(Treatment = factor(Treatment, levels = cereal_levels))

legume_levels <- c(
  "F", "1T:1F", "1T:1F-M", "F-375", "1T:1F-375", "1T:1F-M-375"
)

data_2023B_biomass_legume2_2WW <- data_2023B_biomass_legume_2WW %>%
  mutate(Treatment = factor(Treatment, levels = legume_levels))

cereal_cld_2023B_2WW <- lapply(unique(data_2023B_biomass_cereal2_2WW$CerealGroup), function(g) {
  treatments <- cereal_levels[grepl(g, cereal_levels)]
  df_sub <- data_2023B_biomass_2WW %>%
    filter(Treatment %in% treatments & !is.na(BiomassC))
    cld_df <- get_cld(df_sub, "BiomassC", "CerealGroup", treatments)
    cld_df$CerealGroup <- g
    return(cld_df)
}) %>% bind_rows()

legume_cld_2023B_2WW <- lapply(unique(data_2023B_biomass_legume2_2WW$LegumeGroup), function(g) {
  treatments <- legume_levels[grepl(g, legume_levels)]
  df_sub <- data_2023B_biomass_2WW %>%
    filter(Treatment %in% treatments & !is.na(BiomassL))
  if(nrow(df_sub) > 0) {
    cld_df <- get_cld2(df_sub, "BiomassL", "LegumeGroup", treatments)
    cld_df$LegumeGroup <- g
    return(cld_df)
  }
}) %>% bind_rows()

data_2023B_biomass_cereal_cld_2WW <- data_2023B_biomass_cereal2_2WW %>%
  left_join(cereal_cld_2023B_2WW, by = c("Treatment", "CerealGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2023B_biomass_legume_cld_2WW <- data_2023B_biomass_legume2_2WW %>%
  left_join(legume_cld_2023B_2WW, by = c("Treatment", "LegumeGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2023B_biomass_cereal_cld2_2WW <- data_2023B_biomass_cereal_cld_2WW %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "T", "1T:1F", "1T:1F-M", "T-375", "1T:1F-375", "1T:1F-M-375"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("T", "1T:1F", "1T:1F-M", "T-375", "1T:1F-375", "1T:1F-M-375"))) %>% 
  mutate(CerealGroup = factor(CerealGroup, levels = c("T")))

data_2023B_biomass_legume_cld2_2WW <- data_2023B_biomass_legume_cld_2WW %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "F", "F-375", "1T:1F", "1T:1F-M", "1T:1F-375", "1T:1F-M-375"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("F", "1T:1F", "1T:1F-M", "F-375", "1T:1F-375", "1T:1F-M-375"
                                                  ))) %>% 
  mutate(LegumeGroup = factor(LegumeGroup, levels = c("F")))

plot_cereal_2023B_2WW <- ggplot(data_2023B_biomass_cereal_cld2_2WW, aes(x = Treatment, y = mean_biomass, fill = CerealGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 100), na.rm = TRUE, size = 6) +
  scale_fill_manual(
    name = "Cereal species",
    values = c(
      "T+" = colour_vector[1],
      "T" = colour_vector[1],
      "3T:1F" = colour_vector[1],
      "1T:1F" = colour_vector[1],
      "1T:1F-M" = colour_vector[1],
      "1T:3F" = colour_vector[1]
    ),
    labels = c(
        "T" = "Triticale"
    )
  ) +
  labs(
    subtitle = "Triticale biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

plot_legume_2023B_2WW <- ggplot(data_2023B_biomass_legume_cld2_2WW, aes(x = Treatment, y = mean_biomass, fill = LegumeGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 60), na.rm = TRUE, size = 6) +
  scale_fill_manual(
    name = "Legume species",
    values = c(
      "F+" = colour_vector[2],
      "F" = colour_vector[2],
      "1T:3F" = colour_vector[2],
      "1T:1F" = colour_vector[2],
      "1T:1F-M" = colour_vector[2],
      "3T:1F" = colour_vector[2]
    ),
    labels = c(
        "F" = "Faba bean"
    )
  ) +
  labs(
    subtitle = "Faba bean biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

combined_plot_2023B_2WW <- plot_cereal_2023B_2WW / plot_legume_2023B_2WW +
    plot_layout(ncol = 1
    ) +
    plot_annotation(tag_levels = "a")
combined_plot_2023B_2WW

## 2023B, final harvest, no weeds
data_2023B_biomass_3NW <- data2023_B3NW %>% 
  mutate(Treatment = gsub("TF", "1T:1F", Treatment))


data_2023B_biomass_groups_3NW <- data_2023B_biomass_3NW %>%
  mutate(
    CerealGroup = ifelse(grepl("T", Treatment), "T", NA_character_),
    LegumeGroup = ifelse(grepl("F", Treatment), "F", NA_character_)
  )

data_2023B_biomass_cereal_3NW <- data_2023B_biomass_groups_3NW %>%
  filter(!is.na(CerealGroup)) %>%
  group_by(CerealGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassC, na.rm = TRUE),
    sd_biomass = sd(BiomassC, na.rm = TRUE),
    n = sum(!is.na(BiomassC)),
    .groups = "drop"
  )

data_2023B_biomass_legume_3NW <- data_2023B_biomass_groups_3NW %>%
  filter(!is.na(LegumeGroup)) %>%
  group_by(LegumeGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassL, na.rm = TRUE),
    sd_biomass = sd(BiomassL, na.rm = TRUE),
    n = sum(!is.na(BiomassL)),
    .groups = "drop"
  )

cereal_levels <- c(
  "T", "1T:1F", "1T:1F-M", "T-375", "1T:1F-375", "1T:1F-M-375"
)

data_2023B_biomass_cereal2_3NW <- data_2023B_biomass_cereal_3NW %>%
  mutate(Treatment = factor(Treatment, levels = cereal_levels))

legume_levels <- c(
  "F", "1T:1F", "1T:1F-M", "F-375", "1T:1F-375", "1T:1F-M-375"
)

data_2023B_biomass_legume2_3NW <- data_2023B_biomass_legume_3NW %>%
  mutate(Treatment = factor(Treatment, levels = legume_levels))

cereal_cld_2023B_3NW <- lapply(unique(data_2023B_biomass_cereal2_3NW$CerealGroup), function(g) {
  treatments <- cereal_levels[grepl(g, cereal_levels)]
  df_sub <- data_2023B_biomass_3NW %>%
    filter(Treatment %in% treatments & !is.na(BiomassC))
    cld_df <- get_cld(df_sub, "BiomassC", "CerealGroup", treatments)
    cld_df$CerealGroup <- g
    return(cld_df)
}) %>% bind_rows()

legume_cld_2023B_3NW <- lapply(unique(data_2023B_biomass_legume2_3NW$LegumeGroup), function(g) {
  treatments <- legume_levels[grepl(g, legume_levels)]
  df_sub <- data_2023B_biomass_3NW %>%
    filter(Treatment %in% treatments & !is.na(BiomassL))
  if(nrow(df_sub) > 0) {
    cld_df <- get_cld2(df_sub, "BiomassL", "LegumeGroup", treatments)
    cld_df$LegumeGroup <- g
    return(cld_df)
  }
}) %>% bind_rows()

data_2023B_biomass_cereal_cld_3NW <- data_2023B_biomass_cereal2_3NW %>%
  left_join(cereal_cld_2023B_3NW, by = c("Treatment", "CerealGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2023B_biomass_legume_cld_3NW <- data_2023B_biomass_legume2_3NW %>%
  left_join(legume_cld_2023B_3NW, by = c("Treatment", "LegumeGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2023B_biomass_cereal_cld2_3NW <- data_2023B_biomass_cereal_cld_3NW %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "T", "1T:1F", "1T:1F-M", "T-375", "1T:1F-375", "1T:1F-M-375"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("T", "1T:1F", "1T:1F-M", "T-375", "1T:1F-375", "1T:1F-M-375"))) %>% 
  mutate(CerealGroup = factor(CerealGroup, levels = c("T")))

data_2023B_biomass_legume_cld2_3NW <- data_2023B_biomass_legume_cld_3NW %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "F", "F-375", "1T:1F", "1T:1F-M", "1T:1F-375", "1T:1F-M-375"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("F", "1T:1F", "1T:1F-M", "F-375", "1T:1F-375", "1T:1F-M-375"
                                                  ))) %>% 
  mutate(LegumeGroup = factor(LegumeGroup, levels = c("F")))

plot_cereal_2023B_3NW <- ggplot(data_2023B_biomass_cereal_cld2_3NW, aes(x = Treatment, y = mean_biomass, fill = CerealGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 100), na.rm = TRUE, size = 6) +
  scale_fill_manual(
    name = "Cereal species",
    values = c(
      "T+" = colour_vector[1],
      "T" = colour_vector[1],
      "3T:1F" = colour_vector[1],
      "1T:1F" = colour_vector[1],
      "1T:1F-M" = colour_vector[1],
      "1T:3F" = colour_vector[1]
    ),
    labels = c(
        "T" = "Triticale"
    )
  ) +
  labs(
    subtitle = "Triticale biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

plot_legume_2023B_3NW <- ggplot(data_2023B_biomass_legume_cld2_3NW, aes(x = Treatment, y = mean_biomass, fill = LegumeGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 60), na.rm = TRUE, size = 6) +
  scale_fill_manual(
    name = "Legume species",
    values = c(
      "F+" = colour_vector[2],
      "F" = colour_vector[2],
      "1T:3F" = colour_vector[2],
      "1T:1F" = colour_vector[2],
      "1T:1F-M" = colour_vector[2],
      "3T:1F" = colour_vector[2]
    ),
    labels = c(
        "F" = "Faba bean"
    )
  ) +
  labs(
    subtitle = "Faba bean biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

combined_plot_2023B_3NW <- plot_cereal_2023B_3NW / plot_legume_2023B_3NW +
    plot_layout(ncol = 1
    ) +
    plot_annotation(tag_levels = "a")
combined_plot_2023B_3NW

## 2023B, final harvest, with weeds
data_2023B_biomass_3WW <- data2023_B3WW %>% 
  mutate(Treatment = gsub("TF", "1T:1F", Treatment))


data_2023B_biomass_groups_3WW <- data_2023B_biomass_3WW %>%
  mutate(
    CerealGroup = ifelse(grepl("T", Treatment), "T", NA_character_),
    LegumeGroup = ifelse(grepl("F", Treatment), "F", NA_character_)
  )

data_2023B_biomass_cereal_3WW <- data_2023B_biomass_groups_3WW %>%
  filter(!is.na(CerealGroup)) %>%
  group_by(CerealGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassC, na.rm = TRUE),
    sd_biomass = sd(BiomassC, na.rm = TRUE),
    n = sum(!is.na(BiomassC)),
    .groups = "drop"
  )

data_2023B_biomass_legume_3WW <- data_2023B_biomass_groups_3WW %>%
  filter(!is.na(LegumeGroup)) %>%
  group_by(LegumeGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassL, na.rm = TRUE),
    sd_biomass = sd(BiomassL, na.rm = TRUE),
    n = sum(!is.na(BiomassL)),
    .groups = "drop"
  )

cereal_levels <- c(
  "T", "1T:1F", "1T:1F-M", "T-375", "1T:1F-375", "1T:1F-M-375"
)

data_2023B_biomass_cereal2_3WW <- data_2023B_biomass_cereal_3WW %>%
  mutate(Treatment = factor(Treatment, levels = cereal_levels))

legume_levels <- c(
  "F", "1T:1F", "1T:1F-M", "F-375", "1T:1F-375", "1T:1F-M-375"
)

data_2023B_biomass_legume2_3WW <- data_2023B_biomass_legume_3WW %>%
  mutate(Treatment = factor(Treatment, levels = legume_levels))

cereal_cld_2023B_3WW <- lapply(unique(data_2023B_biomass_cereal2_3WW$CerealGroup), function(g) {
  treatments <- cereal_levels[grepl(g, cereal_levels)]
  df_sub <- data_2023B_biomass_3WW %>%
    filter(Treatment %in% treatments & !is.na(BiomassC))
    cld_df <- get_cld(df_sub, "BiomassC", "CerealGroup", treatments)
    cld_df$CerealGroup <- g
    return(cld_df)
}) %>% bind_rows()

legume_cld_2023B_3WW <- lapply(unique(data_2023B_biomass_legume2_3WW$LegumeGroup), function(g) {
  treatments <- legume_levels[grepl(g, legume_levels)]
  df_sub <- data_2023B_biomass_3WW %>%
    filter(Treatment %in% treatments & !is.na(BiomassL))
  if(nrow(df_sub) > 0) {
    cld_df <- get_cld2(df_sub, "BiomassL", "LegumeGroup", treatments)
    cld_df$LegumeGroup <- g
    return(cld_df)
  }
}) %>% bind_rows()

data_2023B_biomass_cereal_cld_3WW <- data_2023B_biomass_cereal2_3WW %>%
  left_join(cereal_cld_2023B_3WW, by = c("Treatment", "CerealGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2023B_biomass_legume_cld_3WW <- data_2023B_biomass_legume2_3WW %>%
  left_join(legume_cld_2023B_3WW, by = c("Treatment", "LegumeGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2023B_biomass_cereal_cld2_3WW <- data_2023B_biomass_cereal_cld_3WW %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "T", "1T:1F", "1T:1F-M", "T-375", "1T:1F-375", "1T:1F-M-375"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("T", "1T:1F", "1T:1F-M", "T-375", "1T:1F-375", "1T:1F-M-375"))) %>% 
  mutate(CerealGroup = factor(CerealGroup, levels = c("T")))

data_2023B_biomass_legume_cld2_3WW <- data_2023B_biomass_legume_cld_3WW %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "F", "F-375", "1T:1F", "1T:1F-M", "1T:1F-375", "1T:1F-M-375"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("F", "1T:1F", "1T:1F-M", "F-375", "1T:1F-375", "1T:1F-M-375"
                                                  ))) %>% 
  mutate(LegumeGroup = factor(LegumeGroup, levels = c("F")))

plot_cereal_2023B_3WW <- ggplot(data_2023B_biomass_cereal_cld2_3WW, aes(x = Treatment, y = mean_biomass, fill = CerealGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 100), na.rm = TRUE, size = 6) +
  scale_fill_manual(
    name = "Cereal species",
    values = c(
      "T+" = colour_vector[1],
      "T" = colour_vector[1],
      "3T:1F" = colour_vector[1],
      "1T:1F" = colour_vector[1],
      "1T:1F-M" = colour_vector[1],
      "1T:3F" = colour_vector[1]
    ),
    labels = c(
        "T" = "Triticale"
    )
  ) +
  labs(
    subtitle = "Triticale biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

plot_legume_2023B_3WW <- ggplot(data_2023B_biomass_legume_cld2_3WW, aes(x = Treatment, y = mean_biomass, fill = LegumeGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 60), na.rm = TRUE, size = 6) +
  scale_fill_manual(
    name = "Legume species",
    values = c(
      "F+" = colour_vector[2],
      "F" = colour_vector[2],
      "1T:3F" = colour_vector[2],
      "1T:1F" = colour_vector[2],
      "1T:1F-M" = colour_vector[2],
      "3T:1F" = colour_vector[2]
    ),
    labels = c(
        "F" = "Faba bean"
    )
  ) +
  labs(
    subtitle = "Faba bean biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

combined_plot_2023B_3WW <- plot_cereal_2023B_3WW / plot_legume_2023B_3WW +
    plot_layout(ncol = 1
    ) +
    plot_annotation(tag_levels = "a")
combined_plot_2023B_3WW

## 2024, first harvest
data_2024_biomass_1 <- data2024_1 %>% 
  mutate(Treatment = gsub("TF", "1T:1F", Treatment))


data_2024_biomass_groups_1 <- data_2024_biomass_1 %>%
  mutate(
    CerealGroup = ifelse(grepl("T", Treatment), "T", NA_character_),
    LegumeGroup = ifelse(grepl("F", Treatment), "F", NA_character_)
  )

data_2024_biomass_cereal_1 <- data_2024_biomass_groups_1 %>%
  filter(!is.na(CerealGroup)) %>%
  group_by(CerealGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassC, na.rm = TRUE),
    sd_biomass = sd(BiomassC, na.rm = TRUE),
    n = sum(!is.na(BiomassC)),
    .groups = "drop"
  )

data_2024_biomass_legume_1 <- data_2024_biomass_groups_1 %>%
  filter(!is.na(LegumeGroup)) %>%
  group_by(LegumeGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassL, na.rm = TRUE),
    sd_biomass = sd(BiomassL, na.rm = TRUE),
    n = sum(!is.na(BiomassL)),
    .groups = "drop"
  )

cereal_levels <- c(
  "T", "1T:1F"
)

data_2024_biomass_cereal2_1 <- data_2024_biomass_cereal_1 %>%
  mutate(Treatment = factor(Treatment, levels = cereal_levels))

legume_levels <- c(
  "F", "1T:1F"
)

data_2024_biomass_legume2_1 <- data_2024_biomass_legume_1 %>%
  mutate(Treatment = factor(Treatment, levels = legume_levels))

cereal_cld_2024_1 <- lapply(unique(data_2024_biomass_cereal2_1$CerealGroup), function(g) {
  treatments <- cereal_levels[grepl(g, cereal_levels)]
  df_sub <- data_2024_biomass_1 %>%
    filter(Treatment %in% treatments & !is.na(BiomassC))
    cld_df <- get_cld(df_sub, "BiomassC", "CerealGroup", treatments)
    cld_df$CerealGroup <- g
    return(cld_df)
}) %>% bind_rows()

legume_cld_2024_1 <- lapply(unique(data_2024_biomass_legume2_1$LegumeGroup), function(g) {
  treatments <- legume_levels[grepl(g, legume_levels)]
  df_sub <- data_2024_biomass_1 %>%
    filter(Treatment %in% treatments & !is.na(BiomassL))
  if(nrow(df_sub) > 0) {
    cld_df <- get_cld2(df_sub, "BiomassL", "LegumeGroup", treatments)
    cld_df$LegumeGroup <- g
    return(cld_df)
  }
}) %>% bind_rows()

data_2024_biomass_cereal_cld_1 <- data_2024_biomass_cereal2_1 %>%
  left_join(cereal_cld_2024_1, by = c("Treatment", "CerealGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2024_biomass_legume_cld_1 <- data_2024_biomass_legume2_1 %>%
  left_join(legume_cld_2024_1, by = c("Treatment", "LegumeGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2024_biomass_cereal_cld2_1 <- data_2024_biomass_cereal_cld_1 %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "T", "1T:1F"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("T", "1T:1F"))) %>% 
  mutate(CerealGroup = factor(CerealGroup, levels = c("T")))

data_2024_biomass_legume_cld2_1 <- data_2024_biomass_legume_cld_1 %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "F", "1T:1F"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("F", "1T:1F"
                                                  ))) %>% 
  mutate(LegumeGroup = factor(LegumeGroup, levels = c("F")))

plot_cereal_2024_1 <- ggplot(data_2024_biomass_cereal_cld2_1, aes(x = Treatment, y = mean_biomass, fill = CerealGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 6), na.rm = TRUE, size = 6) +
  scale_fill_manual(
    name = "Cereal species",
    values = c(
      "T" = colour_vector[1],
      "1T:1F" = colour_vector[1]
    ),
    labels = c(
        "T" = "Triticale"
    )
  ) +
  labs(
    subtitle = "Triticale biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

plot_legume_2024_1 <- ggplot(data_2024_biomass_legume_cld2_1, aes(x = Treatment, y = mean_biomass, fill = LegumeGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 3), na.rm = TRUE, size = 6) +
  scale_fill_manual(
    name = "Legume species",
    values = c(
      "F" = colour_vector[2],
      "1T:1F" = colour_vector[2]
    ),
    labels = c(
        "F" = "Faba bean"
    )
  ) +
  labs(
    subtitle = "Faba bean biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

combined_plot_2024_1 <- plot_cereal_2024_1 / plot_legume_2024_1 +
    plot_layout(ncol = 1
    ) +
    plot_annotation(tag_levels = "a")
combined_plot_2024_1

## 2024, second harvest, no weeds
data_2024_biomass_2NW <- data2024_2NW %>% 
  mutate(Treatment = gsub("TF", "1T:1F", Treatment))


data_2024_biomass_groups_2NW <- data_2024_biomass_2NW %>%
  mutate(
    CerealGroup = ifelse(grepl("T", Treatment), "T", NA_character_),
    LegumeGroup = ifelse(grepl("F", Treatment), "F", NA_character_)
  )

data_2024_biomass_cereal_2NW <- data_2024_biomass_groups_2NW %>%
  filter(!is.na(CerealGroup)) %>%
  group_by(CerealGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassC, na.rm = TRUE),
    sd_biomass = sd(BiomassC, na.rm = TRUE),
    n = sum(!is.na(BiomassC)),
    .groups = "drop"
  )

data_2024_biomass_legume_2NW <- data_2024_biomass_groups_2NW %>%
  filter(!is.na(LegumeGroup)) %>%
  group_by(LegumeGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassL, na.rm = TRUE),
    sd_biomass = sd(BiomassL, na.rm = TRUE),
    n = sum(!is.na(BiomassL)),
    .groups = "drop"
  )

cereal_levels <- c(
  "T", "1T:1F"
)

data_2024_biomass_cereal2_2NW <- data_2024_biomass_cereal_2NW %>%
  mutate(Treatment = factor(Treatment, levels = cereal_levels))

legume_levels <- c(
  "F", "1T:1F"
)

data_2024_biomass_legume2_2NW <- data_2024_biomass_legume_2NW %>%
  mutate(Treatment = factor(Treatment, levels = legume_levels))

cereal_cld_2024_2NW <- lapply(unique(data_2024_biomass_cereal2_2NW$CerealGroup), function(g) {
  treatments <- cereal_levels[grepl(g, cereal_levels)]
  df_sub <- data_2024_biomass_2NW %>%
    filter(Treatment %in% treatments & !is.na(BiomassC))
    cld_df <- get_cld(df_sub, "BiomassC", "CerealGroup", treatments)
    cld_df$CerealGroup <- g
    return(cld_df)
}) %>% bind_rows()

legume_cld_2024_2NW <- lapply(unique(data_2024_biomass_legume2_2NW$LegumeGroup), function(g) {
  treatments <- legume_levels[grepl(g, legume_levels)]
  df_sub <- data_2024_biomass_2NW %>%
    filter(Treatment %in% treatments & !is.na(BiomassL))
  if(nrow(df_sub) > 0) {
    cld_df <- get_cld2(df_sub, "BiomassL", "LegumeGroup", treatments)
    cld_df$LegumeGroup <- g
    return(cld_df)
  }
}) %>% bind_rows()

data_2024_biomass_cereal_cld_2NW <- data_2024_biomass_cereal2_2NW %>%
  left_join(cereal_cld_2024_2NW, by = c("Treatment", "CerealGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2024_biomass_legume_cld_2NW <- data_2024_biomass_legume2_2NW %>%
  left_join(legume_cld_2024_2NW, by = c("Treatment", "LegumeGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2024_biomass_cereal_cld2_2NW <- data_2024_biomass_cereal_cld_2NW %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "T", "1T:1F"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("T", "1T:1F"))) %>% 
  mutate(CerealGroup = factor(CerealGroup, levels = c("T")))

data_2024_biomass_legume_cld2_2NW <- data_2024_biomass_legume_cld_2NW %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "F", "1T:1F"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("F", "1T:1F"
                                                  ))) %>% 
  mutate(LegumeGroup = factor(LegumeGroup, levels = c("F")))

plot_cereal_2024_2NW <- ggplot(data_2024_biomass_cereal_cld2_2NW, aes(x = Treatment, y = mean_biomass, fill = CerealGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 100), na.rm = TRUE, size = 6) +
  scale_fill_manual(
    name = "Cereal species",
    values = c(
      "T" = colour_vector[1],
      "1T:1F" = colour_vector[1]
    ),
    labels = c(
        "T" = "Triticale"
    )
  ) +
  labs(
    subtitle = "Triticale biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

plot_legume_2024_2NW <- ggplot(data_2024_biomass_legume_cld2_2NW, aes(x = Treatment, y = mean_biomass, fill = LegumeGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 60), na.rm = TRUE, size = 6) +
  scale_fill_manual(
    name = "Legume species",
    values = c(
      "F" = colour_vector[2],
      "1T:1F" = colour_vector[2]
    ),
    labels = c(
        "F" = "Faba bean"
    )
  ) +
  labs(
    subtitle = "Faba bean biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

combined_plot_2024_2NW <- plot_cereal_2024_2NW / plot_legume_2024_2NW +
    plot_layout(ncol = 1
    ) +
    plot_annotation(tag_levels = "a")
combined_plot_2024_2NW

## 2024, second harvest, with weeds
data_2024_biomass_2WW <- data2024_2WW %>% 
  mutate(Treatment = gsub("TF", "1T:1F", Treatment))


data_2024_biomass_groups_2WW <- data_2024_biomass_2WW %>%
  mutate(
    CerealGroup = ifelse(grepl("T", Treatment), "T", NA_character_),
    LegumeGroup = ifelse(grepl("F", Treatment), "F", NA_character_)
  )

data_2024_biomass_cereal_2WW <- data_2024_biomass_groups_2WW %>%
  filter(!is.na(CerealGroup)) %>%
  group_by(CerealGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassC, na.rm = TRUE),
    sd_biomass = sd(BiomassC, na.rm = TRUE),
    n = sum(!is.na(BiomassC)),
    .groups = "drop"
  )

data_2024_biomass_legume_2WW <- data_2024_biomass_groups_2WW %>%
  filter(!is.na(LegumeGroup)) %>%
  group_by(LegumeGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassL, na.rm = TRUE),
    sd_biomass = sd(BiomassL, na.rm = TRUE),
    n = sum(!is.na(BiomassL)),
    .groups = "drop"
  )

cereal_levels <- c(
  "T", "1T:1F"
)

data_2024_biomass_cereal2_2WW <- data_2024_biomass_cereal_2WW %>%
  mutate(Treatment = factor(Treatment, levels = cereal_levels))

legume_levels <- c(
  "F", "1T:1F"
)

data_2024_biomass_legume2_2WW <- data_2024_biomass_legume_2WW %>%
  mutate(Treatment = factor(Treatment, levels = legume_levels))

cereal_cld_2024_2WW <- lapply(unique(data_2024_biomass_cereal2_2WW$CerealGroup), function(g) {
  treatments <- cereal_levels[grepl(g, cereal_levels)]
  df_sub <- data_2024_biomass_2WW %>%
    filter(Treatment %in% treatments & !is.na(BiomassC))
    cld_df <- get_cld(df_sub, "BiomassC", "CerealGroup", treatments)
    cld_df$CerealGroup <- g
    return(cld_df)
}) %>% bind_rows()

legume_cld_2024_2WW <- lapply(unique(data_2024_biomass_legume2_2WW$LegumeGroup), function(g) {
  treatments <- legume_levels[grepl(g, legume_levels)]
  df_sub <- data_2024_biomass_2WW %>%
    filter(Treatment %in% treatments & !is.na(BiomassL))
  if(nrow(df_sub) > 0) {
    cld_df <- get_cld2(df_sub, "BiomassL", "LegumeGroup", treatments)
    cld_df$LegumeGroup <- g
    return(cld_df)
  }
}) %>% bind_rows()

data_2024_biomass_cereal_cld_2WW <- data_2024_biomass_cereal2_2WW %>%
  left_join(cereal_cld_2024_2WW, by = c("Treatment", "CerealGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2024_biomass_legume_cld_2WW <- data_2024_biomass_legume2_2WW %>%
  left_join(legume_cld_2024_2WW, by = c("Treatment", "LegumeGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2024_biomass_cereal_cld2_2WW <- data_2024_biomass_cereal_cld_2WW %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "T", "1T:1F"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("T", "1T:1F"))) %>% 
  mutate(CerealGroup = factor(CerealGroup, levels = c("T")))

data_2024_biomass_legume_cld2_2WW <- data_2024_biomass_legume_cld_2WW %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "F", "1T:1F"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("F", "1T:1F"
                                                  ))) %>% 
  mutate(LegumeGroup = factor(LegumeGroup, levels = c("F")))

plot_cereal_2024_2WW <- ggplot(data_2024_biomass_cereal_cld2_2WW, aes(x = Treatment, y = mean_biomass, fill = CerealGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 100), na.rm = TRUE, size = 6) +
  scale_fill_manual(
    name = "Cereal species",
    values = c(
      "T" = colour_vector[1],
      "1T:1F" = colour_vector[1]
    ),
    labels = c(
        "T" = "Triticale"
    )
  ) +
  labs(
    subtitle = "Triticale biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

plot_legume_2024_2WW <- ggplot(data_2024_biomass_legume_cld2_2WW, aes(x = Treatment, y = mean_biomass, fill = LegumeGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 60), na.rm = TRUE, size = 6) +
  scale_fill_manual(
    name = "Legume species",
    values = c(
      "F" = colour_vector[2],
      "1T:1F" = colour_vector[2]
    ),
    labels = c(
        "F" = "Faba bean"
    )
  ) +
  labs(
    subtitle = "Faba bean biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

combined_plot_2024_2WW <- plot_cereal_2024_2WW / plot_legume_2024_2WW +
    plot_layout(ncol = 1
    ) +
    plot_annotation(tag_levels = "a")
combined_plot_2024_2WW

## 2024, final harvest, no weeds
data_2024_biomass_3NW <- data2024_3NW %>% 
  mutate(Treatment = gsub("TF", "1T:1F", Treatment))


data_2024_biomass_groups_3NW <- data_2024_biomass_3NW %>%
  mutate(
    CerealGroup = ifelse(grepl("T", Treatment), "T", NA_character_),
    LegumeGroup = ifelse(grepl("F", Treatment), "F", NA_character_)
  )

data_2024_biomass_cereal_3NW <- data_2024_biomass_groups_3NW %>%
  filter(!is.na(CerealGroup)) %>%
  group_by(CerealGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassC, na.rm = TRUE),
    sd_biomass = sd(BiomassC, na.rm = TRUE),
    n = sum(!is.na(BiomassC)),
    .groups = "drop"
  )

data_2024_biomass_legume_3NW <- data_2024_biomass_groups_3NW %>%
  filter(!is.na(LegumeGroup)) %>%
  group_by(LegumeGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassL, na.rm = TRUE),
    sd_biomass = sd(BiomassL, na.rm = TRUE),
    n = sum(!is.na(BiomassL)),
    .groups = "drop"
  )

cereal_levels <- c(
  "T", "1T:1F"
)

data_2024_biomass_cereal2_3NW <- data_2024_biomass_cereal_3NW %>%
  mutate(Treatment = factor(Treatment, levels = cereal_levels))

legume_levels <- c(
  "F", "1T:1F"
)

data_2024_biomass_legume2_3NW <- data_2024_biomass_legume_3NW %>%
  mutate(Treatment = factor(Treatment, levels = legume_levels))

cereal_cld_2024_3NW <- lapply(unique(data_2024_biomass_cereal2_3NW$CerealGroup), function(g) {
  treatments <- cereal_levels[grepl(g, cereal_levels)]
  df_sub <- data_2024_biomass_3NW %>%
    filter(Treatment %in% treatments & !is.na(BiomassC))
    cld_df <- get_cld(df_sub, "BiomassC", "CerealGroup", treatments)
    cld_df$CerealGroup <- g
    return(cld_df)
}) %>% bind_rows()

legume_cld_2024_3NW <- lapply(unique(data_2024_biomass_legume2_3NW$LegumeGroup), function(g) {
  treatments <- legume_levels[grepl(g, legume_levels)]
  df_sub <- data_2024_biomass_3NW %>%
    filter(Treatment %in% treatments & !is.na(BiomassL))
  if(nrow(df_sub) > 0) {
    cld_df <- get_cld2(df_sub, "BiomassL", "LegumeGroup", treatments)
    cld_df$LegumeGroup <- g
    return(cld_df)
  }
}) %>% bind_rows()

data_2024_biomass_cereal_cld_3NW <- data_2024_biomass_cereal2_3NW %>%
  left_join(cereal_cld_2024_3NW, by = c("Treatment", "CerealGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2024_biomass_legume_cld_3NW <- data_2024_biomass_legume2_3NW %>%
  left_join(legume_cld_2024_3NW, by = c("Treatment", "LegumeGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2024_biomass_cereal_cld2_3NW <- data_2024_biomass_cereal_cld_3NW %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "T", "1T:1F"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("T", "1T:1F"))) %>% 
  mutate(CerealGroup = factor(CerealGroup, levels = c("T")))

data_2024_biomass_legume_cld2_3NW <- data_2024_biomass_legume_cld_3NW %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "F", "1T:1F"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("F", "1T:1F"
                                                  ))) %>% 
  mutate(LegumeGroup = factor(LegumeGroup, levels = c("F")))

plot_cereal_2024_3NW <- ggplot(data_2024_biomass_cereal_cld2_3NW, aes(x = Treatment, y = mean_biomass, fill = CerealGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 100), na.rm = TRUE, size = 6) +
  scale_fill_manual(
    name = "Cereal species",
    values = c(
      "T" = colour_vector[1],
      "1T:1F" = colour_vector[1]
    ),
    labels = c(
        "T" = "Triticale"
    )
  ) +
  labs(
    subtitle = "Triticale biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

plot_legume_2024_3NW <- ggplot(data_2024_biomass_legume_cld2_3NW, aes(x = Treatment, y = mean_biomass, fill = LegumeGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 60), na.rm = TRUE, size = 6) +
  scale_fill_manual(
    name = "Legume species",
    values = c(
      "F" = colour_vector[2],
      "1T:1F" = colour_vector[2]
    ),
    labels = c(
        "F" = "Faba bean"
    )
  ) +
  labs(
    subtitle = "Faba bean biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

combined_plot_2024_3NW <- plot_cereal_2024_3NW / plot_legume_2024_3NW +
    plot_layout(ncol = 1
    ) +
    plot_annotation(tag_levels = "a")
combined_plot_2024_3NW

## 2024, final harvest, with weeds
data_2024_biomass_3WW <- data2024_3WW %>% 
  mutate(Treatment = gsub("TF", "1T:1F", Treatment))


data_2024_biomass_groups_3WW <- data_2024_biomass_3WW %>%
  mutate(
    CerealGroup = ifelse(grepl("T", Treatment), "T", NA_character_),
    LegumeGroup = ifelse(grepl("F", Treatment), "F", NA_character_)
  )

data_2024_biomass_cereal_3WW <- data_2024_biomass_groups_3WW %>%
  filter(!is.na(CerealGroup)) %>%
  group_by(CerealGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassC, na.rm = TRUE),
    sd_biomass = sd(BiomassC, na.rm = TRUE),
    n = sum(!is.na(BiomassC)),
    .groups = "drop"
  )

data_2024_biomass_legume_3WW <- data_2024_biomass_groups_3WW %>%
  filter(!is.na(LegumeGroup)) %>%
  group_by(LegumeGroup, Treatment) %>%
  summarise(
    mean_biomass = mean(BiomassL, na.rm = TRUE),
    sd_biomass = sd(BiomassL, na.rm = TRUE),
    n = sum(!is.na(BiomassL)),
    .groups = "drop"
  )

cereal_levels <- c(
  "T", "1T:1F"
)

data_2024_biomass_cereal2_3WW <- data_2024_biomass_cereal_3WW %>%
  mutate(Treatment = factor(Treatment, levels = cereal_levels))

legume_levels <- c(
  "F", "1T:1F"
)

data_2024_biomass_legume2_3WW <- data_2024_biomass_legume_3WW %>%
  mutate(Treatment = factor(Treatment, levels = legume_levels))

cereal_cld_2024_3WW <- lapply(unique(data_2024_biomass_cereal2_3WW$CerealGroup), function(g) {
  treatments <- cereal_levels[grepl(g, cereal_levels)]
  df_sub <- data_2024_biomass_3WW %>%
    filter(Treatment %in% treatments & !is.na(BiomassC))
    cld_df <- get_cld(df_sub, "BiomassC", "CerealGroup", treatments)
    cld_df$CerealGroup <- g
    return(cld_df)
}) %>% bind_rows()

legume_cld_2024_3WW <- lapply(unique(data_2024_biomass_legume2_3WW$LegumeGroup), function(g) {
  treatments <- legume_levels[grepl(g, legume_levels)]
  df_sub <- data_2024_biomass_3WW %>%
    filter(Treatment %in% treatments & !is.na(BiomassL))
  if(nrow(df_sub) > 0) {
    cld_df <- get_cld2(df_sub, "BiomassL", "LegumeGroup", treatments)
    cld_df$LegumeGroup <- g
    return(cld_df)
  }
}) %>% bind_rows()

data_2024_biomass_cereal_cld_3WW <- data_2024_biomass_cereal2_3WW %>%
  left_join(cereal_cld_2024_3WW, by = c("Treatment", "CerealGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2024_biomass_legume_cld_3WW <- data_2024_biomass_legume2_3WW %>%
  left_join(legume_cld_2024_3WW, by = c("Treatment", "LegumeGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2024_biomass_cereal_cld2_3WW <- data_2024_biomass_cereal_cld_3WW %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "T", "1T:1F"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("T", "1T:1F"))) %>% 
  mutate(CerealGroup = factor(CerealGroup, levels = c("T")))

data_2024_biomass_legume_cld2_3WW <- data_2024_biomass_legume_cld_3WW %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "F", "1T:1F"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("F", "1T:1F"
                                                  ))) %>% 
  mutate(LegumeGroup = factor(LegumeGroup, levels = c("F")))

plot_cereal_2024_3WW <- ggplot(data_2024_biomass_cereal_cld2_3WW, aes(x = Treatment, y = mean_biomass, fill = CerealGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 100), na.rm = TRUE, size = 6) +
  scale_fill_manual(
    name = "Cereal species",
    values = c(
      "T" = colour_vector[1],
      "1T:1F" = colour_vector[1]
    ),
    labels = c(
        "T" = "Triticale"
    )
  ) +
  labs(
    subtitle = "Triticale biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

plot_legume_2024_3WW <- ggplot(data_2024_biomass_legume_cld2_3WW, aes(x = Treatment, y = mean_biomass, fill = LegumeGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_biomass - sd_biomass, ymax = mean_biomass + sd_biomass), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_biomass + sd_biomass + 60), na.rm = TRUE, size = 6) +
  scale_fill_manual(
    name = "Legume species",
    values = c(
      "F" = colour_vector[2],
      "1T:1F" = colour_vector[2]
    ),
    labels = c(
        "F" = "Faba bean"
    )
  ) +
  labs(
    subtitle = "Faba bean biomass",
    x = "Treatment",
    y = "Biomass (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.subtitle = element_text(hjust = 0.5))

combined_plot_2024_3WW <- plot_cereal_2024_3WW / plot_legume_2024_3WW +
    plot_layout(ncol = 1
    ) +
    plot_annotation(tag_levels = "a")
combined_plot_2024_3WW

dir.create("Figures/biomass")

#########################
## Save plots to files ##
#########################

png("Figures/biomass/plot_biomass_2022_1.png", units = "px", width = 4850, height = 4350, res = 300)
combined_plot_2022_1
dev.off()

png("Figures/biomass/plot_biomass_2022_2_noweeds.png", units = "px", width = 4850, height = 4350, res = 300)
combined_plot_2022_2NW
dev.off()

png("Figures/biomass/plot_biomass_2022_2_weeds.png", units = "px", width = 4850, height = 4350, res = 300)
combined_plot_2022_2WW
dev.off()

png("Figures/biomass/plot_biomass_2022_3_noweeds.png", units = "px", width = 4850, height = 4350, res = 300)
combined_plot_2022_3
dev.off()


png("Figures/biomass/plot_biomass_2023A_1.png", units = "px", width = 3000, height = 4350, res = 300)
combined_plot_2023A_1
dev.off()

png("Figures/biomass/plot_biomass_2023A_2_noweeds.png", units = "px", width = 3000, height = 4350, res = 300)
combined_plot_2023A_2NW
dev.off()

png("Figures/biomass/plot_biomass_2023A_2_weeds.png", units = "px", width = 3000, height = 4350, res = 300)
combined_plot_2023A_2WW
dev.off()

png("Figures/biomass/plot_biomass_2023A_3_noweeds.png", units = "px", width = 3000, height = 4350, res = 300)
combined_plot_2023A_3NW
dev.off()

png("Figures/biomass/plot_biomass_2023A_3_weeds.png", units = "px", width = 3000, height = 4350, res = 300)
combined_plot_2023A_3WW
dev.off()


png("Figures/biomass/plot_biomass_2023B_1.png", units = "px", width = 3000, height = 4350, res = 300)
combined_plot_2023B_1
dev.off()

png("Figures/biomass/plot_biomass_2023B_2_noweeds.png", units = "px", width = 3000, height = 4350, res = 300)
combined_plot_2023B_2NW
dev.off()

png("Figures/biomass/plot_biomass_2023B_2_weeds.png", units = "px", width = 3000, height = 4350, res = 300)
combined_plot_2023B_2WW
dev.off()

png("Figures/biomass/plot_biomass_2023B_3_noweeds.png", units = "px", width = 3000, height = 4350, res = 300)
combined_plot_2023B_3NW
dev.off()

png("Figures/biomass/plot_biomass_2023B_3_weeds.png", units = "px", width = 3000, height = 4350, res = 300)
combined_plot_2023B_3WW
dev.off()


png("Figures/biomass/plot_biomass_2024_1.png", units = "px", width = 2500, height = 4000, res = 300)
combined_plot_2024_1
dev.off()

png("Figures/biomass/plot_biomass_2024_2_noweeds.png", units = "px", width = 2500, height = 4000, res = 300)
combined_plot_2024_2NW
dev.off()

png("Figures/biomass/plot_biomass_2024_2_weeds.png", units = "px", width = 2500, height = 4000, res = 300)
combined_plot_2024_2WW
dev.off()

png("Figures/biomass/plot_biomass_2024_3_noweeds.png", units = "px", width = 2500, height = 4000, res = 300)
combined_plot_2024_3NW
dev.off()

png("Figures/biomass/plot_biomass_2024_3_weeds.png", units = "px", width = 2500, height = 4000, res = 300)
combined_plot_2024_3WW
dev.off()
