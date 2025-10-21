## Author: David Kottelenberg
## Institute: Wageningen University & Research
## Last edited: 21-10-2025
## Summary: Analyses of 2022-2024 field experiments, yield

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

# Set directory
# setwd("")

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

## Load data ##

# Read data
data_2022_yield <- read_xlsx("WCF_2022_data.xlsx", sheet = "FinalHarvest", range = "A1:P77", col_names = TRUE) %>% 
  dplyr::select(Plot, Treatment, "Block", "BiomassSeedsCereal", "1000SeedWeightCereal", "NSpikesCereal", "NSeedsCereal", "BiomassSeedsLegume", "1000SeedWeightLegume", "NPodsLegume", "NSeedsLegume") %>% 
  rename("YieldC" = "BiomassSeedsCereal",
         "YieldL" = "BiomassSeedsLegume",
         "TSWC" = "1000SeedWeightCereal",
         "TSWL" = "1000SeedWeightLegume",
         "NFC" = "NSpikesCereal",
         "NFL" = "NPodsLegume",
         "NSC" = "NSeedsCereal",
         "NSL" = "NSeedsLegume") %>% 
  arrange(Plot) %>% 
  filter(!grepl("Lupine", Treatment)) %>% # Lupine was abandoned in the experiment due to a reaction to the herbicide, and no space left in the no-herbicide parts of the plots
  dplyr::select(Plot, Treatment, Block, YieldC, YieldL) %>% 
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
  ))

#####################
##                 ##
## Experiment 2023 ##
##                 ##
#####################

# Read data
data_2023A_yield <- read_xlsx(
  "WCF_2023_data.xlsx",
  sheet = "FinalHarvestA", 
  range = "A1:AJ81", 
  col_names = TRUE) %>%
  dplyr::select(Plot, Treatment, Block, Weeds,
                "BiomassSeedsTriticaleA",
                "BiomassSeedsTriticaleM",
                "BiomassSeedsTriticale",
                "1000SeedWeightTriticaleA",
                "1000SeedWeightTriticaleM",
                "1000SeedWeightTriticale",
                "NSpikesTriticaleA",
                "NSpikesTriticaleM",
                "NSpikesTriticale",
                "NSeedsTriticaleA",
                "NSeedsTriticaleM",
                "NSeedsTriticale",
                "BiomassSeedsFabaA",
                "BiomassSeedsFabaM",
                "BiomassSeedsFaba",
                "1000SeedWeightFabaA",
                "1000SeedWeightFabaM",
                "1000SeedWeightFaba",
                "NPodsFabaA",
                "NPodsFabaM",
                "NPodsFaba",
                "NSeedsFabaA",
                "NSeedsFabaM",
                "NSeedsFaba") %>%
  rename("YieldTA" = "BiomassSeedsTriticaleA",
         "YieldTM" = "BiomassSeedsTriticaleM",
         "YieldT" = "BiomassSeedsTriticale",
         "YieldFA" = "BiomassSeedsFabaA",
         "YieldFM" = "BiomassSeedsFabaM",
         "YieldF" = "BiomassSeedsFaba",
         "NFTA" = "NSpikesTriticaleA",
         "NFTM" = "NSpikesTriticaleM",
         "NFT" = "NSpikesTriticale",
         "NFFA" = "NPodsFabaA",
         "NFFM" = "NPodsFabaM",
         "NFF" = "NPodsFaba",
         "TSWTA" = "1000SeedWeightTriticaleA",
         "TSWTM" = "1000SeedWeightTriticaleM",
         "TSWT" = "1000SeedWeightTriticale",
         "TSWFA" = "1000SeedWeightFabaA",
         "TSWFM" = "1000SeedWeightFabaM",
         "TSWF" = "1000SeedWeightFaba",
         "NSTA" = "NSeedsTriticaleA",
         "NSTM" = "NSeedsTriticaleM",
         "NST" = "NSeedsTriticale",
         "NSFA" = "NSeedsFabaA",
         "NSFM" = "NSeedsFabaM",
         "NSF" = "NSeedsFaba") %>%
  mutate(YieldTA = as.numeric(YieldTA),
         YieldTM = as.numeric(YieldTM),
         YieldT = as.numeric(YieldT),
         NFTA = as.numeric(NFTA),
         NFTM = as.numeric(NFTM),
         NFT = as.numeric(NFT),
         TSWTA = as.numeric(TSWTA),
         TSWTM = as.numeric(TSWTM),
         TSWT = as.numeric(TSWT),
         NSTA = as.numeric(NSTA),
         NSTM = as.numeric(NSTM),
         NST = as.numeric(NST),
         YieldFA = as.numeric(YieldFA),
         YieldFM = as.numeric(YieldFM),
         YieldF = as.numeric(YieldF),
         NFFA = as.numeric(NFFA),
         NFFM = as.numeric(NFFM),
         NFF = as.numeric(NFF),
         TSWFA = as.numeric(TSWFA),
         TSWFM = as.numeric(TSWFM),
         TSWF = as.numeric(TSWF),
         NSFA = as.numeric(NSFA),
         NSFM = as.numeric(NSFM),
         NSF = as.numeric(NSF)) %>% 
  mutate(
    YieldT = ifelse(
      is.na(YieldT) & grepl("T", Treatment), 
      YieldTA + YieldTM, 
      YieldT
    ),
    YieldF = ifelse(
      is.na(YieldF) & grepl("T", Treatment), 
      YieldFA + YieldFM, 
      YieldF
    ),
    NFT = ifelse(
      is.na(NFT) & grepl("T", Treatment), 
      NFTA + NFTM, 
      NFT
    ),
    NFF = ifelse(
      is.na(NFF) & grepl("T", Treatment), 
      NFFA + NFFM, 
      NFF
    ),
    TSWT = ifelse(
      is.na(TSWT) & grepl("T", Treatment),
      TSWTA + TSWTM,
      TSWT
    ),
    TSWF = ifelse(
      is.na(TSWF) & grepl("T", Treatment),
      TSWFA + TSWFM,
      TSWF
    ),
    NST = ifelse(
      is.na(NST) & grepl("T", Treatment),
      NSTA + NSTM,
      NST
    ),
    NSF = ifelse(
      is.na(NSF) & grepl("T", Treatment),
      NSFA + NSFM,
      NSF
    )) %>%
  arrange(Plot) %>%
  dplyr::select(Plot, Treatment, Block, Weeds, YieldT, YieldF)

data_2023A_yield_T <- data_2023A_yield %>%
  dplyr::select(Plot, Treatment, Block, Weeds, YieldT)

data_2023A_yield_T_W <- data_2023A_yield_T %>%
  mutate(Treatment = paste0(Treatment, "_", Weeds)) %>%
  dplyr::select(!Weeds)

data_2023A_yield_F <- data_2023A_yield %>%
  dplyr::select(Plot, Treatment, Block, Weeds, YieldF)

data_2023A_yield_F_W <- data_2023A_yield_F %>%
  mutate(Treatment = paste0(Treatment, "_", Weeds)) %>%
  dplyr::select(!Weeds)

######################
##                  ##
## Experiment 2023B ##
##                  ##
######################

# Read data
data_2023B_yield <- read_xlsx("WCF_2023_data.xlsx", sheet = "FinalHarvestB", range = "A1:P89", col_names = TRUE) %>% 
  dplyr::select(Plot, Treatment, Block, Weeds, 
                "BiomassSeedsTriticale", 
                "1000SeedWeightTriticale", 
                "NSpikesTriticale", 
                "NSeedsTriticale", 
                "BiomassSeedsFaba", 
                "1000SeedWeightFaba", 
                "NPodsFaba", 
                "NSeedsFaba") %>% 
  rename("YieldT" = "BiomassSeedsTriticale",
         "YieldF" = "BiomassSeedsFaba",
         "NFT" = "NSpikesTriticale",
         "NFF" = "NPodsFaba",
         "TSWT" = "1000SeedWeightTriticale",
         "TSWF" = "1000SeedWeightFaba",
         "NST" = "NSeedsTriticale",
         "NSF" = "NSeedsFaba") %>%
  mutate(YieldT = as.numeric(YieldT),
         NFT = as.numeric(NFT),
         TSWT = as.numeric(TSWT),
         NST = as.numeric(NST),
         YieldF = as.numeric(YieldF),
         NFF = as.numeric(NFF),
         TSWF = as.numeric(TSWF),
         NSF = as.numeric(NSF)) %>% 
  arrange(Plot) %>% 
  filter(!Treatment %in% c("T-25", "TF-M-25", "1T:1F-25", "F-25")) %>% 
  dplyr::select(Plot, Treatment, Block, Weeds, YieldT, YieldF) %>% 
  mutate(Block = Block + 5)

data_2023B_yield_T <- data_2023B_yield %>% 
  dplyr::select(Plot, Treatment, Block, Weeds, YieldT)

data_2023B_yield_T_W <- data_2023B_yield_T %>%
  mutate(Treatment = paste0(Treatment, "_", Weeds)) %>% 
  dplyr::select(!Weeds)

data_2023B_yield_F <- data_2023B_yield %>% 
  dplyr::select(Plot, Treatment, Block, Weeds, YieldF)

data_2023B_yield_F_W <- data_2023B_yield_F %>%
  mutate(Treatment = paste0(Treatment, "_", Weeds)) %>% 
  dplyr::select(!Weeds)

######################
##                  ##
## Experiment 2024 ##
##                  ##
######################

# Read data
data_2024_yield <- read_xlsx("WCF_2024_data.xlsx", sheet = "FinalHarvest", range = "A1:O41", col_names = TRUE) %>% 
  dplyr::select(Plot, Treatment, Block, Weeds, 
                "BiomassSeedsTriticale", 
                "TSWTriticale", 
                "NSpikesTriticale", 
                "NSeedsTriticale", 
                "BiomassSeedsFaba", 
                "TSWFaba", 
                "NPodsFaba", 
                "NSeedsFaba") %>% 
  rename("YieldT" = "BiomassSeedsTriticale",
         "YieldF" = "BiomassSeedsFaba",
         "NFT" = "NSpikesTriticale",
         "NFF" = "NPodsFaba",
         "TSWT" = "TSWTriticale",
         "TSWF" = "TSWFaba",
         "NST" = "NSeedsTriticale",
         "NSF" = "NSeedsFaba") %>%
  mutate(YieldT = as.numeric(YieldT),
         NFT = as.numeric(NFT),
         TSWT = as.numeric(TSWT),
         NST = as.numeric(NST),
         YieldF = as.numeric(YieldF),
         NFF = as.numeric(NFF),
         TSWF = as.numeric(TSWF),
         NSF = as.numeric(NSF)) %>% 
  arrange(Plot) %>% 
  filter(Treatment != "W")

data_2024_yield_T <- data_2024_yield %>% 
  dplyr::select(Plot, Treatment, Block, Weeds, YieldT)

data_2024_yield_T_W <- data_2024_yield_T %>%
  mutate(Treatment = paste0(Treatment, "_", Weeds)) %>% 
  dplyr::select(!Weeds)

data_2024_yield_F <- data_2024_yield %>% 
  dplyr::select(Plot, Treatment, Block, Weeds, YieldF)

data_2024_yield_F_W <- data_2024_yield_F %>%
  mutate(Treatment = paste0(Treatment, "_", Weeds)) %>% 
  dplyr::select(!Weeds)

## Create visualizations ##

# Create yield plots
## Figure 2: 2022
cereals <- c("R", "B", "T", "W")
legumes <- c("P", "F")

data_2022_yield_groups <- data_2022_yield %>%
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

data_2022_yield_cereal <- data_2022_yield_groups %>%
  filter(!is.na(CerealGroup)) %>%
  group_by(CerealGroup, Treatment) %>%
  summarise(
    mean_yield = mean(YieldC, na.rm = TRUE),
    sd_yield = sd(YieldC, na.rm = TRUE),
    n = sum(!is.na(YieldC)),
    .groups = "drop"
  )

data_2022_yield_legume <- data_2022_yield_groups %>%
  filter(!is.na(LegumeGroup)) %>%
  group_by(LegumeGroup, Treatment) %>%
  summarise(
    mean_yield = mean(YieldL, na.rm = TRUE),
    sd_yield = sd(YieldL, na.rm = TRUE),
    n = sum(!is.na(YieldL)),
    .groups = "drop"
  )

cereal_levels <- c(
  "R", "1R:1P", "1R:1F",
  "B", "1B:1P", "1B:1F",
  "T", "1T:1P", "1T:1F",
  "W", "1W:1P", "1W:1F"
)

data_2022_yield_cereal2 <- data_2022_yield_cereal %>%
  mutate(Treatment = factor(Treatment, levels = cereal_levels))

legume_levels <- c(
  "P", "1R:1P", "1B:1P", "1T:1P", "1W:1P",
  "F", "1R:1F", "1B:1F", "1T:1F", "1W:1F"
)

data_2022_yield_legume2 <- data_2022_yield_legume %>%
  mutate(Treatment = factor(Treatment, levels = legume_levels))

cereal_cld_2022 <- lapply(unique(data_2022_yield_cereal2$CerealGroup), function(g) {
  treatments <- cereal_levels[grepl(g, cereal_levels)]
  df_sub <- data_2022_yield %>%
    filter(Treatment %in% treatments & !is.na(YieldC))
    cld_df <- get_cld(df_sub, "YieldC", "CerealGroup", treatments)
    cld_df$CerealGroup <- g
    return(cld_df)
}) %>% bind_rows()

legume_cld_2022 <- lapply(unique(data_2022_yield_legume2$LegumeGroup), function(g) {
  if(g == "P") {
    treatments <- legume_levels[legume_levels %in% c("P", "1R:1P", "1B:1P", "1T:1P", "1W:1P")]
  } else if(g == "F") {
    treatments <- legume_levels[legume_levels %in% c("F", "1R:1F", "1B:1F", "1T:1F", "1W:1F")]
  } else {
    treatments <- character(0)
  }
  df_sub <- data_2022_yield %>%
    filter(Treatment %in% treatments & !is.na(YieldL))
  if(nrow(df_sub) > 0) {
    cld_df <- get_cld2(df_sub, "YieldL", "LegumeGroup", treatments)
    cld_df$LegumeGroup <- g
    return(cld_df)
  }
}) %>% bind_rows()

data_2022_yield_cereal_cld <- data_2022_yield_cereal2 %>%
  left_join(cereal_cld_2022, by = c("Treatment", "CerealGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2022_yield_legume_cld <- data_2022_yield_legume2 %>%
  left_join(legume_cld_2022, by = c("Treatment", "LegumeGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

  
spacer_cereals <- tibble(
  CerealGroup = c("R", "B", "T"),
  Treatment = factor(c("spacer1", "spacer2", "spacer3"), levels = c(levels(data_2022_yield_cereal_cld$Treatment), "spacer1", "spacer2", "spacer3")),
  mean_yield = NA_real_,
  sd_yield = NA_real_,
  n = NA_integer_,
  .group = NA_character_
)

data_2022_yield_cereal_cld2 <- data_2022_yield_cereal_cld %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "R", "1R:1P", "1R:1F", "spacer1",
    "B", "1B:1P", "1B:1F", "spacer2",
    "T", "1T:1P", "1T:1F", "spacer3",
    "W", "1W:1P", "1W:1F"
  ))) %>%
  bind_rows(spacer_cereals) %>% 
  mutate(Treatment = factor(Treatment, levels = c("R", "1R:1F", "1R:1P", "spacer1",
                                                  "B", "1B:1F", "1B:1P", "spacer2",
                                                  "T", "1T:1F", "1T:1P", "spacer3",
                                                  "W", "1W:1F", "1W:1P"))) %>% 
  mutate(CerealGroup = factor(CerealGroup, levels = c("R", "B", "T", "W")))

spacer_legumes <- tibble(
  LegumeGroup = c("P"),
  Treatment = factor("spacer1", levels = c(levels(data_2022_yield_legume_cld$Treatment), "spacer1")),
  mean_yield = NA_real_,
  sd_yield = NA_real_,
  n = NA_integer_,
  .group = NA_character_
)

data_2022_yield_legume_cld2 <- data_2022_yield_legume_cld %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "F", "1R:1F", "1B:1F", "1T:1F", "1W:1F",
    "spacer1",
    "P", "1R:1P", "1B:1P", "1T:1P", "1W:1P"
  ))) %>%
  bind_rows(spacer_legumes) %>% 
  mutate(Treatment = factor(Treatment, levels = c("F", "1R:1F", "1B:1F", "1T:1F", "1W:1F", 
                                                  "spacer1",
                                                  "P", "1R:1P", "1B:1P", "1T:1P", "1W:1P"
                                                  ))) %>% 
  mutate(LegumeGroup = factor(LegumeGroup, levels = c("F", "P")))

plot_cereal_2022 <- ggplot(data_2022_yield_cereal_cld2, aes(x = Treatment, y = mean_yield, fill = CerealGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_yield - sd_yield, ymax = mean_yield + sd_yield), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_yield + sd_yield + 40), na.rm = TRUE, size = 6) +
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
    subtitle = "Cereal yield",
    x = "Treatment",
    y = "Yield (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank())

plot_legume_2022 <- ggplot(data_2022_yield_legume_cld2, aes(x = Treatment, y = mean_yield, fill = LegumeGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_yield - sd_yield, ymax = mean_yield + sd_yield), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_yield + sd_yield + 40), na.rm = TRUE, size = 6) +
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
    subtitle = "Legume yield",
    x = "Treatment",
    y = "Yield (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank())

combined_plot_2022 <- plot_cereal_2022 / plot_legume_2022 +
    plot_layout(ncol = 1) +
    plot_annotation(tag_levels = "a")
combined_plot_2022

# Figure 4: 2023A
data_2023A_yield_NW <- data_2023A_yield %>% 
  mutate(Treatment = gsub("TF", "1T:1F", Treatment)) %>% 
  filter(Weeds == "N")

cereals <- c("T+", "T")
legumes <- c("F+", "F")

data_2023A_yield_groups <- data_2023A_yield_NW %>%
  mutate(
    CerealGroup = ifelse(grepl("T", Treatment), "T", NA_character_),
    LegumeGroup = ifelse(grepl("F", Treatment), "F", NA_character_)
  )

data_2023A_yield_cereal <- data_2023A_yield_groups %>%
  filter(!is.na(CerealGroup)) %>%
  group_by(CerealGroup, Treatment) %>%
  summarise(
    mean_yield = mean(YieldT, na.rm = TRUE),
    sd_yield = sd(YieldT, na.rm = TRUE),
    n = sum(!is.na(YieldT)),
    .groups = "drop"
  )

data_2023A_yield_legume <- data_2023A_yield_groups %>%
  filter(!is.na(LegumeGroup)) %>%
  group_by(LegumeGroup, Treatment) %>%
  summarise(
    mean_yield = mean(YieldF, na.rm = TRUE),
    sd_yield = sd(YieldF, na.rm = TRUE),
    n = sum(!is.na(YieldF)),
    .groups = "drop"
  )

cereal_levels <- c(
  "T+", "T", "3T:1F", "1T:1F", "1T:1F-M", "1T:3F"
)

data_2023A_yield_cereal2 <- data_2023A_yield_cereal %>%
  mutate(Treatment = factor(Treatment, levels = cereal_levels))

legume_levels <- c(
  "F+", "F", "1T:3F", "1T:1F", "1T:1F-M", "3T:1F"
)

data_2023A_yield_legume2 <- data_2023A_yield_legume %>%
  mutate(Treatment = factor(Treatment, levels = legume_levels))

cereal_cld_2023A <- lapply(unique(data_2023A_yield_cereal2$CerealGroup), function(g) {
  treatments <- cereal_levels[grepl(g, cereal_levels)]
  df_sub <- data_2023A_yield_NW %>%
    filter(Treatment %in% treatments & !is.na(YieldT))
    cld_df <- get_cld(df_sub, "YieldT", "CerealGroup", treatments)
    cld_df$CerealGroup <- g
    return(cld_df)
}) %>% bind_rows()

legume_cld_2023A <- lapply(unique(data_2023A_yield_legume2$LegumeGroup), function(g) {
  treatments <- legume_levels[grepl(g, legume_levels)]
  df_sub <- data_2023A_yield_NW %>%
    filter(Treatment %in% treatments & !is.na(YieldF))
  if(nrow(df_sub) > 0) {
    cld_df <- get_cld2(df_sub, "YieldF", "LegumeGroup", treatments)
    cld_df$LegumeGroup <- g
    return(cld_df)
  }
}) %>% bind_rows()

data_2023A_yield_cereal_cld <- data_2023A_yield_cereal2 %>%
  left_join(cereal_cld_2023A, by = c("Treatment", "CerealGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2023A_yield_legume_cld <- data_2023A_yield_legume2 %>%
  left_join(legume_cld_2023A, by = c("Treatment", "LegumeGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2023A_yield_cereal_cld2 <- data_2023A_yield_cereal_cld %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "T+", "T", "3T:1F", "1T:1F", "1T:1F-M", "1T:3F"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("T+", "T", "3T:1F", "1T:1F", "1T:1F-M", "1T:3F"))) %>% 
  mutate(CerealGroup = factor(CerealGroup, levels = c("T")))

data_2023A_yield_legume_cld2 <- data_2023A_yield_legume_cld %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "F+", "F", "1T:3F", "1T:1F", "1T:1F-M", "3T:1F"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("F+", "F", "1T:3F", "1T:1F", "1T:1F-M", "3T:1F"
                                                  ))) %>% 
  mutate(LegumeGroup = factor(LegumeGroup, levels = c("F")))

plot_cereal_2023A <- ggplot(data_2023A_yield_cereal_cld2, aes(x = Treatment, y = mean_yield, fill = CerealGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_yield - sd_yield, ymax = mean_yield + sd_yield), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_yield + sd_yield + 40), na.rm = TRUE, size = 6) +
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
    subtitle = "Cereal yield",
    x = "Treatment",
    y = "Yield (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none")

plot_legume_2023A <- ggplot(data_2023A_yield_legume_cld2, aes(x = Treatment, y = mean_yield, fill = LegumeGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_yield - sd_yield, ymax = mean_yield + sd_yield), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_yield + sd_yield + 40), na.rm = TRUE, size = 6) +
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
    subtitle = "Legume yield",
    x = "Treatment",
    y = "Yield (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none")

combined_plot_2023A <- plot_cereal_2023A / plot_legume_2023A +
    plot_layout(ncol = 1) +
    plot_annotation(tag_levels = "a")
combined_plot_2023A

# Figure 6: 2023B
data_2023B_yield_NW <- data_2023B_yield %>% 
  mutate(Treatment = gsub("TF", "1T:1F", Treatment)) %>% 
  filter(Weeds == "N")

cereals <- c("T", "T-375")
legumes <- c("F", "F-375")

data_2023B_yield_groups <- data_2023B_yield_NW %>%
  mutate(
    CerealGroup = ifelse(grepl("T", Treatment), "T", NA_character_),
    LegumeGroup = ifelse(grepl("F", Treatment), "F", NA_character_)
  )

data_2023B_yield_cereal <- data_2023B_yield_groups %>%
  filter(!is.na(CerealGroup)) %>%
  group_by(CerealGroup, Treatment) %>%
  summarise(
    mean_yield = mean(YieldT, na.rm = TRUE),
    sd_yield = sd(YieldT, na.rm = TRUE),
    n = sum(!is.na(YieldT)),
    .groups = "drop"
  )

data_2023B_yield_legume <- data_2023B_yield_groups %>%
  filter(!is.na(LegumeGroup)) %>%
  group_by(LegumeGroup, Treatment) %>%
  summarise(
    mean_yield = mean(YieldF, na.rm = TRUE),
    sd_yield = sd(YieldF, na.rm = TRUE),
    n = sum(!is.na(YieldF)),
    .groups = "drop"
  )

cereal_levels <- c(
  "T", "1T:1F", "1T:1F-M", "T-375", "1T:1F-375", "1T:1F-M-375"
)

data_2023B_yield_cereal2 <- data_2023B_yield_cereal %>%
  mutate(Treatment = factor(Treatment, levels = cereal_levels))

legume_levels <- c(
  "F", "1T:1F", "1T:1F-M", "F-375", "1T:1F-375", "1T:1F-M-375"
)

data_2023B_yield_legume2 <- data_2023B_yield_legume %>%
  mutate(Treatment = factor(Treatment, levels = legume_levels))

cereal_cld_2023B <- lapply(unique(data_2023B_yield_cereal2$CerealGroup), function(g) {
  treatments <- cereal_levels[grepl(g, cereal_levels)]
  df_sub <- data_2023B_yield_NW %>%
    filter(Treatment %in% treatments & !is.na(YieldT))
    cld_df <- get_cld(df_sub, "YieldT", "CerealGroup", treatments)
    cld_df$CerealGroup <- g
    return(cld_df)
}) %>% bind_rows()

legume_cld_2023B <- lapply(unique(data_2023B_yield_legume2$LegumeGroup), function(g) {
  treatments <- legume_levels[grepl(g, legume_levels)]
  df_sub <- data_2023B_yield_NW %>%
    filter(Treatment %in% treatments & !is.na(YieldF))
  if(nrow(df_sub) > 0) {
    cld_df <- get_cld2(df_sub, "YieldF", "LegumeGroup", treatments)
    cld_df$LegumeGroup <- g
    return(cld_df)
  }
}) %>% bind_rows()

data_2023B_yield_cereal_cld <- data_2023B_yield_cereal2 %>%
  left_join(cereal_cld_2023B, by = c("Treatment", "CerealGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2023B_yield_legume_cld <- data_2023B_yield_legume2 %>%
  left_join(legume_cld_2023B, by = c("Treatment", "LegumeGroup")) %>%
  mutate(.group = ifelse(is.na(.group), "", .group))

data_2023B_yield_cereal_cld2 <- data_2023B_yield_cereal_cld %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "T", "1T:1F", "1T:1F-M", "T-375", "1T:1F-375", "1T:1F-M-375"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("T", "1T:1F", "1T:1F-M", "T-375", "1T:1F-375", "1T:1F-M-375"))) %>% 
  mutate(CerealGroup = factor(CerealGroup, levels = c("T")))

data_2023B_yield_legume_cld2 <- data_2023B_yield_legume_cld %>%
  mutate(Treatment = factor(as.character(Treatment), levels = c(
    "F", "F-375", "1T:1F", "1T:1F-M", "1T:1F-375", "1T:1F-M-375"
  ))) %>%
  mutate(Treatment = factor(Treatment, levels = c("F", "1T:1F", "1T:1F-M", "F-375", "1T:1F-375", "1T:1F-M-375"
                                                  ))) %>% 
  mutate(LegumeGroup = factor(LegumeGroup, levels = c("F")))

plot_cereal_2023B <- ggplot(data_2023B_yield_cereal_cld2, aes(x = Treatment, y = mean_yield, fill = CerealGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_yield - sd_yield, ymax = mean_yield + sd_yield), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_yield + sd_yield + 40), na.rm = TRUE, size = 6) +
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
    subtitle = "Cereal yield",
    x = "Treatment",
    y = "Yield (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none")

plot_legume_2023B <- ggplot(data_2023B_yield_legume_cld2, aes(x = Treatment, y = mean_yield, fill = LegumeGroup)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean_yield - sd_yield, ymax = mean_yield + sd_yield), width = 0.3, size = 1.2, na.rm = TRUE) +
  geom_text(aes(label = gsub(" ", "", .group), y = mean_yield + sd_yield + 40), na.rm = TRUE, size = 6) +
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
    subtitle = "Legume yield",
    x = "Treatment",
    y = "Yield (g m⁻²)"
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.position = "none")

combined_plot_2023B <- plot_cereal_2023B / plot_legume_2023B +
    plot_layout(ncol = 1) +
    plot_annotation(tag_levels = "a")
combined_plot_2023B

## Figure 8: across years including 2024
data_2024_yield_NW <- data_2024_yield %>% 
  mutate(Treatment = gsub("TF", "1T:1F", Treatment)) %>% 
  filter(Weeds == "N")

cld_cereal_2022 <- get_cld(data_2022_yield %>% filter(Treatment %in% c("T", "1T:1F")), "YieldC", "Species", c("T", "1T:1F")) %>%
  mutate(Year = "2022", TreatmentGroup = paste0(Treatment, "-2022"))

cld_cereal_2023A <- get_cld(data_2023A_yield_NW %>% filter(Treatment %in% c("T", "1T:1F")), "YieldT", "Species", c("T", "1T:1F")) %>%
  mutate(Year = "2022", TreatmentGroup = paste0(Treatment, "-2023A"))

cld_cereal_2023B <- get_cld(data_2023B_yield_NW %>% filter(Treatment %in% c("T", "1T:1F")), "YieldT", "Species", c("T", "1T:1F")) %>%
  mutate(Year = "2022", TreatmentGroup = paste0(Treatment, "-2023B"))

cld_cereal_2024 <- get_cld(data_2024_yield_NW %>% filter(Treatment %in% c("T", "1T:1F")), "YieldT", "Species", c("T", "1T:1F")) %>%
  mutate(Year = "2022", TreatmentGroup = paste0(Treatment, "-2024"))


cereal_2022 <- prepare_year_data(data_2022_yield, "2022", "YieldC") %>% filter(Species != "Faba bean")
cereal_2023A <- prepare_year_data(data_2023A_yield_NW, "2023A", "YieldT") %>% filter(Species != "Faba bean")
cereal_2023B <- prepare_year_data(data_2023B_yield_NW, "2023B", "YieldT") %>% filter(Species != "Faba bean")
cereal_2024 <- prepare_year_data(data_2024_yield_NW, "2024", "YieldT") %>% filter(Species != "Faba bean")

cld_legume_2022 <- get_cld(data_2022_yield %>% filter(Treatment %in% c("F", "1T:1F")), "YieldL", "Species", c("F", "1T:1F")) %>%
  mutate(Year = "2022", TreatmentGroup = paste0(Treatment, "-2022"))

cld_legume_2023A <- get_cld(data_2023A_yield_NW %>% filter(Treatment %in% c("F", "1T:1F")), "YieldF", "Species", c("F", "1T:1F")) %>%
  mutate(Year = "2022", TreatmentGroup = paste0(Treatment, "-2023A"))

cld_legume_2023B <- get_cld(data_2023B_yield_NW %>% filter(Treatment %in% c("F", "1T:1F")), "YieldF", "Species", c("F", "1T:1F")) %>%
  mutate(Year = "2022", TreatmentGroup = paste0(Treatment, "-2023B"))

cld_legume_2024 <- get_cld(data_2024_yield_NW %>% filter(Treatment %in% c("F", "1T:1F")), "YieldF", "Species", c("F", "1T:1F")) %>%
  mutate(Year = "2022", TreatmentGroup = paste0(Treatment, "-2024"))

legume_2022 <- prepare_year_data(data_2022_yield, "2022", "YieldL") %>% filter(Species != "Triticale")
legume_2023A <- prepare_year_data(data_2023A_yield_NW, "2023A", "YieldF") %>% filter(Species != "Triticale")
legume_2023B <- prepare_year_data(data_2023B_yield_NW, "2023B", "YieldF") %>% filter(Species != "Triticale")
legume_2024 <- prepare_year_data(data_2024_yield_NW, "2024", "YieldF") %>% filter(Species != "Triticale")

cereal_2022_cld <- cereal_2022 %>%
  left_join(cld_cereal_2022, by = c("TreatmentGroup")) %>%
  mutate(.group = gsub(" ", "", .group))

cereal_2023A_cld <- cereal_2023A %>%
  left_join(cld_cereal_2023A, by = c("TreatmentGroup")) %>%
  mutate(.group = gsub(" ", "", .group))

cereal_2023B_cld <- cereal_2023B %>%
  left_join(cld_cereal_2023B, by = c("TreatmentGroup")) %>%
  mutate(.group = gsub(" ", "", .group))

cereal_2024_cld <- cereal_2024 %>%
  left_join(cld_cereal_2024, by = c("TreatmentGroup")) %>%
  mutate(.group = gsub(" ", "", .group))

legume_2022_cld <- legume_2022 %>%
  left_join(cld_legume_2022, by = c("TreatmentGroup")) %>%
  mutate(.group = gsub(" ", "", .group))

legume_2023A_cld <- legume_2023A %>%
  left_join(cld_legume_2023A, by = c("TreatmentGroup")) %>%
  mutate(.group = gsub(" ", "", .group))

legume_2023B_cld <- legume_2023B %>%
  left_join(cld_legume_2023B, by = c("TreatmentGroup")) %>%
  mutate(.group = gsub(" ", "", .group))

legume_2024_cld <- legume_2024 %>%
  left_join(cld_legume_2024, by = c("TreatmentGroup")) %>%
  mutate(.group = gsub(" ", "", .group))

add_spacers <- function(df, species) {
  year <- unique(df$TreatmentGroup) %>%
    stringr::str_extract("(?<=-).*")
  spacer_rows <- tibble(
    TreatmentGroup = paste0("spacer-", year),
    Species = species,
    mean_yield = NA, sd_yield = NA, n = NA
  )
  bind_rows(df, spacer_rows)
}

cereal_combined <- bind_rows(
  add_spacers(cereal_2022_cld, "Triticale"),
  add_spacers(cereal_2023A_cld, "Triticale"),
  add_spacers(cereal_2023B_cld, "Triticale"),
  cereal_2024_cld
) %>% 
  mutate(TreatmentGroup = gsub("2023A", "2023-SP", TreatmentGroup),
         TreatmentGroup = gsub("2023B", "2023-RD", TreatmentGroup))


legume_combined <- bind_rows(
  add_spacers(legume_2022_cld, "Faba bean"),
  add_spacers(legume_2023A_cld, "Faba bean"),
  add_spacers(legume_2023B_cld, "Faba bean"),
  legume_2024_cld
) %>% 
  mutate(TreatmentGroup = gsub("2023A", "2023-SP", TreatmentGroup),
         TreatmentGroup = gsub("2023B", "2023-RD", TreatmentGroup))

cereal_levels <- c("T-2022", "1T:1F-2022", "spacer-2022",
                   "T-2023-SP", "1T:1F-2023-SP", "spacer-2023-SP",
                   "T-2023-RD", "1T:1F-2023-RD", "spacer-2023-RD",
                   "T-2024", "1T:1F-2024")
legume_levels <- c("F-2022", "1T:1F-2022", "spacer-2022",
                   "F-2023-SP", "1T:1F-2023-SP", "spacer-2023-SP",
                   "F-2023-RD", "1T:1F-2023-RD", "spacer-2023-RD",
                   "F-2024", "1T:1F-2024")

cereal_combined <- cereal_combined %>%
  mutate(TreatmentGroup = factor(TreatmentGroup, levels = cereal_levels))

legume_combined <- legume_combined %>%
  mutate(TreatmentGroup = factor(TreatmentGroup, levels = legume_levels))

cereal_plot_df <- cereal_combined %>%
  filter(! str_starts(as.character(TreatmentGroup), "spacer")) %>%
  mutate(
    Year      = str_extract(TreatmentGroup, "\\d{4}(?:-[A-Z]+)?"),
    Treatment = str_remove(TreatmentGroup, "-[0-9].*$")
  ) %>%
  dplyr::select(Year, Treatment, mean_yield, sd_yield, .group, Species) %>%
  mutate(
    Year      = factor(Year, levels = c("2022","2023-SP","2023-RD","2024")),
    Treatment = factor(Treatment, levels = c("T","1T:1F"))
  )

legume_plot_df <- legume_combined %>%
  filter(! str_starts(as.character(TreatmentGroup), "spacer")) %>%
  mutate(
    Year      = str_extract(TreatmentGroup, "\\d{4}(?:-[A-Z]+)?"),
    Treatment = str_remove(TreatmentGroup, "-[0-9].*$")
  ) %>%
  dplyr::select(Year, Treatment, mean_yield, sd_yield, .group, Species) %>%
  mutate(
    Year      = factor(Year, levels = c("2022","2023-SP","2023-RD","2024")),
    Treatment = factor(Treatment, levels = c("F","1T:1F"))
  )
cereal_plot_df <- cereal_plot_df %>%
  mutate(SpeciesLabel = ifelse(Species == "Intercrop",
                               "Triticale–faba bean",
                               "Triticale"))

legume_plot_df <- legume_plot_df %>%
  mutate(SpeciesLabel = ifelse(Species == "Intercrop",
                               "Triticale–faba bean",
                               "Faba bean"))

base_theme <- theme_classic(base_size = 18) +
  theme(
    strip.background = element_rect(color = "black", fill = "white"),
    strip.placement  = "outside",
    strip.text       = element_text(size = 14),
    axis.title.x     = element_blank(),
    axis.text.x      = element_text(size = 12)
  )

plot_cereal_all <- ggplot(cereal_plot_df,
                          aes(x = Treatment, y = mean_yield, fill = SpeciesLabel)) +
  geom_col(position = position_dodge(.7), width = .7, color = "black") +
  geom_errorbar(aes(ymin = mean_yield - sd_yield,
                    ymax = mean_yield + sd_yield),
                position = position_dodge(.7),
                width = .2) +
  geom_text(aes(label = .group,
                y     = mean_yield + sd_yield + 40),
            position = position_dodge(.7),
            size = 5) +
  facet_grid(. ~ Year, scales = "free_x", space = "free_x") +
  scale_fill_manual(
    name   = "Crop type",
    values = c(
      "Triticale"             = colour_vector[1],
      "Faba bean"             = colour_vector[2],
      "Triticale–faba bean"   = colour_vector[3]
    ),
    breaks = c("Triticale","Triticale–faba bean"),
    labels = c("Triticale","Triticale–faba bean")
  ) +
  labs(y = "Yield (g m⁻²)") +
  base_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

plot_legume_all <- ggplot(legume_plot_df,
                          aes(x = Treatment, y = mean_yield, fill = SpeciesLabel)) +
  geom_col(position = position_dodge(.7), width = .7, color = "black") +
  geom_errorbar(aes(ymin = mean_yield - sd_yield,
                    ymax = mean_yield + sd_yield),
                position = position_dodge(.7),
                width = .2) +
  geom_text(aes(label = .group,
                y     = mean_yield + sd_yield + 20),
            position = position_dodge(.7),
            size = 5) +
  facet_grid(. ~ Year, scales = "free_x", space = "free_x") +
  scale_fill_manual(
    name   = "Crop type",
    values = c(
      "Triticale"             = colour_vector[1],
      "Faba bean"             = colour_vector[2],
      "Triticale–faba bean"   = colour_vector[3]
    ),
    breaks = c("Faba bean","Triticale–faba bean"),
    labels = c("Faba bean","Triticale–faba bean")
  ) +
  labs(y = "Yield (g m⁻²)") +
  base_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

combined_plot_years <- plot_cereal_all / plot_legume_all +
  plot_layout(ncol = 1, guides = "collect") +
  plot_annotation(tag_levels = "a")

combined_plot_years

## LER and competitive ratio analysis ##

# Create tables
## Table 1: Relative performance 2022
sole_means_2022 <- data_2022_yield %>%
  filter(Treatment %in% c("R", "B", "T", "W", "P", "F")) %>%
  summarise(
    R = mean(YieldC[Treatment == "R"]),
    B = mean(YieldC[Treatment == "B"]),
    T = mean(YieldC[Treatment == "T"]),
    W = mean(YieldC[Treatment == "W"]),
    P = mean(YieldL[Treatment == "P"]),
    F = mean(YieldL[Treatment == "F"])
  )

data_2022_yield_LER <- data_2022_yield %>%
  mutate(
    pLERC = case_when(
      str_detect(Treatment, "^[^:]*R") ~ YieldC / sole_means_2022$R,
      str_detect(Treatment, "^[^:]*B") ~ YieldC / sole_means_2022$B,
      str_detect(Treatment, "^[^:]*T") ~ YieldC / sole_means_2022$T,
      str_detect(Treatment, "^[^:]*W") ~ YieldC / sole_means_2022$W,
      TRUE ~ NA_real_
    ),
    pLERL = case_when(
      str_detect(Treatment, "(?<=:)[^:]*P") ~ YieldL / sole_means_2022$P,
      str_detect(Treatment, "(?<=:)[^:]*F") ~ YieldL / sole_means_2022$F,
      TRUE ~ NA_real_
    ),
    NERC = case_when(
      str_detect(Treatment, "^[^:]*R") ~ YieldC / (0.5 * sole_means_2022$R),
      str_detect(Treatment, "^[^:]*B") ~ YieldC / (0.5 * sole_means_2022$B),
      str_detect(Treatment, "^[^:]*T") ~ YieldC / (0.5 * sole_means_2022$T),
      str_detect(Treatment, "^[^:]*W") ~ YieldC / (0.5 * sole_means_2022$W),
      TRUE ~ NA_real_
    ),
    NERL = case_when(
      str_detect(Treatment, "(?<=:)[^:]*P") ~ YieldL / (0.5 * sole_means_2022$P),
      str_detect(Treatment, "(?<=:)[^:]*F") ~ YieldL / (0.5 * sole_means_2022$F),
      TRUE ~ NA_real_
    ),
    LER = pLERC + pLERL,
    CR = (pLERC / 0.5) / (pLERL / 0.5)
  )

weed_biomass_2022 <- c(
  "Rye"       =  24.0,
  "Barley"    =  46.8,
  "Triticale" =  49.4,
  "Wheat"     = 179.4,
  "Pea"       = 171.8,
  "Faba"      = 453.8
  )

data_2022_yield_LER_agg <- data_2022_yield_LER %>%
  group_by(Treatment) %>%
  summarise(
    mPLERC = mean(pLERC, na.rm = TRUE),
    mPLERL = mean(pLERL, na.rm = TRUE),
    mNERC = mean(NERC, na.rm = TRUE),
    mNERL = mean(NERL, na.rm = TRUE),
    mLER = mean(LER, na.rm = TRUE),
    mCR = mean(CR, na.rm = TRUE),
    P_PLERC = tryCatch(t.test(pLERC, mu = 0.5)$p.value, error = function(e) NA_real_),
    P_PLERL = tryCatch(t.test(pLERL, mu = 0.5)$p.value, error = function(e) NA_real_),
    P_NERC = tryCatch(t.test(NERC, mu = 1.0)$p.value, error = function(e) NA_real_),
    P_NERL = tryCatch(t.test(NERL, mu = 1.0)$p.value, error = function(e) NA_real_),
    P_LER = tryCatch(t.test(LER, mu = 1.0)$p.value, error = function(e) NA_real_),
    P_CR = tryCatch(t.test(CR, mu = 1.0)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  na.omit() %>%
  mutate(
    WR = c(
      weed_biomass_2022["Faba"] / weed_biomass_2022["Barley"],
      weed_biomass_2022["Pea"] / weed_biomass_2022["Barley"],
      weed_biomass_2022["Faba"] / weed_biomass_2022["Rye"],
      weed_biomass_2022["Pea"] / weed_biomass_2022["Rye"],
      weed_biomass_2022["Faba"] / weed_biomass_2022["Triticale"],
      weed_biomass_2022["Pea"] / weed_biomass_2022["Triticale"],
      weed_biomass_2022["Faba"] / weed_biomass_2022["Wheat"],
      weed_biomass_2022["Pea"] / weed_biomass_2022["Wheat"]
    )
  )

data_2022_yield_LER_agg

## Table 2: Relative performance 2023-SP
sole_means_2023A <- data_2023A_yield_NW %>%
  filter(Treatment %in% c("T", "F")) %>%
  summarise(
    T = mean(YieldT[Treatment == "T"], na.rm = TRUE),
    F = mean(YieldF[Treatment == "F"], na.rm = TRUE)
  )

weed_biomass_2023A <- c(
  "T+"       =  4.8,
  "T"        = 15.3,
  "3T:1F"    = 12.6,
  "1T:1F-M"  = 15.8,
  "1T:1F"    = 20.3,
  "1T:3F"    = 30.0,
  "F"        = 88.4,
  "F+"       = 76.7
  )

data_2023A_yield_LER <- data_2023A_yield_NW %>%
  mutate(Treatment = gsub("TF", "1T:1F", Treatment)) %>%
  mutate(
    pLERT = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M", "1T:3F", "3T:1F") ~ YieldT / sole_means_2023A$T,
      TRUE ~ NA_real_
    ),
    pLERF = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M", "1T:3F", "3T:1F") ~ YieldF / sole_means_2023A$F,
      TRUE ~ NA_real_
    ),
    NERT = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M") ~ YieldT / (0.5 * sole_means_2023A$T),
      Treatment %in% c("1T:3F") ~ YieldT / (0.25 * sole_means_2023A$T),
      Treatment %in% c("3T:1F") ~ YieldT / (0.75 * sole_means_2023A$T),
      TRUE ~ NA_real_
    ),
    NERF = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M") ~ YieldF / (0.5 * sole_means_2023A$F),
      Treatment %in% c("1T:3F") ~ YieldF / (0.75 * sole_means_2023A$F),
      Treatment %in% c("3T:1F") ~ YieldF / (0.25 * sole_means_2023A$F),
      TRUE ~ NA_real_
    ),
    LER = pLERT + pLERF,
    CR = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M") ~ (pLERT / 0.5) / (pLERF / 0.5),
      Treatment == "1T:3F" ~ (pLERT / 0.25) / (pLERF / 0.75),
      Treatment == "3T:1F" ~ (pLERT / 0.75) / (pLERF / 0.25),
      TRUE ~ NA_real_
    )
  )

data_2023A_yield_LER_agg <- data_2023A_yield_LER %>% 
  group_split(Treatment) %>% 
  map_dfr(function(df) {
    treat <- unique(df$Treatment)
    list(
      Treatment = treat,
      mPLERT = mean(df$pLERT, na.rm = TRUE),
      mPLERF = mean(df$pLERF, na.rm = TRUE),
      mNERT = mean(df$NERT, na.rm = TRUE),
      mNERF = mean(df$NERF, na.rm = TRUE),
      mLER   = mean(df$LER, na.rm = TRUE),
      mCR    = mean(df$CR, na.rm = TRUE),
      P_PLERC = tryCatch(
        t.test(df$pLERT, mu = 1.0)$p.value,
        error = function(e) NA_real_
      ),
      P_PLERL = tryCatch(
        t.test(df$pLERF, mu = 1.0)$p.value,
        error = function(e) NA_real_
      ),
      P_NERC = tryCatch(
        t.test(df$NERT, mu = 1.0)$p.value,
        error = function(e) NA_real_
      ),
      P_NERL = tryCatch(
        t.test(df$NERF, mu = 1.0)$p.value,
        error = function(e) NA_real_
      ),
      P_LER = tryCatch(
        t.test(df$LER, mu = 1.0)$p.value,
        error = function(e) NA_real_
      ),
      P_CR = tryCatch(
        t.test(df$CR, mu = 1.0)$p.value,
        error = function(e) NA_real_
      )
    )
  }) %>% 
  filter(Treatment %in% c("1T:1F", "1T:1F-M", "1T:3F", "3T:1F")) %>% 
  mutate(
    WR = weed_biomass_2023A["F"] / weed_biomass_2023A["T"]
  )

data_2023A_yield_LER_agg

## Table 3: Relative performance 2023-RD
sole_means_2023B <- data_2023B_yield_NW %>%
  filter(Treatment %in% c("T", "F", "T-375", "F-375")) %>%
  summarise(
    T = mean(YieldT[Treatment == "T"], na.rm = TRUE),
    F = mean(YieldF[Treatment == "F"], na.rm = TRUE),
    T375 = mean(YieldT[Treatment == "T-375"], na.rm = TRUE),
    F375 = mean(YieldF[Treatment == "F-375"], na.rm = TRUE)
  )

weed_biomass_2023B <- c(
  "T"           =  7.2,
  "T-375"       = 16.5,
  "1T:1F"       = 13.6,
  "1T:1F-M"     = 14.8,
  "1T:1F-375"   = 31.2,
  "1T:1F-M-375" = 25.3,
  "F"           = 53.9,
  "F-375"       = 66.2
)
  
data_2023B_yield_LER <- data_2023B_yield_NW %>%
  mutate(Treatment = gsub("TF", "1T:1F", Treatment)) %>%
  mutate(
    pLERT = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M") ~ YieldT / sole_means_2023B$T,
      Treatment %in% c("1T:1F-375", "1T:1F-M-375") ~ YieldT / sole_means_2023B$T375,
      TRUE ~ NA_real_
    ),
    pLERF = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M") ~ YieldF / sole_means_2023B$F,
      Treatment %in% c("1T:1F-375", "1T:1F-M-375") ~ YieldF / sole_means_2023B$F375,
      TRUE ~ NA_real_
    ),
    NERT = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M") ~ YieldT / (0.5 * sole_means_2023B$T),
      Treatment %in% c("1T:1F-375", "1T:1F-M-375") ~ YieldT / (0.5 * sole_means_2023B$T375),
      TRUE ~ NA_real_
    ),
    NERF = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M") ~ YieldF / (0.5 * sole_means_2023B$F),
      Treatment %in% c("1T:1F-375", "1T:1F-M-375") ~ YieldF / (0.5 * sole_means_2023B$F375),
      TRUE ~ NA_real_
    ),
    LER = pLERT + pLERF,
    CR = (pLERT / 0.5) / (pLERF / 0.5)
  )

data_2023B_yield_LER_agg <- data_2023B_yield_LER %>% 
  group_by(Treatment) %>% 
  summarise(
    mPLERT = mean(pLERT, na.rm = TRUE),
    mPLERF = mean(pLERF, na.rm = TRUE),
    mNERT = mean(NERT, na.rm = TRUE),
    mNERF = mean(NERF, na.rm = TRUE),
    mLER = mean(LER, na.rm = TRUE),
    mCR = mean(CR, na.rm = TRUE),
    P_PLERC = tryCatch(t.test(pLERT, mu = 0.5)$p.value, error = function(e) NA_real_),
    P_PLERL = tryCatch(t.test(pLERF, mu = 0.5)$p.value, error = function(e) NA_real_),
    P_NERC = tryCatch(t.test(NERT, mu = 1.0)$p.value, error = function(e) NA_real_),
    P_NERL = tryCatch(t.test(NERF, mu = 1.0)$p.value, error = function(e) NA_real_),
    P_LER = tryCatch(t.test(LER, mu = 1.0)$p.value, error = function(e) NA_real_),
    P_CR = tryCatch(t.test(CR, mu = 1.0)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>% 
  filter(Treatment %in% c("1T:1F", "1T:1F-M", "1T:1F-375", "1T:1F-M-375")) %>%
  mutate(
    WR = c(
      weed_biomass_2023B["F"] / weed_biomass_2023B["T"],
      weed_biomass_2023B["F-375"] / weed_biomass_2023B["T-375"],
      weed_biomass_2023B["F"] / weed_biomass_2023B["T"],
      weed_biomass_2023B["F-375"] / weed_biomass_2023B["T-375"]
    )
  )

data_2023B_yield_LER_agg

# Table 4: across years
sole_means_2024 <- data_2024_yield_NW %>%
  filter(Treatment %in% c("T", "F")) %>%
  summarise(
    T = mean(YieldT[Treatment == "T"], na.rm = TRUE),
    F = mean(YieldF[Treatment == "F"], na.rm = TRUE)
  )

weed_biomass_2024 <- c(
  "T" = 49.4,
  "1T:1F" = 99.7,
  "F" = 453.8
)

data_2024_yield_LER <- data_2024_yield_NW %>%
  mutate(Treatment = gsub("TF", "1T:1F", Treatment)) %>%
  mutate(
    pLERT = case_when(
      Treatment %in% c("1T:1F") ~ YieldT / sole_means_2024$T,
      TRUE ~ NA_real_
    ),
    pLERF = case_when(
      Treatment %in% c("1T:1F") ~ YieldF / sole_means_2024$F,
      TRUE ~ NA_real_
    ),
    NERT = case_when(
      Treatment %in% c("1T:1F") ~ YieldT / (0.5 * sole_means_2024$T),
      TRUE ~ NA_real_
    ),
    NERF = case_when(
      Treatment %in% c("1T:1F") ~ YieldF / (0.5 * sole_means_2024$F),
      TRUE ~ NA_real_
    ),
    LER = pLERT + pLERF,
    CR = (pLERT / 0.5) / (pLERF / 0.5)
  )

data_2024_yield_LER_agg <- data_2024_yield_LER %>% 
  group_by(Treatment) %>% 
  summarise(
    mPLERT = mean(pLERT, na.rm = TRUE),
    mPLERF = mean(pLERF, na.rm = TRUE),
    mNERT = mean(NERT, na.rm = TRUE),
    mNERF = mean(NERF, na.rm = TRUE),
    mLER = mean(LER, na.rm = TRUE),
    mCR = mean(CR, na.rm = TRUE),
    P_PLERC = tryCatch(t.test(pLERT, mu = 0.5)$p.value, error = function(e) NA_real_),
    P_PLERL = tryCatch(t.test(pLERF, mu = 0.5)$p.value, error = function(e) NA_real_),
    P_NERC = tryCatch(t.test(NERT, mu = 1.0)$p.value, error = function(e) NA_real_),
    P_NERL = tryCatch(t.test(NERF, mu = 1.0)$p.value, error = function(e) NA_real_),
    P_LER = tryCatch(t.test(LER, mu = 1.0)$p.value, error = function(e) NA_real_),
    P_CR = tryCatch(t.test(CR, mu = 1.0)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>% 
  filter(Treatment %in% c("1T:1F")) %>%
  mutate(
    WR = weed_biomass_2024["F"] / weed_biomass_2024["T"]
  )

data_2024_yield_LER_agg

data_yield_LER_all <- bind_rows(
  data_2022_yield_LER_agg %>% mutate(Year = "2022"),
  data_2023A_yield_LER_agg %>% rename(mPLERC = mPLERT, mPLERL = mPLERF, mNERC = mNERT, mNERL = mNERF) %>% mutate(Year = "2023-SP"), 
  data_2023B_yield_LER_agg %>% rename(mPLERC = mPLERT, mPLERL = mPLERF, mNERC = mNERT, mNERL = mNERF) %>% mutate(Year = "2023-RD"), 
  data_2024_yield_LER_agg %>% rename(mPLERC = mPLERT, mPLERL = mPLERF, mNERC = mNERT, mNERL = mNERF) %>% mutate(Year = "2024"), 
) %>%
 filter(Treatment == "1T:1F")

## Weed effect analysis ##

# Table 5: Across years, compare sole crops with and without herbicide
data_2023A_sole_W <- data_2023A_yield %>% 
    filter(Treatment %in% c("T", "F")) %>%
    pivot_wider(
      id_cols = c(Plot, Treatment, Block),
      names_from = Weeds,
      values_from = c(YieldT, YieldF),
      names_sep = "_"
    ) %>% 
  mutate(
    Yield_Y = case_when(
      Treatment == "T" ~ YieldT_Y,
      Treatment == "F" ~ YieldF_Y,
      TRUE ~ NA_real_
    ),
    Yield_N = case_when(
      Treatment == "T" ~ YieldT_N,
      Treatment == "F" ~ YieldF_N,
      TRUE ~ NA_real_
    )
  ) %>% 
  dplyr::select(-YieldT_Y, -YieldT_N, -YieldF_Y, -YieldF_N) %>% 
  mutate(RY = Yield_Y / Yield_N)

data_2023A_sole_W_agg <- data_2023A_sole_W %>% 
  group_by(Treatment) %>% 
  summarise(
    mYield_Y = mean(Yield_Y, na.rm = TRUE),
    mYield_N = mean(Yield_N, na.rm = TRUE),
    mRY = mean(RY, na.rm = TRUE),
    P_val = t.test(RY, mu = 1)$p.value,
    .groups = "drop"
  ) %>% 
  mutate(Experiment = "2023-SP")

data_2023B_sole_W <- data_2023B_yield %>% 
    filter(Treatment %in% c("T", "F", "T-375", "F-375")) %>%
    pivot_wider(
      id_cols = c(Plot, Treatment, Block),
      names_from = Weeds,
      values_from = c(YieldT, YieldF),
      names_sep = "_"
    ) %>% 
  mutate(
    Yield_Y = case_when(
      Treatment %in% c("T", "T-375") ~ YieldT_Y,
      Treatment %in% c("F", "F-375") ~ YieldF_Y,
      TRUE ~ NA_real_
    ),
    Yield_N = case_when(
      Treatment %in% c("T", "T-375") ~ YieldT_N,
      Treatment %in% c("F", "F-375") ~ YieldF_N,
      TRUE ~ NA_real_
    )
  ) %>% 
  dplyr::select(-YieldT_Y, -YieldT_N, -YieldF_Y, -YieldF_N) %>% 
  mutate(RY = Yield_Y / Yield_N)

data_2023B_sole_W_agg <- data_2023B_sole_W %>% 
  group_by(Treatment) %>% 
  summarise(
    mYield_Y = mean(Yield_Y, na.rm = TRUE),
    mYield_N = mean(Yield_N, na.rm = TRUE),
    mRY = mean(RY, na.rm = TRUE),
    P_val = t.test(RY, mu = 1)$p.value,
    .groups = "drop"
  ) %>% 
  mutate(Experiment = "2023-RD")

data_2024_sole_W <- data_2024_yield %>% 
    filter(Treatment %in% c("T", "F")) %>%
    pivot_wider(
      id_cols = c(Plot, Treatment, Block),
      names_from = Weeds,
      values_from = c(YieldT, YieldF),
      names_sep = "_"
    ) %>% 
  mutate(
    Yield_Y = case_when(
      Treatment == "T" ~ YieldT_Y,
      Treatment == "F" ~ YieldF_Y,
      TRUE ~ NA_real_
    ),
    Yield_N = case_when(
      Treatment == "T" ~ YieldT_N,
      Treatment == "F" ~ YieldF_N,
      TRUE ~ NA_real_
    )
  ) %>% 
  dplyr::select(-YieldT_Y, -YieldT_N, -YieldF_Y, -YieldF_N) %>% 
  mutate(RY = Yield_Y / Yield_N)

data_2024_sole_W_agg <- data_2024_sole_W %>% 
  group_by(Treatment) %>% 
  summarise(
    mYield_Y = mean(Yield_Y, na.rm = TRUE),
    mYield_N = mean(Yield_N, na.rm = TRUE),
    mRY = mean(RY, na.rm = TRUE),
    P_val = t.test(RY, mu = 1)$p.value,
    .groups = "drop"
  ) %>% 
  mutate(Experiment = "2024")

data_sole_W_all <- bind_rows(
  data_2023A_sole_W_agg,
  data_2023B_sole_W_agg,
  data_2024_sole_W_agg
)

# Table 6: Across years, compare intercrops with and without herbicide
sole_means_2023A_W <- data_2023A_yield %>%
  filter(Treatment %in% c("T", "F")) %>%
  group_by(Treatment, Weeds) %>% 
  summarise(
    T = mean(YieldT[Treatment == "T"], na.rm = TRUE),
    F = mean(YieldF[Treatment == "F"], na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  mutate(Treatment = paste0(Treatment, "_", Weeds),
    Yield = ifelse(is.na(T), F, T)) %>% 
  dplyr::select(Treatment, Yield)
sole_means_2023A_W2 <- sole_means_2023A_W$Yield
names(sole_means_2023A_W2) <- sole_means_2023A_W$Treatment

sole_means_2023B_W <- data_2023B_yield %>%
  filter(Treatment %in% c("T", "F", "T-375", "F-375")) %>%
  group_by(Treatment, Weeds) %>% 
  summarise(
    T = mean(YieldT[Treatment == "T"], na.rm = TRUE),
    F = mean(YieldF[Treatment == "F"], na.rm = TRUE),
    T375 = mean(YieldT[Treatment == "T-375"], na.rm = TRUE),
    F375 = mean(YieldF[Treatment == "F-375"], na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  mutate(Treatment = paste0(Treatment, "_", Weeds),
    Yield = case_when(
      grepl("375", Treatment) & is.na(T375) ~ F375,
      grepl("375", Treatment) & is.na(F375) ~ T375,
      !grepl("375", Treatment) & is.na(T) ~ F,
      !grepl("375", Treatment) & is.na(F) ~ T,
      TRUE ~ NA_real_
    )) %>% 
  dplyr::select(Treatment, Yield)
sole_means_2023B_W2 <- sole_means_2023B_W$Yield
names(sole_means_2023B_W2) <- sole_means_2023B_W$Treatment

sole_means_2024_W <- data_2024_yield %>%
  filter(Treatment %in% c("T", "F")) %>%
  group_by(Treatment, Weeds) %>% 
  summarise(
    T = mean(YieldT[Treatment == "T"], na.rm = TRUE),
    F = mean(YieldF[Treatment == "F"], na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  mutate(Treatment = paste0(Treatment, "_", Weeds),
    Yield = ifelse(is.na(T), F, T)) %>% 
  dplyr::select(Treatment, Yield)
sole_means_2024_W2 <- sole_means_2024_W$Yield
names(sole_means_2024_W2) <- sole_means_2024_W$Treatment

data_2023A_yield_LER_W <- data_2023A_yield %>% 
  mutate(Treatment = gsub("TF", "1T:1F", Treatment)) %>%
  filter(Treatment %in% c("3T:1F", "1T:1F", "1T:1F-M", "1T:3F")) %>% 
  pivot_wider(
      id_cols = c(Plot, Treatment, Block),
      names_from = Weeds,
      values_from = c(YieldT, YieldF),
      names_sep = "_"
    ) %>% 
  mutate(Treatment = gsub("TF", "1T:1F", Treatment)) %>%
  mutate(
    pLERT_N = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M", "1T:3F", "3T:1F") ~ YieldT_N / sole_means_2023A_W2["T_N"],
      TRUE ~ NA_real_
    ),
    pLERT_Y = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M", "1T:3F", "3T:1F") ~ YieldT_Y / sole_means_2023A_W2["T_Y"],
      TRUE ~ NA_real_
    ),
    pLERF_N = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M", "1T:3F", "3T:1F") ~ YieldF_N / sole_means_2023A_W2["F_N"],
      TRUE ~ NA_real_
    ),
    pLERF_Y = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M", "1T:3F", "3T:1F") ~ YieldF_Y / sole_means_2023A_W2["F_Y"],
      TRUE ~ NA_real_
    ),
    NERT_N = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M") ~ YieldT_N / (0.5 * sole_means_2023A_W2["T_N"]),
      Treatment %in% c("1T:3F") ~ YieldT_N / (0.25 * sole_means_2023A_W2["T_N"]),
      Treatment %in% c("3T:1F") ~ YieldT_N / (0.75 * sole_means_2023A_W2["T_N"]),
      TRUE ~ NA_real_
    ),
    NERT_Y = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M") ~ YieldT_Y / (0.5 * sole_means_2023A_W2["T_Y"]),
      Treatment %in% c("1T:3F") ~ YieldT_Y / (0.25 * sole_means_2023A_W2["T_Y"]),
      Treatment %in% c("3T:1F") ~ YieldT_Y / (0.75 * sole_means_2023A_W2["T_Y"]),
      TRUE ~ NA_real_
    ),
    NERF_N = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M") ~ YieldF_N / (0.5 * sole_means_2023A_W2["F_N"]),
      Treatment %in% c("1T:3F") ~ YieldF_N / (0.75 * sole_means_2023A_W2["F_N"]),
      Treatment %in% c("3T:1F") ~ YieldF_N / (0.25 * sole_means_2023A_W2["F_N"]),
      TRUE ~ NA_real_
    ),
    NERF_Y = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M") ~ YieldF_Y / (0.5 * sole_means_2023A_W2["F_Y"]),
      Treatment %in% c("1T:3F") ~ YieldF_Y / (0.75 * sole_means_2023A_W2["F_Y"]),
      Treatment %in% c("3T:1F") ~ YieldF_Y / (0.25 * sole_means_2023A_W2["F_Y"]),
      TRUE ~ NA_real_
    ),
    LER_N = pLERT_N + pLERF_N,
    LER_Y = pLERT_Y + pLERF_Y,
    CR_N = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M") ~ (pLERT_N / 0.5) / (pLERF_N / 0.5),
      Treatment == "1T:3F" ~ (pLERT_N / 0.25) / (pLERF_N / 0.75),
      Treatment == "3T:1F" ~ (pLERT_N / 0.75) / (pLERF_N / 0.25),
      TRUE ~ NA_real_
    ),
    CR_Y = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M") ~ (pLERT_Y / 0.5) / (pLERF_Y / 0.5),
      Treatment == "1T:3F" ~ (pLERT_Y / 0.25) / (pLERF_Y / 0.75),
      Treatment == "3T:1F" ~ (pLERT_Y / 0.75) / (pLERF_Y / 0.25),
      TRUE ~ NA_real_
    ),
    RY_C = YieldT_Y / YieldT_N,
    RY_L = YieldF_Y / YieldF_N,
    RY_ICW_C = YieldT_Y / sole_means_2023A_W2["T_N"] * 2,
    RY_ICW_L = YieldF_Y / sole_means_2023A_W2["F_N"] * 2
  )

data_2023B_yield_LER_W <- data_2023B_yield %>%
  mutate(Treatment = gsub("TF", "1T:1F", Treatment)) %>%
  filter(Treatment %in% c("1T:1F", "1T:1F-M", "1T:1F-375", "1T:1F-M-375")) %>% 
  pivot_wider(
      id_cols = c(Plot, Treatment, Block),
      names_from = Weeds,
      values_from = c(YieldT, YieldF),
      names_sep = "_"
    ) %>% 
  mutate(
    pLERT_N = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M") ~ YieldT_N / sole_means_2023B_W2["T_N"],
      Treatment %in% c("1T:1F-375", "1T:1F-M-375") ~ YieldT_N / sole_means_2023B_W2["T-375_N"],
      TRUE ~ NA_real_
    ),
    pLERT_Y = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M") ~ YieldT_Y / sole_means_2023B_W2["T_Y"],
      Treatment %in% c("1T:1F-375", "1T:1F-M-375") ~ YieldT_Y / sole_means_2023B_W2["T-375_Y"],
      TRUE ~ NA_real_
    ),
    pLERF_N = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M") ~ YieldF_N / sole_means_2023B_W2["F_N"],
      Treatment %in% c("1T:1F-375", "1T:1F-M-375") ~ YieldF_N / sole_means_2023B_W2["F-375_N"],
      TRUE ~ NA_real_
    ),
    pLERF_Y = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M") ~ YieldF_Y / sole_means_2023B_W2["F_Y"],
      Treatment %in% c("1T:1F-375", "1T:1F-M-375") ~ YieldF_Y / sole_means_2023B_W2["F-375_Y"],
      TRUE ~ NA_real_
    ),
    NERT_N = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M") ~ YieldT_N / (0.5 * sole_means_2023B_W2["T_N"]),
      Treatment %in% c("1T:1F-375", "1T:1F-M-375") ~ YieldT_N / (0.5 * sole_means_2023B_W2["T-375_N"]),
      TRUE ~ NA_real_
    ),
    NERT_Y = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M") ~ YieldT_Y / (0.5 * sole_means_2023B_W2["T_Y"]),
      Treatment %in% c("1T:1F-375", "1T:1F-M-375") ~ YieldT_Y / (0.5 * sole_means_2023B_W2["T-375_Y"]),
      TRUE ~ NA_real_
    ),
    NERF_N = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M") ~ YieldF_N / (0.5 * sole_means_2023B_W2["F_N"]),
      Treatment %in% c("1T:1F-375", "1T:1F-M-375") ~ YieldF_N / (0.5 * sole_means_2023B_W2["F-375_N"]),
      TRUE ~ NA_real_
    ),
    NERF_Y = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M") ~ YieldF_Y / (0.5 * sole_means_2023B_W2["F_Y"]),
      Treatment %in% c("1T:1F-375", "1T:1F-M-375") ~ YieldF_Y / (0.5 * sole_means_2023B_W2["F-375_Y"]),
      TRUE ~ NA_real_
    ),
    LER_N = pLERT_N + pLERF_N,
    LER_Y = pLERT_Y + pLERF_Y,
    CR_N = (pLERT_N / 0.5) / (pLERF_N / 0.5),
    CR_Y = (pLERT_Y / 0.5) / (pLERF_Y / 0.5),
    RY_C = case_when(
      Treatment %in% c("1T:1F", "1T:1F-M") ~ YieldT_Y / YieldT_N,
      Treatment %in% c("1T:1F-375", "1T:1F-M-375") ~ YieldT_Y / YieldT_N,
      TRUE ~ NA_real_
    ),
    RY_L = YieldF_Y / YieldF_N,
    RY_ICW_C = YieldT_Y / sole_means_2023B_W2["T_N"] * 2,
    RY_ICW_L = YieldF_Y / sole_means_2023B_W2["F_N"] * 2
  )

data_2024_yield_LER_W <- data_2024_yield %>%
  mutate(Treatment = gsub("TF", "1T:1F", Treatment)) %>%
  filter(Treatment %in% c("1T:1F")) %>% 
  pivot_wider(
      id_cols = c(Plot, Treatment, Block),
      names_from = Weeds,
      values_from = c(YieldT, YieldF),
      names_sep = "_"
    ) %>% 
  mutate(
    pLERT_N = case_when(
      Treatment %in% c("1T:1F") ~ YieldT_N / sole_means_2024_W2["T_N"],
      TRUE ~ NA_real_
    ),
    pLERT_Y = case_when(
      Treatment %in% c("1T:1F") ~ YieldT_Y / sole_means_2024_W2["T_Y"],
      TRUE ~ NA_real_
    ),
    pLERF_N = case_when(
      Treatment %in% c("1T:1F") ~ YieldF_N / sole_means_2024_W2["F_N"],
      TRUE ~ NA_real_
    ),
    pLERF_Y = case_when(
      Treatment %in% c("1T:1F") ~ YieldF_Y / sole_means_2024_W2["F_Y"],
      TRUE ~ NA_real_
    ),
    NERT_N = case_when(
      Treatment %in% c("1T:1F") ~ YieldT_N / (0.5 * sole_means_2024_W2["T_N"]),
      TRUE ~ NA_real_
    ),
    NERT_Y = case_when(
      Treatment %in% c("1T:1F") ~ YieldT_Y / (0.5 * sole_means_2024_W2["T_Y"]),
      TRUE ~ NA_real_
    ),
    NERF_N = case_when(
      Treatment %in% c("1T:1F") ~ YieldF_N / (0.5 * sole_means_2024_W2["F_N"]),
      TRUE ~ NA_real_
    ),
    NERF_Y = case_when(
      Treatment %in% c("1T:1F") ~ YieldF_Y / (0.5 * sole_means_2024_W2["F_Y"]),
      TRUE ~ NA_real_
    ),
    LER_N = pLERT_N + pLERF_N,
    LER_Y = pLERT_Y + pLERF_Y,
    CR_N = (pLERT_N / 0.5) / (pLERF_N / 0.5),
    CR_Y = (pLERT_Y / 0.5) / (pLERF_Y / 0.5),
    RY_C = YieldT_Y / YieldT_N,
    RY_L = YieldF_Y / YieldF_N,
    RY_ICW_C = YieldT_Y / sole_means_2024_W2["T_N"] * 2,
    RY_ICW_L = YieldF_Y / sole_means_2024_W2["F_N"] * 2
  )

data_2023A_yield_LER_agg_W <- data_2023A_yield_LER_W %>% 
  group_split(Treatment) %>% 
  map_dfr(function(df) {
    treat <- unique(df$Treatment)
    list(
      Treatment = treat,
      mPLERT_N = mean(df$pLERT_N, na.rm = TRUE),
      mPLERT_Y = mean(df$pLERT_Y, na.rm = TRUE),
      mPLERF_N = mean(df$pLERF_N, na.rm = TRUE),
      mPLERF_Y = mean(df$pLERF_Y, na.rm = TRUE),
      mNERT_N = mean(df$NERT_N, na.rm = TRUE),
      mNERT_Y = mean(df$NERT_Y, na.rm = TRUE),
      mNERF_N = mean(df$NERF_N, na.rm = TRUE),
      mNERF_Y = mean(df$NERF_Y, na.rm = TRUE),
      mLER_N = mean(df$LER_N, na.rm = TRUE),
      mLER_Y = mean(df$LER_Y, na.rm = TRUE),
      mCR_N = mean(df$CR_N, na.rm = TRUE),
      mCR_Y = mean(df$CR_Y, na.rm = TRUE),
      mRY_C = mean(df$RY_C, na.rm = TRUE),
      mRY_L = mean(df$RY_L, na.rm = TRUE),
      mRY_ICW_C = mean(df$RY_ICW_C, na.rm = TRUE),
      mRY_ICW_L = mean(df$RY_ICW_L, na.rm = TRUE),
      P_PLERC_N = tryCatch(
        t.test(df$pLERT_N, mu = ifelse(treat == "1T:3F", 0.25, 0.5))$p.value,
        error = function(e) NA_real_
      ),
      P_PLERC_Y = tryCatch(
        t.test(df$pLERT_Y, mu = ifelse(treat == "1T:3F", 0.25, 0.5))$p.value,
        error = function(e) NA_real_
      ),
      P_PLERL_N = tryCatch(
        t.test(df$pLERF_N, mu = ifelse(treat == "3T:1F", 0.25, 0.5))$p.value,
        error = function(e) NA_real_
      ),
      P_PLERL_Y = tryCatch(
        t.test(df$pLERF_Y, mu = ifelse(treat == "3T:1F", 0.25, 0.5))$p.value,
        error = function(e) NA_real_
      ),
      P_NERC_N = tryCatch(
        t.test(df$NERT_N, mu = 1.0)$p.value,
        error = function(e) NA_real_
      ),
      P_NERC_Y = tryCatch(
        t.test(df$NERT_Y, mu = 1.0)$p.value,
        error = function(e) NA_real_
      ),
      P_NERL_N = tryCatch(
        t.test(df$NERF_N, mu = 1.0)$p.value,
        error = function(e) NA_real_
      ),
      P_NERL_Y = tryCatch(
        t.test(df$NERF_Y, mu = 1.0)$p.value,
        error = function(e) NA_real_
      ),
      P_LER_N = tryCatch(
        t.test(df$LER_N, mu = 1.0)$p.value,
        error = function(e) NA_real_
      ),
      P_LER_Y = tryCatch(
        t.test(df$LER_Y, mu = 1.0)$p.value,
        error = function(e) NA_real_
      ),
      P_CR_N = tryCatch(
        t.test(df$CR_N, mu = 1.0)$p.value,
        error = function(e) NA_real_
      ),
      P_CR_Y = tryCatch(
        t.test(df$CR_Y, mu = 1.0)$p.value,
        error = function(e) NA_real_
      ),
      P_val_mRY_C = tryCatch(t.test(df$RY_C, mu = 1.0)$p.value, error = function(e) NA_real_),
      P_val_mRY_L = tryCatch(t.test(df$RY_L, mu = 1.0)$p.value, error = function(e) NA_real_),
      P_val_mRY_ICW_C = tryCatch(t.test(df$RY_ICW_C, mu = 1.0)$p.value, error = function(e) NA_real_),
      P_val_mRY_ICW_L = tryCatch(t.test(df$RY_ICW_L, mu = 1.0)$p.value, error = function(e) NA_real_)
    )
  }) %>% 
  filter(Treatment %in% c("1T:1F", "1T:1F-M", "1T:3F", "3T:1F"))

data_2023B_yield_LER_agg_W <- data_2023B_yield_LER_W %>% 
  group_by(Treatment) %>% 
  summarise(
    mPLERT_N = mean(pLERT_N, na.rm = TRUE),
    mPLERT_Y = mean(pLERT_Y, na.rm = TRUE),
    mPLERF_N = mean(pLERF_N, na.rm = TRUE),
    mPLERF_Y = mean(pLERF_Y, na.rm = TRUE),
    mNERT_N = mean(NERT_N, na.rm = TRUE),
    mNERT_Y = mean(NERT_Y, na.rm = TRUE),
    mNERF_N = mean(NERF_N, na.rm = TRUE),
    mNERF_Y = mean(NERF_Y, na.rm = TRUE),
    mLER_N = mean(LER_N, na.rm = TRUE),
    mLER_Y = mean(LER_Y, na.rm = TRUE),
    mCR_N = mean(CR_N, na.rm = TRUE),
    mCR_Y = mean(CR_Y, na.rm = TRUE),
    mRY_C = mean(RY_C, na.rm = TRUE),
    mRY_L = mean(RY_L, na.rm = TRUE),
    mRY_ICW_C = mean(RY_ICW_C, na.rm = TRUE),
    mRY_ICW_L = mean(RY_ICW_L, na.rm = TRUE),
    P_PLERC_N = tryCatch(t.test(pLERT_N, mu = 0.5)$p.value, error = function(e) NA_real_),
    P_PLERC_Y = tryCatch(t.test(pLERT_Y, mu = 0.5)$p.value, error = function(e) NA_real_),
    P_PLERL_N = tryCatch(t.test(pLERF_N, mu = 0.5)$p.value, error = function(e) NA_real_),
    P_PLERL_Y = tryCatch(t.test(pLERF_Y, mu = 0.5)$p.value, error = function(e) NA_real_),
    P_NERC_N = tryCatch(t.test(NERT_N, mu = 1.0)$p.value, error = function(e) NA_real_),
    P_NERC_Y = tryCatch(t.test(NERT_Y, mu = 1.0)$p.value, error = function(e) NA_real_),
    P_NERL_N = tryCatch(t.test(NERF_N, mu = 1.0)$p.value, error = function(e) NA_real_),
    P_NERL_Y = tryCatch(t.test(NERF_Y, mu = 1.0)$p.value, error = function(e) NA_real_),
    P_LER_N = tryCatch(t.test(LER_N, mu = 1.0)$p.value, error = function(e) NA_real_),
    P_LER_Y = tryCatch(t.test(LER_Y, mu = 1.0)$p.value, error = function(e) NA_real_),
    P_CR_N = tryCatch(t.test(CR_N, mu = 1.0)$p.value, error = function(e) NA_real_),
    P_CR_Y = tryCatch(t.test(CR_Y, mu = 1.0)$p.value, error = function(e) NA_real_),
    P_val_mRY_C = t.test(RY_C, mu = 1.0)$p.value,
    P_val_mRY_L = t.test(RY_L, mu = 1.0)$p.value,
    P_val_mRY_ICW_C = t.test(RY_ICW_C, mu = 1.0)$p.value,
    P_val_mRY_ICW_L = t.test(RY_ICW_L, mu = 1.0)$p.value,
    .groups = "drop"
  ) %>% 
  filter(Treatment %in% c("1T:1F", "1T:1F-M", "1T:1F-375", "1T:1F-M-375"))

data_2024_yield_LER_agg_W <- data_2024_yield_LER_W %>% 
  group_by(Treatment) %>% 
  summarise(
    mPLERT_N = mean(pLERT_N, na.rm = TRUE),
    mPLERT_Y = mean(pLERT_Y, na.rm = TRUE),
    mPLERF_N = mean(pLERF_N, na.rm = TRUE),
    mPLERF_Y = mean(pLERF_Y, na.rm = TRUE),
    mNERT_N = mean(NERT_N, na.rm = TRUE),
    mNERT_Y = mean(NERT_Y, na.rm = TRUE),
    mNERF_N = mean(NERF_N, na.rm = TRUE),
    mNERF_Y = mean(NERF_Y, na.rm = TRUE),
    mLER_N = mean(LER_N, na.rm = TRUE),
    mLER_Y = mean(LER_Y, na.rm = TRUE),
    mCR_N = mean(CR_N, na.rm = TRUE),
    mCR_Y = mean(CR_Y, na.rm = TRUE),
    mRY_C = mean(RY_C, na.rm = TRUE),
    mRY_L = mean(RY_L, na.rm = TRUE),
    mRY_ICW_C = mean(RY_ICW_C, na.rm = TRUE),
    mRY_ICW_L = mean(RY_ICW_L, na.rm = TRUE),
    P_PLERC_N = tryCatch(t.test(pLERT_N, mu = 0.5)$p.value, error = function(e) NA_real_),
    P_PLERC_Y = tryCatch(t.test(pLERT_Y, mu = 0.5)$p.value, error = function(e) NA_real_),
    P_PLERL_N = tryCatch(t.test(pLERF_N, mu = 0.5)$p.value, error = function(e) NA_real_),
    P_PLERL_Y = tryCatch(t.test(pLERF_Y, mu = 0.5)$p.value, error = function(e) NA_real_),
    P_NERC_N = tryCatch(t.test(NERT_N, mu = 0.5)$p.value, error = function(e) NA_real_),
    P_NERC_Y = tryCatch(t.test(NERT_Y, mu = 0.5)$p.value, error = function(e) NA_real_),
    P_NERL_N = tryCatch(t.test(NERF_N, mu = 0.5)$p.value, error = function(e) NA_real_),
    P_NERL_Y = tryCatch(t.test(NERF_Y, mu = 0.5)$p.value, error = function(e) NA_real_),
    P_LER_N = tryCatch(t.test(LER_N, mu = 1.0)$p.value, error = function(e) NA_real_),
    P_LER_Y = tryCatch(t.test(LER_Y, mu = 1.0)$p.value, error = function(e) NA_real_),
    P_CR_N = tryCatch(t.test(CR_N, mu = 1.0)$p.value, error = function(e) NA_real_),
    P_CR_Y = tryCatch(t.test(CR_Y, mu = 1.0)$p.value, error = function(e) NA_real_),
    P_val_mRY_C = t.test(RY_C, mu = 1.0)$p.value,
    P_val_mRY_L = t.test(RY_L, mu = 1.0)$p.value,
    P_val_mRY_ICW_C = t.test(RY_ICW_C, mu = 1.0)$p.value,
    P_val_mRY_ICW_L = t.test(RY_ICW_L, mu = 1.0)$p.value,
    .groups = "drop"
  ) %>% 
  filter(Treatment %in% c("1T:1F"))

data_yield_W_RY_all <- bind_rows(
  mutate(data_2023A_yield_LER_agg_W, Experiment = "2023-SP"),
  mutate(data_2023B_yield_LER_agg_W, Experiment = "2023-RD"),
  mutate(data_2024_yield_LER_agg_W, Experiment = "2024")
) %>%
  dplyr::select(Experiment, Treatment, mNERT_N, mNERT_Y, mRY_ICW_C, mNERF_N, mNERF_Y, mRY_ICW_L, mCR_N, mCR_Y, P_NERC_N, P_NERC_Y, P_val_mRY_ICW_C, P_NERL_N, P_NERL_Y, P_val_mRY_ICW_L, P_CR_N, P_CR_Y)

## CR-WR relationship analysis ##

# Create combined dataset with individual CR and WR values for all experiments
data_CR_WR_combined <- bind_rows(
  # 2022 data
  data_2022_yield_LER %>%
    filter(!is.na(CR)) %>%
    mutate(
      Experiment = "2022",
      cereal_code = str_extract(Treatment, "^[^:]*"),
      cereal_name = case_when(
        cereal_code == "1R" ~ "Rye",
        cereal_code == "1B" ~ "Barley",
        cereal_code == "1T" ~ "Triticale",
        cereal_code == "1W" ~ "Wheat",
        TRUE ~ NA_character_
      ),
      legume_name = case_when(
        str_detect(Treatment, "P") ~ "Pea",
        str_detect(Treatment, "F") ~ "Faba",
        TRUE ~ NA_character_
      ),
      WR = case_when(
        !is.na(legume_name) & !is.na(cereal_name) ~ weed_biomass_2022[legume_name] / weed_biomass_2022[cereal_name],
        TRUE ~ NA_real_
      )
    ) %>%
    dplyr::select(Experiment, Treatment, CR, WR),

  # 2023A data
  data_2023A_yield_LER %>%
    filter(!is.na(CR)) %>%
    mutate(
      Experiment = "2023-SP",
      WR = weed_biomass_2023A["F"] / weed_biomass_2023A["T"]
    ) %>%
    dplyr::select(Experiment, Treatment, CR, WR),

  # 2023B data
  data_2023B_yield_LER %>%
    filter(!is.na(CR)) %>%
    mutate(
      Experiment = "2023-RD",
      WR = case_when(
        Treatment %in% c("1T:1F", "1T:1F-M") ~ weed_biomass_2023B["F"] / weed_biomass_2023B["T"],
        Treatment %in% c("1T:1F-375", "1T:1F-M-375") ~ weed_biomass_2023B["F-375"] / weed_biomass_2023B["T-375"],
        TRUE ~ NA_real_
      )
    ) %>%
    dplyr::select(Experiment, Treatment, CR, WR),

  # 2024 data
  data_2024_yield_LER %>%
    filter(!is.na(CR)) %>%
    mutate(
      Experiment = "2024",
      WR = weed_biomass_2024["F"] / weed_biomass_2024["T"]
    ) %>%
    dplyr::select(Experiment, Treatment, CR, WR)
) %>%
  filter(!is.na(CR) & !is.na(WR) & CR > 0 & WR > 0) %>%
  mutate(
    ln_CR = log(CR),
    ln_WR = log(WR)
  )

# Perform linear regression analysis
lm_CR_WR <- lm(ln_CR ~ ln_WR, data = data_CR_WR_combined)
summary_lm <- summary(lm_CR_WR)
p_value <- summary_lm$coefficients[2, 4]  # p-value for slope
intercept <- round(summary_lm$coefficients[1, 1], 3)
slope <- round(summary_lm$coefficients[2, 1], 3)

# Format p-value
p_text <- ifelse(p_value < 0.0001, "< 0.0001", paste("=", round(p_value, 4)))

# Create the scatter plot
plot_CR_WR <- ggplot(data_CR_WR_combined, aes(x = ln_WR, y = ln_CR)) +
  geom_point(aes(fill = Treatment, shape = Experiment), size = 4, alpha = 0.8, color = "black") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "solid") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  labs(
    x = expression(ln(WR[LC])),
    y = expression(ln(CR[CL])),
    fill = "Treatment",
    shape = "Experiment"
  ) +
  scale_fill_manual(values = c(
    "1R:1P" = "#E41A1C", "1B:1F" = "#377EB8", "1T:1F" = "#4DAF4A", "1B:1P" = "#984EA3",
    "1W:1P" = "#FF7F00", "1T:1P" = "#FFFF33", "1W:1F" = "#A65628", "1R:1F" = "#F781BF",
    "1T:1F-M" = "#999999", "1T:3F" = "#66C2A5", "3T:1F" = "#FC8D62", "1T:1F-375" = "#8DA0CB",
    "1T:1F-M-375" = "#E78AC3"
  )) +
  scale_shape_manual(values = c("2022" = 21, "2023-SP" = 24, "2023-RD" = 22, "2024" = 23)) +
  guides(fill = guide_legend(override.aes = list(shape = 21, color = "black"))) +
  annotate("text",
           x = min(data_CR_WR_combined$ln_WR) + 0.1 * diff(range(data_CR_WR_combined$ln_WR)) - 0.1,
           y = max(data_CR_WR_combined$ln_CR) - 0.1 * diff(range(data_CR_WR_combined$ln_CR)) + 0.15,
           label = paste0("ln(CR[CL]) == ", intercept, " + ", slope, " %*% ln(WR[LC])"),
           hjust = 0, vjust = 1, size = 6, parse = TRUE) +
  annotate("text",
           x = min(data_CR_WR_combined$ln_WR) + 0.1 * diff(range(data_CR_WR_combined$ln_WR)) - 0.1,
           y = max(data_CR_WR_combined$ln_CR) - 0.15 * diff(range(data_CR_WR_combined$ln_CR)) + 0.15,
           label = paste0("P ", p_text),
           hjust = 0, vjust = 1, size = 6) +
  theme_classic(base_size = 20) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

plot_CR_WR



## Save outputs ##

# Output figures
png("Figures/plot_yield_2022.png", units = "px", width = 4850, height = 4350, res = 300)
combined_plot_2022
dev.off()

png("Figures/plot_yield_2023A.png", units = "px", width = 3000, height = 4000, res = 300)
combined_plot_2023A
dev.off()

png("Figures/plot_yield_2023B.png", units = "px", width = 3000, height = 4000, res = 300)
combined_plot_2023B
dev.off()

png("Figures/plot_yield_all_years.png", units = "px", width = 3300, height = 3000, res = 300)
combined_plot_years
dev.off()

png("Figures/plot_CR_WR_relationship.png", units = "px", width = 3000, height = 2500, res = 300)
plot_CR_WR
dev.off()

# Output tables
data_2022_yield_LER_agg
data_2023A_yield_LER_agg
data_2023B_yield_LER_agg
data_yield_LER_all
data_sole_W_all
data_yield_W_RY_all


## Create comprehensive yield components table ##

# Load 2022 data with yield components
data_2022_components <- read_xlsx("WCF_2022_data.xlsx", sheet = "FinalHarvest", range = "A1:P77", col_names = TRUE) %>%
  dplyr::select(Plot, Treatment, Block,
                BiomassSeedsCereal, `1000SeedWeightCereal`, NSpikesCereal, NSeedsCereal,
                BiomassSeedsLegume, `1000SeedWeightLegume`, NPodsLegume, NSeedsLegume) %>%
  filter(!grepl("Lupine", Treatment)) %>%
  mutate(Year = "2022",
         Experiment = "2022") %>%
  rename(YieldC = BiomassSeedsCereal,
         YieldL = BiomassSeedsLegume,
         TSWC = `1000SeedWeightCereal`,
         TSWL = `1000SeedWeightLegume`,
         NFC = NSpikesCereal,
         NFL = NPodsLegume,
         NSC = NSeedsCereal,
         NSL = NSeedsLegume)

# Load 2023A data with yield components
data_2023A_components <- read_xlsx("WCF_2023_data.xlsx", sheet = "FinalHarvestA", range = "A1:AJ81", col_names = TRUE) %>%
  dplyr::select(Plot, Treatment, Block,
                BiomassSeedsTriticale, `1000SeedWeightTriticale`, NSpikesTriticale, NSeedsTriticale,
                BiomassSeedsFaba, `1000SeedWeightFaba`, NPodsFaba, NSeedsFaba) %>%
  mutate(Year = "2023",
         Experiment = "2023-SP") %>%
  rename(YieldC = BiomassSeedsTriticale,
         YieldL = BiomassSeedsFaba,
         TSWC = `1000SeedWeightTriticale`,
         TSWL = `1000SeedWeightFaba`,
         NFC = NSpikesTriticale,
         NFL = NPodsFaba,
         NSC = NSeedsTriticale,
         NSL = NSeedsFaba)

# Load 2023B data with yield components
data_2023B_components <- read_xlsx("WCF_2023_data.xlsx", sheet = "FinalHarvestB", range = "A1:P89", col_names = TRUE) %>%
  dplyr::select(Plot, Treatment, Block,
                BiomassSeedsTriticale, `1000SeedWeightTriticale`, NSpikesTriticale, NSeedsTriticale,
                BiomassSeedsFaba, `1000SeedWeightFaba`, NPodsFaba, NSeedsFaba) %>%
  mutate(Year = "2023",
         Experiment = "2023-RD") %>%
  rename(YieldC = BiomassSeedsTriticale,
         YieldL = BiomassSeedsFaba,
         TSWC = `1000SeedWeightTriticale`,
         TSWL = `1000SeedWeightFaba`,
         NFC = NSpikesTriticale,
         NFL = NPodsFaba,
         NSC = NSeedsTriticale,
         NSL = NSeedsFaba)

# Load 2024 data with yield components
data_2024_components <- read_xlsx("WCF_2024_data.xlsx", sheet = "FinalHarvest", range = "A1:O41", col_names = TRUE) %>%
  dplyr::select(Plot, Treatment, Block,
                BiomassSeedsTriticale, TSWTriticale, NSpikesTriticale, NSeedsTriticale,
                BiomassSeedsFaba, TSWFaba, NPodsFaba, NSeedsFaba) %>%
  mutate(Year = "2024",
         Experiment = "2024") %>%
  rename(YieldC = BiomassSeedsTriticale,
         YieldL = BiomassSeedsFaba,
         TSWC = TSWTriticale,
         TSWL = TSWFaba,
         NFC = NSpikesTriticale,
         NFL = NPodsFaba,
         NSC = NSeedsTriticale,
         NSL = NSeedsFaba)

# Combine all years
yield_components_all <- bind_rows(
  data_2022_components,
  data_2023A_components,
  data_2023B_components,
  data_2024_components
) %>%
  arrange(Year, Experiment, Treatment, Plot)

# Create summary table with means and standard errors by treatment and year
yield_components_summary <- yield_components_all %>%
  group_by(Year, Experiment, Treatment) %>%
  summarise(
    n = n(),
    YieldC_mean = mean(YieldC, na.rm = TRUE),
    YieldC_se = sd(YieldC, na.rm = TRUE) / sqrt(n()),
    YieldL_mean = mean(YieldL, na.rm = TRUE),
    YieldL_se = sd(YieldL, na.rm = TRUE) / sqrt(n()),
    TSWC_mean = mean(TSWC, na.rm = TRUE),
    TSWC_se = sd(TSWC, na.rm = TRUE) / sqrt(n()),
    TSWL_mean = mean(TSWL, na.rm = TRUE),
    TSWL_se = sd(TSWL, na.rm = TRUE) / sqrt(n()),
    NFC_mean = mean(NFC, na.rm = TRUE),
    NFC_se = sd(NFC, na.rm = TRUE) / sqrt(n()),
    NFL_mean = mean(NFL, na.rm = TRUE),
    NFL_se = sd(NFL, na.rm = TRUE) / sqrt(n()),
    NSC_mean = mean(NSC, na.rm = TRUE),
    NSC_se = sd(NSC, na.rm = TRUE) / sqrt(n()),
    NSL_mean = mean(NSL, na.rm = TRUE),
    NSL_se = sd(NSL, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Display tables
cat("\n=== Raw Yield Components Data (all observations) ===\n")
print(yield_components_all, n = Inf)

cat("\n=== Yield Components Summary (means ± SE by treatment and year) ===\n")
print(yield_components_summary, n = Inf)

# Optionally save to CSV
write.csv(yield_components_all, "yield_components_raw.csv", row.names = FALSE)
write.csv(yield_components_summary, "yield_components_summary.csv", row.names = FALSE)
