## Author: David Kottelenberg
## Institute: Wageningen University & Research
## Last edited: 21-10-2025
## Summary: Analyses of 2022-2024 field experiments, light interception

rm(list = ls())

library(tidyverse)
library(scales)
library(readxl)
library(ggpubr)
library(lubridate)
library(emdbook)
library(emmeans)
library(RColorBrewer)
library(xtable)
library(patchwork)

# Set directory
# setwd("")
dir.create("light_interception")

colourVector <- c(brewer.pal(n=8,"Set1"), "#FDC086")

modelColors <- c("Cereal" = colourVector[1],
                 "Legume" = colourVector[2],
                 "Faba" = colourVector[2],
                 "Triticale" = colourVector[1],
                 "Wheat" = colourVector[1],
                 "Wheat_Faba" = colourVector[3],
                 "Weed" = colourVector[4],
                 "T" = colourVector[1],
                 "F" = colourVector[2],
                 "T-375" = colourVector[1],
                 "F-375" = colourVector[2],
                 "TA" = colourVector[5],
                 "TM" = colourVector[9],
                 "FA" = colourVector[7],
                 "FM" =  colourVector[8],
                 "Intercrop" = colourVector[3],
                 "Sole cereal" = colourVector[1],
                 "Sole legume" = colourVector[2],
                 "1T:1F" = colourVector[3],
                 "1T:3F" = colourVector[3],
                 "3T:1F" = colourVector[3],
                 "TF-M" = colourVector[3],
                 "Sole triticale" = colourVector[1],
                 "Sole faba bean" = colourVector[2],
                 "Triticale-faba bean" = colourVector[3],
                 "Triticale_Faba" = colourVector[3])


##################################################################
##################################################################
######################                      ######################
######################  Light interception  ######################
######################                      ######################
##################################################################
##################################################################

######################
##                  ##
## Experiment 2022  ##
##                  ##
######################

# Set sowing date
sowingDate2022 <- as.Date("2022-04-19")

# Load graph generation and model fitting code
source("WCFModelFitting.R")
source("WCFExploratoryGraphs.R")

# Load temperature data
tempData2022 <- read_delim("TempData2022.txt") %>% 
  mutate(Time = as.integer(as.Date(Date) - sowingDate2022))

cumulativeTemp2022 <- ((tempData2022$TMin + tempData2022$TMax) / 2)
cumulativeTemp2022[which(cumulativeTemp2022 < 0.0)] <- 0
cumulativeTemp2022 <- cumsum(cumulativeTemp2022)

tempData2022 <- tempData2022 %>% 
  mutate(CumulTemp = cumulativeTemp2022)

# Read data
data2022Full <- read_xlsx("WCF_2022_data.xlsx", sheet = "LightInterception", range = "A1:G2008", col_names = TRUE) %>% 
  rename("LightInterception" = "PLightInterception") %>% 
  mutate(Time = tempData2022$CumulTemp[match(DAS, tempData2022$Time)]) %>% 
  dplyr::select(Plot, Treatment, TreatmentN, Block, Time, Sample, LightInterception) %>% 
  mutate(LightInterception = as.numeric(LightInterception)) %>% 
  filter(!grepl("Lupine", Treatment))

data2022 <- as_tibble(aggregate(LightInterception ~ Plot + Treatment + TreatmentN + Block + Time, data = data2022Full, FUN = function(x)mean(x, na.rm = TRUE)))
data2022$LightInterception[which(data2022$LightInterception < 0.001)] <- 0.001
data2022$LightInterception[which(data2022$LightInterception > 0.999)] <- 0.999

plotFit <- function(data, par, fun){
  plot <- ggplot() +
    geom_point(data = data, aes(x = Time, y = LightInterception), size = 3, shape = 1) +
    geom_function(fun = function(x)allDeterministicFunctions[[fun]](par, x),
                  size = 1.5)
  return(plot)
}

## Manual fits

dataTreatment2022 <- split(data2022, data2022$Treatment)
dataTreatment2022 <- lapply(dataTreatment2022, function(dat){
  datNew <- dat %>% 
    dplyr::select(Time, LightInterception)
  return(datNew)
})




##################################
###                            ###
### Fit logistic normal curves ###
###                            ###
##################################     

fitsLogisticNormal2022 <- rep(list(NA), 14)

# Fit of Barley
fitsLogisticNormal2022[[1]] <- fitDataToModel(data = dataTreatment2022[[1]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.01))
plotFit(dataTreatment2022[[1]], fitsLogisticNormal2022[[1]]$par, "logistic")

# Fit of Barley_Faba
fitsLogisticNormal2022[[2]] <- fitDataToModel(data = dataTreatment2022[[2]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[2]], fitsLogisticNormal2022[[2]]$par, "logistic")

# Fit of Barley_Pea
fitsLogisticNormal2022[[3]] <- fitDataToModel(data = dataTreatment2022[[3]], model = logisticNormal, startParameters = c("r" = 0.005, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[3]], fitsLogisticNormal2022[[3]]$par, "logistic")

# Fit of Faba
fitsLogisticNormal2022[[4]] <- fitDataToModel(data = dataTreatment2022[[4]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[4]], fitsLogisticNormal2022[[4]]$par, "logistic")

# Fit of Pea
fitsLogisticNormal2022[[5]] <- fitDataToModel(data = dataTreatment2022[[5]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[5]], fitsLogisticNormal2022[[5]]$par, "logistic")

# Fit of Rye
fitsLogisticNormal2022[[6]] <- fitDataToModel(data = dataTreatment2022[[6]], model = logisticNormal, startParameters = c("r" = 0.012, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[6]], fitsLogisticNormal2022[[6]]$par, "logistic")

# Fit Rye_Faba
fitsLogisticNormal2022[[7]] <- fitDataToModel(data = dataTreatment2022[[7]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[7]], fitsLogisticNormal2022[[7]]$par, "logistic")

# Fit Rye_Pea
fitsLogisticNormal2022[[8]] <- fitDataToModel(data = dataTreatment2022[[8]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[8]], fitsLogisticNormal2022[[8]]$par, "logistic")

# Fit Triticale
fitsLogisticNormal2022[[9]] <- fitDataToModel(data = dataTreatment2022[[9]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[9]], fitsLogisticNormal2022[[9]]$par, "logistic")

# Fit Triticale_Faba
fitsLogisticNormal2022[[10]] <- fitDataToModel(data = dataTreatment2022[[10]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[10]], fitsLogisticNormal2022[[10]]$par, "logistic")

# Fit Triticale_Pea
fitsLogisticNormal2022[[11]] <- fitDataToModel(data = dataTreatment2022[[11]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[11]], fitsLogisticNormal2022[[11]]$par, "logistic")

# Fit Wheat
fitsLogisticNormal2022[[12]] <- fitDataToModel(data = dataTreatment2022[[12]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[12]], fitsLogisticNormal2022[[12]]$par, "logistic")

# Fit Wheat_Faba
fitsLogisticNormal2022[[13]] <- fitDataToModel(data = dataTreatment2022[[13]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[13]], fitsLogisticNormal2022[[13]]$par, "logistic")

# Fit Wheat_Pea
fitsLogisticNormal2022[[14]] <- fitDataToModel(data = dataTreatment2022[[14]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[14]], fitsLogisticNormal2022[[14]]$par, "logistic")

AICLogisticNormal2022 <- sum(unlist(lapply(fitsLogisticNormal2022, function(fits)2 * 3 + 2 * fits$value)))


##############################################
###                                        ###
### Fit transformed logistic normal curves ###
###                                        ###
##############################################      

fitsTLogisticNormal2022 <- rep(list(NA), 14)

# Fit of Barley
fitsTLogisticNormal2022[[1]] <- fitDataToModel(data = dataTreatment2022[[1]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 300, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[1]], fitsTLogisticNormal2022[[1]]$par, "logistic")

# Fit of Barley_Faba
fitsTLogisticNormal2022[[2]] <- fitDataToModel(data = dataTreatment2022[[2]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[2]], fitsTLogisticNormal2022[[2]]$par, "logistic")

# Fit of Barley_Pea
fitsTLogisticNormal2022[[3]] <- fitDataToModel(data = dataTreatment2022[[3]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[3]], fitsTLogisticNormal2022[[3]]$par, "logistic")

# Fit of Faba
fitsTLogisticNormal2022[[4]] <- fitDataToModel(data = dataTreatment2022[[4]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 300, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[4]], fitsTLogisticNormal2022[[4]]$par, "logistic")

# Fit of Pea
fitsTLogisticNormal2022[[5]] <- fitDataToModel(data = dataTreatment2022[[5]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 300, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[5]], fitsTLogisticNormal2022[[5]]$par, "logistic")

# Fit of Rye
fitsTLogisticNormal2022[[6]] <- fitDataToModel(data = dataTreatment2022[[6]], model = tLogisticNormal, startParameters = c("r" = 0.03, "h" = 280, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[6]], fitsTLogisticNormal2022[[6]]$par, "logistic")

# Fit Rye_Faba
fitsTLogisticNormal2022[[7]] <- fitDataToModel(data = dataTreatment2022[[7]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 300, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[7]], fitsTLogisticNormal2022[[7]]$par, "logistic")

# Fit Rye_Pea
fitsTLogisticNormal2022[[8]] <- fitDataToModel(data = dataTreatment2022[[8]], model = tLogisticNormal, startParameters = c("r" = 0.02, "h" = 280, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[8]], fitsTLogisticNormal2022[[8]]$par, "logistic")

# Fit Triticale
fitsTLogisticNormal2022[[9]] <- fitDataToModel(data = dataTreatment2022[[9]], model = tLogisticNormal, startParameters = c("r" = 0.02, "h" = 280, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[9]], fitsTLogisticNormal2022[[9]]$par, "logistic")

# Fit Triticale_Faba
fitsTLogisticNormal2022[[10]] <- fitDataToModel(data = dataTreatment2022[[10]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 300, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[10]], fitsTLogisticNormal2022[[10]]$par, "logistic")

# Fit Triticale_Pea
fitsTLogisticNormal2022[[11]] <- fitDataToModel(data = dataTreatment2022[[11]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 300, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[11]], fitsTLogisticNormal2022[[11]]$par, "logistic")

# Fit Wheat
fitsTLogisticNormal2022[[12]] <- fitDataToModel(data = dataTreatment2022[[12]], model = tLogisticNormal, startParameters = c("r" = 0.02, "h" = 280, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[12]], fitsTLogisticNormal2022[[12]]$par, "logistic")

# Fit Wheat_Faba
fitsTLogisticNormal2022[[13]] <- fitDataToModel(data = dataTreatment2022[[13]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 300, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[13]], fitsTLogisticNormal2022[[13]]$par, "logistic")

# Fit Wheat_Pea
fitsTLogisticNormal2022[[14]] <- fitDataToModel(data = dataTreatment2022[[14]], model = tLogisticNormal, startParameters = c("r" = 0.02, "h" = 280, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[14]], fitsTLogisticNormal2022[[14]]$par, "logistic")

AICTLogisticNormal2022 <- sum(unlist(lapply(fitsTLogisticNormal2022, function(fits)2 * 4 + 2 * fits$value)))

#################################
###                           ###
### Fit logistic gamma curves ###
###                           ###
#################################     

fitsLogisticGamma2022 <- rep(list(NA), 14)

# Fit of Barley
fitsLogisticGamma2022[[1]] <- fitDataToModel(data = dataTreatment2022[[1]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0))
plotFit(dataTreatment2022[[1]], fitsLogisticGamma2022[[1]]$par, "logistic")

# Fit of Barley_Faba
fitsLogisticGamma2022[[2]] <- fitDataToModel(data = dataTreatment2022[[2]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 350, "shape" = 1.0))
plotFit(dataTreatment2022[[2]], fitsLogisticGamma2022[[2]]$par, "logistic")

# Fit of Barley_Pea
fitsLogisticGamma2022[[3]] <- fitDataToModel(data = dataTreatment2022[[3]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 350, "shape" = 1.0))
plotFit(dataTreatment2022[[3]], fitsLogisticGamma2022[[3]]$par, "logistic")

# Fit of Faba
fitsLogisticGamma2022[[4]] <- fitDataToModel(data = dataTreatment2022[[4]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0))
plotFit(dataTreatment2022[[4]], fitsLogisticGamma2022[[4]]$par, "logistic")

# Fit of Pea
fitsLogisticGamma2022[[5]] <- fitDataToModel(data = dataTreatment2022[[5]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0))
plotFit(dataTreatment2022[[5]], fitsLogisticGamma2022[[5]]$par, "logistic")

# Fit of Rye
fitsLogisticGamma2022[[6]] <- fitDataToModel(data = dataTreatment2022[[6]], model = logisticGamma, startParameters = c("r" = 0.03, "h" = 280, "shape" = 1.0))
plotFit(dataTreatment2022[[6]], fitsLogisticGamma2022[[6]]$par, "logistic")

# Fit Rye_Faba
fitsLogisticGamma2022[[7]] <- fitDataToModel(data = dataTreatment2022[[7]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0))
plotFit(dataTreatment2022[[7]], fitsLogisticGamma2022[[7]]$par, "logistic")

# Fit Rye_Pea
fitsLogisticGamma2022[[8]] <- fitDataToModel(data = dataTreatment2022[[8]], model = logisticGamma, startParameters = c("r" = 0.02, "h" = 280, "shape" = 1.0))
plotFit(dataTreatment2022[[8]], fitsLogisticGamma2022[[8]]$par, "logistic")

# Fit Triticale
fitsLogisticGamma2022[[9]] <- fitDataToModel(data = dataTreatment2022[[9]], model = logisticGamma, startParameters = c("r" = 0.02, "h" = 280, "shape" = 1.0))
plotFit(dataTreatment2022[[9]], fitsLogisticGamma2022[[9]]$par, "logistic")

# Fit Triticale_Faba
fitsLogisticGamma2022[[10]] <- fitDataToModel(data = dataTreatment2022[[10]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0))
plotFit(dataTreatment2022[[10]], fitsLogisticGamma2022[[10]]$par, "logistic")

# Fit Triticale_Pea
fitsLogisticGamma2022[[11]] <- fitDataToModel(data = dataTreatment2022[[11]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0))
plotFit(dataTreatment2022[[11]], fitsLogisticGamma2022[[11]]$par, "logistic")

# Fit Wheat
fitsLogisticGamma2022[[12]] <- fitDataToModel(data = dataTreatment2022[[12]], model = logisticGamma, startParameters = c("r" = 0.02, "h" = 280, "shape" = 1.0))
plotFit(dataTreatment2022[[12]], fitsLogisticGamma2022[[12]]$par, "logistic")

# Fit Wheat_Faba
fitsLogisticGamma2022[[13]] <- fitDataToModel(data = dataTreatment2022[[13]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0))
plotFit(dataTreatment2022[[13]], fitsLogisticGamma2022[[13]]$par, "logistic")

# Fit Wheat_Pea
fitsLogisticGamma2022[[14]] <- fitDataToModel(data = dataTreatment2022[[14]], model = logisticGamma, startParameters = c("r" = 0.02, "h" = 280, "shape" = 1.0))
plotFit(dataTreatment2022[[14]], fitsLogisticGamma2022[[14]]$par, "logistic")

AICLogisticGamma2022 <- sum(unlist(lapply(fitsLogisticGamma2022, function(fits)2 * 3 + 2 * fits$value)))


#############################################
###                                       ###
### Fit transformed logistic gamma curves ###
###                                       ###
#############################################      

fitsTLogisticGamma2022 <- rep(list(NA), 14)

# Fit of Barley
fitsTLogisticGamma2022[[1]] <- fitDataToModel(data = dataTreatment2022[[1]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[1]], fitsTLogisticGamma2022[[1]]$par, "logistic")

# Fit of Barley_Faba
fitsTLogisticGamma2022[[2]] <- fitDataToModel(data = dataTreatment2022[[2]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 350, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[2]], fitsTLogisticGamma2022[[2]]$par, "logistic")

# Fit of Barley_Pea
fitsTLogisticGamma2022[[3]] <- fitDataToModel(data = dataTreatment2022[[3]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 350, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[3]], fitsTLogisticGamma2022[[3]]$par, "logistic")

# Fit of Faba
fitsTLogisticGamma2022[[4]] <- fitDataToModel(data = dataTreatment2022[[4]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[4]], fitsTLogisticGamma2022[[4]]$par, "logistic")

# Fit of Pea
fitsTLogisticGamma2022[[5]] <- fitDataToModel(data = dataTreatment2022[[5]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[5]], fitsTLogisticGamma2022[[5]]$par, "logistic")

# Fit of Rye
fitsTLogisticGamma2022[[6]] <- fitDataToModel(data = dataTreatment2022[[6]], model = tLogisticGamma, startParameters = c("r" = 0.03, "h" = 280, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[6]], fitsTLogisticGamma2022[[6]]$par, "logistic")

# Fit Rye_Faba
fitsTLogisticGamma2022[[7]] <- fitDataToModel(data = dataTreatment2022[[7]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[7]], fitsTLogisticGamma2022[[7]]$par, "logistic")

# Fit Rye_Pea
fitsTLogisticGamma2022[[8]] <- fitDataToModel(data = dataTreatment2022[[8]], model = tLogisticGamma, startParameters = c("r" = 0.02, "h" = 280, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[8]], fitsTLogisticGamma2022[[8]]$par, "logistic")

# Fit Triticale
fitsTLogisticGamma2022[[9]] <- fitDataToModel(data = dataTreatment2022[[9]], model = tLogisticGamma, startParameters = c("r" = 0.02, "h" = 250, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[9]], fitsTLogisticGamma2022[[9]]$par, "logistic")

# Fit Triticale_Faba
fitsTLogisticGamma2022[[10]] <- fitDataToModel(data = dataTreatment2022[[10]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[10]], fitsTLogisticGamma2022[[10]]$par, "logistic")

# Fit Triticale_Pea
fitsTLogisticGamma2022[[11]] <- fitDataToModel(data = dataTreatment2022[[11]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[11]], fitsTLogisticGamma2022[[11]]$par, "logistic")

# Fit Wheat
fitsTLogisticGamma2022[[12]] <- fitDataToModel(data = dataTreatment2022[[12]], model = tLogisticGamma, startParameters = c("r" = 0.02, "h" = 280, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[12]], fitsTLogisticGamma2022[[12]]$par, "logistic")

# Fit Wheat_Faba
fitsTLogisticGamma2022[[13]] <- fitDataToModel(data = dataTreatment2022[[13]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[13]], fitsTLogisticGamma2022[[13]]$par, "logistic")

# Fit Wheat_Pea
fitsTLogisticGamma2022[[14]] <- fitDataToModel(data = dataTreatment2022[[14]], model = tLogisticGamma, startParameters = c("r" = 0.02, "h" = 280, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[14]], fitsTLogisticGamma2022[[14]]$par, "logistic")

AICTLogisticGamma2022 <- sum(unlist(lapply(fitsTLogisticGamma2022, function(fits)2 * 4 + 2 * fits$value)))

AICLogisticNormal2022
AICTLogisticNormal2022
AICLogisticGamma2022
AICTLogisticGamma2022

# Transformed logistic normal is the best fitting model

# Calculate 95% CIs
parameters2022 <- names(fitsTLogisticNormal2022[[1]]$par)
parameterCIs2022 <- rep(NA, times = 3 * length(parameters2022))
names(parameterCIs2022) <- unlist(lapply(parameters2022, function(par)c(par, paste(par, "Lower", sep = ""), paste(par, "Upper", sep = ""))))
parameterCIs2022 <- rep(list(parameterCIs2022), length(fitsTLogisticNormal2022))

cI2022 <- function(par, x, z, opt, dMax){
  pars <- par
  pars[names(opt)] <- opt
  nll <- tLogisticNormal(pars, x, z, dMax)
  return(nll)
}

for(i in 1:length(fitsTLogisticNormal2022)){
  parameterVecs <- lapply(fitsTLogisticNormal2022, function(fit){
    parStart <- fit$par / 100
    parEnd <- fit$par * 20
    parVec <- list(seq(parStart[1], parEnd[1], length.out = 400), seq(parStart[2], parEnd[2], length.out = 400), seq(parStart[3], parEnd[3], length.out = 400), seq(parStart[4], parEnd[4], length.out = 400))
    names(parVec) <- names(fit$par)
    return(parVec)
  })
  optPar <- fitsTLogisticNormal2022[[i]]$par
  parameterCIs2022[[i]][names(optPar)] <- optPar
  parameterCIs2022[[i]][unlist((lapply(names(optPar), function(n)c(paste(n, "Lower", sep = ""), paste(n, "Upper", sep = "")))))] <-
    unlist(lapply(1:length(parameterVecs[[i]]), function(j){
      t(tibble(calculate95CI(parStartVec = optPar[-which(names(optPar) == names(fitsTLogisticNormal2022[[i]]$par)[[j]])],
                             optVec = parameterVecs[[i]][j],
                             model = cI2022,
                             optNLL = fitsTLogisticNormal2022[[i]]$value,
                             x = dataTreatment2022[[i]][,1][[1]],
                             z = dataTreatment2022[[i]][,2][[1]],
                             dMax = 1.0)))
    }))
  
}
names(parameterCIs2022) <- names(dataTreatment2022)

# Write to table
CIOut2022 <- tibble(Treatment = names(parameterCIs2022)) %>% 
  bind_cols(bind_rows(parameterCIs2022))
write.table(CIOut2022, "light_interception/light_interception_tLogisticNormal_parameters_2022.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

CIOut2022Table <- CIOut2022 %>% 
  mutate(r = paste0(round(r, 4), " (", round(rLower, 4), "-", round(rUpper, 4), ")"),
         h = paste0(round(h, 2), " (", round(hLower, 2), "-", round(hUpper, 2), ")"),
         sd = paste0(round(sd, 4), " (", round(sdLower, 4), "-", round(sdUpper, 4), ")"),
         d = paste0(round(d, 3), " (", round(dLower, 3), "-", round(dUpper, 3), ")")) %>% 
  dplyr::select(Treatment, r, h, sd, d)

print(xtable(CIOut2022Table, type = "latex"), include.rownames = FALSE)

dataTreatment2022 <- lapply(1:length(dataTreatment2022), function(i){
  dat <- dataTreatment2022[[i]]
  treatment <- names(dataTreatment2022)[i]
  dat <- dat %>% 
    mutate(Treatment = treatment)
})

dataTreatment2022 <- lapply(dataTreatment2022, function(Data){
  Average <- aggregate(Data$LightInterception, list(Data$Time), function(x)mean(x, na.rm = TRUE))
  Data$Average <- Average$x[match(Data$Time, Average$Group.1)]
  return(Data)
})

# Plot triple combination graphs
geomPoints2022 <- lapply(dataTreatment2022, function(d){
  return(geom_point(data = d, aes(x = Time, y = Average, col = as.factor(Treatment[1])), size = 3, shape = 1))
})

geomFunctions2022 <- lapply(1:length(fitsTLogisticNormal2022), function(i){
  o <- fitsTLogisticNormal2022[[i]]
  return(geom_function(fun = function(x)o$par["d"] * (exp(o$par["r"]*(x - o$par["h"])))/(1 + exp(o$par["r"]*(x - o$par["h"]))),
                       aes(col = as.factor(dataTreatment2022[[i]]$Treatment[1])),
                       size = 1.5))
})

### Get all triple combinations
treatmentsToCompare <- unique(data2022$Treatment)

triples <- unlist(lapply(c("Barley", "Rye","Triticale", "Wheat"), function(c){
  lapply(c("Pea", "Faba"), function(l){
    return(c(paste(c, l, sep = "_"), c, l))
  })
}), recursive = FALSE)

names(dataTreatment2022) <- lapply(dataTreatment2022, function(l)l$Treatment[1])
names(geomPoints2022) <- names(dataTreatment2022)
names(geomFunctions2022) <- names(dataTreatment2022)

triplePlots2022 <- lapply(triples, function(triple){
  plotMultipleGraphs(geomPoints2022[triple], geomFunctions2022[triple], triple, colourVector[c(3, 1, 2)], XLab = bquote(T[c]), YLab = bquote(p[li]), Title = "")
})

# Create a combined plot for intercrops and normal density sole crops
combinedPlot <- plotMultipleGraphs(geomPoints2022[treatmentsToCompare], geomFunctions2022[treatmentsToCompare], treatmentsToCompare, c(colourVector[1], colourVector[2], colourVector[5], colourVector[4], colourVector[3], colourVector[6]), XLab = bquote(T[c]), YLab = bquote(p[li]), Title = "")

### Combine data for AIC comparisons
combinations2022 <- lapply(triples, function(t)combn(t, 2))


# Fit of Barley-Pea
startParametersBarleyPea <- list(c("r" = 0.008, "h" = 400, "sd" = 0.1, "d" = 0.9),
                                 c("r" = 0.01, "h" = 380, "sd" = 0.1, "d" = 0.8),
                                 c("r" = 0.008, "h" = 380, "sd" = 0.1, "d" = 0.9),
                                 c("r" = 0.009, "h" = 420, "sd" = 0.1, "d" = 0.9),
                                 c("r" = 0.009, "h" = 410, "sd" = 0.1, "d" = 0.9))
combinationAICBarleyPea <- round(fitCombinations(data = dataTreatment2022, triple = triples[[1]], model = tLogisticNormal, startParameters = startParametersBarleyPea), 2)

# Fit of Barley-Faba
startParametersBarleyFaba <- list(c("r" = 0.008, "h" = 400, "sd" = 0.1, "d" = 0.9),
                                 c("r" = 0.01, "h" = 380, "sd" = 0.1, "d" = 0.8),
                                 c("r" = 0.008, "h" = 380, "sd" = 0.1, "d" = 0.9),
                                 c("r" = 0.009, "h" = 420, "sd" = 0.1, "d" = 0.9),
                                 c("r" = 0.009, "h" = 410, "sd" = 0.1, "d" = 0.9))
combinationAICBarleyFaba <- round(fitCombinations(data = dataTreatment2022, triple = triples[[2]], model = tLogisticNormal, startParameters = startParametersBarleyFaba), 2)

# Fit of Rye-Pea
startParametersRyePea <- list(c("r" = 0.01, "h" = 280, "sd" = 0.1, "d" = 0.95),
                                  c("r" = 0.008, "h" = 320, "sd" = 0.1, "d" = 0.90),
                                  c("r" = 0.008, "h" = 310, "sd" = 0.1, "d" = 0.9),
                                  c("r" = 0.015, "h" = 320, "sd" = 0.1, "d" = 0.9),
                                  c("r" = 0.01, "h" = 320, "sd" = 0.1, "d" = 0.9))
combinationAICRyePea <- round(fitCombinations(data = dataTreatment2022, triple = triples[[3]], model = tLogisticNormal, startParameters = startParametersRyePea), 2)

# Fit of Rye-Faba
startParametersRyeFaba <- list(c("r" = 0.01, "h" = 280, "sd" = 0.1, "d" = 0.95),
                               c("r" = 0.008, "h" = 320, "sd" = 0.1, "d" = 0.90),
                               c("r" = 0.008, "h" = 380, "sd" = 0.1, "d" = 0.9),
                               c("r" = 0.01, "h" = 320, "sd" = 0.1, "d" = 0.9),
                               c("r" = 0.01, "h" = 320, "sd" = 0.1, "d" = 0.9))
combinationAICRyeFaba <- round(fitCombinations(data = dataTreatment2022, triple = triples[[4]], model = tLogisticNormal, startParameters = startParametersRyeFaba), 2)

# Fit of Triticale-Pea
startParametersTriticalePea <- list(c("r" = 0.01, "h" = 280, "sd" = 0.1, "d" = 0.95),
                               c("r" = 0.008, "h" = 320, "sd" = 0.1, "d" = 0.90),
                               c("r" = 0.008, "h" = 380, "sd" = 0.1, "d" = 0.9),
                               c("r" = 0.01, "h" = 320, "sd" = 0.1, "d" = 0.9),
                               c("r" = 0.01, "h" = 320, "sd" = 0.1, "d" = 0.9))
combinationAICTriticalePea <- round(fitCombinations(data = dataTreatment2022, triple = triples[[5]], model = tLogisticNormal, startParameters = startParametersTriticalePea), 2)

# Fit of Triticale-Faba
startParametersTriticaleFaba <- list(c("r" = 0.01, "h" = 280, "sd" = 0.1, "d" = 0.95),
                               c("r" = 0.008, "h" = 320, "sd" = 0.1, "d" = 0.90),
                               c("r" = 0.008, "h" = 380, "sd" = 0.1, "d" = 0.9),
                               c("r" = 0.01, "h" = 320, "sd" = 0.1, "d" = 0.9),
                               c("r" = 0.01, "h" = 320, "sd" = 0.1, "d" = 0.9))
combinationAICTriticaleFaba <- round(fitCombinations(data = dataTreatment2022, triple = triples[[6]], model = tLogisticNormal, startParameters = startParametersTriticaleFaba), 2)

# Fit of Wheat-Pea
startParametersWheatPea <- list(c("r" = 0.01, "h" = 350, "sd" = 0.1, "d" = 0.95),
                                    c("r" = 0.008, "h" = 320, "sd" = 0.1, "d" = 0.90),
                                    c("r" = 0.008, "h" = 380, "sd" = 0.1, "d" = 0.9),
                                    c("r" = 0.01, "h" = 320, "sd" = 0.1, "d" = 0.9),
                                    c("r" = 0.01, "h" = 320, "sd" = 0.1, "d" = 0.9))
combinationAICWheatPea <- round(fitCombinations(data = dataTreatment2022, triple = triples[[7]], model = tLogisticNormal, startParameters = startParametersWheatPea), 2)

# Fit of Wheat-Faba
startParametersWheatFaba <- list(c("r" = 0.01, "h" = 420, "sd" = 0.1, "d" = 0.95),
                                     c("r" = 0.008, "h" = 320, "sd" = 0.1, "d" = 0.90),
                                     c("r" = 0.008, "h" = 380, "sd" = 0.1, "d" = 0.9),
                                     c("r" = 0.01, "h" = 320, "sd" = 0.1, "d" = 0.9),
                                     c("r" = 0.01, "h" = 320, "sd" = 0.1, "d" = 0.9))
combinationAICWheatFaba <- round(fitCombinations(data = dataTreatment2022, triple = triples[[8]], model = tLogisticNormal, startParameters = startParametersWheatFaba), 2)

combinationAIC2022 <- tibble(Treatment = rep(unique(data2022$Treatment)[1:8], each = 3),
                             Combination = rep(c("IC & C & L", "IC + C & L", "IC + L & C"), times = 8),
                             AIC = c(combinationAICBarleyPea,
                                     combinationAICBarleyFaba,
                                     combinationAICRyePea,
                                     combinationAICRyeFaba,
                                     combinationAICTriticalePea,
                                     combinationAICTriticaleFaba,
                                     combinationAICWheatPea,
                                     combinationAICWheatFaba),
                             Best = c("",
                                      "*",
                                      "",
                                      "",
                                      "*",
                                      "",
                                      "*",
                                      "",
                                      "",
                                      "*",
                                      "",
                                      "",
                                      "*",
                                      "*",
                                      "",
                                      "",
                                      "*",
                                      "",
                                      "*",
                                      "",
                                      "",
                                      "*",
                                      "",
                                      ""))

print(xtable(combinationAIC2022, type = "latex"), include.rownames = FALSE)

########################
##                    ##
## Experiment 2023 A  ##
##                    ##
########################

sowingDate2023 <- as.Date("2023-03-02")

# Read data
data2023AFull <- read_xlsx("WCF_2023_data.xlsx", sheet = "LightInterceptionA", range = "A1:H3201", col_names = TRUE) %>% 
  rename(Time = DAS, Val = PLightInterception) %>% 
  filter(Time != 90) %>% 
  dplyr::select(Plot, Treatment, Weeds, Sample, Time, Val)

data2023A <- as_tibble(aggregate(Val ~ Plot + Treatment + Time + Weeds, data = data2023AFull, FUN = function(x)mean(x, na.rm = TRUE)))
data2023A$Val[which(data2023A$Val < 0.001)] <- 0.001
data2023A$Val[which(data2023A$Val > 0.999)] <- 0.999

tempData2023 <- read_delim("TempData2023.txt", delim = "\t") %>% 
  mutate(Time = as.integer(as.Date(Date) - sowingDate2023))

cumulativeTemp2023 <- ((tempData2023$TMin + tempData2023$TMax) / 2)
cumulativeTemp2023[which(cumulativeTemp2023 < 0.0)] <- 0
cumulativeTemp2023 <- cumsum(cumulativeTemp2023)

tempData2023 <- tempData2023 %>% 
  mutate(cumulativeTemp = cumulativeTemp2023)

dataWeeds2023A <- data2023A %>% 
  filter(Weeds == "Y") %>% 
  dplyr::select(-Weeds) %>% 
  filter(Time != 91) %>% 
  mutate(Time = tempData2023$cumulativeTemp[match(Time, tempData2023$Time)])
dataNoWeeds2023A <- data2023A %>% 
  filter(Weeds == "N") %>% 
  dplyr::select(-Weeds) %>% 
  filter(Time != 91) %>% 
  mutate(Time = tempData2023$cumulativeTemp[match(Time, tempData2023$Time)])

## Manual fits
plotFit <- function(data, par, fun){
  plot <- ggplot() +
    geom_point(data = data, aes(x = Time, y = Val), size = 3, shape = 1) +
    geom_function(fun = function(x)allDeterministicFunctions[[fun]](par, x),
                  size = 1)
  return(plot)
}


##################
##################
####          ####
#### No weeds ####
####          ####
##################
##################

dataTreatment2023ANW <- split(dataNoWeeds2023A, dataNoWeeds2023A$Treatment)
dataTreatment2023ANW <- lapply(dataTreatment2023ANW, function(dat){
  datNew <- dat %>% 
    dplyr::select(Time, Val)
  return(datNew)
})
dataTreatment2023ANW$`F` <- dataTreatment2023ANW$`F` %>% 
  filter(Time != 801.65)




##################################
###                            ###
### Fit logistic normal curves ###
###                            ###
##################################     

fitsNWLogisticNormal2023A <- rep(list(NA), 8)

# Fit of 1T:1F
fitsNWLogisticNormal2023A[[1]] <- fitDataToModel(data = dataTreatment2023ANW[[1]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 520, "sd" = 0.1))
plotFit(dataTreatment2023ANW[[1]], fitsNWLogisticNormal2023A[[1]]$par, "logistic")

# Fit of 1T:3F
fitsNWLogisticNormal2023A[[2]] <- fitDataToModel(data = dataTreatment2023ANW[[2]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 520, "sd" = 0.1))
plotFit(dataTreatment2023ANW[[2]], fitsNWLogisticNormal2023A[[2]]$par, "logistic")

# Fit of 3T:1F
fitsNWLogisticNormal2023A[[3]] <- fitDataToModel(data = dataTreatment2023ANW[[3]], model = logisticNormal, startParameters = c("r" = 0.012, "h" = 475, "sd" = 0.1))
plotFit(dataTreatment2023ANW[[3]], fitsNWLogisticNormal2023A[[3]]$par, "logistic")

# Fit of F
fitsNWLogisticNormal2023A[[4]] <- fitDataToModel(data = dataTreatment2023ANW[[4]], model = logisticNormal, startParameters = c("r" = 0.005, "h" = 650, "sd" = 0.1))
plotFit(dataTreatment2023ANW[[4]], fitsNWLogisticNormal2023A[[4]]$par, "logistic")

# Fit of F+
fitsNWLogisticNormal2023A[[5]] <- fitDataToModel(data = dataTreatment2023ANW[[5]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 575, "sd" = 0.1))
plotFit(dataTreatment2023ANW[[5]], fitsNWLogisticNormal2023A[[5]]$par, "logistic")

# Fit of T
fitsNWLogisticNormal2023A[[6]] <- fitDataToModel(data = dataTreatment2023ANW[[6]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 400, "sd" = 0.1))
plotFit(dataTreatment2023ANW[[6]], fitsNWLogisticNormal2023A[[6]]$par, "logistic")

# Fit of T+
fitsNWLogisticNormal2023A[[7]] <- fitDataToModel(data = dataTreatment2023ANW[[7]], model = logisticNormal, startParameters = c("r" = 0.012, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2023ANW[[7]], fitsNWLogisticNormal2023A[[7]]$par, "logistic")

# Fit TF-M
fitsNWLogisticNormal2023A[[8]] <- fitDataToModel(data = dataTreatment2023ANW[[8]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 520, "sd" = 0.1))
plotFit(dataTreatment2023ANW[[8]], fitsNWLogisticNormal2023A[[8]]$par, "logistic")

AICLogisticNormal2023ANW <- sum(unlist(lapply(fitsNWLogisticNormal2023A, function(fits)2 * 3 + 2 * fits$value)))


##############################################
###                                        ###
### Fit transformed logistic normal curves ###
###                                        ###
##############################################      

fitsNWTLogisticNormal2023A <- rep(list(NA), 8)

# Fit of 1T:1F
fitsNWTLogisticNormal2023A[[1]] <- fitDataToModel(data = dataTreatment2023ANW[[1]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 480, "sd" = 0.1, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2023ANW[[1]], fitsNWTLogisticNormal2023A[[1]]$par, "tLogistic")

# Fit of 1T:3F
fitsNWTLogisticNormal2023A[[2]] <- fitDataToModel(data = dataTreatment2023ANW[[2]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 550, "sd" = 0.1, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2023ANW[[2]], fitsNWTLogisticNormal2023A[[2]]$par, "tLogistic")

# Fit of 3T:1F
fitsNWTLogisticNormal2023A[[3]] <- fitDataToModel(data = dataTreatment2023ANW[[3]], model = tLogisticNormal, startParameters = c("r" = 0.012, "h" = 475, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023ANW[[3]], fitsNWTLogisticNormal2023A[[3]]$par, "tLogistic")

# Fit of F
fitsNWTLogisticNormal2023A[[4]] <- fitDataToModel(data = dataTreatment2023ANW[[4]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 640, "sd" = 0.1, "d" = 0.95), dMax = 1.0)
plotFit(dataTreatment2023ANW[[4]], fitsNWTLogisticNormal2023A[[4]]$par, "tLogistic")

# Fit of F+
fitsNWTLogisticNormal2023A[[5]] <- fitDataToModel(data = dataTreatment2023ANW[[5]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 575, "sd" = 0.1, "d" = 0.95), dMax = 1.0)
plotFit(dataTreatment2023ANW[[5]], fitsNWTLogisticNormal2023A[[5]]$par, "tLogistic")

# Fit of T
fitsNWTLogisticNormal2023A[[6]] <- fitDataToModel(data = dataTreatment2023ANW[[6]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 400, "sd" = 0.1, "d" = 0.90), dMax = 1.0)
plotFit(dataTreatment2023ANW[[6]], fitsNWTLogisticNormal2023A[[6]]$par, "tLogistic")

# Fit of T+
fitsNWTLogisticNormal2023A[[7]] <- fitDataToModel(data = dataTreatment2023ANW[[7]], model = tLogisticNormal, startParameters = c("r" = 0.015, "h" = 450, "sd" = 0.1, "d" = 0.95), dMax = 1.0)
plotFit(dataTreatment2023ANW[[7]], fitsNWTLogisticNormal2023A[[7]]$par, "tLogistic")

# Fit TF-M
fitsNWTLogisticNormal2023A[[8]] <- fitDataToModel(data = dataTreatment2023ANW[[8]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 520, "sd" = 0.1, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2023ANW[[8]], fitsNWTLogisticNormal2023A[[8]]$par, "tLogistic")

AICTLogisticNormal2023ANW <- sum(unlist(lapply(fitsNWTLogisticNormal2023A, function(fits)2 * 4 + 2 * fits$value)))


#################################
###                           ###
### Fit logistic gamma curves ###
###                           ###
#################################     

fitsNWLogisticGamma2023A <- rep(list(NA), 8)

# Fit of 1T:1F
fitsNWLogisticGamma2023A[[1]] <- fitDataToModel(data = dataTreatment2023ANW[[1]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 520, "shape" = 1.0))
plotFit(dataTreatment2023ANW[[1]], fitsNWLogisticGamma2023A[[1]]$par, "logistic")

# Fit of 1T:3F
fitsNWLogisticGamma2023A[[2]] <- fitDataToModel(data = dataTreatment2023ANW[[2]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 520, "shape" = 1.0))
plotFit(dataTreatment2023ANW[[2]], fitsNWLogisticGamma2023A[[2]]$par, "logistic")

# Fit of 3T:1F
fitsNWLogisticGamma2023A[[3]] <- fitDataToModel(data = dataTreatment2023ANW[[3]], model = logisticGamma, startParameters = c("r" = 0.012, "h" = 475, "shape" = 1.0))
plotFit(dataTreatment2023ANW[[3]], fitsNWLogisticGamma2023A[[3]]$par, "logistic")

# Fit of F
fitsNWLogisticGamma2023A[[4]] <- fitDataToModel(data = dataTreatment2023ANW[[4]], model = logisticGamma, startParameters = c("r" = 0.005, "h" = 650, "shape" = 1.0))
plotFit(dataTreatment2023ANW[[4]], fitsNWLogisticGamma2023A[[4]]$par, "logistic")

# Fit of F+
fitsNWLogisticGamma2023A[[5]] <- fitDataToModel(data = dataTreatment2023ANW[[5]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 575, "shape" = 1.0))
plotFit(dataTreatment2023ANW[[5]], fitsNWLogisticGamma2023A[[5]]$par, "logistic")

# Fit of T
fitsNWLogisticGamma2023A[[6]] <- fitDataToModel(data = dataTreatment2023ANW[[6]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 400, "shape" = 1.0))
plotFit(dataTreatment2023ANW[[6]], fitsNWLogisticGamma2023A[[6]]$par, "logistic")

# Fit of T+
fitsNWLogisticGamma2023A[[7]] <- fitDataToModel(data = dataTreatment2023ANW[[7]], model = logisticGamma, startParameters = c("r" = 0.012, "h" = 350, "shape" = 1.0))
plotFit(dataTreatment2023ANW[[7]], fitsNWLogisticGamma2023A[[7]]$par, "logistic")

# Fit TF-M
fitsNWLogisticGamma2023A[[8]] <- fitDataToModel(data = dataTreatment2023ANW[[8]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 520, "shape" = 1.0))
plotFit(dataTreatment2023ANW[[8]], fitsNWLogisticGamma2023A[[8]]$par, "logistic")

AICLogisticGamma2023ANW <- sum(unlist(lapply(fitsNWLogisticGamma2023A, function(fits)2 * 3 + 2 * fits$value)))


#############################################
###                                       ###
### Fit transformed logistic gamma curves ###
###                                       ###
#############################################      

fitsNWTLogisticGamma2023A <- rep(list(NA), 8)

# Fit of 1T:1F
fitsNWTLogisticGamma2023A[[1]] <- fitDataToModel(data = dataTreatment2023ANW[[1]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 520, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023ANW[[1]], fitsNWTLogisticGamma2023A[[1]]$par, "tLogistic")

# Fit of 1T:3F
fitsNWTLogisticGamma2023A[[2]] <- fitDataToModel(data = dataTreatment2023ANW[[2]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 520, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023ANW[[2]], fitsNWTLogisticGamma2023A[[2]]$par, "tLogistic")

# Fit of 3T:1F
fitsNWTLogisticGamma2023A[[3]] <- fitDataToModel(data = dataTreatment2023ANW[[3]], model = tLogisticGamma, startParameters = c("r" = 0.012, "h" = 475, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023ANW[[3]], fitsNWTLogisticGamma2023A[[3]]$par, "tLogistic")

# Fit of F
fitsNWTLogisticGamma2023A[[4]] <- fitDataToModel(data = dataTreatment2023ANW[[4]], model = tLogisticGamma, startParameters = c("r" = 0.012, "h" = 650, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023ANW[[4]], fitsNWTLogisticGamma2023A[[4]]$par, "tLogistic")

# Fit of F+
fitsNWTLogisticGamma2023A[[5]] <- fitDataToModel(data = dataTreatment2023ANW[[5]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 575, "shape" = 1.0, "d" = 1.0), dMax = 1.0)
plotFit(dataTreatment2023ANW[[5]], fitsNWTLogisticGamma2023A[[5]]$par, "tLogistic")

# Fit of T
fitsNWTLogisticGamma2023A[[6]] <- fitDataToModel(data = dataTreatment2023ANW[[6]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 400, "shape" = 1.0, "d" = 0.90), dMax = 1.0)
plotFit(dataTreatment2023ANW[[6]], fitsNWTLogisticGamma2023A[[6]]$par, "tLogistic")

# Fit of T+
fitsNWTLogisticGamma2023A[[7]] <- fitDataToModel(data = dataTreatment2023ANW[[7]], model = tLogisticGamma, startParameters = c("r" = 0.012, "h" = 350, "shape" = 1.0, "d" = 0.95), dMax = 1.0)
plotFit(dataTreatment2023ANW[[7]], fitsNWTLogisticGamma2023A[[7]]$par, "tLogistic")

# Fit TF-M
fitsNWTLogisticGamma2023A[[8]] <- fitDataToModel(data = dataTreatment2023ANW[[8]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 520, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2023ANW[[8]], fitsNWTLogisticGamma2023A[[8]]$par, "tLogistic")

AICTLogisticGamma2023ANW <- sum(unlist(lapply(fitsNWTLogisticGamma2023A, function(fits)2 * 4 + 2 * fits$value)))

AICLogisticNormal2023ANW
AICTLogisticNormal2023ANW
AICLogisticGamma2023ANW
AICTLogisticGamma2023ANW

# Transformed logistic normal is the best fitting model

# Calculate 95% CIs
parameters2023A <- names(fitsNWTLogisticNormal2023A[[1]]$par)
parameterCIs2023A <- rep(NA, times = 3 * length(parameters2023A))
names(parameterCIs2023A) <- unlist(lapply(parameters2023A, function(par)c(par, paste(par, "Lower", sep = ""), paste(par, "Upper", sep = ""))))
parameterCIs2023A <- rep(list(parameterCIs2023A), length(fitsNWTLogisticNormal2023A))

cI2023A <- function(par, x, z, opt, dMax){
  pars <- par
  pars[names(opt)] <- opt
  nll <- tLogisticNormal(pars, x, z, dMax)
  return(nll)
}

for(i in 1:length(fitsNWTLogisticNormal2023A)){
  parameterVecs <- lapply(fitsNWTLogisticNormal2023A, function(fit){
    parStart <- fit$par / 100
    parEnd <- fit$par * 20
    parVec <- list(seq(parStart[1], parEnd[1], length.out = 400), seq(parStart[2], parEnd[2], length.out = 400), seq(parStart[3], parEnd[3], length.out = 400), seq(parStart[4], parEnd[4], length.out = 400))
    names(parVec) <- names(fit$par)
    return(parVec)
  })
  optPar <- fitsNWTLogisticNormal2023A[[i]]$par
  parameterCIs2023A[[i]][names(optPar)] <- optPar
  parameterCIs2023A[[i]][unlist((lapply(names(optPar), function(n)c(paste(n, "Lower", sep = ""), paste(n, "Upper", sep = "")))))] <-
    unlist(lapply(1:length(parameterVecs[[i]]), function(j){
      t(tibble(calculate95CI(parStartVec = optPar[-which(names(optPar) == names(fitsNWTLogisticNormal2023A[[i]]$par)[[j]])],
                             optVec = parameterVecs[[i]][j],
                             model = cI2023A,
                             optNLL = fitsNWTLogisticNormal2023A[[i]]$value,
                             x = dataTreatment2023ANW[[i]][,1][[1]],
                             z = dataTreatment2023ANW[[i]][,2][[1]],
                             dMax = 1.0)))
    }))
  
}
names(parameterCIs2023A) <- names(dataTreatment2023ANW)

# Write to table
CIOut2023ANW <- tibble(Treatment = names(parameterCIs2023A)) %>% 
  bind_cols(bind_rows(parameterCIs2023A))
write.table(CIOut2023ANW, "light_interception/light_interception_NW_tLogisticNormal_parameters_2023A.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

CIOut2023ATable <- CIOut2023ANW %>% 
  mutate(r = paste0(round(r, 4), " (", round(rLower, 4), "-", round(rUpper, 4), ")"),
         h = paste0(round(h, 2), " (", round(hLower, 2), "-", round(hUpper, 2), ")"),
         sd = paste0(round(sd, 4), " (", round(sdLower, 4), "-", round(sdUpper, 4), ")"),
         d = paste0(round(d, 3), " (", round(dLower, 3), "-", round(dUpper, 3), ")")) %>% 
  dplyr::select(Treatment, r, h, sd, d)

print(xtable(CIOut2023ATable, type = "latex"), include.rownames = FALSE)

dataTreatment2023ANWM <- lapply(1:length(dataTreatment2023ANW), function(i){
  dat <- dataTreatment2023ANW[[i]]
  treatment <- names(dataTreatment2023ANW)[i]
  dat <- dat %>% 
    mutate(Treatment = treatment)
})

dataTreatment2023ANWM <- lapply(dataTreatment2023ANWM, function(Data){
  Average <- aggregate(Data$Val, list(Data$Time), function(x)mean(x, na.rm = TRUE))
  Data$Average <- Average$x[match(Data$Time, Average$Group.1)]
  return(Data)
})

# Plot triple combination graphs
geomPoints2023ANW <- lapply(dataTreatment2023ANWM, function(d){
  return(geom_point(data = d, aes(x = Time, y = Average, col = as.factor(Treatment[1])), size = 3, shape = 1))
})

geomFunctions2023ANW <- lapply(1:length(fitsNWTLogisticNormal2023A), function(i){
  o <- fitsNWTLogisticNormal2023A[[i]]
  return(geom_function(fun = function(x)o$par["d"] * (exp(o$par["r"]*(x - o$par["h"])))/(1 + exp(o$par["r"]*(x - o$par["h"]))),
                       aes(col = as.factor(dataTreatment2023ANWM[[i]]$Treatment[1])),
                       size = 1.5))
})

### Get all triple combinations
treatmentsToCompare <- unique(data2023A$Treatment)
treatmentsToCompare <- treatmentsToCompare[-which(treatmentsToCompare %in% c("T+", "F+"))] # The increased density sole crops are not compared here

triples <- lapply(treatmentsToCompare, function(treatment){
  if(!treatment %in% c("T", "F")){
    return(c(treatment, "T", "F"))
  }else{
    return(NA)
  }
})
triples <- triples[!is.na(triples)]
names(dataTreatment2023ANWM) <- lapply(dataTreatment2023ANWM, function(l)l$Treatment[1])
names(geomPoints2023ANW) <- names(dataTreatment2023ANWM)
names(geomFunctions2023ANW) <- names(dataTreatment2023ANWM)

triplePlots2023ANW <- lapply(triples, function(triple){
  plotMultipleGraphs(geomPoints2023ANW[triple], geomFunctions2023ANW[triple], triple, colourVector[c(3, 1, 2)], XLab = bquote(T[c]), YLab = bquote(p[li]), Title = "")
})

# Create double comparison plots for the sole crops
doubles <- list(c("T", "T+"), c("F", "F+"))

doublePlotsNW <- lapply(doubles, function(double){
  plotMultipleGraphs(geomPoints2023ANW[double], geomFunctions2023ANW[double], double, c(colourVector[2], colourVector[1]), XLab = bquote(T[c]), YLab = bquote(p[li]), Title = "")
})

# Create a combined plot for intercrops and normal density sole crops
combinedPlotNW <- plotMultipleGraphs(geomPoints2023ANW[treatmentsToCompare], geomFunctions2023ANW[treatmentsToCompare], treatmentsToCompare, c(colourVector[1], colourVector[2], colourVector[5], colourVector[4], colourVector[3], colourVector[6]), XLab = bquote(T[c]), YLab = bquote(p[li]), Title = "")

### Combine data for AIC comparisons
# Fit of 1T:1F
startParameters1T1F <- list(c("r" = 0.008, "h" = 480, "sd" = 0.1, "d" = 0.75),
                                 c("r" = 0.008, "h" = 500, "sd" = 0.1, "d" = 0.75),
                                 c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75),
                                 c("r" = 0.008, "h" = 550, "sd" = 0.1, "d" = 0.75),
                                 c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75))
combinationAIC1T1F <- round(fitCombinations(data = dataTreatment2023ANW, triple = triples[[1]], model = tLogisticNormal, startParameters = startParameters1T1F), 2)

# Fit of 1T:3F
startParameters1T3F <- list(c("r" = 0.008, "h" = 480, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 500, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 550, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75))
combinationAIC1T3F <- round(fitCombinations(data = dataTreatment2023ANW, triple = triples[[2]], model = tLogisticNormal, startParameters = startParameters1T3F), 2)

# Fit of 3T:1F
startParameters3T1F <- list(c("r" = 0.008, "h" = 480, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 500, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 550, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75))
combinationAIC3T1F <- round(fitCombinations(data = dataTreatment2023ANW, triple = triples[[3]], model = tLogisticNormal, startParameters = startParameters3T1F), 2)

# Fit of TF-M
startParametersTFM <- list(c("r" = 0.008, "h" = 480, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 500, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 550, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75))
combinationAICTFM <- round(fitCombinations(data = dataTreatment2023ANW, triple = triples[[4]], model = tLogisticNormal, startParameters = startParametersTFM), 2)

combinationAIC2023ANW <- tibble(Treatment = rep(unique(data2023A$Treatment)[c(1:3, 8)], each = 3),
                             Combination = rep(c("IC & C & L", "IC + C & L", "IC + L & C"), times = 4),
                             AIC = c(combinationAIC1T1F,
                                     combinationAIC1T3F,
                                     combinationAIC3T1F,
                                     combinationAICTFM),
                             Best = c("*",
                                      "",
                                      "",
                                      "*",
                                      "",
                                      "",
                                      "*",
                                      "",
                                      "",
                                      "*",
                                      "",
                                      ""))

print(xtable(combinationAIC2023ANW, type = "latex"), include.rownames = FALSE)

###################
####################
####            ####
#### With weeds ####
####            ####
####################
####################

dataTreatment2023AWW <- split(dataWeeds2023A, dataWeeds2023A$Treatment)
dataTreatment2023AWW <- lapply(dataTreatment2023AWW, function(dat){
  datNew <- dat %>% 
    dplyr::select(Time, Val)
  return(datNew)
})



##################################
###                            ###
### Fit logistic normal curves ###
###                            ###
##################################     

fitsWWLogisticNormal2023A <- rep(list(NA), 8)

# Fit of 1T:1F
fitsWWLogisticNormal2023A[[1]] <- fitDataToModel(data = dataTreatment2023AWW[[1]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 520, "sd" = 0.1))
plotFit(dataTreatment2023AWW[[1]], fitsWWLogisticNormal2023A[[1]]$par, "logistic")

# Fit of 1T:3F
fitsWWLogisticNormal2023A[[2]] <- fitDataToModel(data = dataTreatment2023AWW[[2]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 520, "sd" = 0.1))
plotFit(dataTreatment2023AWW[[2]], fitsWWLogisticNormal2023A[[2]]$par, "logistic")

# Fit of 3T:1F
fitsWWLogisticNormal2023A[[3]] <- fitDataToModel(data = dataTreatment2023AWW[[3]], model = logisticNormal, startParameters = c("r" = 0.012, "h" = 475, "sd" = 0.1))
plotFit(dataTreatment2023AWW[[3]], fitsWWLogisticNormal2023A[[3]]$par, "logistic")

# Fit of F
fitsWWLogisticNormal2023A[[4]] <- fitDataToModel(data = dataTreatment2023AWW[[4]], model = logisticNormal, startParameters = c("r" = 0.005, "h" = 650, "sd" = 0.1))
plotFit(dataTreatment2023AWW[[4]], fitsWWLogisticNormal2023A[[4]]$par, "logistic")

# Fit of F+
fitsWWLogisticNormal2023A[[5]] <- fitDataToModel(data = dataTreatment2023AWW[[5]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 575, "sd" = 0.1))
plotFit(dataTreatment2023AWW[[5]], fitsWWLogisticNormal2023A[[5]]$par, "logistic")

# Fit of T
fitsWWLogisticNormal2023A[[6]] <- fitDataToModel(data = dataTreatment2023AWW[[6]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 400, "sd" = 0.1))
plotFit(dataTreatment2023AWW[[6]], fitsWWLogisticNormal2023A[[6]]$par, "logistic")

# Fit of T+
fitsWWLogisticNormal2023A[[7]] <- fitDataToModel(data = dataTreatment2023AWW[[7]], model = logisticNormal, startParameters = c("r" = 0.012, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2023AWW[[7]], fitsWWLogisticNormal2023A[[7]]$par, "logistic")

# Fit TF-M
fitsWWLogisticNormal2023A[[8]] <- fitDataToModel(data = dataTreatment2023AWW[[8]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 520, "sd" = 0.1))
plotFit(dataTreatment2023AWW[[8]], fitsWWLogisticNormal2023A[[8]]$par, "logistic")

AICLogisticNormal2023AWW <- sum(unlist(lapply(fitsWWLogisticNormal2023A, function(fits)2 * 3 + 2 * fits$value)))


##############################################
###                                        ###
### Fit transformed logistic normal curves ###
###                                        ###
##############################################      

fitsWWTLogisticNormal2023A <- rep(list(NA), 8)

# Fit of 1T:1F
fitsWWTLogisticNormal2023A[[1]] <- fitDataToModel(data = dataTreatment2023AWW[[1]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 520, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023AWW[[1]], fitsWWTLogisticNormal2023A[[1]]$par, "tLogistic")

# Fit of 1T:3F
fitsWWTLogisticNormal2023A[[2]] <- fitDataToModel(data = dataTreatment2023AWW[[2]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 520, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023AWW[[2]], fitsWWTLogisticNormal2023A[[2]]$par, "tLogistic")

# Fit of 3T:1F
fitsWWTLogisticNormal2023A[[3]] <- fitDataToModel(data = dataTreatment2023AWW[[3]], model = tLogisticNormal, startParameters = c("r" = 0.012, "h" = 475, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023AWW[[3]], fitsWWTLogisticNormal2023A[[3]]$par, "tLogistic")

# Fit of F
fitsWWTLogisticNormal2023A[[4]] <- fitDataToModel(data = dataTreatment2023AWW[[4]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 600, "sd" = 0.1, "d" = 0.95), dMax = 1.0)
plotFit(dataTreatment2023AWW[[4]], fitsWWTLogisticNormal2023A[[4]]$par, "tLogistic")

# Fit of F+
fitsWWTLogisticNormal2023A[[5]] <- fitDataToModel(data = dataTreatment2023AWW[[5]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 575, "sd" = 0.1, "d" = 0.95), dMax = 1.0)
plotFit(dataTreatment2023AWW[[5]], fitsWWTLogisticNormal2023A[[5]]$par, "tLogistic")

# Fit of T
fitsWWTLogisticNormal2023A[[6]] <- fitDataToModel(data = dataTreatment2023AWW[[6]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 400, "sd" = 0.1, "d" = 0.90), dMax = 1.0)
plotFit(dataTreatment2023AWW[[6]], fitsWWTLogisticNormal2023A[[6]]$par, "tLogistic")

# Fit of T+
fitsWWTLogisticNormal2023A[[7]] <- fitDataToModel(data = dataTreatment2023AWW[[7]], model = tLogisticNormal, startParameters = c("r" = 0.012, "h" = 350, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2023AWW[[7]], fitsWWTLogisticNormal2023A[[7]]$par, "tLogistic")

# Fit TF-M
fitsWWTLogisticNormal2023A[[8]] <- fitDataToModel(data = dataTreatment2023AWW[[8]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 520, "sd" = 0.1, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2023AWW[[8]], fitsWWTLogisticNormal2023A[[8]]$par, "tLogistic")

AICTLogisticNormal2023AWW <- sum(unlist(lapply(fitsWWTLogisticNormal2023A, function(fits)2 * 4 + 2 * fits$value)))


#################################
###                           ###
### Fit logistic gamma curves ###
###                           ###
#################################     

fitsWWLogisticGamma2023A <- rep(list(NA), 8)

# Fit of 1T:1F
fitsWWLogisticGamma2023A[[1]] <- fitDataToModel(data = dataTreatment2023AWW[[1]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 520, "shape" = 1.0))
plotFit(dataTreatment2023AWW[[1]], fitsWWLogisticGamma2023A[[1]]$par, "logistic")

# Fit of 1T:3F
fitsWWLogisticGamma2023A[[2]] <- fitDataToModel(data = dataTreatment2023AWW[[2]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 520, "shape" = 1.0))
plotFit(dataTreatment2023AWW[[2]], fitsWWLogisticGamma2023A[[2]]$par, "logistic")

# Fit of 3T:1F
fitsWWLogisticGamma2023A[[3]] <- fitDataToModel(data = dataTreatment2023AWW[[3]], model = logisticGamma, startParameters = c("r" = 0.012, "h" = 475, "shape" = 1.0))
plotFit(dataTreatment2023AWW[[3]], fitsWWLogisticGamma2023A[[3]]$par, "logistic")

# Fit of F
fitsWWLogisticGamma2023A[[4]] <- fitDataToModel(data = dataTreatment2023AWW[[4]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 600, "shape" = 1.0))
plotFit(dataTreatment2023AWW[[4]], fitsWWLogisticGamma2023A[[4]]$par, "logistic")

# Fit of F+
fitsWWLogisticGamma2023A[[5]] <- fitDataToModel(data = dataTreatment2023AWW[[5]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 575, "shape" = 1.0))
plotFit(dataTreatment2023AWW[[5]], fitsWWLogisticGamma2023A[[5]]$par, "logistic")

# Fit of T
fitsWWLogisticGamma2023A[[6]] <- fitDataToModel(data = dataTreatment2023AWW[[6]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 400, "shape" = 1.0))
plotFit(dataTreatment2023AWW[[6]], fitsWWLogisticGamma2023A[[6]]$par, "logistic")

# Fit of T+
fitsWWLogisticGamma2023A[[7]] <- fitDataToModel(data = dataTreatment2023AWW[[7]], model = logisticGamma, startParameters = c("r" = 0.012, "h" = 350, "shape" = 1.0))
plotFit(dataTreatment2023AWW[[7]], fitsWWLogisticGamma2023A[[7]]$par, "logistic")

# Fit TF-M
fitsWWLogisticGamma2023A[[8]] <- fitDataToModel(data = dataTreatment2023AWW[[8]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 520, "shape" = 1.0))
plotFit(dataTreatment2023AWW[[8]], fitsWWLogisticGamma2023A[[8]]$par, "logistic")

AICLogisticGamma2023AWW <- sum(unlist(lapply(fitsWWLogisticGamma2023A, function(fits)2 * 3 + 2 * fits$value)))


#############################################
###                                       ###
### Fit transformed logistic gamma curves ###
###                                       ###
#############################################      

fitsWWTLogisticGamma2023A <- rep(list(NA), 8)

# Fit of 1T:1F
fitsWWTLogisticGamma2023A[[1]] <- fitDataToModel(data = dataTreatment2023AWW[[1]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 520, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023AWW[[1]], fitsWWTLogisticGamma2023A[[1]]$par, "tLogistic")

# Fit of 1T:3F
fitsWWTLogisticGamma2023A[[2]] <- fitDataToModel(data = dataTreatment2023AWW[[2]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 520, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023AWW[[2]], fitsWWTLogisticGamma2023A[[2]]$par, "tLogistic")

# Fit of 3T:1F
fitsWWTLogisticGamma2023A[[3]] <- fitDataToModel(data = dataTreatment2023AWW[[3]], model = tLogisticGamma, startParameters = c("r" = 0.012, "h" = 475, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023AWW[[3]], fitsWWTLogisticGamma2023A[[3]]$par, "tLogistic")

# Fit of F
fitsWWTLogisticGamma2023A[[4]] <- fitDataToModel(data = dataTreatment2023AWW[[4]], model = tLogisticGamma, startParameters = c("r" = 0.012, "h" = 600, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023AWW[[4]], fitsWWTLogisticGamma2023A[[4]]$par, "tLogistic")

# Fit of F+
fitsWWTLogisticGamma2023A[[5]] <- fitDataToModel(data = dataTreatment2023AWW[[5]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 575, "shape" = 1.0, "d" = 1.0), dMax = 1.0)
plotFit(dataTreatment2023AWW[[5]], fitsWWTLogisticGamma2023A[[5]]$par, "tLogistic")

# Fit of T
fitsWWTLogisticGamma2023A[[6]] <- fitDataToModel(data = dataTreatment2023AWW[[6]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 400, "shape" = 1.0, "d" = 0.90), dMax = 1.0)
plotFit(dataTreatment2023AWW[[6]], fitsWWTLogisticGamma2023A[[6]]$par, "tLogistic")

# Fit of T+
fitsWWTLogisticGamma2023A[[7]] <- fitDataToModel(data = dataTreatment2023AWW[[7]], model = tLogisticGamma, startParameters = c("r" = 0.012, "h" = 350, "shape" = 1.0, "d" = 0.95), dMax = 1.0)
plotFit(dataTreatment2023AWW[[7]], fitsWWTLogisticGamma2023A[[7]]$par, "tLogistic")

# Fit TF-M
fitsWWTLogisticGamma2023A[[8]] <- fitDataToModel(data = dataTreatment2023AWW[[8]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 520, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2023AWW[[8]], fitsWWTLogisticGamma2023A[[8]]$par, "tLogistic")

AICTLogisticGamma2023AWW <- sum(unlist(lapply(fitsWWTLogisticGamma2023A, function(fits)2 * 4 + 2 * fits$value)))

AICLogisticNormal2023AWW
AICTLogisticNormal2023AWW
AICLogisticGamma2023AWW
AICTLogisticGamma2023AWW

# Transformed logistic normal is the best fitting model

# Calculate 95% CIs
parameters2023A <- names(fitsWWTLogisticNormal2023A[[1]]$par)
parameterCIs2023A <- rep(NA, times = 3 * length(parameters2023A))
names(parameterCIs2023A) <- unlist(lapply(parameters2023A, function(par)c(par, paste(par, "Lower", sep = ""), paste(par, "Upper", sep = ""))))
parameterCIs2023A <- rep(list(parameterCIs2023A), length(fitsWWTLogisticNormal2023A))

cI2023A <- function(par, x, z, opt, dMax){
  pars <- par
  pars[names(opt)] <- opt
  nll <- tLogisticNormal(pars, x, z, dMax)
  return(nll)
}

for(i in 1:length(fitsWWTLogisticNormal2023A)){
  parameterVecs <- lapply(fitsWWTLogisticNormal2023A, function(fit){
    parStart <- fit$par / 100
    parEnd <- fit$par * 20
    parVec <- list(seq(parStart[1], parEnd[1], length.out = 400), seq(parStart[2], parEnd[2], length.out = 400), seq(parStart[3], parEnd[3], length.out = 400), seq(parStart[4], parEnd[4], length.out = 400))
    names(parVec) <- names(fit$par)
    return(parVec)
  })
  optPar <- fitsWWTLogisticNormal2023A[[i]]$par
  parameterCIs2023A[[i]][names(optPar)] <- optPar
  parameterCIs2023A[[i]][unlist((lapply(names(optPar), function(n)c(paste(n, "Lower", sep = ""), paste(n, "Upper", sep = "")))))] <-
    unlist(lapply(1:length(parameterVecs[[i]]), function(j){
      t(tibble(calculate95CI(parStartVec = optPar[-which(names(optPar) == names(fitsWWTLogisticNormal2023A[[i]]$par)[[j]])],
                             optVec = parameterVecs[[i]][j],
                             model = cI2023A,
                             optNLL = fitsWWTLogisticNormal2023A[[i]]$value,
                             x = dataTreatment2023AWW[[i]][,1][[1]],
                             z = dataTreatment2023AWW[[i]][,2][[1]],
                             dMax = 1.0)))
    }))
  
}
names(parameterCIs2023A) <- names(dataTreatment2023AWW)

# Write to table
CIOut2023AWW <- tibble(Treatment = names(parameterCIs2023A)) %>% 
  bind_cols(bind_rows(parameterCIs2023A))
write.table(CIOut2023AWW, "light_interception/light_interception_WW_tLogisticNormal_parameters_2023A.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

CIOut2023ATable <- CIOut2023AWW %>% 
  mutate(r = paste0(round(r, 4), " (", round(rLower, 4), "-", round(rUpper, 4), ")"),
         h = paste0(round(h, 2), " (", round(hLower, 2), "-", round(hUpper, 2), ")"),
         sd = paste0(round(sd, 4), " (", round(sdLower, 4), "-", round(sdUpper, 4), ")"),
         d = paste0(round(d, 3), " (", round(dLower, 3), "-", round(dUpper, 3), ")")) %>% 
  dplyr::select(Treatment, r, h, sd, d)

print(xtable(CIOut2023ATable, type = "latex"), include.rownames = FALSE)

CIOut2023ANW <- CIOut2023ANW %>% 
  mutate(Treatment = paste0(Treatment, "_H"))
CIOut2023AWW <- CIOut2023AWW %>% 
  mutate(Treatment = paste0(Treatment, "_W"))
  
CIOut2023ATable <- bind_rows(CIOut2023ANW, CIOut2023AWW) %>%
  mutate(r = paste0(round(r, 4), " (", round(rLower, 4), "-", round(rUpper, 4), ")"),
         h = paste0(round(h, 2), " (", round(hLower, 2), "-", round(hUpper, 2), ")"),
         sd = paste0(round(sd, 4), " (", round(sdLower, 4), "-", round(sdUpper, 4), ")"),
         d = paste0(round(d, 3), " (", round(dLower, 3), "-", round(dUpper, 3), ")")) %>% 
  dplyr::select(Treatment, r, h, sd, d)

print(xtable(CIOut2023ATable, type = "latex"), include.rownames = FALSE)

dataTreatment2023AWWM <- lapply(1:length(dataTreatment2023AWW), function(i){
  dat <- dataTreatment2023AWW[[i]]
  treatment <- names(dataTreatment2023AWW)[i]
  dat <- dat %>% 
    mutate(Treatment = treatment)
})

dataTreatment2023AWWM <- lapply(dataTreatment2023AWWM, function(Data){
  Average <- aggregate(Data$Val, list(Data$Time), function(x)mean(x, na.rm = TRUE))
  Data$Average <- Average$x[match(Data$Time, Average$Group.1)]
  return(Data)
})

# Plot triple combination graphs
geomPoints2023AWW <- lapply(dataTreatment2023AWWM, function(d){
  return(geom_point(data = d, aes(x = Time, y = Average, col = as.factor(Treatment[1])), size = 3, shape = 1))
})

geomFunctions2023AWW <- lapply(1:length(fitsWWTLogisticNormal2023A), function(i){
  o <- fitsWWTLogisticNormal2023A[[i]]
  return(geom_function(fun = function(x)o$par["d"] * (exp(o$par["r"]*(x - o$par["h"])))/(1 + exp(o$par["r"]*(x - o$par["h"]))),
                       aes(col = as.factor(dataTreatment2023AWWM[[i]]$Treatment[1])),
                       size = 1.5))
})

### Get all triple combinations
treatmentsToCompare <- unique(data2023A$Treatment)
treatmentsToCompare <- treatmentsToCompare[-which(treatmentsToCompare %in% c("T+", "F+"))] # The increased density sole crops are not compared here

triples <- lapply(treatmentsToCompare, function(treatment){
  if(!treatment %in% c("T", "F")){
    return(c(treatment, "T", "F"))
  }else{
    return(NA)
  }
})
triples <- triples[!is.na(triples)]
names(dataTreatment2023AWWM) <- lapply(dataTreatment2023AWWM, function(l)l$Treatment[1])
names(geomPoints2023AWW) <- names(dataTreatment2023AWWM)
names(geomFunctions2023AWW) <- names(dataTreatment2023AWWM)

triplePlots2023AWW <- lapply(triples, function(triple){
  plotMultipleGraphs(geomPoints2023AWW[triple], geomFunctions2023AWW[triple], triple, colourVector[c(3, 1, 2)], XLab = bquote(T[c]), YLab = bquote(p[li]), Title = "")
})

# Create double comparison plots for the sole crops
doubles <- list(c("T", "T+"), c("F", "F+"))

doublePlots2023AWW <- lapply(doubles, function(double){
  plotMultipleGraphs(geomPoints2023AWW[double], geomFunctions2023AWW[double], double, c(colourVector[2], colourVector[1]), XLab = bquote(T[c]), YLab = bquote(p[li]), Title = "")
})

# Create a combined plot for intercrops and normal density sole crops
combinedPlot2023AWW <- plotMultipleGraphs(geomPoints2023AWW[treatmentsToCompare], geomFunctions2023AWW[treatmentsToCompare], treatmentsToCompare, c(colourVector[1], colourVector[2], colourVector[5], colourVector[4], colourVector[3], colourVector[6]), XLab = bquote(T[c]), YLab = bquote(p[li]), Title = "")

### Combine data for AIC comparisons
# Fit of 1T:1F
startParameters1T1F <- list(c("r" = 0.008, "h" = 480, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 500, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 550, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75))
combinationAIC1T1F <- round(fitCombinations(data = dataTreatment2023AWW, triple = triples[[1]], model = tLogisticNormal, startParameters = startParameters1T1F), 2)

# Fit of 1T:3F
startParameters1T3F <- list(c("r" = 0.008, "h" = 550, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 500, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 550, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75))

combinationAIC1T3F <- round(fitCombinations(data = dataTreatment2023AWW, triple = triples[[2]], model = tLogisticNormal, startParameters = startParameters1T3F), 2)

# Fit of 3T:1F
startParameters3T1F <- list(c("r" = 0.008, "h" = 480, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 500, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 550, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75))
combinationAIC3T1F <- round(fitCombinations(data = dataTreatment2023AWW, triple = triples[[3]], model = tLogisticNormal, startParameters = startParameters3T1F), 2)

# Fit of TF-M
startParametersTFM <- list(c("r" = 0.008, "h" = 480, "sd" = 0.1, "d" = 0.75),
                           c("r" = 0.008, "h" = 500, "sd" = 0.1, "d" = 0.75),
                           c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75),
                           c("r" = 0.008, "h" = 550, "sd" = 0.1, "d" = 0.75),
                           c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75))
combinationAICTFM <- round(fitCombinations(data = dataTreatment2023AWW, triple = triples[[4]], model = tLogisticNormal, startParameters = startParametersTFM), 2)

combinationAIC2023AWW <- tibble(Treatment = rep(unique(data2023A$Treatment)[c(1:3, 8)], each = 3),
                              Combination = rep(c("IC & C & L", "IC + C & L", "IC + L & C"), times = 4),
                              AIC = c(combinationAIC1T1F,
                                      combinationAIC1T3F,
                                      combinationAIC3T1F,
                                      combinationAICTFM),
                              Best = c("*",
                                       "",
                                       "",
                                       "*",
                                       "",
                                       "",
                                       "",
                                       "*",
                                       "",
                                       "*",
                                       "",
                                       ""))

print(xtable(combinationAIC2023AWW, type = "latex"), include.rownames = FALSE)

########################
##                    ##
## Experiment 2023 B  ##
##                    ##
########################

# Read data
data2023BFull <- read_xlsx("WCF_2023_data.xlsx", sheet = "LightInterceptionB", range = "A1:G2641", col_names = TRUE) %>% 
  dplyr::select(Plot, Treatment, Weeds, DAS, Sample, PLightInterception) %>% 
  mutate(PLightInterception = as.numeric(PLightInterception)) %>% 
  rename(Val = PLightInterception,
         Time = DAS) %>% 
  filter(!(Treatment %in% c("T-25", "F-25", "1T:1F-25", "TF-M-25"))) %>% 
  filter(!Time == 35) # Remove first measurement, which was not correctly taken



data2023B <- as_tibble(aggregate(Val ~ Plot + Treatment + Time + Weeds, data = data2023BFull, FUN = function(x)mean(x, na.rm = TRUE)))
data2023B$Val[which(data2023B$Val < 0.001)] <- 0.001
data2023B$Val[which(data2023B$Val > 0.999)] <- 0.999

data2023B <- data2023B %>% 
  mutate(Treatment = factor(Treatment, levels = c("F", "F-375", "TF-M", "TF-M-375", "T", "T-375", "1T:1F", "1T:1F-375"))) %>% 
  arrange(Treatment)

dataWeeds2023B <- data2023B %>% 
  filter(Weeds == "Y") %>% 
  dplyr::select(-Weeds) %>% 
  filter(Time != 91) %>%
  mutate(Time = tempData2023$cumulativeTemp[match(Time, tempData2023$Time)]) %>% 
  dplyr::select(Treatment, Time, Val)
dataNoWeeds2023B <- data2023B %>% 
  filter(Weeds == "N") %>% 
  dplyr::select(-Weeds) %>% 
  filter(Time != 91) %>%
  mutate(Time = tempData2023$cumulativeTemp[match(Time, tempData2023$Time)]) %>% 
  dplyr::select(Treatment, Time, Val)

## Manual fits

##################
##################
####          ####
#### No Weeds ####
####          ####
##################
##################

dataTreatment2023BNW <- split(dataNoWeeds2023B, dataNoWeeds2023B$Treatment)
dataTreatment2023BNW <- lapply(dataTreatment2023BNW, function(dat){
  datNew <- dat %>% 
    dplyr::select(Time, Val)
  return(datNew)
})

##################################
###                            ###
### Fit logistic normal curves ###
###                            ###
##################################     

fitsNWLogisticNormal2023B <- rep(list(NA), 8)

# Fit of F
fitsNWLogisticNormal2023B[[1]] <- fitDataToModel(data = dataTreatment2023BNW[[1]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 700, "sd" = 0.1))
plotFit(dataTreatment2023BNW[[1]], fitsNWLogisticNormal2023B[[1]]$par, "logistic")

# Fit of F-375
fitsNWLogisticNormal2023B[[2]] <- fitDataToModel(data = dataTreatment2023BNW[[2]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 700, "sd" = 0.1))
plotFit(dataTreatment2023BNW[[2]], fitsNWLogisticNormal2023B[[2]]$par, "logistic")

# Fit of TF-M
fitsNWLogisticNormal2023B[[3]] <- fitDataToModel(data = dataTreatment2023BNW[[3]], model = logisticNormal, startParameters = c("r" = 0.02, "h" = 550, "sd" = 0.1))
plotFit(dataTreatment2023BNW[[3]], fitsNWLogisticNormal2023B[[3]]$par, "logistic")

# Fit of TF-M-375
fitsNWLogisticNormal2023B[[4]] <- fitDataToModel(data = dataTreatment2023BNW[[4]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 575, "sd" = 0.1))
plotFit(dataTreatment2023BNW[[4]], fitsNWLogisticNormal2023B[[4]]$par, "logistic")

# Fit of T
fitsNWLogisticNormal2023B[[5]] <- fitDataToModel(data = dataTreatment2023BNW[[5]], model = logisticNormal, startParameters = c("r" = 0.005, "h" = 500, "sd" = 0.1))
plotFit(dataTreatment2023BNW[[5]], fitsNWLogisticNormal2023B[[5]]$par, "logistic")

# Fit T-375
fitsNWLogisticNormal2023B[[6]] <- fitDataToModel(data = dataTreatment2023BNW[[6]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 550, "sd" = 0.1))
plotFit(dataTreatment2023BNW[[6]], fitsNWLogisticNormal2023B[[6]]$par, "logistic")

# Fit 1T:1F
fitsNWLogisticNormal2023B[[7]] <- fitDataToModel(data = dataTreatment2023BNW[[7]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 550, "sd" = 0.1))
plotFit(dataTreatment2023BNW[[7]], fitsNWLogisticNormal2023B[[7]]$par, "logistic")

# Fit 1T:1F-375
fitsNWLogisticNormal2023B[[8]] <- fitDataToModel(data = dataTreatment2023BNW[[8]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 600, "sd" = 0.1))
plotFit(dataTreatment2023BNW[[8]], fitsNWLogisticNormal2023B[[8]]$par, "logistic")


AICLogisticNormal2023BNW <- sum(unlist(lapply(fitsNWLogisticNormal2023B, function(fits)2 * 3 + 2 * fits$value)))


##############################################
###                                        ###
### Fit transformed logistic normal curves ###
###                                        ###
##############################################      

fitsNWTLogisticNormal2023B <- rep(list(NA), 8)

# Fit of F
fitsNWTLogisticNormal2023B[[1]] <- fitDataToModel(data = dataTreatment2023BNW[[1]], model = tLogisticNormal, startParameters = c("r" = 0.005, "h" = 700, "sd" = 0.1, "d" = 0.95), dMax = 1.0)
plotFit(dataTreatment2023BNW[[1]], fitsNWTLogisticNormal2023B[[1]]$par, "tLogistic")

# Fit of F-375
fitsNWTLogisticNormal2023B[[2]] <- fitDataToModel(data = dataTreatment2023BNW[[2]], model = tLogisticNormal, startParameters = c("r" = 0.03, "h" = 700, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BNW[[2]], fitsNWTLogisticNormal2023B[[2]]$par, "tLogistic")

# Fit of TF-M
fitsNWTLogisticNormal2023B[[3]] <- fitDataToModel(data = dataTreatment2023BNW[[3]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 550, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BNW[[3]], fitsNWTLogisticNormal2023B[[3]]$par, "tLogistic")

# Fit of TF-M-375
fitsNWTLogisticNormal2023B[[4]] <- fitDataToModel(data = dataTreatment2023BNW[[4]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 575, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BNW[[4]], fitsNWTLogisticNormal2023B[[4]]$par, "tLogistic")

# Fit of T
fitsNWTLogisticNormal2023B[[5]] <- fitDataToModel(data = dataTreatment2023BNW[[5]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 450, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BNW[[5]], fitsNWTLogisticNormal2023B[[5]]$par, "tLogistic")

# Fit T-375
fitsNWTLogisticNormal2023B[[6]] <- fitDataToModel(data = dataTreatment2023BNW[[6]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 550, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BNW[[6]], fitsNWTLogisticNormal2023B[[6]]$par, "tLogistic")

# Fit 1T:1F
fitsNWTLogisticNormal2023B[[7]] <- fitDataToModel(data = dataTreatment2023BNW[[7]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 550, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BNW[[7]], fitsNWTLogisticNormal2023B[[7]]$par, "tLogistic")

# Fit 1T:1F-375
fitsNWTLogisticNormal2023B[[8]] <- fitDataToModel(data = dataTreatment2023BNW[[8]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 600, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BNW[[8]], fitsNWTLogisticNormal2023B[[8]]$par, "tLogistic")

AICTLogisticNormal2023BNW <- sum(unlist(lapply(fitsNWTLogisticNormal2023B, function(fits)2 * 4 + 2 * fits$value)))


#################################
###                           ###
### Fit logistic gamma curves ###
###                           ###
#################################     

fitsNWLogisticGamma2023B <- rep(list(NA), 8)

# Fit of F
fitsNWLogisticGamma2023B[[1]] <- fitDataToModel(data = dataTreatment2023BNW[[1]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 700, "shape" = 1.0))
plotFit(dataTreatment2023BNW[[1]], fitsNWLogisticGamma2023B[[1]]$par, "logistic")

# Fit of F-375
fitsNWLogisticGamma2023B[[2]] <- fitDataToModel(data = dataTreatment2023BNW[[2]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 700, "shape" = 1.0))
plotFit(dataTreatment2023BNW[[2]], fitsNWLogisticGamma2023B[[2]]$par, "logistic")

# Fit of TF-M
fitsNWLogisticGamma2023B[[3]] <- fitDataToModel(data = dataTreatment2023BNW[[3]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 550, "shape" = 1.0))
plotFit(dataTreatment2023BNW[[3]], fitsNWLogisticGamma2023B[[3]]$par, "logistic")

# Fit of TF-M-375
fitsNWLogisticGamma2023B[[4]] <- fitDataToModel(data = dataTreatment2023BNW[[4]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 575, "shape" = 1.0))
plotFit(dataTreatment2023BNW[[4]], fitsNWLogisticGamma2023B[[4]]$par, "logistic")

# Fit of T
fitsNWLogisticGamma2023B[[5]] <- fitDataToModel(data = dataTreatment2023BNW[[5]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 450, "shape" = 1.0))
plotFit(dataTreatment2023BNW[[5]], fitsNWLogisticGamma2023B[[5]]$par, "logistic")

# Fit T-375
fitsNWLogisticGamma2023B[[6]] <- fitDataToModel(data = dataTreatment2023BNW[[6]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 550, "shape" = 1.0))
plotFit(dataTreatment2023BNW[[6]], fitsNWLogisticGamma2023B[[6]]$par, "logistic")

# Fit 1T:1F
fitsNWLogisticGamma2023B[[7]] <- fitDataToModel(data = dataTreatment2023BNW[[7]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 550, "shape" = 1.0))
plotFit(dataTreatment2023BNW[[7]], fitsNWLogisticGamma2023B[[7]]$par, "logistic")

# Fit 1T:1F-375
fitsNWLogisticGamma2023B[[8]] <- fitDataToModel(data = dataTreatment2023BNW[[8]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 600, "shape" = 1.0))
plotFit(dataTreatment2023BNW[[8]], fitsNWLogisticGamma2023B[[8]]$par, "logistic")

AICLogisticGamma2023BNW <- sum(unlist(lapply(fitsNWLogisticGamma2023B, function(fits)2 * 3 + 2 * fits$value)))


#############################################
###                                       ###
### Fit transformed logistic gamma curves ###
###                                       ###
#############################################      

fitsNWTLogisticGamma2023B <- rep(list(NA), 8)

# Fit of F
fitsNWTLogisticGamma2023B[[1]] <- fitDataToModel(data = dataTreatment2023BNW[[1]], model = tLogisticGamma, startParameters = c("r" = 0.02, "h" = 700, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BNW[[1]], fitsNWTLogisticGamma2023B[[1]]$par, "tLogistic")

# Fit of F-375
fitsNWTLogisticGamma2023B[[2]] <- fitDataToModel(data = dataTreatment2023BNW[[2]], model = tLogisticGamma, startParameters = c("r" = 0.03, "h" = 700, "shape" = 1.0, "d" = 0.90), dMax = 1.0)
plotFit(dataTreatment2023BNW[[2]], fitsNWTLogisticGamma2023B[[2]]$par, "tLogistic")

# Fit of TF-M
fitsNWTLogisticGamma2023B[[3]] <- fitDataToModel(data = dataTreatment2023BNW[[3]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 550, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BNW[[3]], fitsNWTLogisticGamma2023B[[3]]$par, "tLogistic")

# Fit of TF-M-375
fitsNWTLogisticGamma2023B[[4]] <- fitDataToModel(data = dataTreatment2023BNW[[4]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 575, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BNW[[4]], fitsNWTLogisticGamma2023B[[4]]$par, "tLogistic")

# Fit of T
fitsNWTLogisticGamma2023B[[5]] <- fitDataToModel(data = dataTreatment2023BNW[[5]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 450, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BNW[[5]], fitsNWTLogisticGamma2023B[[5]]$par, "tLogistic")

# Fit T-375
fitsNWTLogisticGamma2023B[[6]] <- fitDataToModel(data = dataTreatment2023BNW[[6]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 550, "shape" = 1.0, "d" = 0.75), dMax = 1.0)
plotFit(dataTreatment2023BNW[[6]], fitsNWTLogisticGamma2023B[[6]]$par, "tLogistic")

# Fit 1T:1F
fitsNWTLogisticGamma2023B[[7]] <- fitDataToModel(data = dataTreatment2023BNW[[7]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 550, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BNW[[7]], fitsNWTLogisticGamma2023B[[7]]$par, "tLogistic")

# Fit 1T:1F-375
fitsNWTLogisticGamma2023B[[8]] <- fitDataToModel(data = dataTreatment2023BNW[[8]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 600, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BNW[[8]], fitsNWTLogisticGamma2023B[[8]]$par, "tLogistic")


AICTLogisticGamma2023BNW <- sum(unlist(lapply(fitsNWTLogisticGamma2023B, function(fits)2 * 4 + 2 * fits$value)))

AICLogisticNormal2023BNW
AICTLogisticNormal2023BNW
AICLogisticGamma2023BNW
AICTLogisticGamma2023BNW

# Transformed logistic normal is the best fitting model

# Calculate 95% CIs
parameters2023BNW <- names(fitsNWTLogisticNormal2023B[[1]]$par)
parameterCIs2023BNW <- rep(NA, times = 3 * length(parameters2023BNW))
names(parameterCIs2023BNW) <- unlist(lapply(parameters2023BNW, function(par)c(par, paste(par, "Lower", sep = ""), paste(par, "Upper", sep = ""))))
parameterCIs2023BNW <- rep(list(parameterCIs2023BNW), length(fitsNWTLogisticNormal2023B))

cI2023BNW <- function(par, x, z, opt, dMax){
  pars <- par
  pars[names(opt)] <- opt
  nll <- tLogisticNormal(pars, x, z, dMax)
  return(nll)
}

for(i in 1:length(fitsNWTLogisticNormal2023B)){
  parameterVecs <- lapply(fitsNWTLogisticNormal2023B, function(fit){
    parStart <- fit$par / 100
    parEnd <- fit$par * 20
    parVec <- list(seq(parStart[1], parEnd[1], length.out = 400), seq(parStart[2], parEnd[2], length.out = 400), seq(parStart[3], parEnd[3], length.out = 400), seq(parStart[4], parEnd[4], length.out = 400))
    names(parVec) <- names(fit$par)
    return(parVec)
  })
  optPar <- fitsNWTLogisticNormal2023B[[i]]$par
  parameterCIs2023BNW[[i]][names(optPar)] <- optPar
  parameterCIs2023BNW[[i]][unlist((lapply(names(optPar), function(n)c(paste(n, "Lower", sep = ""), paste(n, "Upper", sep = "")))))] <-
    unlist(lapply(1:length(parameterVecs[[i]]), function(j){
      t(tibble(calculate95CI(parStartVec = optPar[-which(names(optPar) == names(fitsNWTLogisticNormal2023B[[i]]$par)[[j]])],
                             optVec = parameterVecs[[i]][j],
                             model = cI2023BNW,
                             optNLL = fitsNWTLogisticNormal2023B[[i]]$value,
                             x = dataTreatment2023BNW[[i]][,1][[1]],
                             z = dataTreatment2023BNW[[i]][,2][[1]],
                             dMax = 1.0)))
    }))
  
}
names(parameterCIs2023BNW) <- names(dataTreatment2023BNW)

# Write to table
CIOut2023BNW <- tibble(Treatment = names(parameterCIs2023BNW)) %>% 
  bind_cols(bind_rows(parameterCIs2023BNW))
write.table(CIOut2023BNW, "light_interception/light_interception_tLogisticNormal_parameters_2023BNW.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

CIOut2023BNWTable <- CIOut2023BNW %>% 
  mutate(r = paste0(round(r, 4), " (", round(rLower, 4), "-", round(rUpper, 4), ")"),
         h = paste0(round(h, 2), " (", round(hLower, 2), "-", round(hUpper, 2), ")"),
         sd = paste0(round(sd, 4), " (", round(sdLower, 4), "-", round(sdUpper, 4), ")"),
         d = paste0(round(d, 3), " (", round(dLower, 3), "-", round(dUpper, 3), ")")) %>% 
  dplyr::select(Treatment, r, h, sd, d)

print(xtable(CIOut2023BNWTable, type = "latex"), include.rownames = FALSE)


dataTreatment2023BNWM <- lapply(1:length(dataTreatment2023BNW), function(i){
  dat <- dataTreatment2023BNW[[i]]
  treatment <- names(dataTreatment2023BNW)[i]
  dat <- dat %>% 
    mutate(Treatment = treatment)
})

names(dataTreatment2023BNWM) <- names(dataTreatment2023BNW)

dataTreatment2023BNWM <- lapply(dataTreatment2023BNWM, function(Data){
  Average <- aggregate(Data$Val, list(Data$Time), function(x)mean(x, na.rm = TRUE))
  Data$Average <- Average$x[match(Data$Time, Average$Group.1)]
  return(Data)
})

# Plot triple combination graphs
geomPoints2023BNW <- lapply(dataTreatment2023BNWM, function(d){
  return(geom_point(data = d, aes(x = Time, y = Average, col = as.factor(Treatment[1])), size = 3, shape = 1))
})

names(geomPoints2023BNW) <- names(dataTreatment2023BNWM)

geomFunctions2023BNW <- lapply(1:length(fitsNWTLogisticNormal2023B), function(i){
  o <- fitsNWTLogisticNormal2023B[[i]]
  return(geom_function(fun = function(x)o$par["d"] * (exp(o$par["r"]*(x - o$par["h"])))/(1 + exp(o$par["r"]*(x - o$par["h"]))),
                       aes(col = as.factor(dataTreatment2023BNWM[[i]]$Treatment[1])),
                       size = 1.5))
})

names(geomFunctions2023BNW) <- names(dataTreatment2023BNWM)

values2023BNW <- colourVector[1:length(fitsNWTLogisticNormal2023B)]
names(values2023BNW) <- unname(sapply(dataTreatment2023BNWM, function(d)d$Treatment[1]))
labels2023BNW <- unname(sapply(dataTreatment2023BNWM, function(d)d$Treatment[1]))

plots2023BNW <- lapply(1:length(geomPoints2023BNW), function(i){
  plot <- ggplot() +
    geomPoints2023BNW[[i]] +
    geomFunctions2023BNW[[i]] +
    theme_classic(base_size = 25) +
    scale_color_manual(labels = labels2023BNW[i],
                       values = "#7FC97F") +
    labs(title = "",
         y = bquote(p[li]),
         x = bquote(T[c]),
         color = "Treatment") 
})

# Create a combined plot for intercrops and normal density sole crops
combinedPlot2023BNW <- plotMultipleGraphs(geomPoints2023BNW, geomFunctions2023BNW, names(dataTreatment2023BNW), c(colourVector[c(2, 9, 8, 4, 1, 5, 3, 7)]), XLab = bquote(T[c]), YLab = bquote(p[li]), Title = "")

combinedFaba2023BNW <- plotMultipleGraphs(geomPoints2023BNW[c(1, 2)], geomFunctions2023BNW[c(1, 2)], c("F", "F-375"), c(colourVector[c(2, 9)]), XLab = bquote(T[c]), YLab = bquote(p[li]), Title = "")
combinedTriticale2023BNW <- plotMultipleGraphs(geomPoints2023BNW[c(5, 6)], geomFunctions2023BNW[c(5, 6)], c("T", "T-375"), c(colourVector[c(1, 5)]), XLab = bquote(T[c]), YLab = bquote(p[li]), Title = "")
combinedIC2023BNW <- plotMultipleGraphs(geomPoints2023BNW[c(3, 4, 7, 8)], geomFunctions2023BNW[c(3, 4, 7, 8)], c("TF-M", "TF-M-375", "1T:1F", "1T:1F-375"), c(colourVector[c(8, 4, 3, 7)]), XLab = bquote(T[c]), YLab = bquote(p[li]), Title = "")
combined1252023BNW <- plotMultipleGraphs(geomPoints2023BNW[c(1, 3, 5, 7)], geomFunctions2023BNW[c(1, 3, 5, 7)], c("F", "TF-M", "T", "1T:1F"), c(colourVector[c(2, 8, 1, 3)]), XLab = bquote(T[c]), YLab = bquote(p[li]), Title = "")
combined3752023BNW <- plotMultipleGraphs(geomPoints2023BNW[c(2, 4, 6, 8)], geomFunctions2023BNW[c(2, 4, 6, 8)], c("F-375", "TF-M-375", "T-375", "1T:1F-375"), c(colourVector[c(9, 4, 5, 7)]), XLab = bquote(T[c]), YLab = bquote(p[li]), Title = "")



### Combine data for AIC comparisons
# Fit of 1T:1F
startParameters1T1F <- list(c("r" = 0.008, "h" = 520, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 500, "sd" = 0.1, "d" = 0.8),
                            c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 550, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75))
combinationAIC1T1F <- round(fitCombinations(data = dataTreatment2023BNW, triple = c("1T:1F", "T", "F"), model = tLogisticNormal, startParameters = startParameters1T1F), 2)

# Fit of 1T:1F-375
startParameters1T1F375 <- list(c("r" = 0.008, "h" = 480, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 500, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 550, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75))
combinationAIC1T1F375 <- round(fitCombinations(data = dataTreatment2023BNW, triple = c("1T:1F-375", "T-375", "F-375"), model = tLogisticNormal, startParameters = startParameters1T1F375), 2)

# Fit of TF-M
startParametersTFM <- list(c("r" = 0.008, "h" = 520, "sd" = 0.1, "d" = 0.75),
                             c("r" = 0.008, "h" = 500, "sd" = 0.1, "d" = 0.8),
                             c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75),
                             c("r" = 0.008, "h" = 550, "sd" = 0.1, "d" = 0.75),
                             c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75))
combinationAICTFM <- round(fitCombinations(data = dataTreatment2023BNW, triple = c("TF-M", "T", "F"), model = tLogisticNormal, startParameters = startParametersTFM), 2)

# Fit of TF-M-375
startParametersTFM375 <- list(c("r" = 0.008, "h" = 480, "sd" = 0.1, "d" = 0.75),
                             c("r" = 0.008, "h" = 500, "sd" = 0.1, "d" = 0.75),
                             c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75),
                             c("r" = 0.008, "h" = 550, "sd" = 0.1, "d" = 0.75),
                             c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75))
combinationAICTFM375 <- round(fitCombinations(data = dataTreatment2023BNW, triple = c("TF-M-375", "T-375", "F-375"), model = tLogisticNormal, startParameters = startParametersTFM375), 2)

combinationAIC2023BNW <- tibble(Treatment = rep(c("1T:1F", "1T:1F-375", "TF-M", "TF-M-375"), each = 3),
                              Combination = rep(c("IC & C & L", "IC + C & L", "IC + L & C"), times = 4),
                              AIC = c(combinationAIC1T1F,
                                      combinationAIC1T1F375,
                                      combinationAICTFM,
                                      combinationAICTFM375),
                              Best = c("*",
                                       "*",
                                       "",
                                       "*",
                                       "",
                                       "",
                                       "*",
                                       "",
                                       "",
                                       "*",
                                       "",
                                       ""))

print(xtable(combinationAIC2023BNW, type = "latex"), include.rownames = FALSE)

###################
###################
####           ####
#### WithWeeds ####
####           ####
###################
###################

dataTreatment2023BWW <- split(dataWeeds2023B, dataWeeds2023B$Treatment)
dataTreatment2023BWW <- lapply(dataTreatment2023BWW, function(dat){
  datNew <- dat %>% 
    dplyr::select(Time, Val)
  return(datNew)
})


##################################
###                            ###
### Fit logistic normal curves ###
###                            ###
##################################     

fitsWWLogisticNormal2023B <- rep(list(NA), 8)

# Fit of F
fitsWWLogisticNormal2023B[[1]] <- fitDataToModel(data = dataTreatment2023BWW[[1]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 700, "sd" = 0.1))
plotFit(dataTreatment2023BWW[[1]], fitsWWLogisticNormal2023B[[1]]$par, "logistic")

# Fit of F-25
# fitsWWLogisticNormal2023B[[2]] <- fitDataToModel(data = dataTreatment2023BWW[[2]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 700, "sd" = 0.1))
# plotFit(dataTreatment2023BWW[[2]], fitsWWLogisticNormal2023B[[2]]$par, "logistic")

# Fit of F-375
fitsWWLogisticNormal2023B[[2]] <- fitDataToModel(data = dataTreatment2023BWW[[2]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 700, "sd" = 0.1))
plotFit(dataTreatment2023BWW[[2]], fitsWWLogisticNormal2023B[[2]]$par, "logistic")

# Fit of TF-M
fitsWWLogisticNormal2023B[[3]] <- fitDataToModel(data = dataTreatment2023BWW[[3]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 550, "sd" = 0.1))
plotFit(dataTreatment2023BWW[[3]], fitsWWLogisticNormal2023B[[3]]$par, "logistic")

# Fit of TF-M-375
fitsWWLogisticNormal2023B[[4]] <- fitDataToModel(data = dataTreatment2023BWW[[4]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 575, "sd" = 0.1))
plotFit(dataTreatment2023BWW[[4]], fitsWWLogisticNormal2023B[[4]]$par, "logistic")

# Fit of T
fitsWWLogisticNormal2023B[[5]] <- fitDataToModel(data = dataTreatment2023BWW[[5]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 450, "sd" = 0.1))
plotFit(dataTreatment2023BWW[[5]], fitsWWLogisticNormal2023B[[5]]$par, "logistic")

# Fit of T-25
fitsWWLogisticNormal2023B[[7]] <- fitDataToModel(data = dataTreatment2023BWW[[7]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 450, "sd" = 0.1))
plotFit(dataTreatment2023BWW[[7]], fitsWWLogisticNormal2023B[[7]]$par, "logistic")

# Fit T-375
fitsWWLogisticNormal2023B[[6]] <- fitDataToModel(data = dataTreatment2023BWW[[6]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 480, "sd" = 0.1))
plotFit(dataTreatment2023BWW[[6]], fitsWWLogisticNormal2023B[[6]]$par, "logistic")

# Fit 1T:1F
fitsWWLogisticNormal2023B[[7]] <- fitDataToModel(data = dataTreatment2023BWW[[7]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 550, "sd" = 0.1))
plotFit(dataTreatment2023BWW[[7]], fitsWWLogisticNormal2023B[[7]]$par, "logistic")

# Fit 1T:1F-25
# fitsWWLogisticNormal2023B[[10]] <- fitDataToModel(data = dataTreatment2023BWW[[10]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 550, "sd" = 0.1))
# plotFit(dataTreatment2023BWW[[10]], fitsWWLogisticNormal2023B[[10]]$par, "logistic")

# Fit 1T:1F-375
fitsWWLogisticNormal2023B[[8]] <- fitDataToModel(data = dataTreatment2023BWW[[8]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 600, "sd" = 0.1))
plotFit(dataTreatment2023BWW[[8]], fitsWWLogisticNormal2023B[[8]]$par, "logistic")


AICLogisticNormal2023BWW <- sum(unlist(lapply(fitsWWLogisticNormal2023B, function(fits)2 * 3 + 2 * fits$value)))


##############################################
###                                        ###
### Fit transformed logistic normal curves ###
###                                        ###
##############################################      

fitsWWTLogisticNormal2023B <- rep(list(NA), 8)

# Fit of F
fitsWWTLogisticNormal2023B[[1]] <- fitDataToModel(data = dataTreatment2023BWW[[1]], model = tLogisticNormal, startParameters = c("r" = 0.005, "h" = 700, "sd" = 0.1, "d" = 0.95), dMax = 1.0)
plotFit(dataTreatment2023BWW[[1]], fitsWWTLogisticNormal2023B[[1]]$par, "tLogistic")

# Fit of F-25
# fitsWWTLogisticNormal2023B[[2]] <- fitDataToModel(data = dataTreatment2023BWW[[2]], model = tLogisticNormal, startParameters = c("r" = 0.005, "h" = 700, "sd" = 0.1, "d" = 0.9), dMax = 1.0)
# plotFit(dataTreatment2023BWW[[2]], fitsWWTLogisticNormal2023B[[2]]$par, "tLogistic")

# Fit of F-375
fitsWWTLogisticNormal2023B[[2]] <- fitDataToModel(data = dataTreatment2023BWW[[2]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 700, "sd" = 0.1, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2023BWW[[2]], fitsWWTLogisticNormal2023B[[2]]$par, "tLogistic")

# Fit of TF-M
fitsWWTLogisticNormal2023B[[3]] <- fitDataToModel(data = dataTreatment2023BWW[[3]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 550, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BWW[[3]], fitsWWTLogisticNormal2023B[[3]]$par, "tLogistic")

# Fit of TF-M-375
fitsWWTLogisticNormal2023B[[4]] <- fitDataToModel(data = dataTreatment2023BWW[[4]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 575, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BWW[[4]], fitsWWTLogisticNormal2023B[[4]]$par, "tLogistic")

# Fit of T
fitsWWTLogisticNormal2023B[[5]] <- fitDataToModel(data = dataTreatment2023BWW[[5]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 450, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BWW[[5]], fitsWWTLogisticNormal2023B[[5]]$par, "tLogistic")

# Fit of T-25
# fitsWWTLogisticNormal2023B[[7]] <- fitDataToModel(data = dataTreatment2023BWW[[7]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 450, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
# plotFit(dataTreatment2023BWW[[7]], fitsWWTLogisticNormal2023B[[7]]$par, "tLogistic")

# Fit T-375
fitsWWTLogisticNormal2023B[[6]] <- fitDataToModel(data = dataTreatment2023BWW[[6]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 550, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BWW[[6]], fitsWWTLogisticNormal2023B[[6]]$par, "tLogistic")

# Fit 1T:1F
fitsWWTLogisticNormal2023B[[7]] <- fitDataToModel(data = dataTreatment2023BWW[[7]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 550, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BWW[[7]], fitsWWTLogisticNormal2023B[[7]]$par, "tLogistic")

# Fit 1T:1F-25
# fitsWWTLogisticNormal2023B[[10]] <- fitDataToModel(data = dataTreatment2023BWW[[10]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 550, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
# plotFit(dataTreatment2023BWW[[10]], fitsWWTLogisticNormal2023B[[10]]$par, "tLogistic")

# Fit 1T:1F-375
fitsWWTLogisticNormal2023B[[8]] <- fitDataToModel(data = dataTreatment2023BWW[[8]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 600, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BWW[[8]], fitsWWTLogisticNormal2023B[[8]]$par, "tLogistic")

AICTLogisticNormal2023BWW <- sum(unlist(lapply(fitsWWTLogisticNormal2023B, function(fits)2 * 4 + 2 * fits$value)))


#################################
###                           ###
### Fit logistic gamma curves ###
###                           ###
#################################     

fitsWWLogisticGamma2023B <- rep(list(NA), 8)

# Fit of F
fitsWWLogisticGamma2023B[[1]] <- fitDataToModel(data = dataTreatment2023BWW[[1]], model = logisticGamma, startParameters = c("r" = 0.03, "h" = 650, "shape" = 1.0))
plotFit(dataTreatment2023BWW[[1]], fitsWWLogisticGamma2023B[[1]]$par, "logistic")

# Fit of F-25
# fitsWWLogisticGamma2023B[[2]] <- fitDataToModel(data = dataTreatment2023BWW[[2]], model = logisticGamma, startParameters = c("r" = 0.02, "h" = 600, "shape" = 1.0))
# plotFit(dataTreatment2023BWW[[2]], fitsWWLogisticGamma2023B[[2]]$par, "logistic")

# Fit of F-375
fitsWWLogisticGamma2023B[[2]] <- fitDataToModel(data = dataTreatment2023BWW[[2]], model = logisticGamma, startParameters = c("r" = 0.005, "h" = 600, "shape" = 1.0))
plotFit(dataTreatment2023BWW[[2]], fitsWWLogisticGamma2023B[[2]]$par, "logistic")

# Fit of TF-M
fitsWWLogisticGamma2023B[[3]] <- fitDataToModel(data = dataTreatment2023BWW[[3]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 600, "shape" = 1.0))
plotFit(dataTreatment2023BWW[[3]], fitsWWLogisticGamma2023B[[3]]$par, "logistic")

# Fit of TF-M-375
fitsWWLogisticGamma2023B[[4]] <- fitDataToModel(data = dataTreatment2023BWW[[4]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 575, "shape" = 1.0))
plotFit(dataTreatment2023BWW[[4]], fitsWWLogisticGamma2023B[[4]]$par, "logistic")

# Fit of T
fitsWWLogisticGamma2023B[[5]] <- fitDataToModel(data = dataTreatment2023BWW[[5]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 450, "shape" = 1.0))
plotFit(dataTreatment2023BWW[[5]], fitsWWLogisticGamma2023B[[5]]$par, "logistic")

# Fit of T-25
# fitsWWLogisticGamma2023B[[7]] <- fitDataToModel(data = dataTreatment2023BWW[[7]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 450, "shape" = 1.0))
# plotFit(dataTreatment2023BWW[[7]], fitsWWLogisticGamma2023B[[7]]$par, "logistic")

# Fit T-375
fitsWWLogisticGamma2023B[[6]] <- fitDataToModel(data = dataTreatment2023BWW[[6]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 550, "shape" = 1.0))
plotFit(dataTreatment2023BWW[[6]], fitsWWLogisticGamma2023B[[6]]$par, "logistic")

# Fit 1T:1F
fitsWWLogisticGamma2023B[[7]] <- fitDataToModel(data = dataTreatment2023BWW[[7]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 550, "shape" = 1.0))
plotFit(dataTreatment2023BWW[[7]], fitsWWLogisticGamma2023B[[7]]$par, "logistic")

# Fit 1T:1F-25
# fitsWWLogisticGamma2023B[[10]] <- fitDataToModel(data = dataTreatment2023BWW[[10]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 550, "shape" = 1.0))
# plotFit(dataTreatment2023BWW[[10]], fitsWWLogisticGamma2023B[[10]]$par, "logistic")

# Fit 1T:1F-375
fitsWWLogisticGamma2023B[[8]] <- fitDataToModel(data = dataTreatment2023BWW[[8]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 600, "shape" = 1.0))
plotFit(dataTreatment2023BWW[[8]], fitsWWLogisticGamma2023B[[8]]$par, "logistic")

AICLogisticGamma2023BWW <- sum(unlist(lapply(fitsWWLogisticGamma2023B, function(fits)2 * 3 + 2 * fits$value)))


#############################################
###                                       ###
### Fit transformed logistic gamma curves ###
###                                       ###
#############################################      

fitsWWTLogisticGamma2023B <- rep(list(NA), 8)

# Fit of F
fitsWWTLogisticGamma2023B[[1]] <- fitDataToModel(data = dataTreatment2023BWW[[1]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 650, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BWW[[1]], fitsWWTLogisticGamma2023B[[1]]$par, "tLogistic")

# Fit of F-25
# fitsWWTLogisticGamma2023B[[2]] <- fitDataToModel(data = dataTreatment2023BWW[[2]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 650, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
# plotFit(dataTreatment2023BWW[[2]], fitsWWTLogisticGamma2023B[[2]]$par, "tLogistic")

# Fit of F-375
fitsWWTLogisticGamma2023B[[2]] <- fitDataToModel(data = dataTreatment2023BWW[[2]], model = tLogisticGamma, startParameters = c("r" = 0.005, "h" = 600, "shape" = 1.0, "d" = 0.90), dMax = 1.0)
plotFit(dataTreatment2023BWW[[2]], fitsWWTLogisticGamma2023B[[2]]$par, "tLogistic")

# Fit of TF-M
fitsWWTLogisticGamma2023B[[3]] <- fitDataToModel(data = dataTreatment2023BWW[[3]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 650, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BWW[[3]], fitsWWTLogisticGamma2023B[[3]]$par, "tLogistic")

# Fit of TF-M-375
fitsWWTLogisticGamma2023B[[4]] <- fitDataToModel(data = dataTreatment2023BWW[[4]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 575, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BWW[[4]], fitsWWTLogisticGamma2023B[[4]]$par, "tLogistic")

# Fit of T
fitsWWTLogisticGamma2023B[[5]] <- fitDataToModel(data = dataTreatment2023BWW[[5]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 450, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BWW[[5]], fitsWWTLogisticGamma2023B[[5]]$par, "tLogistic")

# Fit of T-25
# fitsWWTLogisticGamma2023B[[7]] <- fitDataToModel(data = dataTreatment2023BWW[[7]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 450, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
# plotFit(dataTreatment2023BWW[[7]], fitsWWTLogisticGamma2023B[[7]]$par, "tLogistic")

# Fit T-375
fitsWWTLogisticGamma2023B[[6]] <- fitDataToModel(data = dataTreatment2023BWW[[6]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 550, "shape" = 1.0, "d" = 0.75), dMax = 1.0)
plotFit(dataTreatment2023BWW[[6]], fitsWWTLogisticGamma2023B[[6]]$par, "tLogistic")

# Fit 1T:1F
fitsWWTLogisticGamma2023B[[7]] <- fitDataToModel(data = dataTreatment2023BWW[[7]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 550, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BWW[[7]], fitsWWTLogisticGamma2023B[[7]]$par, "tLogistic")

# Fit 1T:1F-25
# fitsWWTLogisticGamma2023B[[10]] <- fitDataToModel(data = dataTreatment2023BWW[[10]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 550, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
# plotFit(dataTreatment2023BWW[[10]], fitsWWTLogisticGamma2023B[[10]]$par, "tLogistic")

# Fit 1T:1F-375
fitsWWTLogisticGamma2023B[[8]] <- fitDataToModel(data = dataTreatment2023BWW[[8]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 600, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BWW[[8]], fitsWWTLogisticGamma2023B[[8]]$par, "tLogistic")


AICTLogisticGamma2023BWW <- sum(unlist(lapply(fitsWWTLogisticGamma2023B, function(fits)2 * 4 + 2 * fits$value)))

AICLogisticNormal2023BWW
AICTLogisticNormal2023BWW
AICLogisticGamma2023BWW
AICTLogisticGamma2023BWW

# Transformed logistic normal is the best fitting model

# Calculate 95% CIs
parameters2023BWW <- names(fitsWWTLogisticNormal2023B[[1]]$par)
parameterCIs2023BWW <- rep(NA, times = 3 * length(parameters2023BWW))
names(parameterCIs2023BWW) <- unlist(lapply(parameters2023BWW, function(par)c(par, paste(par, "Lower", sep = ""), paste(par, "Upper", sep = ""))))
parameterCIs2023BWW <- rep(list(parameterCIs2023BWW), length(fitsWWTLogisticNormal2023B))

cI2023BWW <- function(par, x, z, opt, dMax){
  pars <- par
  pars[names(opt)] <- opt
  nll <- tLogisticNormal(pars, x, z, dMax)
  return(nll)
}

for(i in 1:length(fitsWWTLogisticNormal2023B)){
  parameterVecs <- lapply(fitsWWTLogisticNormal2023B, function(fit){
    parStart <- fit$par / 100
    parEnd <- fit$par * 20
    parVec <- list(seq(parStart[1], parEnd[1], length.out = 400), seq(parStart[2], parEnd[2], length.out = 400), seq(parStart[3], parEnd[3], length.out = 400), seq(parStart[4], parEnd[4], length.out = 400))
    names(parVec) <- names(fit$par)
    return(parVec)
  })
  optPar <- fitsWWTLogisticNormal2023B[[i]]$par
  parameterCIs2023BWW[[i]][names(optPar)] <- optPar
  parameterCIs2023BWW[[i]][unlist((lapply(names(optPar), function(n)c(paste(n, "Lower", sep = ""), paste(n, "Upper", sep = "")))))] <-
    unlist(lapply(1:length(parameterVecs[[i]]), function(j){
      t(tibble(calculate95CI(parStartVec = optPar[-which(names(optPar) == names(fitsWWTLogisticNormal2023B[[i]]$par)[[j]])],
                             optVec = parameterVecs[[i]][j],
                             model = cI2023BWW,
                             optNLL = fitsWWTLogisticNormal2023B[[i]]$value,
                             x = dataTreatment2023BWW[[i]][,1][[1]],
                             z = dataTreatment2023BWW[[i]][,2][[1]],
                             dMax = 1.0)))
    }))
  
}
names(parameterCIs2023BWW) <- names(dataTreatment2023BWW)

# Write to table
CIOut2023BWW <- tibble(Treatment = names(parameterCIs2023BWW)) %>% 
  bind_cols(bind_rows(parameterCIs2023BWW))
write.table(CIOut2023BWW, "light_interception/light_interception_tLogisticNormal_parameters_2023BWW.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

CIOut2023BWWTable <- CIOut2023BWW %>% 
  mutate(r = paste0(round(r, 4), " (", round(rLower, 4), "-", round(rUpper, 4), ")"),
         h = paste0(round(h, 2), " (", round(hLower, 2), "-", round(hUpper, 2), ")"),
         sd = paste0(round(sd, 4), " (", round(sdLower, 4), "-", round(sdUpper, 4), ")"),
         d = paste0(round(d, 3), " (", round(dLower, 3), "-", round(dUpper, 3), ")")) %>% 
  dplyr::select(Treatment, r, h, sd, d)

print(xtable(CIOut2023BWWTable, type = "latex"), include.rownames = FALSE)

CIOut2023BNW <- CIOut2023BNW %>% 
  mutate(Treatment = paste0(Treatment, "_H"))
CIOut2023BWW <- CIOut2023BWW %>% 
  mutate(Treatment = paste0(Treatment, "_W"))
  
CIOut2023BTable <- bind_rows(CIOut2023BNW, CIOut2023BWW) %>%
  mutate(r = paste0(round(r, 4), " (", round(rLower, 4), "-", round(rUpper, 4), ")"),
         h = paste0(round(h, 2), " (", round(hLower, 2), "-", round(hUpper, 2), ")"),
         sd = paste0(round(sd, 4), " (", round(sdLower, 4), "-", round(sdUpper, 4), ")"),
         d = paste0(round(d, 3), " (", round(dLower, 3), "-", round(dUpper, 3), ")")) %>% 
  dplyr::select(Treatment, r, h, sd, d)

print(xtable(CIOut2023BTable, type = "latex"), include.rownames = FALSE)

dataTreatment2023BWWM <- lapply(1:length(dataTreatment2023BWW), function(i){
  dat <- dataTreatment2023BWW[[i]]
  treatment <- names(dataTreatment2023BWW)[i]
  dat <- dat %>% 
    mutate(Treatment = treatment)
})

dataTreatment2023BWWM <- lapply(dataTreatment2023BWWM, function(Data){
  Average <- aggregate(Data$Val, list(Data$Time), function(x)mean(x, na.rm = TRUE))
  Data$Average <- Average$x[match(Data$Time, Average$Group.1)]
  return(Data)
})

# Plot triple combination graphs
geomPoints2023BWW <- lapply(dataTreatment2023BWWM, function(d){
  return(geom_point(data = d, aes(x = Time, y = Average, col = as.factor(Treatment[1])), size = 3, shape = 1))
})

geomFunctions2023BWW <- lapply(1:length(fitsWWTLogisticNormal2023B), function(i){
  o <- fitsWWTLogisticNormal2023B[[i]]
  return(geom_function(fun = function(x)o$par["d"] * (exp(o$par["r"]*(x - o$par["h"])))/(1 + exp(o$par["r"]*(x - o$par["h"]))),
                       aes(col = as.factor(dataTreatment2023BWWM[[i]]$Treatment[1])),
                       size = 1.5))
})

values2023BWW <- colourVector[1:length(fitsWWTLogisticNormal2023B)]
names(values2023BWW) <- unname(sapply(dataTreatment2023BWWM, function(d)d$Treatment[1]))
labels2023BWW <- unname(sapply(dataTreatment2023BWWM, function(d)d$Treatment[1]))

plots2023BWW <- lapply(1:length(geomPoints2023BWW), function(i){
  plot <- ggplot() +
    geomPoints2023BWW[[i]] +
    geomFunctions2023BWW[[i]] +
    theme_classic(base_size = 25) +
    scale_color_manual(labels = labels2023BWW[i],
                       values = "#7FC97F") +
    labs(title = "",
         y = bquote(p[li]),
         x = bquote(T[c]),
         color = "Treatment") 
})

# Create a combined plot for intercrops and normal density sole crops
combinedPlot2023BWW <- plotMultipleGraphs(geomPoints2023BWW, geomFunctions2023BWW, names(dataTreatment2023BWW), c(colourVector[c(2, 9, 8, 4, 1, 5, 3, 7)]), XLab = bquote(T[c]), YLab = bquote(p[li]), Title = "")

combinedFaba2023BWW <- plotMultipleGraphs(geomPoints2023BWW[c(1, 2)], geomFunctions2023BWW[c(1, 2)], c("F", "F-375"), c(colourVector[c(2, 9)]), XLab = bquote(T[c]), YLab = bquote(p[li]), Title = "")
combinedTriticale2023BWW <- plotMultipleGraphs(geomPoints2023BWW[c(5, 6)], geomFunctions2023BWW[c(5, 6)], c("T", "T-375"), c(colourVector[c(1, 5)]), XLab = bquote(T[c]), YLab = bquote(p[li]), Title = "")
combinedIC2023BWW <- plotMultipleGraphs(geomPoints2023BWW[c(3, 4, 7, 8)], geomFunctions2023BWW[c(3, 4, 7, 8)], c("TF-M", "TF-M-375", "1T:1F", "1T:1F-375"), c(colourVector[c(8, 4, 3, 7)]), XLab = bquote(T[c]), YLab = bquote(p[li]), Title = "")
combined1252023BWW <- plotMultipleGraphs(geomPoints2023BWW[c(1, 3, 5, 7)], geomFunctions2023BWW[c(1, 3, 5, 7)], c("F", "TF-M", "T", "1T:1F"), c(colourVector[c(2, 8, 1, 3)]), XLab = bquote(T[c]), YLab = bquote(p[li]), Title = "")
combined3752023BWW <- plotMultipleGraphs(geomPoints2023BWW[c(2, 4, 6, 8)], geomFunctions2023BWW[c(2, 4, 6, 8)], c("F-375", "TF-M-375", "T-375", "1T:1F-375"), c(colourVector[c(9, 8, 5, 7)]), XLab = bquote(T[c]), YLab = bquote(p[li]), Title = "")

### Combine data for AIC comparisons
# Fit of 1T:1F
startParameters1T1F <- list(c("r" = 0.008, "h" = 520, "sd" = 0.1, "d" = 0.75),
                             c("r" = 0.008, "h" = 500, "sd" = 0.1, "d" = 0.8),
                             c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75),
                             c("r" = 0.008, "h" = 550, "sd" = 0.1, "d" = 0.75),
                             c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75))
combinationAIC1T1F <- round(fitCombinations(data = dataTreatment2023BWW, triple = c("1T:1F", "T", "F"), model = tLogisticNormal, startParameters = startParameters1T1F), 2)

# Fit of 1T:1F-375
startParameters1T1F375 <- list(c("r" = 0.008, "h" = 480, "sd" = 0.1, "d" = 0.75),
                             c("r" = 0.008, "h" = 500, "sd" = 0.1, "d" = 0.75),
                             c("r" = 0.008, "h" = 650, "sd" = 0.1, "d" = 0.8),
                             c("r" = 0.008, "h" = 550, "sd" = 0.1, "d" = 0.75),
                             c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75))
combinationAIC1T1F375 <- round(fitCombinations(data = dataTreatment2023BWW, triple = c("1T:1F-375", "T-375", "F-375"), model = tLogisticNormal, startParameters = startParameters1T1F375), 2)

# Fit of TF-M
startParametersTFM <- list(c("r" = 0.008, "h" = 520, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 500, "sd" = 0.1, "d" = 0.8),
                            c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 550, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75))
combinationAICTFM <- round(fitCombinations(data = dataTreatment2023BWW, triple = c("TF-M", "T", "F"), model = tLogisticNormal, startParameters = startParametersTFM), 2)

# Fit of TF-M-375
startParametersTFM375 <- list(c("r" = 0.01, "h" = 600, "sd" = 0.1, "d" = 0.8),
                            c("r" = 0.01, "h" = 600, "sd" = 0.1, "d" = 0.85),
                            c("r" = 0.01, "h" = 650, "sd" = 0.1, "d" = 0.9),
                            c("r" = 0.008, "h" = 550, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.01, "h" = 650, "sd" = 0.1, "d" = 0.8))
combinationAICTFM375 <- round(fitCombinations(data = dataTreatment2023BWW, triple = c("TF-M-375", "T-375", "F-375"), model = tLogisticNormal, startParameters = startParametersTFM375), 2)

combinationAIC2023BWW <- tibble(Treatment = rep(c("1T:1F", "1T:1F-375", "TF-M", "TF-M-375"), each = 3),
                                Combination = rep(c("IC & C & L", "IC + C & L", "IC + L & C"), times = 4),
                                AIC = c(combinationAIC1T1F,
                                        combinationAIC1T1F375,
                                        combinationAICTFM,
                                        combinationAICTFM375),
                                Best = c("*",
                                         "*",
                                         "",
                                         "*",
                                         "",
                                         "",
                                         "*",
                                         "*",
                                         "",
                                         "*",
                                         "",
                                         ""))

print(xtable(combinationAIC2023BWW, type = "latex"), include.rownames = FALSE)

######################
##                  ##
## Experiment 2024  ##
##                  ##
######################

sowingDate2024 <- as.Date("2023-04-12")

# Read data
data2024 <- read_xlsx("WCF_2024_data.xlsx", sheet = "LightInterception", range = "A1:E151", col_names = TRUE) %>% 
  mutate(DAS = as.numeric(as.Date(Date) - sowingDate2024)) %>% 
  rename(Time = DAS, Val = PLightInterception) %>% 
  dplyr::select(Plot, Treatment, Time, Weeds, Val)

data2024$Val[which(data2024$Val < 0.001)] <- 0.001
data2024$Val[which(data2024$Val > 0.999)] <- 0.999

tempData2024 <- read_delim("TempData2024.txt", delim = "\t") %>% 
  mutate(Time = as.integer(as.Date(Date) - sowingDate2024))

cumulativeTemp2024 <- ((tempData2024$TMin + tempData2024$TMax) / 2)
cumulativeTemp2024[which(cumulativeTemp2024 < 0.0)] <- 0
cumulativeTemp2024 <- cumsum(cumulativeTemp2024)

tempData2024 <- tempData2024 %>% 
  mutate(cumulativeTemp = cumulativeTemp2024)

dataWeeds2024 <- data2024 %>% 
  filter(Weeds == "Y") %>% 
  dplyr::select(-Weeds) %>% 
  mutate(Time = tempData2024$cumulativeTemp[match(Time, tempData2024$Time)])
dataNoWeeds2024 <- data2024 %>% 
  filter(Weeds == "N") %>% 
  dplyr::select(-Weeds) %>% 
  mutate(Time = tempData2024$cumulativeTemp[match(Time, tempData2024$Time)])

## Manual fits
plotFit <- function(data, par, fun){
  plot <- ggplot() +
    geom_point(data = data, aes(x = Time, y = Val), size = 3, shape = 1) +
    geom_function(fun = function(x)allDeterministicFunctions[[fun]](par, x),
                  size = 1.5)
  return(plot)
}


##################
##################
####          ####
#### No weeds ####
####          ####
##################
##################

dataTreatment2024NW <- split(dataNoWeeds2024, dataNoWeeds2024$Treatment)
dataTreatment2024NW <- lapply(dataTreatment2024NW, function(dat){
  datNew <- dat %>% 
    dplyr::select(Time, Val)
  return(datNew)
})




##################################
###                            ###
### Fit logistic normal curves ###
###                            ###
##################################     

fitsNWLogisticNormal2024 <- rep(list(NA), 3)

# Fit of 1T:1F
fitsNWLogisticNormal2024[[1]] <- fitDataToModel(data = dataTreatment2024NW[[1]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 450, "sd" = 0.1))
plotFit(dataTreatment2024NW[[1]], fitsNWLogisticNormal2024[[1]]$par, "logistic")

# Fit of T
fitsNWLogisticNormal2024[[2]] <- fitDataToModel(data = dataTreatment2024NW[[2]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 520, "sd" = 0.1))
plotFit(dataTreatment2024NW[[2]], fitsNWLogisticNormal2024[[2]]$par, "logistic")

# Fit of F
fitsNWLogisticNormal2024[[3]] <- fitDataToModel(data = dataTreatment2024NW[[3]], model = logisticNormal, startParameters = c("r" = 0.012, "h" = 475, "sd" = 0.1))
plotFit(dataTreatment2024NW[[3]], fitsNWLogisticNormal2024[[3]]$par, "logistic")

AICLogisticNormal2024NW <- sum(unlist(lapply(fitsNWLogisticNormal2024, function(fits)2 * 3 + 2 * fits$value)))


##############################################
###                                        ###
### Fit transformed logistic normal curves ###
###                                        ###
##############################################      

fitsNWTLogisticNormal2024 <- rep(list(NA), 3)

# Fit of 1T:1F
fitsNWTLogisticNormal2024[[1]] <- fitDataToModel(data = dataTreatment2024NW[[1]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 480, "sd" = 0.1, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2024NW[[1]], fitsNWTLogisticNormal2024[[1]]$par, "tLogistic")

# Fit of T
fitsNWTLogisticNormal2024[[2]] <- fitDataToModel(data = dataTreatment2024NW[[2]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 550, "sd" = 0.1, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2024NW[[2]], fitsNWTLogisticNormal2024[[2]]$par, "tLogistic")

# Fit of F
fitsNWTLogisticNormal2024[[3]] <- fitDataToModel(data = dataTreatment2024NW[[3]], model = tLogisticNormal, startParameters = c("r" = 0.012, "h" = 350, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2024NW[[3]], fitsNWTLogisticNormal2024[[3]]$par, "tLogistic")

AICTLogisticNormal2024NW <- sum(unlist(lapply(fitsNWTLogisticNormal2024, function(fits)2 * 4 + 2 * fits$value)))


#################################
###                           ###
### Fit logistic gamma curves ###
###                           ###
#################################     

fitsNWLogisticGamma2024 <- rep(list(NA), 3)

# Fit of 1T:1F
fitsNWLogisticGamma2024[[1]] <- fitDataToModel(data = dataTreatment2024NW[[1]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 520, "shape" = 1.0))
plotFit(dataTreatment2024NW[[1]], fitsNWLogisticGamma2024[[1]]$par, "logistic")

# Fit of T
fitsNWLogisticGamma2024[[2]] <- fitDataToModel(data = dataTreatment2024NW[[2]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 520, "shape" = 1.0))
plotFit(dataTreatment2024NW[[2]], fitsNWLogisticGamma2024[[2]]$par, "logistic")

# Fit of F
fitsNWLogisticGamma2024[[3]] <- fitDataToModel(data = dataTreatment2024NW[[3]], model = logisticGamma, startParameters = c("r" = 0.012, "h" = 475, "shape" = 1.0))
plotFit(dataTreatment2024NW[[3]], fitsNWLogisticGamma2024[[3]]$par, "logistic")

AICLogisticGamma2024NW <- sum(unlist(lapply(fitsNWLogisticGamma2024, function(fits)2 * 3 + 2 * fits$value)))


#############################################
###                                       ###
### Fit transformed logistic gamma curves ###
###                                       ###
#############################################      

fitsNWTLogisticGamma2024 <- rep(list(NA), 3)

# Fit of 1T:1F
fitsNWTLogisticGamma2024[[1]] <- fitDataToModel(data = dataTreatment2024NW[[1]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 420, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2024NW[[1]], fitsNWTLogisticGamma2024[[1]]$par, "tLogistic")

# Fit of T
fitsNWTLogisticGamma2024[[2]] <- fitDataToModel(data = dataTreatment2024NW[[2]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 520, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2024NW[[2]], fitsNWTLogisticGamma2024[[2]]$par, "tLogistic")

# Fit of F
fitsNWTLogisticGamma2024[[3]] <- fitDataToModel(data = dataTreatment2024NW[[3]], model = tLogisticGamma, startParameters = c("r" = 0.012, "h" = 350, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2024NW[[3]], fitsNWTLogisticGamma2024[[3]]$par, "tLogistic")

AICTLogisticGamma2024NW <- sum(unlist(lapply(fitsNWTLogisticGamma2024, function(fits)2 * 4 + 2 * fits$value)))

AICLogisticNormal2024NW
AICTLogisticNormal2024NW
AICLogisticGamma2024NW
AICTLogisticGamma2024NW

# Transformed logistic normal is the best fitting model

# Calculate 95% CIs
parameters2024 <- names(fitsNWTLogisticNormal2024[[1]]$par)
parameterCIs2024 <- rep(NA, times = 3 * length(parameters2024))
names(parameterCIs2024) <- unlist(lapply(parameters2024, function(par)c(par, paste(par, "Lower", sep = ""), paste(par, "Upper", sep = ""))))
parameterCIs2024 <- rep(list(parameterCIs2024), length(fitsNWTLogisticNormal2024))

cI2024 <- function(par, x, z, opt, dMax){
  pars <- par
  pars[names(opt)] <- opt
  nll <- tLogisticNormal(pars, x, z, dMax)
  return(nll)
}

for(i in 1:length(fitsNWTLogisticNormal2024)){
  parameterVecs <- lapply(fitsNWTLogisticNormal2024, function(fit){
    parStart <- fit$par / 100
    parEnd <- fit$par * 20
    parVec <- list(seq(parStart[1], parEnd[1], length.out = 400), seq(parStart[2], parEnd[2], length.out = 400), seq(parStart[3], parEnd[3], length.out = 400), seq(parStart[4], parEnd[4], length.out = 400))
    names(parVec) <- names(fit$par)
    return(parVec)
  })
  optPar <- fitsNWTLogisticNormal2024[[i]]$par
  parameterCIs2024[[i]][names(optPar)] <- optPar
  parameterCIs2024[[i]][unlist((lapply(names(optPar), function(n)c(paste(n, "Lower", sep = ""), paste(n, "Upper", sep = "")))))] <-
    unlist(lapply(1:length(parameterVecs[[i]]), function(j){
      t(tibble(calculate95CI(parStartVec = optPar[-which(names(optPar) == names(fitsNWTLogisticNormal2024[[i]]$par)[[j]])],
                             optVec = parameterVecs[[i]][j],
                             model = cI2024,
                             optNLL = fitsNWTLogisticNormal2024[[i]]$value,
                             x = dataTreatment2024NW[[i]][,1][[1]],
                             z = dataTreatment2024NW[[i]][,2][[1]],
                             dMax = 1.0)))
    }))
  
}
names(parameterCIs2024) <- names(dataTreatment2024NW)

# Write to table
CIOut2024NW <- tibble(Treatment = names(parameterCIs2024)) %>% 
  bind_cols(bind_rows(parameterCIs2024))
write.table(CIOut2024NW, "light_interception/light_interception_NW_tLogisticNormal_parameters_2024.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

CIOut2024Table <- CIOut2024NW %>% 
  mutate(r = paste0(round(r, 4), " (", round(rLower, 4), "-", round(rUpper, 4), ")"),
         h = paste0(round(h, 2), " (", round(hLower, 2), "-", round(hUpper, 2), ")"),
         sd = paste0(round(sd, 4), " (", round(sdLower, 4), "-", round(sdUpper, 4), ")"),
         d = paste0(round(d, 3), " (", round(dLower, 3), "-", round(dUpper, 3), ")")) %>% 
  dplyr::select(Treatment, r, h, sd, d)

print(xtable(CIOut2024Table, type = "latex"), include.rownames = FALSE)

dataTreatment2024NWM <- lapply(1:length(dataTreatment2024NW), function(i){
  dat <- dataTreatment2024NW[[i]]
  treatment <- names(dataTreatment2024NW)[i]
  dat <- dat %>% 
    mutate(Treatment = treatment)
})

dataTreatment2024NWM <- lapply(dataTreatment2024NWM, function(Data){
  Average <- aggregate(Data$Val, list(Data$Time), function(x)mean(x, na.rm = TRUE))
  Data$Average <- Average$x[match(Data$Time, Average$Group.1)]
  return(Data)
})

# Plot triple combination graphs
geomPoints2024NW <- lapply(dataTreatment2024NWM, function(d){
  return(geom_point(data = d, aes(x = Time, y = Average, col = as.factor(Treatment[1])), size = 3, shape = 1))
})

geomFunctions2024NW <- lapply(1:length(fitsNWTLogisticNormal2024), function(i){
  o <- fitsNWTLogisticNormal2024[[i]]
  return(geom_function(fun = function(x)o$par["d"] * (exp(o$par["r"]*(x - o$par["h"])))/(1 + exp(o$par["r"]*(x - o$par["h"]))),
                       aes(col = as.factor(dataTreatment2024NWM[[i]]$Treatment[1])),
                       size = 1.5))
})


### Get all triple combinations
treatmentsToCompare <- unique(data2024$Treatment)
triples <- list(c("1T:1F", "T", "F"))
names(dataTreatment2024NWM) <- lapply(dataTreatment2024NWM, function(l)l$Treatment[1])
names(geomPoints2024NW) <- names(dataTreatment2024NWM)
names(geomFunctions2024NW) <- names(dataTreatment2024NWM)

triplePlots2024NW <- lapply(triples, function(triple){
  plotMultipleGraphs(geomPoints2024NW[triple], geomFunctions2024NW[triple], triple, colourVector[c(3, 1, 2)], XLab = bquote(T[c]), YLab = bquote(p[li]), Title = "")
})

### Combine data for AIC comparisons
# Fit of 1T:1F
startParameters1T1F <- list(c("r" = 0.01, "h" = 480, "sd" = 0.1, "d" = 0.9),
                            c("r" = 0.01, "h" = 350, "sd" = 0.1, "d" = 0.9),
                            c("r" = 0.012, "h" = 550, "sd" = 0.1, "d" = 0.85),
                            c("r" = 0.01, "h" = 420, "sd" = 0.1, "d" = 0.9),
                            c("r" = 0.01, "h" = 480, "sd" = 0.1, "d" = 0.9))
combinationAIC1T1F <- round(fitCombinations(data = dataTreatment2024NW, triple = triples[[1]], model = tLogisticNormal, startParameters = startParameters1T1F), 2)

combinationAIC2024NW <- tibble(Treatment = rep(unique(data2024$Treatment)[c(1)], each = 3),
                                Combination = rep(c("IC & C & L", "IC + C & L", "IC + L & C"), times = 1),
                                AIC = c(combinationAIC1T1F),
                                Best = c("*",
                                         "",
                                         ""))

print(xtable(combinationAIC2024NW, type = "latex"), include.rownames = FALSE)

# # Print all plots to PDF
# pdf("light_interception/light_interception_NW_tLogisticNormal_graphs_2024_NoWeeds.pdf")
# print(triplePlotsNW)
# dev.off()

# png("figures/light_interception/plotLightInterception2024Triple1T1FNW.png", units = "px", width = 4000, height = 3200, res = 300)
# triplePlotsNW[[1]] + labs(subtitle = "2024, herbicide") + ylim(c(0, 1.0))
# dev.off()

###################
####################
####            ####
#### With weeds ####
####            ####
####################
####################

dataTreatment2024WW <- split(dataWeeds2024, dataWeeds2024$Treatment)
dataTreatment2024WW <- lapply(dataTreatment2024WW, function(dat){
  datNew <- dat %>% 
    dplyr::select(Time, Val)
  return(datNew)
})




##################################
###                            ###
### Fit logistic normal curves ###
###                            ###
##################################     

fitsWWLogisticNormal2024 <- rep(list(NA), 3)

# Fit of 1T:1F
fitsWWLogisticNormal2024[[1]] <- fitDataToModel(data = dataTreatment2024WW[[1]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 450, "sd" = 0.1))
plotFit(dataTreatment2024WW[[1]], fitsWWLogisticNormal2024[[1]]$par, "logistic")

# Fit of T
fitsWWLogisticNormal2024[[2]] <- fitDataToModel(data = dataTreatment2024WW[[2]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 520, "sd" = 0.1))
plotFit(dataTreatment2024WW[[2]], fitsWWLogisticNormal2024[[2]]$par, "logistic")

# Fit of F
fitsWWLogisticNormal2024[[3]] <- fitDataToModel(data = dataTreatment2024WW[[3]], model = logisticNormal, startParameters = c("r" = 0.012, "h" = 475, "sd" = 0.1))
plotFit(dataTreatment2024WW[[3]], fitsWWLogisticNormal2024[[3]]$par, "logistic")

AICLogisticNormal2024WW <- sum(unlist(lapply(fitsWWLogisticNormal2024, function(fits)2 * 3 + 2 * fits$value)))


##############################################
###                                        ###
### Fit transformed logistic normal curves ###
###                                        ###
##############################################      

fitsWWTLogisticNormal2024 <- rep(list(NA), 3)

# Fit of 1T:1F
fitsWWTLogisticNormal2024[[1]] <- fitDataToModel(data = dataTreatment2024WW[[1]], model = tLogisticNormal, startParameters = c("r" = 0.02, "h" = 450, "sd" = 0.1, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2024WW[[1]], fitsWWTLogisticNormal2024[[1]]$par, "tLogistic")

# Fit of T
fitsWWTLogisticNormal2024[[2]] <- fitDataToModel(data = dataTreatment2024WW[[2]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 550, "sd" = 0.1, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2024WW[[2]], fitsWWTLogisticNormal2024[[2]]$par, "tLogistic")

# Fit of F
fitsWWTLogisticNormal2024[[3]] <- fitDataToModel(data = dataTreatment2024WW[[3]], model = tLogisticNormal, startParameters = c("r" = 0.012, "h" = 350, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2024WW[[3]], fitsWWTLogisticNormal2024[[3]]$par, "tLogistic")

AICTLogisticNormal2024WW <- sum(unlist(lapply(fitsWWTLogisticNormal2024, function(fits)2 * 4 + 2 * fits$value)))


#################################
###                           ###
### Fit logistic gamma curves ###
###                           ###
#################################     

fitsWWLogisticGamma2024 <- rep(list(NA), 3)

# Fit of 1T:1F
fitsWWLogisticGamma2024[[1]] <- fitDataToModel(data = dataTreatment2024WW[[1]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 520, "shape" = 1.0))
plotFit(dataTreatment2024WW[[1]], fitsWWLogisticGamma2024[[1]]$par, "logistic")

# Fit of T
fitsWWLogisticGamma2024[[2]] <- fitDataToModel(data = dataTreatment2024WW[[2]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 520, "shape" = 1.0))
plotFit(dataTreatment2024WW[[2]], fitsWWLogisticGamma2024[[2]]$par, "logistic")

# Fit of F
fitsWWLogisticGamma2024[[3]] <- fitDataToModel(data = dataTreatment2024WW[[3]], model = logisticGamma, startParameters = c("r" = 0.012, "h" = 475, "shape" = 1.0))
plotFit(dataTreatment2024WW[[3]], fitsWWLogisticGamma2024[[3]]$par, "logistic")

AICLogisticGamma2024WW <- sum(unlist(lapply(fitsWWLogisticGamma2024, function(fits)2 * 3 + 2 * fits$value)))


#############################################
###                                       ###
### Fit transformed logistic gamma curves ###
###                                       ###
#############################################      

fitsWWTLogisticGamma2024 <- rep(list(NA), 3)

# Fit of 1T:1F
fitsWWTLogisticGamma2024[[1]] <- fitDataToModel(data = dataTreatment2024WW[[1]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 420, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2024WW[[1]], fitsWWTLogisticGamma2024[[1]]$par, "tLogistic")

# Fit of T
fitsWWTLogisticGamma2024[[2]] <- fitDataToModel(data = dataTreatment2024WW[[2]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 520, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2024WW[[2]], fitsWWTLogisticGamma2024[[2]]$par, "tLogistic")

# Fit of F
fitsWWTLogisticGamma2024[[3]] <- fitDataToModel(data = dataTreatment2024WW[[3]], model = tLogisticGamma, startParameters = c("r" = 0.012, "h" = 350, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2024WW[[3]], fitsWWTLogisticGamma2024[[3]]$par, "tLogistic")

AICTLogisticGamma2024WW <- sum(unlist(lapply(fitsWWTLogisticGamma2024, function(fits)2 * 4 + 2 * fits$value)))

AICLogisticNormal2024WW
AICTLogisticNormal2024WW
AICLogisticGamma2024WW
AICTLogisticGamma2024WW

# Transformed logistic normal is the best fitting model

# Calculate 95% CIs
parameters2024 <- names(fitsWWTLogisticNormal2024[[1]]$par)
parameterCIs2024 <- rep(NA, times = 3 * length(parameters2024))
names(parameterCIs2024) <- unlist(lapply(parameters2024, function(par)c(par, paste(par, "Lower", sep = ""), paste(par, "Upper", sep = ""))))
parameterCIs2024 <- rep(list(parameterCIs2024), length(fitsWWTLogisticNormal2024))

cI2024 <- function(par, x, z, opt, dMax){
  pars <- par
  pars[names(opt)] <- opt
  nll <- tLogisticNormal(pars, x, z, dMax)
  return(nll)
}

for(i in 1:length(fitsWWTLogisticNormal2024)){
  parameterVecs <- lapply(fitsWWTLogisticNormal2024, function(fit){
    parStart <- fit$par / 100
    parEnd <- fit$par * 20
    parVec <- list(seq(parStart[1], parEnd[1], length.out = 400), seq(parStart[2], parEnd[2], length.out = 400), seq(parStart[3], parEnd[3], length.out = 400), seq(parStart[4], parEnd[4], length.out = 400))
    names(parVec) <- names(fit$par)
    return(parVec)
  })
  optPar <- fitsWWTLogisticNormal2024[[i]]$par
  parameterCIs2024[[i]][names(optPar)] <- optPar
  parameterCIs2024[[i]][unlist((lapply(names(optPar), function(n)c(paste(n, "Lower", sep = ""), paste(n, "Upper", sep = "")))))] <-
    unlist(lapply(1:length(parameterVecs[[i]]), function(j){
      t(tibble(calculate95CI(parStartVec = optPar[-which(names(optPar) == names(fitsWWTLogisticNormal2024[[i]]$par)[[j]])],
                             optVec = parameterVecs[[i]][j],
                             model = cI2024,
                             optNLL = fitsWWTLogisticNormal2024[[i]]$value,
                             x = dataTreatment2024WW[[i]][,1][[1]],
                             z = dataTreatment2024WW[[i]][,2][[1]],
                             dMax = 1.0)))
    }))
  
}
names(parameterCIs2024) <- names(dataTreatment2024WW)

# Write to table
CIOut2024WW <- tibble(Treatment = names(parameterCIs2024)) %>% 
  bind_cols(bind_rows(parameterCIs2024))
write.table(CIOut2024WW, "light_interception/light_interception_WW_tLogisticNormal_parameters_2024.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

CIOut2024Table <- CIOut2024WW %>% 
  mutate(r = paste0(round(r, 4), " (", round(rLower, 4), "-", round(rUpper, 4), ")"),
         h = paste0(round(h, 2), " (", round(hLower, 2), "-", round(hUpper, 2), ")"),
         sd = paste0(round(sd, 4), " (", round(sdLower, 4), "-", round(sdUpper, 4), ")"),
         d = paste0(round(d, 3), " (", round(dLower, 3), "-", round(dUpper, 3), ")")) %>% 
  dplyr::select(Treatment, r, h, sd, d)

print(xtable(CIOut2024Table, type = "latex"), include.rownames = FALSE)

CIOut2024NW <- CIOut2024NW %>% 
  mutate(Treatment = paste0(Treatment, "_H"))
CIOut2024WW <- CIOut2024WW %>% 
  mutate(Treatment = paste0(Treatment, "_W"))
  
CIOut2024Table <- bind_rows(CIOut2024NW, CIOut2024WW) %>%
  mutate(r = paste0(round(r, 4), " (", round(rLower, 4), "-", round(rUpper, 4), ")"),
         h = paste0(round(h, 2), " (", round(hLower, 2), "-", round(hUpper, 2), ")"),
         sd = paste0(round(sd, 4), " (", round(sdLower, 4), "-", round(sdUpper, 4), ")"),
         d = paste0(round(d, 3), " (", round(dLower, 3), "-", round(dUpper, 3), ")")) %>% 
  dplyr::select(Treatment, r, h, sd, d)

print(xtable(CIOut2024Table, type = "latex"), include.rownames = FALSE)


dataTreatment2024WWM <- lapply(1:length(dataTreatment2024WW), function(i){
  dat <- dataTreatment2024WW[[i]]
  treatment <- names(dataTreatment2024WW)[i]
  dat <- dat %>% 
    mutate(Treatment = treatment)
})

dataTreatment2024WWM <- lapply(dataTreatment2024WWM, function(Data){
  Average <- aggregate(Data$Val, list(Data$Time), function(x)mean(x, na.rm = TRUE))
  Data$Average <- Average$x[match(Data$Time, Average$Group.1)]
  return(Data)
})

# Plot triple combination graphs
geomPoints2024WW <- lapply(dataTreatment2024WWM, function(d){
  return(geom_point(data = d, aes(x = Time, y = Average, col = as.factor(Treatment[1])), size = 3, shape = 1))
})

geomFunctions2024WW <- lapply(1:length(fitsWWTLogisticNormal2024), function(i){
  o <- fitsWWTLogisticNormal2024[[i]]
  return(geom_function(fun = function(x)o$par["d"] * (exp(o$par["r"]*(x - o$par["h"])))/(1 + exp(o$par["r"]*(x - o$par["h"]))),
                       aes(col = as.factor(dataTreatment2024WWM[[i]]$Treatment[1])),
                       size = 1.5))
})

### Get all triple combinations
treatmentsToCompare <- unique(data2024$Treatment)
triples <- list(c("1T:1F", "T", "F"))
names(dataTreatment2024WWM) <- lapply(dataTreatment2024WWM, function(l)l$Treatment[1])
names(geomPoints2024WW) <- names(dataTreatment2024WWM)
names(geomFunctions2024WW) <- names(dataTreatment2024WWM)

triplePlots2024WW <- lapply(triples, function(triple){
  plotMultipleGraphs(geomPoints2024WW[triple], geomFunctions2024WW[triple], triple, colourVector[c(3, 1, 2)], XLab = bquote(T[c]), YLab = bquote(p[li]), Title = "")
})

### Combine data for AIC comparisons
# Fit of 1T:1F
startParameters1T1F <- list(c("r" = 0.01, "h" = 480, "sd" = 0.1, "d" = 0.9),
                            c("r" = 0.01, "h" = 350, "sd" = 0.1, "d" = 0.9),
                            c("r" = 0.012, "h" = 550, "sd" = 0.1, "d" = 0.85),
                            c("r" = 0.01, "h" = 420, "sd" = 0.1, "d" = 0.9),
                            c("r" = 0.01, "h" = 480, "sd" = 0.1, "d" = 0.9))
combinationAIC1T1F <- round(fitCombinations(data = dataTreatment2024WW, triple = triples[[1]], model = tLogisticNormal, startParameters = startParameters1T1F), 2)

combinationAIC2024WW <- tibble(Treatment = rep(unique(data2024$Treatment)[c(1)], each = 3),
                               Combination = rep(c("IC & C & L", "IC + C & L", "IC + L & C"), times = 1),
                               AIC = c(combinationAIC1T1F),
                               Best = c("*",
                                        "",
                                        ""))

print(xtable(combinationAIC2024WW, type = "latex"), include.rownames = FALSE)

png("Figures/plot_light_interception_2023-RD_2024.png", units = "px", width = 5000, height = 8000, res = 300)
combined1252023BWW + labs(subtitle = "") + 
  ylim(c(0, 1.0)) + 
  theme(axis.title.y = element_text(size = 35),
        axis.title.x = element_text(size = 35), 
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30),
        legend.position = "none") + 
  labs(subtitle = "2023-RD, 12.5 cm RD, weed-infested") + 
combined1252023BNW + labs(subtitle = "") + 
  ylim(c(0, 1.0)) + 
  theme(axis.title.y = element_text(size = 35),
        axis.title.x = element_text(size = 35),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), 
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30)) + 
  labs(subtitle = "2023-RD, 12.5 cm RD, herbicide-treated") + 
combined3752023BWW + labs(subtitle = "") + 
  ylim(c(0, 1.0)) + 
  labs(subtitle = "2023-RD, 37.5 cm RD, weed-infested") + 
  theme(axis.title.y = element_text(size = 35),
        axis.title.x = element_text(size = 35), 
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30),
        legend.position = "none") +
combined3752023BNW + labs(subtitle = "") + 
  ylim(c(0, 1.0)) + 
  theme(axis.title.y = element_text(size = 35),
        axis.title.x = element_text(size = 35),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), 
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30)) + 
  labs(subtitle = "2023-RD, 37.5 cm RD, herbicide-treated") + 
triplePlotsWW[[1]] + labs(subtitle = "2024, weed-infested", x = bquote(T[c])) + theme(legend.position = "none") + ylim(c(0, 1.0)) +
  theme(axis.title.y = element_text(size = 35),
        axis.title.x = element_text(size = 35), 
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30)) + 
triplePlotsNW[[1]] + labs(subtitle = "2024, herbicide-treated", x = bquote(T[c])) + 
  theme(axis.title.y = element_text(size = 35),
        axis.title.x = element_text(size = 35),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank(), 
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30)) + 
plot_layout(ncol = 2, nrow = 3, axis_titles = "collect") +
plot_annotation(tag_levels = "a",
                theme = theme(plot.title = element_text(size = 35)))
dev.off()






AICsTable <- tibble(Year = c(rep("2022", 4), rep("2023A", 4), rep("2023B", 4), rep("2024", 4)),
                   Model = rep(c("Logistic normal",
                             "Transformed logistic normal",
                             "Logistic gamma",
                             "Transformed logistic gamma"), times = 4),
                   AIC_H = c(AICLogisticNormal2022,
                           AICTLogisticNormal2022,
                           AICLogisticGamma2022,
                           AICTLogisticGamma2022,
                           AICLogisticNormal2023ANW,
                           AICTLogisticNormal2023ANW,
                           AICLogisticGamma2023ANW,
                           AICTLogisticGamma2023AWW,
                           AICLogisticNormal2023BNW,
                           AICTLogisticNormal2023BNW,
                           AICLogisticGamma2023BNW,
                           AICTLogisticGamma2023BNW,
                           AICLogisticNormal2024NW,
                           AICTLogisticNormal2024NW,
                           AICLogisticGamma2024NW,
                           AICTLogisticGamma2024NW),
                   AIC_W = c(NA,
                             NA,
                             NA,
                             NA,
                             AICLogisticNormal2023AWW,
                             AICTLogisticNormal2023AWW,
                             AICLogisticGamma2023AWW,
                             AICTLogisticGamma2023AWW,
                             AICLogisticNormal2023BWW,
                             AICTLogisticNormal2023BWW,
                             AICLogisticGamma2023BWW,
                             AICTLogisticGamma2023BWW,
                             AICLogisticNormal2024WW,
                             AICTLogisticNormal2024WW,
                             AICLogisticGamma2024WW,
                             AICTLogisticGamma2024WW)) %>% 
  mutate(AIC_H = round(AIC_H, 2),
         AIC_W = round(AIC_W, 2))
print(xtable(AICsTable, type = "latex"), include.rownames = FALSE)


combinationAICAll <- bind_rows(combinationAIC2022,
                               combinationAIC2023AWW,
                               combinationAIC2023ANW,
                               combinationAIC2023BWW,
                               combinationAIC2023BNW,
                               combinationAIC2024WW,
                               combinationAIC2024NW) %>% 
  mutate(Experiment = c(rep("2022", nrow(combinationAIC2022)),
                      rep("2023A", nrow(combinationAIC2023AWW) * 2),
                      rep("2023B", nrow(combinationAIC2023BWW) * 2),
                      rep("2024", nrow(combinationAIC2024WW) * 2)),
         HT = c(rep("H", nrow(combinationAIC2022)),
                rep("W", nrow(combinationAIC2023AWW)),
                rep("H", nrow(combinationAIC2023ANW)),
                rep("W", nrow(combinationAIC2023BWW)),
                rep("H", nrow(combinationAIC2023BNW)),
                rep("W", nrow(combinationAIC2024WW)),
                rep("H", nrow(combinationAIC2024NW)))) %>% 
  dplyr::select(Experiment, Treatment, HT, Combination, AIC) %>% 
  pivot_wider(values_from = AIC, names_from = HT) %>% 
  rename(AIC_H = H, AIC_W = W)

print(xtable(combinationAICAll, type = "latex"), include.rownames = FALSE)

CICDiff <- tibble(Experiment = c(rep("2023A", 4),
                                 rep("2023B", 4),
                                 "2024"),
                  Treatment = c("1T:1F",
                                "1T:3F",
                                "3T:1F",
                                "TF-M",
                                "1T:1F",
                                "1T:1F-375",
                                "TF-M",
                                "TF-M-375",
                                "1T:1F"),
                  Difference_H = c(filter(combinationAICAll, Combination == "IC & C & L")$AIC_H - filter(combinationAICAll, Combination == "IC + C & L")$AIC_H)[9:17],
                  Difference_W = c(filter(combinationAICAll, Combination == "IC & C & L")$AIC_W - filter(combinationAICAll, Combination == "IC + C & L")$AIC_W)[9:17],
)
print(xtable(CICDiff, type = "latex"), include.rownames = FALSE)


### New figure combinations ###
windowsFonts(LiberationMono = "Liberation Mono")
## Fig. 1: a) 2022 cereal and legume sole crop light interception; b) combination with intercrops
species_order <- c("Rye","Barley","Triticale","Wheat","Pea","Faba")
points_df_2022 <- bind_rows(
  lapply(dataTreatment2022, function(df) {
    df %>% mutate(Species = Treatment[1])
  })
) %>% 
  filter(Species %in% species_order) %>% 
  mutate(Species = factor(Species, levels = species_order))

time_seq <- seq(min(points_df_2022$Time), filter(tempData2022, Date == as.Date("2022-06-21"))$CumulTemp + 10, length.out = 200)

functions_df_2022 <- map2_dfr(fitsTLogisticNormal2022, dataTreatment2022, function(fit, df) {
  sp <- df$Treatment[1]
  d <- fit$par["d"]; r <- fit$par["r"]; h <- fit$par["h"]
  tibble(
    Time      = time_seq,
    Predicted = d * exp(r * (Time - h)) / (1 + exp(r * (Time - h))),
    Species   = sp
  )
}) %>% 
  filter(Species %in% species_order) %>% 
  mutate(Species = factor(Species, levels = species_order))

colours_2022 <- setNames(colourVector[c(1,2,3,4,5,7)],
                         c("Rye","Barley","Triticale","Wheat","Pea","Faba"))
shapes_2022  <- c(Rye=16, Barley=17, Triticale=15, Wheat=18, Pea=8, Faba=3)

weed_biomass <- c(
  "Rye       :  24.0",
  "Barley    :  46.8",
  "Triticale :  49.4",
  "Wheat     : 179.4",
  "Pea       : 171.8",
  "Faba      : 453.8"
  )

custom_labels <- paste0(names(weed_biomass), ": ", weed_biomass)

filter(tempData2022, Date == as.Date("2022-06-21"))$CumulTemp

plot_sole_2022 <- ggplot() +
  geom_point(data = points_df_2022,
             aes(Time, Average, colour = Species, shape = Species),
             size = 4) +
  geom_line(data = functions_df_2022,
            aes(Time, Predicted, colour = Species),
            size = 1.2) +
  scale_colour_manual(
    name   = "Crop species",
    values = colours_2022,
    labels = custom_labels
  ) +
  scale_shape_manual(
    name   = "Crop species",
    values = shapes_2022,
    labels = custom_labels
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(shape    = unname(shapes_2022),
                          linetype = 0)
    ),
    shape = "none"
  ) +
  theme_classic(base_size = 25) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.05),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = alpha("white", 1.0), color = NA),
    legend.key.size = unit(2.0, "lines"),
    legend.text = element_text(family = "LiberationMono", size = 20),
    legend.title = element_text(size = 22)
  ) +
  labs(x = "Cumulative daily average temperature (°Cd)",
       y = "Proportion PAR intercepted",
       title = "") +
  geom_vline(xintercept = filter(tempData2022, Date == as.Date("2022-06-21"))$CumulTemp,
             linetype = "dashed", colour = "black", size = 1)

plot_sole_2022

plot_intercrops_2022 <- (
  (triplePlots2022[[3]] + 
     labs(subtitle = "Rye-Pea", 
       x = "Cumulative daily average temperature (°Cd)",
       y = "Proportion PAR intercepted") + 
     ylim(0, 1) +
     theme(legend.position = "none",
           axis.title.x = element_blank(),
           axis.text.x  = element_blank(),
           axis.ticks.x = element_blank())) +

  (triplePlots2022[[4]] + 
     labs(subtitle = "Rye-Faba", 
       x = "Cumulative daily average temperature (°Cd)",
       y = "Proportion PAR intercepted") + 
     ylim(0, 1) +
     theme(legend.position = "none",
           axis.title.x = element_blank(),
           axis.text.x  = element_blank(),
           axis.ticks.x = element_blank(),
           axis.title.y = element_blank(),
           axis.text.y  = element_blank(),
           axis.ticks.y = element_blank())) +

  (triplePlots2022[[1]] + 
     labs(subtitle = "Barley-Pea", 
       x = "Cumulative daily average temperature (°Cd)",
       y = "Proportion PAR intercepted") + 
     ylim(0, 1) +
     theme(legend.position = "none",
           axis.title.x = element_blank(),
           axis.text.x  = element_blank(),
           axis.ticks.x = element_blank())) +

  (triplePlots2022[[2]] + 
     labs(subtitle = "Barley-Faba", 
       x = "Cumulative daily average temperature (°Cd)",
       y = "Proportion PAR intercepted") + 
     ylim(0, 1) +
     theme(legend.position = "none",
           axis.title.x = element_blank(),
           axis.text.x  = element_blank(),
           axis.ticks.x = element_blank(),
           axis.title.y = element_blank(),
           axis.text.y  = element_blank(),
           axis.ticks.y = element_blank())) +

  (triplePlots2022[[5]] + 
     labs(subtitle = "Triticale-Pea", 
       x = "Cumulative daily average temperature (°Cd)",
       y = "Proportion PAR intercepted") + 
     ylim(0, 1) +
     theme(legend.position = "none",
           axis.title.x = element_blank(),
           axis.text.x  = element_blank(),
           axis.ticks.x = element_blank())) +

  (triplePlots2022[[6]] + 
     labs(subtitle = "Triticale-Faba", 
       x = "Cumulative daily average temperature (°Cd)",
       y = "Proportion PAR intercepted") + 
     ylim(0, 1) +
     theme(legend.position = "none",
           axis.title.x = element_blank(),
           axis.text.x  = element_blank(),
           axis.ticks.x = element_blank(),
           axis.title.y = element_blank(),
           axis.text.y  = element_blank(),
           axis.ticks.y = element_blank())) +

  (triplePlots2022[[7]] + 
     labs(subtitle = "Wheat-Pea", 
       x = "Cumulative daily average temperature (°Cd)",
       y = "Proportion PAR intercepted") + 
     ylim(0, 1) +
     theme(legend.position = "none")) +

  (triplePlots2022[[8]] + 
     labs(subtitle = "Wheat-Faba", 
       x = "Cumulative daily average temperature (°Cd)",
       y = "Proportion PAR intercepted") +
     ylim(0, 1) +
     scale_color_manual(values = modelColors,
                        labels = c("Intercrop", "Sole cereal", "Sole legume")) +
     theme(axis.title.y = element_blank(),
           axis.text.y  = element_blank(),
           axis.ticks.y = element_blank()))
) +
  plot_layout(ncol = 2, axis_titles = "collect")

combined_plot_2022 <- (
  plot_sole_2022 | plot_intercrops_2022
) +
  plot_layout(widths = c(1, 0.6),
              axis_titles = "collect") +
  plot_annotation(
    tag_levels = "a"
  )

combined_plot_2022

## Fig. 3: 2023-SP intercrop light interception
## a: sole crops
treatment_order_2023A_sole <- c("T+", "T", "F", "F+")
dataTreatment2023ANWM2 <- lapply(dataTreatment2023ANWM, function(df){
  return(df %>% mutate(Treatment = gsub("TF", "1T:1F", Treatment)))
  })

points_df_2023A_sole <- bind_rows(
  lapply(dataTreatment2023ANWM2, function(df) {
    df %>% mutate(Species = Treatment[1])
  })
) %>% 
  filter(Species %in% treatment_order_2023A_sole) %>% 
  mutate(Species = factor(Species, levels = treatment_order_2023A_sole))

time_seq <- seq(min(points_df_2023A_sole$Time), filter(tempData2023, Date == as.Date("2023-05-23"))$cumulativeTemp + 30, length.out = 200)

functions_df_2023A_sole <- map2_dfr(fitsNWTLogisticNormal2023A[c(4, 5, 6, 7)], dataTreatment2023ANWM2[c(4, 5, 6, 7)], function(fit, df) {
  sp <- df$Treatment[1]
  d <- fit$par["d"]; r <- fit$par["r"]; h <- fit$par["h"]
  tibble(
    Time      = time_seq,
    Predicted = d * exp(r * (Time - h)) / (1 + exp(r * (Time - h))),
    Species   = sp
  )
}) %>% 
  filter(Species %in% treatment_order_2023A_sole) %>% 
  mutate(Species = factor(Species, levels = treatment_order_2023A_sole))

colours_2023A_sole <- setNames(colourVector[c(1, 2, 8, 9)],
                         c("T+", "T", "F", "F+"))
shapes_2023A_sole  <- c("T+"=16, "T"=17, "F"=4, "F+"=5)

weed_biomass_sole <- c(
  "T+       :  4.8",
  "T        : 15.3",
  "F        : 88.4",
  "F+       : 76.7"
  )

custom_labels <- paste0(names(weed_biomass_sole), ": ", weed_biomass_sole)

plot_2023A_sole <- ggplot() +
  geom_point(data = points_df_2023A_sole,
             aes(Time, Average, colour = Species, shape = Species),
             size = 4) +
  geom_line(data =functions_df_2023A_sole,
            aes(Time, Predicted, colour = Species),
            size = 1.2) +
  scale_colour_manual(
    name   = "Crop species",
    values = colours_2023A_sole,
    labels = custom_labels
  ) +
  scale_shape_manual(
    name   = "Crop species",
    values = shapes_2023A_sole,
    labels = custom_labels
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(shape    = unname(shapes_2023A_sole),
                          linetype = 0)
    ),
    shape = "none"
  ) +
  theme_classic(base_size = 25) +
    theme(
    legend.position = "inside",
    legend.position.inside = c(1.0, 0.05),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = alpha("white", 1.0), color = NA),
    legend.key.size = unit(2.0, "lines"),
    legend.text = element_text(family = "LiberationMono", size = 20),
    legend.title = element_text(size = 22)
  ) +
  labs(x = "Cumulative daily average temperature (°Cd)",
       y = "Proportion PAR intercepted",
       title = "") +
  geom_vline(xintercept = filter(tempData2023, Date == as.Date("2023-05-23"))$cumulativeTemp,
             linetype = "dashed", colour = "black", size = 1) +
  ylim(c(0, 1))

plot_2023A_sole

## b: intercrops
treatment_order_2023A_intercrop <- c("3T:1F", "1T:1F-M", "1T:1F", "1T:3F")
dataTreatment2023ANWM2 <- lapply(dataTreatment2023ANWM, function(df){
  return(df %>% mutate(Treatment = gsub("TF", "1T:1F", Treatment)))
  })

points_df_2023A_intercrop <- bind_rows(
  lapply(dataTreatment2023ANWM2, function(df) {
    df %>% mutate(Species = Treatment[1])
  })
) %>% 
  filter(Species %in% treatment_order_2023A_intercrop) %>% 
  mutate(Species = factor(Species, levels = treatment_order_2023A_intercrop))

time_seq <- seq(min(points_df_2023A_intercrop$Time), filter(tempData2023, Date == as.Date("2023-05-23"))$cumulativeTemp + 30, length.out = 200)

functions_df_2023A_intercrop <- map2_dfr(fitsNWTLogisticNormal2023A[c(1, 2, 3, 8)], dataTreatment2023ANWM2[c(1, 2, 3, 8)], function(fit, df) {
  sp <- df$Treatment[1]
  d <- fit$par["d"]; r <- fit$par["r"]; h <- fit$par["h"]
  tibble(
    Time      = time_seq,
    Predicted = d * exp(r * (Time - h)) / (1 + exp(r * (Time - h))),
    Species   = sp
  )
}) %>% 
  filter(Species %in% treatment_order_2023A_intercrop) %>% 
  mutate(Species = factor(Species, levels = treatment_order_2023A_intercrop))

colours_2023A_intercrop <- setNames(colourVector[c(3, 4, 5, 7)],
                         c("3T:1F", "1T:1F-M", "1T:1F", "1T:3F"))
shapes_2023A_intercrop  <- c("3T:1F"=15, "1T:1F-M"=18, "1T:1F"=8, "1T:3F"=3)

weed_biomass_intercrop <- c(
  "3T:1F    : 12.6",
  "1T:1F-M  : 15.8",
  "1T:1F    : 20.3",
  "1T:3F    : 30.0"
  )

custom_labels <- paste0(names(weed_biomass_intercrop), ": ", weed_biomass_intercrop)

# Extract T and F sole crop data for grey reference lines
functions_df_2023A_sole_TF <- functions_df_2023A_sole %>%
  filter(Species %in% c("T", "F"))

plot_2023A_intercrop <- ggplot() +
  geom_line(data = functions_df_2023A_sole_TF,
            aes(Time, Predicted, group = Species),
            color = "grey75", size = 1.2, linetype = "solid") +
  geom_point(data = points_df_2023A_intercrop,
             aes(Time, Average, colour = Species, shape = Species),
             size = 4) +
  geom_line(data = functions_df_2023A_intercrop,
            aes(Time, Predicted, colour = Species),
            size = 1.2) +
  scale_colour_manual(
    name   = "Crop species",
    values = c(colours_2023A_intercrop, colours_2023A_sole),
    labels = custom_labels
  ) +
  scale_shape_manual(
    name   = "Crop species",
    values = shapes_2023A_intercrop,
    labels = custom_labels
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(shape    = unname(shapes_2023A_intercrop),
                          linetype = 0)
    ),
    shape = "none"
  ) +
  theme_classic(base_size = 25) +
    theme(
    legend.position = "inside",
    legend.position.inside = c(1.0, 0.05),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = alpha("white", 1.0), color = NA),
    legend.key.size = unit(2.0, "lines"),
    legend.text = element_text(family = "LiberationMono", size = 20),
    legend.title = element_text(size = 22)
  ) +
  labs(x = "Cumulative daily average temperature (°Cd)",
       y = "Proportion PAR intercepted",
       title = "") +
  geom_vline(xintercept = filter(tempData2023, Date == as.Date("2023-05-23"))$cumulativeTemp,
             linetype = "dashed", colour = "black", size = 1) +
  ylim(c(0, 1))

plot_2023A_intercrop

plot_2023A <- plot_2023A_sole + plot_2023A_intercrop +
  plot_layout(ncol = 2, axis_titles = "collect") +
  plot_annotation(tag_levels = "a")

plot_2023A

## Fig. 5: 2023-RD intercrop light interception
### a: sole crops
treatment_order_2023B_sole <- c("T", "T-375", "F", "F-375")
dataTreatment2023BNWM2 <- lapply(dataTreatment2023BNWM, function(df){
  return(df %>% mutate(Treatment = gsub("TF", "1T:1F", Treatment)))
  })

points_df_2023B_sole <- bind_rows(
  lapply(dataTreatment2023BNWM2, function(df) {
    df %>% mutate(Species = Treatment[1])
  })
) %>% 
  filter(Species %in% treatment_order_2023B_sole) %>% 
  mutate(Species = factor(Species, levels = treatment_order_2023B_sole))

time_seq <- seq(min(points_df_2023B_sole$Time), filter(tempData2023, Date == as.Date("2023-05-23"))$cumulativeTemp + 200, length.out = 200)

functions_df_2023B_sole <- map2_dfr(fitsNWTLogisticNormal2023B, dataTreatment2023BNWM2, function(fit, df) {
  sp <- df$Treatment[1]
  d <- fit$par["d"]; r <- fit$par["r"]; h <- fit$par["h"]
  tibble(
    Time      = time_seq,
    Predicted = d * exp(r * (Time - h)) / (1 + exp(r * (Time - h))),
    Species   = sp
  )
}) %>% 
  filter(Species %in% treatment_order_2023B_sole) %>% 
  mutate(Species = factor(Species, levels = treatment_order_2023B_sole))

colours_2023B_sole <- setNames(colourVector[c(1, 2, 8, 9)],
                         c("T", "T-375", "F", "F-375"))
shapes_2023B_sole  <- c("T"=16, "T-375"=17, "F"=4, "F-375"=5)

custom_labels_sole <- c(
  "T           :  7.2",
  "T-375       : 16.5",
  "F           : 53.9",
  "F-375       : 66.2"
)

plot_2023B_sole <- ggplot() +
  geom_point(data = points_df_2023B_sole,
             aes(Time, Average, colour = Species, shape = Species),
             size = 4) +
  geom_line(data = functions_df_2023B_sole,
            aes(Time, Predicted, colour = Species),
            size = 1.2) +
  scale_colour_manual(
    name   = "Crop species",
    values = colours_2023B_sole,
    labels = custom_labels_sole
  ) +
  scale_shape_manual(
    name   = "Crop species",
    values = shapes_2023B_sole,
    labels = custom_labels_sole
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(shape    = unname(shapes_2023B_sole),
                          linetype = 0)
    ),
    shape = "none"
  ) +
  theme_classic(base_size = 25) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.05),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = alpha("white", 1.0), color = NA),
    legend.key.size = unit(2.0, "lines"),
    legend.text = element_text(family = "LiberationMono", size = 20),
    legend.title = element_text(size = 22)
  ) +
  labs(x = "Cumulative daily average temperature (°Cd)",
       y = "Proportion PAR intercepted",
       title = "") +
  geom_vline(xintercept = filter(tempData2023, Date == as.Date("2023-05-23"))$cumulativeTemp,
             linetype = "dashed", colour = "black", size = 1) +
  ylim(c(0, 1))

plot_2023B_sole

### b: 12.5 cm
treatment_order_2023B_125 <- c("T", "1T:1F", "1T:1F-M", "F")
dataTreatment2023BNWM2 <- lapply(dataTreatment2023BNWM, function(df){
  return(df %>% mutate(Treatment = gsub("TF", "1T:1F", Treatment)))
  })

points_df_2023B_125 <- bind_rows(
  lapply(dataTreatment2023BNWM2, function(df) {
    df %>% mutate(Species = Treatment[1])
  })
) %>% 
  filter(Species %in% treatment_order_2023B_125) %>% 
  mutate(Species = factor(Species, levels = treatment_order_2023B_125))

time_seq <- seq(min(points_df_2023B_125$Time), filter(tempData2023, Date == as.Date("2023-05-23"))$cumulativeTemp + 200, length.out = 200)

functions_df_2023B_125 <- map2_dfr(fitsNWTLogisticNormal2023B, dataTreatment2023BNWM2, function(fit, df) {
  sp <- df$Treatment[1]
  d <- fit$par["d"]; r <- fit$par["r"]; h <- fit$par["h"]
  tibble(
    Time      = time_seq,
    Predicted = d * exp(r * (Time - h)) / (1 + exp(r * (Time - h))),
    Species   = sp
  )
}) %>% 
  filter(Species %in% treatment_order_2023B_125) %>% 
  mutate(Species = factor(Species, levels = treatment_order_2023B_125))

colours_2023B_125 <- setNames(colourVector[c(1, 3, 4, 8)],
                         c("T", "1T:1F", "1T:1F-M", "F"))
shapes_2023B_125  <- c("T"=16, "1T:1F"=15, "1T:1F-M"=18, "F"=4)

custom_labels_125 <- c(
  "T           :  7.2",
  "1T:1F       : 13.6",
  "1T:1F-M     : 14.8",
  "F           : 53.9"
)

plot_2023B_125 <- ggplot() +
  geom_point(data = points_df_2023B_125,
             aes(Time, Average, colour = Species, shape = Species),
             size = 4) +
  geom_line(data = functions_df_2023B_125,
            aes(Time, Predicted, colour = Species),
            size = 1.2) +
  scale_colour_manual(
    name   = "Crop species",
    values = colours_2023B_125,
    labels = custom_labels_125
  ) +
  scale_shape_manual(
    name   = "Crop species",
    values = shapes_2023B_125,
    labels = custom_labels_125
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(shape    = unname(shapes_2023B_125),
                          linetype = 0)
    ),
    shape = "none"
  ) +
  theme_classic(base_size = 25) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.05),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = alpha("white", 1.0), color = NA),
    legend.key.size = unit(2.0, "lines"),
    legend.text = element_text(family = "LiberationMono", size = 20),
    legend.title = element_text(size = 22)
  ) +
  labs(x = "Cumulative daily average temperature (°Cd)",
       y = "Proportion PAR intercepted",
       title = "") +
  geom_vline(xintercept = filter(tempData2023, Date == as.Date("2023-05-23"))$cumulativeTemp,
             linetype = "dashed", colour = "black", size = 1) +
  ylim(c(0, 1))

plot_2023B_125

### c: 37.5 cm
treatment_order_2023B_375 <- c("T-375", "1T:1F-375", "1T:1F-M-375", "F-375")
dataTreatment2023BNWM2 <- lapply(dataTreatment2023BNWM, function(df){
  return(df %>% mutate(Treatment = gsub("TF", "1T:1F", Treatment)))
  })

points_df_2023B_375 <- bind_rows(
  lapply(dataTreatment2023BNWM2, function(df) {
    df %>% mutate(Species = Treatment[1])
  })
) %>% 
  filter(Species %in% treatment_order_2023B_375) %>% 
  mutate(Species = factor(Species, levels = treatment_order_2023B_375))

time_seq <- seq(min(points_df_2023B_375$Time), filter(tempData2023, Date == as.Date("2023-05-23"))$cumulativeTemp + 200, length.out = 200)

functions_df_2023B_375 <- map2_dfr(fitsNWTLogisticNormal2023B, dataTreatment2023BNWM2, function(fit, df) {
  sp <- df$Treatment[1]
  d <- fit$par["d"]; r <- fit$par["r"]; h <- fit$par["h"]
  tibble(
    Time      = time_seq,
    Predicted = d * exp(r * (Time - h)) / (1 + exp(r * (Time - h))),
    Species   = sp
  )
}) %>% 
  filter(Species %in% treatment_order_2023B_375) %>% 
  mutate(Species = factor(Species, levels = treatment_order_2023B_375))

colours_2023B_375 <- setNames(colourVector[c(2, 5, 7, 9)],
                         c("T-375", "1T:1F-375", "1T:1F-M-375", "F-375"))
shapes_2023B_375  <- c("T-375"=17, "1T:1F-375"=8, "1T:1F-M-375"=3, "F-375"=5)

custom_labels_375 <- c(
  "T-375       : 16.5",
  "1T:1F-375   : 31.2",
  "1T:1F-M-375 : 25.3",
  "F-375       : 66.2"
)

plot_2023B_375 <- ggplot() +
  geom_point(data = points_df_2023B_375,
             aes(Time, Average, colour = Species, shape = Species),
             size = 4) +
  geom_line(data = functions_df_2023B_375,
            aes(Time, Predicted, colour = Species),
            size = 1.2) +
  scale_colour_manual(
    name   = "Crop species",
    values = colours_2023B_375,
    labels = custom_labels_375
  ) +
  scale_shape_manual(
    name   = "Crop species",
    values = shapes_2023B_375,
    labels = custom_labels_375
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(shape    = unname(shapes_2023B_375),
                          linetype = 0)
    ),
    shape = "none"
  ) +
  theme_classic(base_size = 25) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.05),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = alpha("white", 1.0), color = NA),
    legend.key.size = unit(2.0, "lines"),
    legend.text = element_text(family = "LiberationMono", size = 20),
    legend.title = element_text(size = 22)
  ) +
  labs(x = "Cumulative daily average temperature (°Cd)",
       y = "Proportion PAR intercepted",
       title = "") +
  geom_vline(xintercept = filter(tempData2023, Date == as.Date("2023-05-23"))$cumulativeTemp,
             linetype = "dashed", colour = "black", size = 1) +
  ylim(c(0, 1))

plot_2023B_375

plot_2023B <- plot_2023B_sole + plot_2023B_125 + plot_2023B_375 +
  plot_layout(design = "
    12
    3#
  ",
  axis_titles = "collect") +
  plot_annotation(tag_levels = "a")

plot_2023B

# Figure 7: all T, 1T:1F, F treatments across years.
custom_labels_2022 <- c(
  "T     :  49.4",
  "1T:1F :  99.7",
  "F     : 453.8"
)

(
plot_1t1f_2022 <- ggplot() +
  geomPoints2022[c("Triticale", "Triticale_Faba", "Faba")] +
  geomFunctions2022[c("Triticale", "Triticale_Faba", "Faba")] +
  scale_colour_manual(
    name   = "Treatment",
    values = modelColors,
    labels = c(
      "Triticale" = custom_labels_2022[1],
      "Triticale_Faba" = custom_labels_2022[2],
      "Faba" = custom_labels_2022[3]
    )
  ) +
  theme_classic(base_size = 25) +
  labs(
    x = "Cumulative daily average temperature (°Cd)",
    y = "Proportion PAR intercepted",
    subtitle = "2022"
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.9, 0.05),
    legend.justification   = c("right", "bottom"),
    legend.background      = element_rect(fill = alpha("white", 1.0), color = NA),
    legend.key.size = unit(2.0, "lines"),
    legend.text = element_text(family = "LiberationMono", size = 16),
    legend.title = element_text(size = 18)
  ) +
  geom_vline(xintercept = filter(tempData2022, Date == as.Date("2022-06-21"))$CumulTemp,
             linetype = "dashed", colour = "black", size = 1) +
  xlim(c(100, filter(tempData2022, Date == as.Date("2022-06-21"))$CumulTemp + 10))
)

custom_labels_2023A <- c(
  "T     : 15.3",
  "1T:1F : 20.3",
  "F     : 88.4"
)

(
plot_1t1f_2023A <- ggplot() +
  geomPoints2023ANW[c("T", "1T:1F", "F")] +
  geomFunctions2023ANW[c("T", "1T:1F", "F")] +
  scale_colour_manual(
    name   = "Treatment",
    values = modelColors,
    labels = c(
      "T" = custom_labels_2023A[1],
      "1T:1F" = custom_labels_2023A[2],
      "F" = custom_labels_2023A[3]
    )
  ) +
  theme_classic(base_size = 25) +
  labs(
    x = "Cumulative daily average temperature (°Cd)",
    y = "Proportion PAR intercepted",
    subtitle = "2023-SP"
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.985, 0.05),
    legend.justification   = c("right", "bottom"),
    legend.background      = element_rect(fill = alpha("white", 1.0), color = NA),
    legend.key.size = unit(2.0, "lines"),
    legend.text = element_text(family = "LiberationMono", size = 16),
    legend.title = element_text(size = 18)
  ) +
  geom_vline(xintercept = filter(tempData2023, Date == as.Date("2023-05-23"))$cumulativeTemp,
             linetype = "dashed", colour = "black", size = 1) +
  xlim(c(200, filter(tempData2023, Date == as.Date("2023-05-23"))$cumulativeTemp + 50))
)

custom_labels_2023B <- c(
  "T     :  7.2",
  "1T:1F : 13.6",
  "F     : 53.9"
)

(
plot_1t1f_2023B <- ggplot() +
  geomPoints2023BNW[c("T", "1T:1F", "F")] +
  geomFunctions2023BNW[c("T", "1T:1F", "F")] +
  scale_colour_manual(
    name   = "Treatment",
    values = modelColors,
    labels = c(
      "T" = custom_labels_2023B[1],
      "1T:1F" = custom_labels_2023B[2],
      "F" = custom_labels_2023B[3]
    )
  ) +
  theme_classic(base_size = 25) +
  labs(
    x = "Cumulative daily average temperature (°Cd)",
    y = "Proportion PAR intercepted",
    subtitle = "2023-RD"
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.9, 0.05),
    legend.justification   = c("right", "bottom"),
    legend.background      = element_rect(fill = alpha("white", 1.0), color = NA),
    legend.key.size = unit(2.0, "lines"),
    legend.text = element_text(family = "LiberationMono", size = 16),
    legend.title = element_text(size = 18)
  ) +
  geom_vline(xintercept = filter(tempData2023, Date == as.Date("2023-05-23"))$cumulativeTemp,
             linetype = "dashed", colour = "black", size = 1) +
  xlim(c(250, filter(tempData2023, Date == as.Date("2023-05-23"))$cumulativeTemp + 200))
)

custom_labels_2024 <- c(
  "T     :  60.2",
  "1T:1F :  63.6",
  "F     : 129.9"
)

(
plot_1t1f_2024 <- ggplot() +
  geomPoints2024NW[c("T", "1T:1F", "F")] +
  geomFunctions2024NW[c("T", "1T:1F", "F")] +
  scale_colour_manual(
    name   = "Treatment",
    values = modelColors,
    labels = c(
      "T" = custom_labels_2024[1],
      "1T:1F" = custom_labels_2024[2],
      "F" = custom_labels_2024[3]
    )
  ) +
  theme_classic(base_size = 25) +
  labs(
    x = "Cumulative daily average temperature (°Cd)",
    y = "Proportion PAR intercepted",
    subtitle = "2024"
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.9, 0.05),
    legend.justification   = c("right", "bottom"),
    legend.background      = element_rect(fill = alpha("white", 1.0), color = NA),
    legend.key.size = unit(2.0, "lines"),
    legend.text = element_text(family = "LiberationMono", size = 16),
    legend.title = element_text(size = 18)
  ) +
  geom_vline(xintercept = filter(tempData2024, Date == as.Date("2024-06-26"))$cumulativeTemp,
             linetype = "dashed", colour = "black", size = 1) +
  xlim(c(150, filter(tempData2024, Date == as.Date("2024-06-26"))$cumulativeTemp + 10))
)

(
combined_plot_1t1f <- 
  plot_1t1f_2022 + ylim(c(0, 1)) + 
  plot_1t1f_2023A + ylim(c(0, 1)) + 
  plot_1t1f_2023B + ylim(c(0, 1)) + 
  plot_1t1f_2024 + ylim(c(0, 1)) + 
  plot_layout(ncol = 2, nrow = 2, axis_titles = "collect") +
  plot_annotation(tag_levels = "a")
)

png("EJA/Draft/Code and data/Figures/plot_light_2022.png", units = "px", width = 7000, height = 5000, res = 300)
combined_plot_2022
dev.off()

png("EJA/Draft/Code and data/Figures/plot_light_2023ANW.png", units = "px", width = 6200, height = 3100, res = 300)
plot_2023A
dev.off()

png("EJA/Draft/Code and data/Figures/plot_light_2023BNW.png", units = "px", width = 6000, height = 6000, res = 300)
plot_2023B
dev.off()

png("EJA/Draft/Code and data/Figures/plot_light_1t1f.png", units = "px", width = 4850, height = 4350, res = 300)
combined_plot_1t1f
dev.off()

