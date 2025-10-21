library(tidyverse)
library(emmeans)
library(multcomp)
library(stringr)

source("WCFTreatments.R")

getBiomass2022 <- function(data){
  treatmentBiomass = sapply(treatments2022, function(t){
    splitType = strsplit(t, split = "_")[[1]]
    if(length(splitType) > 2){
      treatment = paste(splitType[1], splitType[2], sep = "_")
      type = splitType[3]
    }else{
      treatment = splitType[1]
      type = splitType[2]
    }
    column = ifelse(is.na(type), ifelse(treatment %in% c("Pea", "Lupine", "Faba"), "L", "C"), type) # if type is NA, treatment is a T or F (+) SC
    typeData = data[which(data$Treatment == treatment), column]
    return(typeData)
  })
  empty <- which(lapply(treatmentBiomass, function(t)length(t)) == 0)
  if(length(empty) > 0){
    treatmentBiomass[empty] <- NULL
    
    dataTypes <- tibble(Type = as.factor(rep(treatments2022[-empty], each = 4)), # repeat row type for every replicate
                        Block = rep(1:4, times = length(treatments2022[-empty])),
                        Biomass = unlist(treatmentBiomass))
  }else{
    dataTypes <- tibble(Type = as.factor(rep(treatments2022, each = 4)), # repeat row type for every replicate
                        Block = rep(1:4, times = length(treatments2022)),
                        Biomass = unlist(treatmentBiomass))
  }
  
  
  return(dataTypes)
}

getBiomass2023A <- function(data, separateStrips = FALSE){
  types <- if(separateStrips){rowTypes2023A}else{treatments2023A}
  typeBiomass = sapply(types, function(type){
    splitType = strsplit(type, split = "_")[[1]]
    treatment = splitType[1]
    type = splitType[2]
    column = ifelse(is.na(type), ifelse(grepl("T", treatment), "T", "F"), type) # if type is NA, treatment is a T or F (+) SC
    typeData = data[which(data$Treatment == treatment), column]
    return(typeData)
  })
  
  dataTypes <- tibble(Type = as.factor(rep(types, each = 5)), # repeat row type for every replicate
                      Block = rep(1:5, times = length(types)),
                      Biomass = unlist(typeBiomass))
  
  return(dataTypes)
}

getBiomass2023B <- function(data){
  treatmentBiomass = sapply(treatments2023B, function(treatment){
    splitTreatment = strsplit(treatment, split = "_")[[1]]
    treatment = splitTreatment[1]
    type = splitTreatment[2]
    column = ifelse(is.na(type), ifelse(grepl("T", treatment), "T", "F"), type) # if type is NA, treatment is a T or F (+) SC
    typeData = data[which(data$Treatment == treatment), column]
    return(typeData)
  })
  
  dataTreatments <- tibble(Type = as.factor(rep(treatments2023B, each = 4)), # repeat row type for every replicate
                           Block = rep(1:4, times = length(treatments2023B)),
                           Biomass = unlist(treatmentBiomass))
  
  return(dataTreatments)
}

getBiomass2024 <- function(data, NW = FALSE){
  types <- treatments2024
  typeBiomass = sapply(types, function(type){
    splitType = strsplit(type, split = "_")[[1]]
    treatment = splitType[1]
    type = splitType[2]
    column = ifelse(is.na(type), ifelse(grepl("T", treatment), "T", "F"), type) # if type is NA, treatment is a T or F (+) SC
    typeData = data[which(data$Treatment == treatment), column]
    return(typeData)
  })
  
  if(NW){
    types <- types[which(types != "W")]
    dataTypes <- tibble(Type = as.factor(rep(types, each = 5)), # repeat row type for every replicate
                        Block = rep(1:5, times = length(types)),
                        Biomass = unlist(typeBiomass))
    
  }else{
    dataTypes <- tibble(Type = as.factor(rep(types, each = 5)), # repeat row type for every replicate
                        Block = rep(1:5, times = length(types)),
                        Biomass = unlist(typeBiomass))
    
  }
  
  return(dataTypes)
}

calcCorrectedBiomass2022 <- function(data){
  dataTypes <- getBiomass2022(data)
  multiplication <- sapply(treatments2022, function(t)ifelse(grepl("_", t), 2, 1))
  correctedDataTypes <- dataTypes %>% 
    mutate(Biomass = Biomass * multiplication[match(Type, names(multiplication))])
  
  return(correctedDataTypes)
}

calcCorrectedBiomass2023A <- function(data, separateStrips = FALSE){
  dataTypes <- getBiomass2023A(data, separateStrips)
  multiplication <- c("1T:1F_T" = 2,
                      "1T:1F_F" = 2,
                      "1T:3F_T" = 4,
                      "1T:3F_F" = 4/3,
                      "3T:1F_T" = 4/3,
                      "3T:1F_F" = 4,
                      "TF-M_T" = 2,
                      "TF-M_F" = 2,
                      "T" = 1,
                      "F" = 1,
                      "T+" = 4/6,
                      "F+" = 4/6)
  correctedDataTypes <- dataTypes %>% 
    mutate(Biomass = Biomass * multiplication[match(Type, names(multiplication))])
  
  return(correctedDataTypes)
}

calcCorrectedBiomass2023B <- function(data){
  dataTreatments <- getBiomass2023B(data)
  multiplication <- c("T" = 1,
                      "1T:1F_T" = 2,
                      "1T:1F_F" = 2,
                      "TF-M_T" = 2,
                      "TF-M_F" = 2,
                      "F" = 1,
                      #"T-25" = 1,
                      #"1T:1F-25_T" = 2,
                      #"1T:1F-25_F" = 2,
                      #"F-25" = 1,
                      "T-375" = 1,
                      "1T:1F-375_T" = 2,
                      "1T:1F-375_F" = 2,
                      "TF-M-375_T" = 2,
                      "TF-M-375_F" = 2,
                      "F-375" = 1)
  correctedDataTreatments <- dataTreatments %>% 
    mutate(Biomass = Biomass * multiplication[match(Type, names(multiplication))])
  
  return(correctedDataTreatments)
}

calcCorrectedBiomass2024 <- function(data, NW = FALSE){
  dataTreatments <- getBiomass2024(data, NW)
  multiplication <- c("T" = 1,
                      "F" = 1,
                      "1T:1F_T" = 2,
                      "1T:1F_F" = 2,
                      "W" = 1)
  correctedDataTreatments <- dataTreatments %>% 
    mutate(Biomass = Biomass * multiplication[match(Type, names(multiplication))])
  
  return(correctedDataTreatments)
}

calcRelativeBiomass <- function(data, experiment){
  getDataFunctions <- list("2022" = getBiomass2022,
                           "2023A" = getBiomass2023A,
                           "2023B" = getBiomass2023B,
                           "2024" = getBiomass2024)
  dataTypes <- getDataFunctions[[experiment]](data) 
  crop <- sapply(dataTypes$Type, function(t){
    crop <- strsplit(as.character(t), "_")[[1]][2]
    crop <- ifelse(is.na(crop), strsplit(as.character(t), "\\+")[[1]][1], crop)
    return(crop)
  })
  dataTypesCrop <- dataTypes %>% 
    mutate(Crop = crop)
  
  meanSoleCropCerealBiomass <- mean(filter(dataTypesCrop, Crop %in% c("T", "C"))$Biomass, na.rm = TRUE)
  meanSoleCropLegumeBiomass <- mean(filter(dataTypesCrop, Crop %in% c("F", "L"))$Biomass, na.rm = TRUE)
  relativeBiomass <- sapply(1:nrow(dataTypesCrop), function(i){
    return(ifelse(dataTypesCrop$Crop[i] %in% c("T", "C"), 
                  dataTypesCrop$Biomass[i] / meanSoleCropCerealBiomass, 
                  dataTypesCrop$Biomass[i] / meanSoleCropLegumeBiomass))
  })
  newData <- dataTypes %>% 
    mutate(RelativeBiomass = relativeBiomass)
  
  return(newData)
}

getBiomassCLD <- function(data, relative = FALSE){
  data <- na.omit(data)
  variable <- ifelse(relative, "RelativeBiomass", "Biomass")
  # Check for normality and heteroscedasticity
  ## T
  modbiomass <- lmer(data[,variable][[1]] ~ (1|Block) + Treatment, data = data)
  
  nonNormal <- shapiro.test(residuals(modbiomass))$p.value < 0.05
  unequalVariance <- bartlett.test(residuals(modbiomass) ~ data$Treatment)$p.value < 0.05
  
  newDat <- data
  if(nonNormal | unequalVariance){
    bc <- boxcox(data[,variable][[1]] ~ data$Treatment, lambda = seq(-5, 5, 1/10))
    lambda <- bc$x[which.max(bc$y)]
    if(lambda > -0.3 & lambda < 0.3){
      newDat[,variable] <- log(newDat[,variable])
    }else{
      newDat[,variable] <- newDat[,variable]^(bc$x[which.max(bc$y)])  
    }
    modbiomass <- lmer(newDat[,variable][[1]] ~ (1|Block) + Treatment, data = newDat)
    nonNormal <- shapiro.test(residuals(modbiomass))$p.value < 0.05
    unequalVariance <- bartlett.test(residuals(modbiomass) ~ data$Treatment)$p.value < 0.05
    if(nonNormal | unequalVariance){
      cat("Cannot guarantee normal distribution and/or equal variances. Please manually check for these assumptions. \n")
    }
  }
  cat(paste("p-value = ", summary(aov(newDat[,variable][[1]] ~ (1|Block) + Treatment, data = newDat))[[1]]$`Pr(>F)`[1], "\n", sep = ""))
  if(relative){
    modbiomass <- lmer(RelativeBiomass ~ (1|Block) + Treatment, data = newDat)
    PHBiomass <- emmeans(modbiomass, list(pairwise ~ Treatment), adjust = "tukey")
    CLDBiomass <- cld(PHBiomass$emmeans,
                      Letters = letters,
                      descending = TRUE)
  }else{
    modbiomass <- lmer(Biomass ~ (1|Block) + Treatment, data = newDat)
    PHBiomass <- emmeans(modbiomass, list(pairwise ~ Treatment), adjust = "tukey")
    CLDBiomass <- cld(PHBiomass$emmeans,
                      Letters = letters,
                      descending = TRUE)
  }
  
  return(CLDBiomass)
}

calcPLERData2022 <- function(data){
  soleCrops <- c("Barley", "Rye", "Triticale", "Wheat", "Pea", "Lupine", "Faba")
  dataTypes <- getBiomass2022(data)
  
  dataTypesIC <- dataTypes %>% 
    filter(!(Type %in% soleCrops))
  
  # Get data of the sole crops at normal density
  dataSC <- data[data$Treatment %in% soleCrops,]
  averageSCC <- aggregate(C ~ Treatment, data = dataSC, FUN = function(x)mean(x, na.rm = TRUE)) %>% 
    filter(C != 0)
  averageSCL <- aggregate(L ~ Treatment, data = dataSC, FUN = function(x)mean(x, na.rm = TRUE)) %>% 
    filter(L != 0)
  averageSC <- c(averageSCC$C, averageSCL$L)
  names(averageSC) <- c(averageSCC$Treatment, averageSCL$Treatment)
  
  dataTypesICPLER <- sapply(split(dataTypesIC, seq(nrow(dataTypesIC))), function(dat){
    cropType = match(substr(strsplit(as.character(dat$Type), split = "_")[[1]][3], 1, 1), c("cereal" = "C", "legume" = "L"))
    crop = ifelse(cropType == 1, strsplit(as.character(dat$Type), "_")[[1]][1], strsplit(as.character(dat$Type), "_")[[1]][2])
    pLER = dat$Biomass / averageSC[crop]
    return(pLER)
  })
  
  dataTypesIC <- dataTypesIC %>% 
    mutate(pLER = dataTypesICPLER)
  
  return(dataTypesIC)
}

calcPLERData2023A <- function(data, separateStrips = FALSE){
  dataTypes <- getBiomass2023A(data, separateStrips)
  
  dataTypesIC <- dataTypes %>% 
    filter(!(Type %in% c("T", "F", "T+", "F+")))
  
  # Get data of the sole crops at normal density
  dataSC <- data[data$Treatment %in% c("T", "F"),]
  averageSC <- c(triticale = mean(dataSC$T, na.rm = TRUE), faba = mean(dataSC$F, na.rm = TRUE))
  
  dataTypesICPLER <- sapply(split(dataTypesIC, seq(nrow(dataTypesIC))), function(dat){
    crop = match(substr(strsplit(as.character(dat$Type), split = "_")[[1]][2], 1, 1), c("triticale" = "T", "faba" = "F"))
    pLER = dat$Biomass / averageSC[crop]
    return(pLER)
  })
  
  dataTypesIC <- dataTypesIC %>% 
    mutate(pLER = dataTypesICPLER)
  
  return(dataTypesIC)
}

calcPLERData2023B <- function(data){
  dataTreatments <- getBiomass2023B(data)
  
  dataTreatmentsIC <- dataTreatments %>% 
    filter(!(Type %in% c("T", "T-375", "F", "F-375")))
  
  # Get data of the sole crops at normal density
  dataSC <- data[data$Treatment %in% c("T", "T-375", "F", "F-375"),]
  averageSC <- c(triticale125 = mean(filter(dataSC, Treatment == "T")$T, na.rm = TRUE), 
                 faba125 = mean(filter(dataSC, Treatment == "F")$F, na.rm = TRUE),
                 triticale375 = mean(filter(dataSC, Treatment == "T-375")$T, na.rm = TRUE), 
                 faba375 = mean(filter(dataSC, Treatment == "F-375")$F, na.rm = TRUE))
  
  dataTreatmentsICPLER <- sapply(split(dataTreatmentsIC, seq(nrow(dataTreatmentsIC))), function(dat){
    crop = ifelse(strsplit(as.character(dat$Type), split = "_")[[1]][2], ifelse(grepl("375", dat$Type), "triticale375", "triticale125"), ifelse(grepl("375", dat$Type), "faba375", "faba125"))
    pLER = dat$Biomass / averageSC[crop]
    return(pLER)
  })
  
  dataTreatmentsIC <- dataTreatmentsIC %>% 
    mutate(pLER = dataTreatmentsICPLER)
  
  return(dataTreatmentsIC)
}

calcPLERData2024 <- function(data, NW = FALSE){
  dataTreatments <- getBiomass2024(data, NW)
  
  dataTypesIC <- dataTreatments %>% 
    filter(!(Type %in% c("T", "F")))
  
  # Get data of the sole crops at normal density
  dataSC <- data[data$Treatment %in% c("T", "F"),]
  averageSC <- c(triticale = mean(filter(dataSC, Treatment == "T")$T, na.rm = TRUE), 
                 faba = mean(filter(dataSC, Treatment == "F")$F, na.rm = TRUE))
  
  dataTypesICPLER <- sapply(split(dataTypesIC, seq(nrow(dataTypesIC))), function(dat){
    crop = match(substr(strsplit(as.character(dat$Type), split = "_")[[1]][2], 1, 1), c("triticale" = "T", "faba" = "F"))
    pLER = dat$Biomass / averageSC[crop]
    return(pLER)
  })
  
  dataTypesIC <- dataTypesIC %>% 
    mutate(pLER = dataTypesICPLER)
  
  return(dataTypesIC)
}

separateCL <- function(data, cTreatments, lTreatments, cSubscript, lSubscript){
  dataC <- data %>% 
    filter(grepl(cSubscript, Type) | Type %in% cTreatments) %>% 
    rename(Treatment = Type)
  treatmentsC <- sapply(as.character(dataC$Treatment), function(t){
    spl <- strsplit(t, split = "_")[[1]]
    if(length(spl) > 2){
      return(paste(spl[1], spl[2], sep = "_"))
    }else{
      return(spl[1])
    }
  })
  dataC <- mutate(dataC, Treatment = as.factor(treatmentsC))
  
  dataL <- data %>% 
    filter(grepl(lSubscript, Type) | Type %in% lTreatments) %>% 
    rename(Treatment = Type)
  treatmentsL <- sapply(as.character(dataL$Treatment), function(t){
    spl <- strsplit(t, split = "_")[[1]]
    if(length(spl) > 2){
      return(paste(spl[1], spl[2], sep = "_"))
    }else{
      return(spl[1])
    }
  })
  dataL <- mutate(dataL, Treatment = as.factor(treatmentsL))
  
  return(list(dataC, dataL))
}

calcPLERAvg <- function(data){  
  # Check normal distribution
  shapiroPValues <- sapply(data, function(d)shapiro.test(d$pLER)$p.value)
  nonNormal <- 1:length(shapiroPValues) %in% c(which(shapiroPValues < 0.05))

  toTransform = which(nonNormal)
  
  pValues <- sapply(1:length(data), function(i){
    cat(paste("Comparing pLER means of " , names(data)[i], "\n", sep = ""))
    if(i %in% toTransform){
      newDat <- data[[i]] %>%
        mutate(pLER = log(pLER))
      cat("Non-normal data or unequal variance detected. Applying log transformation.\n")
      if(shapiro.test(newDat$pLER)$p.value < 0.05){
        cat("Data is still non-normally distributed. Resorting to non-parametric test.\n")
        dat <- data[[i]]
        if(grepl("TM", dat$Type[1]) | grepl("FM", dat$Type[1]) | (strsplit(as.character(dat$Type[1]), split = "_")[[1]][2] == "F" & strsplit(as.character(dat$Type[1]), split = "_")[[1]][1] == "3T:1F") | (strsplit(as.character(dat$Type[1]), split = "_")[[1]][2] == "T" & strsplit(as.character(dat$Type[1]), split = "_")[[1]][1] == "1T:3F")){
          return(wilcox.test(dat$pLER, mu = 0.25, alternative = "two.sided")$p.value)
        }else if(dat$Type[1] == "3T:1F_T" | dat$Type[1] == "1T:3F_F"){
          return(wilcox.test(dat$pLER, mu = 0.75, alternative = "two.sided")$p.value)
        }else{
          return(wilcox.test(dat$pLER, mu = 0.5, alternative = "two.sided")$p.value)
        }
      }else{
        if(grepl("TM", newDat$Type[1]) | grepl("FM", newDat$Type[1]) | (strsplit(as.character(newDat$Type[1]), split = "_")[[1]][2] == "F" & strsplit(as.character(newDat$Type[1]), split = "_")[[1]][1] == "3T:1F") | (strsplit(as.character(newDat$Type[1]), split = "_")[[1]][2] == "T" & strsplit(as.character(newDat$Type[1]), split = "_")[[1]][1] == "1T:3F")){
          return(t.test(newDat$pLER, mu = log(0.25), alternative = "two.sided")$p.value)
        }else if(newDat$Type[1] == "3T:1F_T" | newDat$Type[1] == "1T:3F_F"){
          return(t.test(newDat$pLER, mu = log(0.75), alternative = "two.sided")$p.value)
        }else{
          return(t.test(newDat$pLER, mu = log(0.5), alternative = "two.sided")$p.value)
        }
      }
    }else{
      dat <- data[[i]]
      if(grepl("TM", dat$Type[1]) | grepl("FM", dat$Type[1]) | (strsplit(as.character(dat$Type[1]), split = "_")[[1]][2] == "F" & strsplit(as.character(dat$Type[1]), split = "_")[[1]][1] == "3T:1F") | (strsplit(as.character(dat$Type[1]), split = "_")[[1]][2] == "T" & strsplit(as.character(dat$Type[1]), split = "_")[[1]][1] == "1T:3F")){
        return(t.test(dat$pLER, mu = 0.25, alternative = "two.sided")$p.value)
      }else if(dat$Type[1] == "3T:1F_T" | dat$Type[1] == "1T:3F_F"){
        return(t.test(dat$pLER, mu = 0.75, alternative = "two.sided")$p.value)
      }else{
        return(t.test(dat$pLER, mu = 0.5, alternative = "two.sided")$p.value)
      }
    }
  })
  
  dataAvg <- tibble(Type = sapply(data, function(dat)dat$Type[1]),
                    pLER = sapply(data, function(dat)mean(dat$pLER, na.rm = TRUE)),
                    pLERSD = sapply(data, function(dat)sd(dat$pLER, na.rm = TRUE)),
                    PValue = pValues)
  
  treatments <- sapply(1:nrow(dataAvg), function(i){
    spl <- strsplit(as.character(dataAvg$Type[i]), split = "_")[[1]]
    if(length(spl) > 2){
      return(paste(spl[1], spl[2], sep = "_"))
    }else{
      return(spl[1])
    }
  })
  
  types <- sapply(1:nrow(dataAvg), function(i){
    spl <- strsplit(as.character(dataAvg$Type[i]), split = "_")[[1]]
    if(length(spl) > 2){
      return(spl[3])
    }else{
      return(spl[2])
    }
  })
  
  dataAvg <- mutate(dataAvg, Treatment = treatments) %>% 
    mutate(Type = types)
  
  return(dataAvg)
}

calcLER2022 <- function(data){
  soleCrops <- c("Barley", "Rye", "Triticale", "Wheat", "Pea", "Lupine", "Faba")
  
  dataIC <- data %>% 
    filter(!(Treatment %in% soleCrops))
  
  # Get data of the sole crops at normal density
  dataSC <- data[data$Treatment %in% soleCrops,]
  averageSCC <- aggregate(C ~ Treatment, data = dataSC, FUN = function(x)mean(x, na.rm = TRUE)) %>% 
    filter(C != 0)
  averageSCL <- aggregate(L ~ Treatment, data = dataSC, FUN = function(x)mean(x, na.rm = TRUE)) %>% 
    filter(L != 0)
  averageSC <- c(averageSCC$C, averageSCL$L)
  names(averageSC) <- c(averageSCC$Treatment, averageSCL$Treatment)
  
  dataTypesICLER <- sapply(split(dataIC, seq(nrow(dataIC))), function(dat){
    LER = dat$C / averageSC[strsplit(dat$Treatment, "_")[[1]][1]] + dat$L / averageSC[strsplit(dat$Treatment, "_")[[1]][2]]
    return(LER)
  })
  
  newData <- data %>% 
    dplyr::select(Treatment, Block, C, L) %>% 
    filter(!Treatment %in% c("Barley", "Rye", "Triticale", "Wheat", "Pea", "Lupine", "Faba")) %>% 
    mutate(LER = dataTypesICLER)
  
  return(newData)
}

calcLER2023A <- function(data){
  newData <- data %>% 
    dplyr::select(Treatment, Block, T, F) %>% 
    mutate(LER = T / mean(data[data$Treatment == "T","T"][[1]], na.rm = TRUE) + F / mean(data[data$Treatment == "F","F"][[1]], na.rm = TRUE)) %>% 
    filter(!Treatment %in% c("T", "F", "T+", "F+"))
  return(newData)
}

calcLER2023B <- function(data){
  LER <- sapply(1:nrow(data), function(i){
    L <- ifelse(data$Treatment[i] %in% c("1T:1F", "TF-M"), data$T[i] / mean(filter(data, Treatment == "T")$T, na.rm = TRUE) + data$F[i] / mean(filter(data, Treatment == "F")$F, na.rm = TRUE),
                ifelse(data$Treatment[i] %in% c("1T:1F-375", "TF-M-375"), data$T[i] / mean(filter(data, Treatment == "T-375")$T, na.rm = TRUE) + data$F[i] / mean(filter(data, Treatment == "F-375")$F, na.rm = TRUE), NA))
    return(L)
  })
  newData <- data %>% 
    dplyr::select(Treatment, Block, T, F) %>% 
    mutate(LER = LER) %>% 
    filter(!Treatment %in% c("T", "F", "T-375", "F-375"))
  return(newData)
}

calcLER2024 <- function(data){
  newData <- data %>% 
    dplyr::select(Treatment, Block, T, F) %>% 
    mutate(LER = T / mean(data[data$Treatment == "T","T"][[1]], na.rm = TRUE) + F / mean(data[data$Treatment == "F","F"][[1]], na.rm = TRUE)) %>% 
    filter(!Treatment %in% c("T", "F"))
  return(newData)
}

calcLERAvg <- function(data){
  data <- na.omit(data)
  variable <- "LER"
  modLER <- lm(LER ~ Treatment, data = data)
  nonNormal <- shapiro.test(residuals(modLER))$p.value < 0.05
  unequalVariance <- bartlett.test(residuals(modLER) ~ data$Treatment)$p.value < 0.05
  newDat <- data
  if(nonNormal | unequalVariance){
    bc <- boxcox(data[,variable][[1]] ~ data$Treatment, lambda = seq(-5, 5, 1/10))
    newDat[,variable] <- newDat[,variable]^(bc$x[which.max(bc$y)])
    modLER <- lmer(newDat[,variable][[1]] ~ (1|Block) + Treatment, data = newDat)
    nonNormal <- shapiro.test(residuals(modLER))$p.value < 0.05
    unequalVariance <- bartlett.test(residuals(modLER) ~ data$Treatment)$p.value < 0.05
    if(nonNormal | unequalVariance){
      cat("Cannot guarantee normal distribution and/or equal variances. Please manually check for these assumptions. \n")
    }
  }
  
  dataAvg <- aggregate(LER ~ Treatment, data = data, function(x)mean(x, na.rm = TRUE)) %>% 
    mutate(SD = aggregate(LER ~ Treatment, data = data, function(x)sd(x, na.rm = TRUE))$LER) %>% 
    mutate(PValue = sapply(split(newDat, newDat$Treatment), function(dat)t.test(dat$LER, mu = 1.0, alternative = "two.sided")$p.value))
  
  return(dataAvg)
}

calcArithmeticHarmonic <- function(data, dataSC){
  soleC2022 <- c("Barley", "Rye", "Triticale", "Wheat")
  soleL2022 <- c("Pea", "Lupine", "Faba")
  soleC2023B <- c("T", "T-25", "T-375")
  soleL2023B <- c("F", "F-25", "F-375")
  inter2023B <- c("1T:1F", "1T:1F-25", "1T:1F-375", "TF-M", "TF-M-25", "TF-M-375")
  dataComb <- lapply(data, function(dat){
    datSCCSel <- dataSC[which(dataSC$Treatment %in% c(soleC2022, soleC2023B, "T")),]
    datSCCSelW <- datSCCSel[which(datSCCSel$Treatment %in% c(strsplit(dat$Treatment[1], "_")[[1]][1], "T", paste("T-", strsplit(dat$Treatment[1], "-")[[1]][length(strsplit(dat$Treatment[1], "-")[[1]])], sep = ""))),]
    if(nrow(datSCCSelW) > nrow(dat)){
      datSCCSelWW <- filter(datSCCSelW, Treatment == "T-375")$W
    }else{
      datSCCSelWW <- filter(datSCCSelW, Treatment %in% c("T", datSCCSelW$Treatment[1]))$W
    }
    dat$SCC <- datSCCSelWW
    dat$SCC[which(is.na(dat$SCC))] <- mean(dat$SCC, na.rm = TRUE)
    datSCLSel <- dataSC[which(dataSC$Treatment %in% c(soleL2022, soleL2023B, "F")),]
    datSCLSelW <- datSCLSel[which(datSCLSel$Treatment %in% c(strsplit(dat$Treatment[1], "_")[[1]][2], "F", paste("F-", strsplit(dat$Treatment[1], "-")[[1]][length(strsplit(dat$Treatment[1], "-")[[1]])], sep = ""))),]
    if(nrow(datSCLSelW) > nrow(dat)){
      datSCLSelWW <- filter(datSCLSelW, Treatment == "F-375")$W
    }else{
      datSCLSelWW <- filter(datSCLSelW, Treatment %in% c("F", datSCLSelW$Treatment[1]))$W
    }
    dat$SCL <- datSCLSelWW
    dat$SCL[which(is.na(dat$SCL))] <- mean(dat$SCL, na.rm = TRUE)
    if(dat$Treatment[1] %in% c("1T:1F", "TF-M", inter2023B) | grepl("_", dat$Treatment[1])){
      dat <- mutate(dat, Arithmetic = (dat$SCC + dat$SCL) / 2)
      dat$Harmonic <- 1 / (0.5 * (1 / dat$SCC) + 0.5 * (1 / dat$SCL))
    }else if(dat$Treatment[1] == "1T:3F"){
      dat$Arithmetic <- (0.5 * dat$SCC + 1.5 * dat$SCL) / 2
      dat$Harmonic <- 1 / (0.25 * (1 / dat$SCC) + 0.75 * (1 / dat$SCL))
    }else{
      dat$Arithmetic <- (0.5 * dat$SCL + 1.5 * dat$SCC) / 2
      dat$Harmonic <- 1 / (0.25 * (1 / dat$SCL) + 0.75 * (1 / dat$SCC))
    }
    return(dat)
  })
  return(dataComb)
}


calcWeeds <- function(data){
  # modWeedsA2 <- glmer(W ~ Treatment + (1|Block), data = dataWeedsA2, family = "Gamma")
  modWeeds <- lmer(W ~ Treatment + (1|Block), data = data)
  
  # variable <- ifelse(relative, "RelativeBiomass", "Biomass")
  # Check for normality and heteroscedasticity
  ## T
  # modbiomass <- lmer(W ~ (1|Block) + Treatment, data = data)
  
  nonNormal <- shapiro.test(residuals(modWeeds))$p.value < 0.05
  unequalVariance <- bartlett.test(residuals(modWeeds) ~ data$Treatment)$p.value < 0.05
  
  newDat <- data
  if(nonNormal | unequalVariance){
    bc <- boxcox(data$W ~ data$Treatment, lambda = seq(-5, 5, 1/10))
    lambda <- bc$x[which.max(bc$y)]
    if(lambda > -0.3 & lambda < 0.3){
      newDat$W <- log(newDat$W)
    }else{
      newDat$W <- newDat$W^(lambda)  
    }
    
    modWeeds <- lmer(newDat$W ~ (1|Block) + Treatment, data = newDat)
    nonNormal <- shapiro.test(residuals(modWeeds))$p.value < 0.05
    unequalVariance <- bartlett.test(residuals(modWeeds) ~ newDat$Treatment)$p.value < 0.05
    if(nonNormal | unequalVariance){
      cat("Cannot guarantee normal distribution and/or equal variances. Please manually check for these assumptions. \n")
    }
  }
  
  cat(paste("p-value = ", summary(aov(newDat$W ~ (1|Block) + Treatment, data = newDat))[[1]]$`Pr(>F)`[1], "\n", sep = ""))
  
  PHWeeds <- emmeans(modWeeds, list(pairwise ~ Treatment), adjust = "tukey")
  CLDWeeds <- cld(PHWeeds$emmeans,
                    Letters = letters,
                    decreasing = TRUE)
  
  return(CLDWeeds)
  
  # dataWeedsA2 <- mutate(dataWeedsA2, group = sapply(dataWeedsA2$Treatment, function(t)ifelse(t == "T" | t == "T+", "Sole triticale", ifelse(t == "F" | t == "F+", "Sole faba", "Intercrop"))))
  
}

calcWeedsGrouped <- function(data){
  modWeeds <- lmer(W ~ Group + (1|Block), data = data)
  nonNormal <- shapiro.test(residuals(modWeeds))$p.value < 0.05
  unequalVariance <- bartlett.test(residuals(modWeeds) ~ data$Group)$p.value < 0.05
  
  newDat <- data
  if(nonNormal | unequalVariance){
    bc <- boxcox(data$W ~ data$Group, lambda = seq(-5, 5, 1/10))
    newDat$WeedBiomass <- newDat$W^(bc$x[which.max(bc$y)])
    modWeeds <- lmer(newDat$W ~ (1|Block) + Group, data = newDat)
    nonNormal <- shapiro.test(residuals(modWeeds))$p.value < 0.05
    unequalVariance <- bartlett.test(residuals(modWeeds) ~ newDat$Group)$p.value < 0.05
    if(nonNormal | unequalVariance){
      cat("Cannot guarantee normal distribution and/or equal variances. Please manually check for these assumptions. \n")
    }
  }
  
  cat(paste("p-value = ", summary(aov(newDat$W ~ (1|Block) + Group, data = newDat))[[1]]$`Pr(>F)`[1], "\n", sep = ""))
  
  PHWeeds <- emmeans(modWeeds, list(pairwise ~ Group), adjust = "tukey")
  CLDWeeds <- cld(PHWeeds$emmeans,
                  Letters = letters,
                  decreasing = TRUE)
  return(CLDWeeds)
}

predictWeeds <- function(dataWeeds){
  sole2022 <- c("Barley", "Rye", "Triticale", "Wheat", "Lupine", "Pea", "Faba")
  sole2023A <- c("T", "T+", "F", "F+")
  sole2023B <- c("T", "T-25", "T-375", "F", "F-25", "F-375")
  sole2024 <- c("T", "F")
  dataWeedsSC <- dataWeeds %>% 
    filter(Treatment %in% c(sole2022, sole2023A, sole2023B, sole2024))
  dataWeedsIC <- dataWeeds %>% 
    filter(!Treatment %in% c(sole2022, sole2023A, sole2023B, sole2024))

  dataWeedsICSplit <- split(dataWeedsIC, dataWeedsIC$Treatment)
  dataWeedsICSplitPred <- calcArithmeticHarmonic(dataWeedsICSplit, dataWeedsSC)
  dataWeedsICPred <- bind_rows(dataWeedsICSplitPred)

  lmArithmetic <- lm(W ~ Arithmetic, data = dataWeedsICPred)
  lmHarmonic <- lm(W ~ Harmonic, data = dataWeedsICPred)
  
  dataWeedsICPredLong <- dataWeedsICPred %>% 
  pivot_longer(cols = c(Arithmetic, Harmonic), names_to = "Model", values_to = "Predicted")
  dataWeedsICPredLong <- mutate(dataWeedsICPredLong, Treatment = factor(Treatment))
  plotWeedPred <- ggplot(data = dataWeedsICPredLong, aes(x = Predicted, y = W, group_by = Model, color = Model)) +
    geom_point(aes(shape = Treatment), size = 4) + 
    geom_abline(intercept = 0, slope = 1,
                size = 2, 
                linetype = "dashed") +
    geom_function(fun = function(x)return(summary(lmArithmetic)$coefficients[1,1] + summary(lmArithmetic)$coefficients[2,1] * x),
                  col = "#F8766D",
                  size = 2) +
    geom_function(fun = function(x)return(summary(lmHarmonic)$coefficients[1,1] + summary(lmHarmonic)$coefficients[2,1] * x),
                  col = "#00BFC4",
                  size = 2) +
    theme_bw(base_size = 30) +
    labs(y = bquote("Observed weed biomass (g "~m^-2~")"),
         x = bquote("Predicted weed biomass (g "~m^-2~")"),
         title = "Predicted and observed weed biomass") +
    scale_shape_manual(values = 1:nlevels(dataWeedsICPredLong$Treatment))
  
  predWeeds2022Long <- dataWeedsICPred %>% 
    pivot_longer(cols = c("W", "Arithmetic", "Harmonic"), names_to = "Type", values_to = "WeedBiomass")
  
  modWeedPred <- lmer(WeedBiomass ~ Type + (1|Block), data = predWeeds2022Long)
  
  nonNormal <- shapiro.test(residuals(modWeedPred))$p.value < 0.05
  unequalVariance <- bartlett.test(residuals(modWeedPred) ~ na.omit(predWeeds2022Long)$Type)$p.value < 0.05
  
  newDat <- predWeeds2022Long
  if(nonNormal | unequalVariance){
    bc <- boxcox(predWeeds2022Long$WeedBiomass ~ predWeeds2022Long$Type, lambda = seq(-5, 5, 1/10))
    modVal <- bc$x[which.max(bc$y)]
    if(modVal < 0.2 & modVal > -0.2){
      newDat$WeedBiomass <- log(newDat$WeedBiomass)
    }else{
      newDat$WeedBiomass <- newDat$WeedBiomass^(bc$x[which.max(bc$y)])
    }
    modWeedPred <- lmer(newDat$WeedBiomass ~ (1|Block) + Type, data = newDat)
    nonNormal <- shapiro.test(residuals(modWeedPred))$p.value < 0.05
    unequalVariance <- bartlett.test(residuals(modWeedPred) ~ na.omit(newDat)$Type)$p.value < 0.05
    if(nonNormal | unequalVariance){
      cat("Cannot guarantee normal distribution and/or equal variances. Please manually check for these assumptions. \n")
    }
  }

  #cat(paste("p-value = ", summary(aov(newDat$WeedBiomass ~ (1|Block) + Treatment, data = newDat))[[1]]$`Pr(>F)`[1], "\n", sep = ""))

  #PHWeeds <- emmeans(modWeedPred, list(pairwise ~ Type), adjust = "tukey")
  #CLDWeeds <- cld(PHWeeds$emmeans,
  #                Letters = letters,
  #                decreasing = TRUE)
  
    return(list(dataWeedsICPred, plotWeedPred))
}
