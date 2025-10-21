library(tidyverse)
library(emmeans)
library(multcomp)
library(multcompView)
library(stringr)
library(xtable)

source("WCFTreatments.R")
source("WCFBiomassFunctions.R")

getCorrectedYield <- function(data, year){
  if(year == "2022"){
    multiplication <- c("Barley_Pea_C" = 2, 
                        "Barley_Lupine_C" = 2, 
                        "Barley_Faba_C" = 2, 
                        "Rye_Pea_C" = 2,
                        "Rye_Lupine_C" = 2,
                        "Rye_Faba_C" = 2,
                        "Triticale_Pea_C" = 2,
                        "Triticale_Lupine_C" = 2,
                        "Triticale_Faba_C" = 2,
                        "Wheat_Pea_C" = 2,
                        "Wheat_Lupine_C" = 2,
                        "Wheat_Faba_C" = 2,
                        "Barley_Pea_L" = 2, 
                        "Barley_Lupine_L" = 2, 
                        "Barley_Faba_L" = 2, 
                        "Rye_Pea_L" = 2,
                        "Rye_Lupine_L" = 2,
                        "Rye_Faba_L" = 2,
                        "Triticale_Pea_L" = 2,
                        "Triticale_Lupine_L" = 2,
                        "Triticale_Faba_L" = 2,
                        "Wheat_Pea_L" = 2,
                        "Wheat_Lupine_L" = 2,
                        "Wheat_Faba_L" = 2,
                        "Barley_C" = 1,
                        "Rye_C" = 1,
                        "Triticale_C" = 1,
                        "Wheat_C" = 1,
                        "Pea_L" = 1,
                        "Lupine_L" = 1,
                        "Faba_L" = 1)
  }else if(year == "2023A"){
    multiplication <- c("1T:1F_T" = 2,
                        "1T:1F_F" = 2,
                        "1T:3F_T" = 4,
                        "1T:3F_F" = 4/3,
                        "3T:1F_T" = 4/3,
                        "3T:1F_F" = 4,
                        "TF-M_T" = 2,
                        "TF-M_F" = 2,
                        "T_T" = 1,
                        "F_F" = 1,
                        "T+_T" = 4/6,
                        "F+_F" = 4/6)
  }else if(year == "2023B"){
    multiplication <- c("T" = 1,
                        "1T:1F" = 2,
                        "1T:1F" = 2,
                        "TF-M" = 2,
                        "TF-M" = 2,
                        "F" = 1,
                        "T-375" = 1,
                        "1T:1F-375" = 2,
                        "1T:1F-375" = 2,
                        "TF-M-375" = 2,
                        "TF-M-375" = 2,
                        "F-375" = 1)
  }else if(year == "2024"){
    multiplication <- c("1T:1F" = 2,
                        "T" = 1,
                        "F" = 1)
  }
  
  cnames <- colnames(data)
  cnamesCorr <- cnames[which(!cnames %in% c("Plot", "Treatment", "TreatmentN", "Block", "Weeds", "Date", "Time"))]
  
  correctedData <- data
  if(year == "2023A"){
    correctedData$YieldTA <- correctedData$YieldTA * multiplication[match(paste(correctedData$Treatment, "_T", sep = ""), names(multiplication))]
    correctedData$YieldTM <- correctedData$YieldTM * multiplication[match(paste(correctedData$Treatment, "_T", sep = ""), names(multiplication))]
    correctedData$YieldT <- correctedData$YieldT * multiplication[match(paste(correctedData$Treatment, "_T", sep = ""), names(multiplication))]
    correctedData$TSWTA <- correctedData$TSWTA * multiplication[match(paste(correctedData$Treatment, "_T", sep = ""), names(multiplication))]
    correctedData$TSWTM <- correctedData$TSWTM * multiplication[match(paste(correctedData$Treatment, "_T", sep = ""), names(multiplication))]
    correctedData$TSWT <- correctedData$TSWT * multiplication[match(paste(correctedData$Treatment, "_T", sep = ""), names(multiplication))]
    correctedData$NFTA <- correctedData$NFTA * multiplication[match(paste(correctedData$Treatment, "_T", sep = ""), names(multiplication))]
    correctedData$NFTM <- correctedData$NFTM * multiplication[match(paste(correctedData$Treatment, "_T", sep = ""), names(multiplication))]
    correctedData$NFT <- correctedData$NFT * multiplication[match(paste(correctedData$Treatment, "_T", sep = ""), names(multiplication))]
    correctedData$NSTA <- correctedData$NSTA * multiplication[match(paste(correctedData$Treatment, "_T", sep = ""), names(multiplication))]
    correctedData$NSTM <- correctedData$NSTM * multiplication[match(paste(correctedData$Treatment, "_T", sep = ""), names(multiplication))]
    correctedData$NST <- correctedData$NST * multiplication[match(paste(correctedData$Treatment, "_T", sep = ""), names(multiplication))]
    correctedData$YieldFA <- correctedData$YieldFA * multiplication[match(paste(correctedData$Treatment, "_F", sep = ""), names(multiplication))]
    correctedData$YieldFM <- correctedData$YieldFM * multiplication[match(paste(correctedData$Treatment, "_F", sep = ""), names(multiplication))]
    correctedData$YieldF <- correctedData$YieldF * multiplication[match(paste(correctedData$Treatment, "_F", sep = ""), names(multiplication))]
    correctedData$TSWFA <- correctedData$TSWFA * multiplication[match(paste(correctedData$Treatment, "_F", sep = ""), names(multiplication))]
    correctedData$TSWFM <- correctedData$TSWFM * multiplication[match(paste(correctedData$Treatment, "_F", sep = ""), names(multiplication))]
    correctedData$TSWF <- correctedData$TSWF * multiplication[match(paste(correctedData$Treatment, "_F", sep = ""), names(multiplication))]
    correctedData$NFFA <- correctedData$NFFA * multiplication[match(paste(correctedData$Treatment, "_F", sep = ""), names(multiplication))]
    correctedData$NFFM <- correctedData$NFFM * multiplication[match(paste(correctedData$Treatment, "_F", sep = ""), names(multiplication))]
    correctedData$NFF <- correctedData$NFF * multiplication[match(paste(correctedData$Treatment, "_F", sep = ""), names(multiplication))]
    correctedData$NSFA <- correctedData$NSFA * multiplication[match(paste(correctedData$Treatment, "_F", sep = ""), names(multiplication))]
    correctedData$NSFM <- correctedData$NSFM * multiplication[match(paste(correctedData$Treatment, "_F", sep = ""), names(multiplication))]
    correctedData$NSF <- correctedData$NSF * multiplication[match(paste(correctedData$Treatment, "_F", sep = ""), names(multiplication))]
  }else if(year %in% c("2023B", "2024")){
    correctedData$YieldT <- correctedData$YieldT * multiplication[match(correctedData$Treatment, names(multiplication))]
    correctedData$TSWT <- correctedData$TSWT * multiplication[match(correctedData$Treatment, names(multiplication))]
    correctedData$NFT <- correctedData$NFT * multiplication[match(correctedData$Treatment, names(multiplication))]
    correctedData$NST <- correctedData$NST * multiplication[match(correctedData$Treatment, names(multiplication))]
    correctedData$YieldF <- correctedData$YieldF * multiplication[match(correctedData$Treatment, names(multiplication))]
    correctedData$TSWF <- correctedData$TSWF * multiplication[match(correctedData$Treatment, names(multiplication))]
    correctedData$NFF <- correctedData$NFF * multiplication[match(correctedData$Treatment, names(multiplication))]
    correctedData$NSF <- correctedData$NSF * multiplication[match(correctedData$Treatment, names(multiplication))]
  }else{
    correctedData$YieldC <- correctedData$YieldC * multiplication[match(paste(correctedData$Treatment, "_C", sep = ""), names(multiplication))]
    correctedData$TSWC <- correctedData$TSWC * multiplication[match(paste(correctedData$Treatment, "_C", sep = ""), names(multiplication))]
    correctedData$NFC <- correctedData$NFC * multiplication[match(paste(correctedData$Treatment, "_C", sep = ""), names(multiplication))]
    correctedData$NSC <- correctedData$NSC * multiplication[match(paste(correctedData$Treatment, "_C", sep = ""), names(multiplication))]
    correctedData$YieldL <- correctedData$YieldL * multiplication[match(paste(correctedData$Treatment, "_L", sep = ""), names(multiplication))]
    correctedData$TSWL <- correctedData$TSWL * multiplication[match(paste(correctedData$Treatment, "_L", sep = ""), names(multiplication))]
    correctedData$NFL <- correctedData$NFL * multiplication[match(paste(correctedData$Treatment, "_L", sep = ""), names(multiplication))]
    correctedData$NSL <- correctedData$NSL * multiplication[match(paste(correctedData$Treatment, "_L", sep = ""), names(multiplication))]
  }
  return(correctedData)
}

getYieldCLD <- function(data, variable){
  data <- na.omit(data)
  # Check for normality and heteroscedasticity
  ## T
  modyield <- lmer(data[,variable][[1]] ~ (1|Block) + Treatment, data = data)
  
  nonNormal <- shapiro.test(residuals(modyield))$p.value < 0.05
  unequalVariance <- bartlett.test(residuals(modyield) ~ data$Treatment)$p.value < 0.05
  
  newDat <- data
  if(nonNormal | unequalVariance){
    bc <- boxcox(data[,variable][[1]] ~ data$Treatment, lambda = seq(-5, 5, 1/10))
    lambda <- bc$x[which.max(bc$y)]
    if(lambda > -0.3 & lambda < 0.3){
      newDat[,variable] <- log(newDat[,variable])
    }else{
      newDat[,variable] <- newDat[,variable]^(bc$x[which.max(bc$y)])  
    }
    modyield <- lmer(newDat[,variable][[1]] ~ (1|Block) + Treatment, data = newDat)
    nonNormal <- shapiro.test(residuals(modyield))$p.value < 0.05
    unequalVariance <- bartlett.test(residuals(modyield) ~ data$Treatment)$p.value < 0.05
    if(nonNormal | unequalVariance){
      cat("Cannot guarantee normal distribution and/or equal variances. Please manually check for these assumptions. \n")
    }
  }
  cat(paste("p-value = ", summary(aov(newDat[,variable][[1]] ~ (1|Block) + Treatment, data = newDat))[[1]]$`Pr(>F)`[1], "\n", sep = ""))
  modyield <- lmer(paste(variable, "~ (1|Block) + Treatment", sep = ""), data = newDat)
  PHYield <- emmeans(modyield, list(pairwise ~ Treatment), adjust = "tukey")
  CLDYield <- cld(PHYield$emmeans,
                    Letters = letters,
                    descending = TRUE)
  
  return(CLDYield)
}

get_yield_cld <- function(data, variable){
  data <- na.omit(data)
  # Check for normality and heteroscedasticity
  ## T
  mod_yield <- lmer(data[,variable][[1]] ~ (1|Block) + Treatment, data = data)
  
  non_normal <- shapiro.test(residuals(mod_yield))$p.value < 0.05
  unequal_variance <- bartlett.test(residuals(mod_yield) ~ data$Treatment)$p.value < 0.05
  
  new_dat <- data
  if(non_normal | unequal_variance){
    bc <- boxcox(data[,variable][[1]] ~ data$Treatment, lambda = seq(-5, 5, 1/10))
    lambda <- bc$x[which.max(bc$y)]
    print(paste0("lambda: ", lambda))
    if(lambda > -0.3 & lambda < 0.3){
      new_dat[,variable] <- log(new_dat[,variable])
    }else{
      new_dat[,variable] <- new_dat[,variable]^(bc$x[which.max(bc$y)])  
    }
    modyield <- lmer(new_dat[,variable][[1]] ~ (1|Block) + Treatment, data = new_dat)
    nonNormal <- shapiro.test(residuals(mod_yield))$p.value < 0.05
    unequalVariance <- bartlett.test(residuals(mod_yield) ~ data$Treatment)$p.value < 0.05
    if(non_normal | unequal_variance){
      cat("Cannot guarantee normal distribution and/or equal variances. Please manually check for these assumptions. \n")
    }
  }
  cat(paste("p-value = ", summary(aov(new_dat[,variable][[1]] ~ (1|Block) + Treatment, data = new_dat))[[1]]$`Pr(>F)`[1], "\n", sep = ""))
  #browser()
  mod_yield <- lmer(paste(variable, "~ (1|Block) + Treatment", sep = ""), data = new_dat)
  ph_yield <- emmeans(mod_yield, list(pairwise ~ Treatment), adjust = "tukey")
  cld_yield <- cld(ph_yield,
                    Letters = letters,
                    sort = FALSE)
  
  return(cld_yield)
}

calc_yield_pLER_data <- function(data, vars){
  sole_crops <- c("Barley", "Rye", "Triticale", "Wheat", "Pea", "Faba", "T", "T+", "F", "F+", "T", "F", "T-375", "F-375")
  data_types <- data %>% 
    dplyr::select(Treatment, Block, all_of(vars))
  
  data_types_IC <- data_types %>% 
    filter(!(Treatment %in% sole_crops))
  
  # Get data of the sole crops at normal density
  data_SC <- data[data$Treatment %in% sole_crops,]
  data_SC <- data_SC[,colSums(is.na(data_SC))<nrow(data_SC)]
  
  agg_list <- lapply(vars[which(vars %in% colnames(data_SC))], function(var){
    return(aggregate(data_SC[,var][[1]] ~ Treatment, data = data_SC, FUN = function(x)mean(x, na.rm = TRUE)) %>% 
             filter(var != 0))
  })
  average_SC <- unlist(lapply(agg_list, function(l)l[,2]))
  names(average_SC) <- unlist(lapply(agg_list, function(l)l[,1]))

  pLERs <- lapply(vars, function(var){
    data_types_IC[,var][[1]] / sapply(1:nrow(data_types_IC), function(i){
      if(var %in% c("YieldC", "YieldT", "YieldTA", "YieldTM", "NFTA", "NFTM", "NFT", "NSTA", "NSTM", "NST", "TSWTA", "TSWTM", "TSWT")){
        name <- strsplit(data_types_IC[i,]$Treatment, split = "_")[[1]][1]
        if(name %in% sole_crops){
          return(average_SC[name])
        }else if(paste("T", substr(data_types_IC[i,]$Treatment, nchar(data_types_IC[i,]$Treatment) - 2, nchar(data_types_IC[i,]$Treatment)), sep = "") %in% sole_crops){
          return(average_SC[paste("T", substr(data_types_IC[i,]$Treatment, nchar(data_types_IC[i,]$Treatment) - 2, nchar(data_types_IC[i,]$Treatment)), sep = "")])
        }else{
          return(average_SC["T"])
        }
      }else{
        name <- strsplit(data_types_IC[i,]$Treatment, split = "_")[[1]][2]
        if(name %in% sole_crops){
          return(average_SC[name])
        }else if(paste("F", substr(data_types_IC[i,]$Treatment, nchar(data_types_IC[i,]$Treatment) - 2, nchar(data_types_IC[i,]$Treatment)), sep = "") %in% sole_crops){
          return(average_SC[paste("F", substr(data_types_IC[i,]$Treatment, nchar(data_types_IC[i,]$Treatment) - 2, nchar(data_types_IC[i,]$Treatment)), sep = "")])
        }else{
          return(average_SC["F"])
        }
      }
    })
  })
  
  pLERs_df <- bind_cols(pLERs)
  colnames(pLERs_df) <- sapply(1:length(pLERs), function(i)paste("pLER", i, sep = ""))
  
  data_types_IC <- bind_cols(data_types_IC, pLERs_df)
  
  return(data_types_IC)
}


# calc_pLER_avg_yield <- function(data, crop = "T"){ 
#   # Check normal distribution
#   non_normal <- shapiro.test(data$pLER)$p.value < 0.05
  
#   cat(paste("Comparing pLER means of " , data$Treatment[1], "\n", sep = ""))
#   if(non_normal){
#     cat("Potential non-normally distributed data. Applying Box-Cox transformation.\n")
#     bc <- boxcox(data$pLER ~ data$Block, lambda = seq(-5, 5, 1/10))
#     lambda <- bc$x[which.max(bc$y)]
#     new_dat <- data
#     if(lambda > -0.3 & lambda < 0.3){
#       new_dat$pLER <- log(new_dat$pLER)
#     }else{
#       new_dat$pLER <- new_dat$pLER^(bc$x[which.max(bc$y)])  
#     }
#     if(shapiro.test(new_dat$pLER)$p.value < 0.05){
#       cat("Data is still potentially non-normally distributed. Resorting to non-parametric test.\n")
#       if(!(TRUE & !(grepl("TM", data$Treatment[1]) | grepl("FM", data$Treatment[1]) | (crop == "F" & data$Treatment[1] == "3T:1F") | (crop == "T" & data$Treatment[1] == "1T:3F") | data$Treatment[1] == "1T:3F_T" | data$Treatment[1] == "3T:1F_F"))){
#         pvalue = wilcox.test(data$pLER, mu = 0.25, alternative = "two.sided")$p.value
#       }else if(data$Treatment[1] == "3T:1F" | data$Treatment[1] == "1T:3F" | data$Treatment[1] == "3T:1F_T" | data$Treatment[1] == "1T:3F_F"){
#         pvalue = wilcox.test(data$pLER, mu = 0.75, alternative = "two.sided")$p.value
#       }else{
#         pvalue = wilcox.test(data$pLER, mu = 0.5, alternative = "two.sided")$p.value
#       }
#     }else{
#       if(!(TRUE & !(grepl("TM", new_dat$Treatment[1]) | grepl("FM", new_dat$Treatment[1]) | (crop == "F" & new_dat$Treatment[1] == "3T:1F") | (crop == "T" & new_dat$Treatment[1] == "1T:3F") | new_dat$Treatment[1] == "1T:3F_T" | new_dat$Treatment[1] == "3T:1F_F"))){
#         pvalue = t.test(new_dat$pLER, mu = 0.25^(bc$x[which.max(bc$y)]) , alternative = "two.sided")$p.value
#       }else if(new_dat$Treatment[1] == "3T:1F" | new_dat$Treatment[1] == "1T:3F" | new_dat$Treatment[1] == "3T:1F_T" | new_dat$Treatment[1] == "1T:3F_F"){
#         pvalue = t.test(new_dat$pLER, mu = 0.75^(bc$x[which.max(bc$y)]) , alternative = "two.sided")$p.value
#       }else{
#         pvalue = t.test(new_dat$pLER, mu = 0.5^(bc$x[which.max(bc$y)]) , alternative = "two.sided")$p.value
#       }
#     }
#   }else{
#     if(!(TRUE & !(grepl("TM", data$Treatment[1]) | grepl("FM", data$Treatment[1]) | (crop == "F" & data$Treatment[1] == "3T:1F") | (crop == "T" & data$Treatment[1] == "1T:3F" | data$Treatment[1] == "1T:3F_T" | data$Treatment[1] == "3T:1F_F")))){
#       pvalue = t.test(data$pLER, mu = 0.25, alternative = "two.sided")$p.value
#     }else if(data$Treatment[1] == "3T:1F" | data$Treatment[1] == "1T:3F" | data$Treatment[1] == "3T:1F_T" | data$Treatment[1] == "1T:3F_F"){
#       pvalue = t.test(data$pLER, mu = 0.75, alternative = "two.sided")$p.value
#     }else{
#       pvalue = t.test(data$pLER, mu = 0.5, alternative = "two.sided")$p.value
#     }
#   }
  
#   data_avg <- tibble(Treatment = data$Treatment[1],
#                     pLER = mean(data$pLER, na.rm = TRUE),
#                     pLERSD = sd(data$pLER, na.rm = TRUE),
#                     PValue = pvalue)
  
#   return(data_avg)
# }

# calc_LER <- function(data){
#   new_data <- data %>% 
#     dplyr::select(Treatment, Block, YieldT, YieldF) %>% 
#     mutate(LER = YieldT / mean(data[data$Treatment == "T","YieldT"][[1]], na.rm = TRUE) + YieldF / mean(data[data$Treatment == "F","YieldF"][[1]], na.rm = TRUE)) %>% 
#     filter(Treatment %in% c("1T:1F"))
  
#   new_data2 <- na.omit(new_data)
#   non_normal <- shapiro.test(new_data2$LER)$p.value < 0.05
#   new_data3 <- new_data2
#   if(non_normal){
#     bc <- boxcox(new_data2$LER ~ new_data2$Block, lambda = seq(-5, 5, 1/10))
#     new_data3$LER <- new_data3$LER^(bc$x[which.max(bc$y)])
#     non_normal <- shapiro.test(new_data3$LER)$p.value < 0.05
#     if(non_normal){
#       cat("Cannot guarantee normal distribution and/or equal variances. Please manually check for these assumptions. \n")
#     }
#   }
  
#   data_avg <- aggregate(LER ~ Treatment, data = new_data3, function(x)mean(x, na.rm = TRUE)) %>% 
#     mutate(SD = aggregate(LER ~ Treatment, data = new_data3, function(x)sd(x, na.rm = TRUE))$LER) %>% 
#     mutate(PValue = t.test(new_data3$LER, mu = 1.0, alternative = "two.sided")$p.value)
#   return(data_avg)
# }


# calcYieldPLERData2022 <- function(data, varC, varL){
#   soleCrops <- c("Barley", "Rye", "Triticale", "Wheat", "Pea", "Faba")
#   dataTypes <- data %>% 
#     dplyr::select(Treatment, Block, all_of(varC), all_of(varL))
  
#   dataTypesIC <- dataTypes %>% 
#     filter(!(Treatment %in% soleCrops))
  
#   # Get data of the sole crops at normal density
#   dataSC <- data[data$Treatment %in% soleCrops,]
#   averageSCC <- aggregate(dataSC[,varC][[1]] ~ Treatment, data = dataSC, FUN = function(x)mean(x, na.rm = TRUE)) %>% 
#     filter(varC != 0)
#   averageSCL <- aggregate(dataSC[,varL][[1]] ~ Treatment, data = dataSC, FUN = function(x)mean(x, na.rm = TRUE)) %>% 
#     filter(varL != 0)
#   averageSC <- c(averageSCC[,2], averageSCL[,2])
#   names(averageSC) <- c(averageSCC$Treatment, averageSCL$Treatment)
  
#   pLERC <- dataTypesIC[,varC][[1]] / sapply(1:nrow(dataTypesIC), function(i)averageSC[strsplit(dataTypesIC[i,]$Treatment, split = "_")[[1]][1]])
#   pLERL <- dataTypesIC[,varL][[1]] / sapply(1:nrow(dataTypesIC), function(i)averageSC[strsplit(dataTypesIC[i,]$Treatment, split = "_")[[1]][2]])
  
#   dataTypesIC <- dataTypesIC %>% 
#     mutate(pLERC = pLERC,
#            pLERL = pLERL)
  
#   return(dataTypesIC)
# }

calcYieldPLERData <- function(data, vars){
  soleCrops <- c("Barley", "Rye", "Triticale", "Wheat", "Pea", "Faba", "T", "T+", "F", "F+", "T", "F", "T-375", "F-375")
  dataTypes <- data %>% 
    dplyr::select(Treatment, Block, all_of(vars))
  
  dataTypesIC <- dataTypes %>% 
    filter(!(Treatment %in% soleCrops))
  
  # Get data of the sole crops at normal density
  dataSC <- data[data$Treatment %in% soleCrops,]
  dataSC <- dataSC[,colSums(is.na(dataSC))<nrow(dataSC)]
  
  aggList <- lapply(vars[which(vars %in% colnames(dataSC))], function(var){
    return(aggregate(dataSC[,var][[1]] ~ Treatment, data = dataSC, FUN = function(x)mean(x, na.rm = TRUE)) %>% 
             filter(var != 0))
  })
  averageSC <- unlist(lapply(aggList, function(l)l[,2]))
  names(averageSC) <- unlist(lapply(aggList, function(l)l[,1]))

  pLERs <- lapply(vars, function(var){
    dataTypesIC[,var][[1]] / sapply(1:nrow(dataTypesIC), function(i){
      if(var %in% c("YieldC", "YieldT", "YieldTA", "YieldTM", "NFTA", "NFTM", "NFT", "NSTA", "NSTM", "NST", "TSWTA", "TSWTM", "TSWT")){
        name <- strsplit(dataTypesIC[i,]$Treatment, split = "_")[[1]][1]
        if(name %in% soleCrops){
          return(averageSC[name])
        }else if(paste("T", substr(dataTypesIC[i,]$Treatment, nchar(dataTypesIC[i,]$Treatment) - 2, nchar(dataTypesIC[i,]$Treatment)), sep = "") %in% soleCrops){
          return(averageSC[paste("T", substr(dataTypesIC[i,]$Treatment, nchar(dataTypesIC[i,]$Treatment) - 2, nchar(dataTypesIC[i,]$Treatment)), sep = "")])
        }else{
          return(averageSC["T"])
        }
      }else{
        name <- strsplit(dataTypesIC[i,]$Treatment, split = "_")[[1]][2]
        if(name %in% soleCrops){
          return(averageSC[name])
        }else if(paste("F", substr(dataTypesIC[i,]$Treatment, nchar(dataTypesIC[i,]$Treatment) - 2, nchar(dataTypesIC[i,]$Treatment)), sep = "") %in% soleCrops){
          return(averageSC[paste("F", substr(dataTypesIC[i,]$Treatment, nchar(dataTypesIC[i,]$Treatment) - 2, nchar(dataTypesIC[i,]$Treatment)), sep = "")])
        }else{
          return(averageSC["F"])
        }
      }
    })
  })
  
  pLERsTibble <- bind_cols(pLERs)
  colnames(pLERsTibble) <- sapply(1:length(pLERs), function(i)paste("pLER", i, sep = ""))
  
  dataTypesIC <- bind_cols(dataTypesIC, pLERsTibble)
  
  return(dataTypesIC)
}

calcPLERAvgYield <- function(data, crop = "T"){ 
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
        if(!(TRUE & !(grepl("TM", dat$Treatment[1]) | grepl("FM", dat$Treatment[1]) | (crop == "F" & dat$Treatment[1] == "3T:1F") | (crop == "T" & dat$Treatment[1] == "1T:3F") | dat$Treatment[1] == "1T:3F_T" | dat$Treatment[1] == "3T:1F_F"))){
          return(wilcox.test(dat$pLER, mu = 0.25, alternative = "two.sided")$p.value)
        }else if(dat$Treatment[1] == "3T:1F" | dat$Treatment[1] == "1T:3F" | dat$Treatment[1] == "3T:1F_T" | dat$Treatment[1] == "1T:3F_F"){
          return(wilcox.test(dat$pLER, mu = 0.75, alternative = "two.sided")$p.value)
        }else{
          return(wilcox.test(dat$pLER, mu = 0.5, alternative = "two.sided")$p.value)
        }
      }else{
        if(!(TRUE & !(grepl("TM", newDat$Treatment[1]) | grepl("FM", newDat$Treatment[1]) | (crop == "F" & newDat$Treatment[1] == "3T:1F") | (crop == "T" & newDat$Treatment[1] == "1T:3F") | newDat$Treatment[1] == "1T:3F_T" | newDat$Treatment[1] == "3T:1F_F"))){
          return(t.test(newDat$pLER, mu = log(0.25), alternative = "two.sided")$p.value)
        }else if(newDat$Treatment[1] == "3T:1F" | newDat$Treatment[1] == "1T:3F" | newDat$Treatment[1] == "3T:1F_T" | newDat$Treatment[1] == "1T:3F_F"){
          return(t.test(newDat$pLER, mu = log(0.75), alternative = "two.sided")$p.value)
        }else{
          return(t.test(newDat$pLER, mu = log(0.5), alternative = "two.sided")$p.value)
        }
      }
    }else{
      dat <- data[[i]]
      if(!(TRUE & !(grepl("TM", dat$Treatment[1]) | grepl("FM", dat$Treatment[1]) | (crop == "F" & dat$Treatment[1] == "3T:1F") | (crop == "T" & dat$Treatment[1] == "1T:3F" | dat$Treatment[1] == "1T:3F_T" | dat$Treatment[1] == "3T:1F_F")))){
        return(t.test(dat$pLER, mu = 0.25, alternative = "two.sided")$p.value)
      }else if(dat$Treatment[1] == "3T:1F" | dat$Treatment[1] == "1T:3F" | dat$Treatment[1] == "3T:1F_T" | dat$Treatment[1] == "1T:3F_F"){
        return(t.test(dat$pLER, mu = 0.75, alternative = "two.sided")$p.value)
      }else{
        return(t.test(dat$pLER, mu = 0.5, alternative = "two.sided")$p.value)
      }
    }
  })
  
  dataAvg <- tibble(Treatment = sapply(data, function(dat)dat$Treatment[1]),
                    pLER = sapply(data, function(dat)mean(dat$pLER, na.rm = TRUE)),
                    pLERSD = sapply(data, function(dat)sd(dat$pLER, na.rm = TRUE)),
                    PValue = pValues)
  
  return(dataAvg)
}




### new ###

get_cld <- function(df, yield_col, group_var, treatments) {
  df_sub <- df %>% filter(Treatment %in% treatments)
  form <- as.formula(paste0(yield_col, " ~ Treatment + (1|Block)"))
  model <- lmer(form, data = df_sub, REML = FALSE)
  emm <- emmeans(model, list(pairwise ~ Treatment), adjust = "tukey")
  cld_df <- cld(emm,
                Letters = letters,
                decreasing = TRUE)
  return(cld_df)
}

get_cld2 <- function(df, yield_col, group_var, treatments) {
  df_sub <- df %>% filter(Treatment %in% treatments)
  form <- as.formula(paste0(yield_col, " ~ Treatment + (1|Block)"))
  model <- lmer(form, data = df_sub, REML = FALSE)
  emm <- emmeans(model, list(pairwise ~ Treatment), adjust = "tukey")
  cld_df <- cld(emm,
                Letters = letters,
                decreasing = TRUE)
  return(cld_df)
}

prepare_year_data <- function(data, year_label, yield_col) {
  data %>%
    filter(Treatment %in% c("T", "1T:1F", "F")) %>%
    mutate(
      Year = year_label,
      TreatmentGroup = paste0(Treatment, "-", year_label),
      Species = case_when(
        Treatment == "T" ~ "Triticale",
        Treatment == "F" ~ "Faba bean",
        Treatment == "1T:1F" ~ "Intercrop"
      )
    ) %>%
    group_by(TreatmentGroup, Species) %>%
    summarise(
      mean_yield = mean(!!sym(yield_col), na.rm = TRUE),
      sd_yield = sd(!!sym(yield_col), na.rm = TRUE),
      n = sum(!is.na(!!sym(yield_col))),
      .groups = "drop"
    )
}
