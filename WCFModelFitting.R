# Logistic function with binomial distribution to fit on emergence data
logisticBinomial <- function(par, x, z, n){
  mu <- (exp(par["r"]*(x - par["h"])))/(1 + exp(par["r"]*(x - par["h"])))
  nll <- -sum(dbinom(z, size = n, prob = mu, log = TRUE))
  return(nll)
}

# Logistic with poisson distribution
logisticPoisson <- function(par, x, z){
  mu <- (exp(par["r"]*x - par["h"]))/(1 + exp(par["r"]*(x - par["h"])))
  nll <- -sum(dpois(z, lambda = mu, log = TRUE))
  return(nll)
}

# Logistic function with beta-binomial distribution
logisticBetaBinomial <- function(par, x, z, n){
  mu <- (exp(par["r"]*(x - par["h"])))/(1 + exp(par["r"]*(x - par["h"])))
  nll <- -sum(dbetabinom(z, size = n, prob = mu, theta = par["theta"], log = TRUE))
  return(nll)
}

# Logistic function with normal distribution
logisticNormal <- function(par, x, z){
  mu <- (exp(par["r"]*(x - par["h"])))/(1 + exp(par["r"]*(x - par["h"])))
  nll <- -sum(dnorm(z, mean = mu, sd = par["sd"], log = TRUE))
  return(nll)
}

# Logistic function with log-normal distribution
logisticLogNormal <- function(par, x, z){
  mu <- (exp(par["r"]*(x - par["h"])))/(1 + exp(par["r"]*(x - par["h"])))
  nll <- -sum(dlnorm(z, meanlog = log(mu), sdlog = log(par["sd"]), log = TRUE))
  return(nll)
}

# Logistic function with normal distribution
tLogisticNormal <- function(par, x, z, dMax = NA){
  mu <- par["d"] * ((exp(par["r"]*(x - par["h"])))/(1 + exp(par["r"]*(x - par["h"]))))
  nll <- -sum(dnorm(z, mean = mu, sd = par["sd"], log = TRUE))
  if(!is.na(dMax)){
    if(par["d"] > dMax){
      return(999999999)
    }else{
      return(nll)
    }
  }
}

# Logistic function with log-normal distribution
tLogisticLogNormal <- function(par, x, z, dMax = NA){
  mu <- par["d"] * ((exp(par["r"]*(x - par["h"])))/(1 + exp(par["r"]*(x - par["h"]))))
  nll <- -sum(dlnorm(z, meanlog = log(mu), sdlog = log(par["sd"]), log = TRUE))
  if(!is.na(dMax)){
    if(par["d"] > dMax){
      return(999999999)
    }else{
      return(nll)
    }
  }
}

# Logistic function with gamma distribution
logisticGamma <- function(par, x, z){
  mu <- (exp(par["r"]*(x - par["h"])))/(1 + exp(par["r"]*(x - par["h"])))
  nll <- -sum(dgamma(z, shape = par["shape"], scale = mu / par["shape"], log = TRUE))
  return(nll)
}

# Logistic function with beta distribution
logisticBeta <- function(par, x, z){
  mu <- (exp(par["r"]*(x - par["h"])))/(1 + exp(par["r"]*(x - par["h"])))
  nll <- -sum(dbeta(z, shape1 = par["shape1"], shape2 = par["shape2"], log = TRUE))
  return(nll)
}

# Logistic function with gamma distribution
tLogisticGamma <- function(par, x, z, dMax = NA){
  mu <- par["d"] * ((exp(par["r"]*(x - par["h"])))/(1 + exp(par["r"]*(x - par["h"]))))
  nll <- -sum(dgamma(z, shape = par["shape"], scale = mu / par["shape"], log = TRUE))
  if(!is.na(dMax)){
    if(par["d"] > dMax){
      return(999999999)
    }else{
      return(nll)
    }
  }
}

# Logistic function with beta distribution
tLogisticBeta <- function(par, x, z, dMax = NA){
  mu <- par["d"] * ((exp(par["r"]*(x - par["h"])))/(1 + exp(par["r"]*(x - par["h"]))))
  nll <- -sum(dbeta(z, shape1 = par["shape1"], shape2 = par["shape2"], log = TRUE))
  if(!is.na(dMax)){
    if(par["d"] > dMax){
      return(999999999)
    }else{
      return(nll)
    }
  }
}

# Zero-inflated gamma distribution
dzigamma <- function(x, mu, shape, scale, zprob, log = FALSE){
  ifelse(x == 0, zprob + (1 - zprob) * dbinom(0, prob = mu, size = 1, log = log), (1 - zprob) * dgamma(x, shape = shape, scale = scale, log = log))
}

# Logistic function with zero-inflated gamma distribution
logisticZIGamma <- function(par, x, z){
  mu <- (exp(par["r"]*(x - par["h"])))/(1 + exp(par["r"]*(x - par["h"])))
  muzprob <- (exp(-par["r"]*(x - par["h"])))/(1 + exp(-par["r"]*(x - par["h"])))
  nll <- -sum(dzigamma(z, mu = mu, shape = as.numeric(par["shape"]), scale = mu / par["shape"], zprob = muzprob, log = TRUE))
  return(nll)
}

# monomolecular with normal distribution
monoNormal <- function(par, x, z){
  mu <- 1-exp(-par["r"]*x)
  nll <- -sum(dnorm(z, mean = mu, sd = par["sd"], log = TRUE))
  return(nll)
}

# monomolecular with log-normal distribution
monoLogNormal <- function(par, x, z){
  mu <- 1-exp(-par["r"]*x)
  nll <- -sum(dlnorm(z, meanlog = log(mu), sdlog = log(par["sd"]), log = TRUE))
  return(nll)
}

# exponential model with normal distribution
expNormal <- function(par, x, z){
  mu <- par["r"] * x ^ par["p"]
  nll <- -sum(dnorm(z, mean = mu, sd = par["sd"], log = TRUE))
  return(nll)
}

allNLLModels <- list("logisticNormal" = logisticNormal,
                     "logisticLogNormal" = logisticLogNormal,
                     "monoNormal" = monoNormal,
                     "monoLogNormal" = monoLogNormal,
                     "logisticPoisson" = logisticPoisson,
                     "logisticBinomial" = logisticBinomial,
                     "logisticBetaBinomial" = logisticBetaBinomial,
                     "logisticGamma" = logisticGamma,
                     "logisticBeta" = logisticBeta,
                     "logisticZIGamma" = logisticZIGamma,
                     "expNormal" = expNormal)

# Set functions
logistic <- function(par, x){
  return((exp(par["r"]*(x - par["h"])))/(1 + exp(par["r"]*(x - par["h"]))))
}
tLogistic <- function(par, x){
  return(par["d"] * ((exp(par["r"]*(x - par["h"])))/(1 + exp(par["r"]*(x - par["h"])))))
}
monomolecular <- function(par, x){
  return(1-exp(-par["r"]*x))
}
exponential <- function(par, x){
  return(par["r"] * x ^ par["p"])
}

allDeterministicFunctions <- list("logistic" = logistic,
                                  "monomolecular" = monomolecular,
                                  "exponential" = exponential,
                                  "tLogistic" = tLogistic)

# Create functions to fit models
fitDataToModel <- function(data, model, startParameters, method = "Nelder-Mead", maxit = 1000, dMax = NA){
  # browser()
  tryCatch({
    X <- data[,1][[1]]
    Z <- data[,2][[1]]
    if(any(is.na(Z))){
      X <- X[-which(is.na(Z))]
      Z <- Z[-which(is.na(Z))]
    }
    if("n" %in% names(formals(model))){
      N <- data[,3][[1]]
      Opt <- optim(par = startParameters, fn = model, x = X, z = Z, n = N,
                   method = method,
                   control = list(maxit = maxit))
    }else{
      if(!is.na(dMax)){
        Opt <- optim(par = startParameters, fn = model, x = X, z = Z, dMax = dMax,
                     method = method,
                     control = list(maxit = maxit))
      }else{
        Opt <- optim(par = startParameters, fn = model, x = X, z = Z,
                     method = method,
                     control = list(maxit = maxit))
      }
    }
    
    return(Opt)
  },
  error=function(cond){
    message(paste("Error in treatment ", data[1,1][[1]], "\n", sep = ""))
    message("here is the original error message: ")
    message(cond)
    message("\n")
    message(".")
    message("\n")
    return(NA)
  })
}

# Create function to calculate 95% CIs
calculate95CI <- function(parStartVec, optVec, model, optNLL, x, z, dMax = NA){
  tryCatch({
    if(any(is.na(z))){
      x <- x[-which(is.na(z))]
      z <- z[-which(is.na(z))]
    }
    par <- c(parStartVec, optVec[[1]][1] * 100)
    names(par)[length(par)] <- names(optVec)
    prof <- numeric(length = length(optVec[[1]]))
    coefs <- list(rep(list(), times = length(optVec[[1]])))
    for(i in 1:length(optVec[[1]])){
      tryCatch({
        optVal <- optVec[[1]][[i]]
        names(optVal) <- names(optVec)[length(optVec)]
        if(is.na(dMax)){
          opt <- optim(par = parStartVec, fn = model,
                       x = x,
                       z = z,
                       opt = optVal)
        }else{
          opt <- optim(par = parStartVec, fn = model, dMax = dMax,
                       x = x,
                       z = z,
                       opt = optVal)
        }

      },
      error = function(cond){
        opt <- list(par = NA, value = NA)
      })
      
      coefs[[i]] <- opt$par
      prof[i] <- opt$value
    }
    profLower <- prof[1:which.min(prof)]
    vecLower <- optVec[[1]][1:which.min(prof)]
    profHigher <- prof[which.min(prof):length(prof)]
    vecHigher <- optVec[[1]][which.min(prof):length(prof)]
    CI95Low <- approx(profLower, vecLower, xout = optNLL + qchisq(0.95, length(coefs[[1]]))/2)
    CI95High <- approx(profHigher, vecHigher, xout = optNLL + qchisq(0.95, length(coefs[[1]]))/2)
    
    return(c(CI95Low$y, CI95High$y))
  },
  error=function(cond){
    return(NA)
  })
}

fitAndPlotData <- function(data, dataType, weedSuffix, xLab, yLab, dir, fitOnAverage = FALSE){
  # browser()
  dataTreatment <- split(data, data$Treatment)
  dataTreatment <- lapply(dataTreatment, function(dat){
    datNew <- dat %>% 
      select(Time, Val)
    return(datNew)
  })
  
  if(fitOnAverage){
    dataTreatment <- lapply(dataTreatment, function(Data){
      Average <- aggregate(Data$Val, list(Data$Time), function(x)mean(x, na.rm = TRUE))
      Data$Val <- Average$x[match(Data$Time, Average$Group.1)]
      return(Data)
    })
  }
  
  
  fitsLogisticNormal <- lapply(dataTreatment, function(dat){
    fit <- fitDataToModel(data = dat, model = logisticNormal, startParameters = c("r" = 0.005, "h" = 650, "sd" = 0.1))
    # print(fit$convergence)
    return(fit)
  })
  AICLogisticNormal <- sum(unlist(lapply(fitsLogisticNormal, function(fits)fits$value)))
  
  fitsLogisticLogNormal <- lapply(dataTreatment, function(dat){
    fit <- fitDataToModel(data = dat, model = logisticLogNormal, startParameters = c("r" = 0.005, "h" = 650, "sd" = exp(0.1)))
    # print(fit$convergence)
    return(fit)
  })
  AICLogisticLogNormal <- sum(unlist(lapply(fitsLogisticLogNormal, function(fits)fits$value)))
  
  fitsTLogisticNormal <- lapply(dataTreatment, function(dat){
    fit <- fitDataToModel(data = dat, model = tLogisticNormal, startParameters = c("r" = 0.005, "h" = 650, "sd" = 0.1, "d" = 0.85), dMax = 1)
    # print(fit$convergence)
    return(fit)
  })
  AICTLogisticNormal <- sum(unlist(lapply(fitsTLogisticNormal, function(fits)fits$value)))
  
  fitsTLogisticLogNormal <- lapply(dataTreatment, function(dat){
    fit <- fitDataToModel(data = dat, model = tLogisticLogNormal, startParameters = c("r" = 0.005, "h" = 650, "sd" = exp(0.1), "d" = 0.85), dMax = 1)
    # print(fit$convergence)
    return(fit)
  })
  AICTLogisticLogNormal <- sum(unlist(lapply(fitsTLogisticLogNormal, function(fits)fits$value)))
  
  fitsLogisticGamma <- lapply(dataTreatment, function(dat){
    fit <- fitDataToModel(data = dat, model = logisticGamma, startParameters = c("r" = 0.005, "h" = 650, "shape" = 1))
    # print(fit$convergence)
    return(fit)
  })
  AICLogisticGamma <- sum(unlist(lapply(fitsLogisticGamma, function(fits)fits$value)))
  
  fitsLogisticBeta <- lapply(dataTreatment, function(dat){
    fit <- fitDataToModel(data = dat, model = logisticBeta, startParameters = c("r" = 0.005, "h" = 650, "shape1" = 0.5, "shape2" = 0.5))
    # print(fit$convergence)
    return(fit)
  })
  AICLogisticBeta <- sum(unlist(lapply(fitsLogisticBeta, function(fits)fits$value)))
  
  fitsTLogisticGamma <- lapply(dataTreatment, function(dat){
    fit <- fitDataToModel(data = dat, model = tLogisticGamma, startParameters = c("r" = 0.005, "h" = 650, "shape" = 1, "d" = 0.85), dMax = 1)
    # print(fit$convergence)
    return(fit)
  })
  AICTLogisticGamma <- sum(unlist(lapply(fitsTLogisticGamma, function(fits)fits$value)))
  
  fitsTLogisticBeta <- lapply(dataTreatment, function(dat){
    fit <- fitDataToModel(data = dat, model = tLogisticBeta, startParameters = c("r" = 0.005, "h" = 650, "shape1" = 0.5, "shape2" = 0.5, "d" = 0.85), dMax = 1)
    # print(fit$convergence)
    return(fit)
  })
  AICTLogisticBeta <- sum(unlist(lapply(fitsTLogisticBeta, function(fits)fits$value)))
  
  allAICs <- c("logisticNormal" = AICLogisticNormal,
               "logisticLogNormal" = AICLogisticLogNormal,
               "tLogisticNormal" = AICTLogisticNormal,
               "tLogisticLogNormal" = AICTLogisticLogNormal,
               "logisticGamma" = AICLogisticGamma,
               "logisticBeta" = AICLogisticBeta,
               "tLogisticGamma" = AICTLogisticGamma,
               "tLogisticBeta" = AICTLogisticBeta)
  
  # The lowest AIC belongs to the tLogisticGamma function
  sort(allAICs)
  
  bestModel <- names(sort(allAICs))[1]
  
  # Check all plots
  parLogisticNormal     <- lapply(fitsLogisticNormal, function(fits)fits$par)
  parLogisticLogNormal  <- lapply(fitsLogisticLogNormal, function(fits)fits$par)
  parTLogisticNormal    <- lapply(fitsTLogisticNormal, function(fits)fits$par)
  parTLogisticLogNormal <- lapply(fitsTLogisticLogNormal, function(fits)fits$par)
  parLogisticGamma      <- lapply(fitsLogisticGamma, function(fits)fits$par)
  parLogisticBeta       <- lapply(fitsLogisticBeta, function(fits)fits$par)
  parTLogisticGamma     <- lapply(fitsTLogisticGamma, function(fits)fits$par)
  parTLogisticBeta      <- lapply(fitsTLogisticBeta, function(fits)fits$par)
  
  plotFit <- function(data, par, fun){
    plot <- ggplot() +
      geom_point(data = data, aes(x = Time, y = Val)) +
      geom_function(fun = function(x)allDeterministicFunctions[[fun]](par, x),
                    size = 1) + 
      labs(title = par["d"])
    return(plot)
  }
  
  ## Use below lines to check the fits
  # lapply(1:length(parLogisticNormal), function(i){
  #   plotFit(dataTreatment[[i]], parLogisticNormal[[i]], "logistic")
  # })
  # lapply(1:length(parLogisticLogNormal), function(i){
  #   plotFit(dataTreatment[[i]], parLogisticLogNormal[[i]], "logistic")
  # })
  # lapply(1:length(parTLogisticNormal), function(i){
  #   plotFit(dataTreatment[[i]], parTLogisticNormal[[i]], "tLogistic")
  # })
  # lapply(1:length(parTLogisticLogNormal), function(i){
  #   plotFit(dataTreatment[[i]], parTLogisticLogNormal[[i]], "tLogistic")
  # })
  # lapply(1:length(parLogisticGamma), function(i){
  #   plotFit(dataTreatment[[i]], parLogisticGamma[[i]], "logistic")
  # })
  # lapply(1:length(parLogisticBeta), function(i){
  #   plotFit(dataTreatment[[i]], parLogisticBeta[[i]], "logistic")
  # })
  # lapply(1:length(parTLogisticGamma), function(i){
  #   plotFit(dataTreatment[[i]], parTLogisticGamma[[i]], "tLogistic")
  # })
  # lapply(1:length(parTLogisticBeta), function(i){
  #   plotFit(dataTreatment[[i]], parTLogisticBeta[[i]], "tLogistic")
  # })
  
  allFits <- list("logisticNormal" = fitsLogisticNormal, 
                  "logisticLogNormal" = fitsLogisticLogNormal, 
                  "tLogisticNormal" = fitsTLogisticNormal, 
                  "tLogisticLogNormal" = fitsTLogisticLogNormal,
                  "logisticGamma" = fitsLogisticGamma, 
                  "tLogisticGamma" = fitsTLogisticGamma, 
                  "logisticBeta" = fitsLogisticBeta, 
                  "tLogisticBeta" = fitsTLogisticBeta)
  
  parameters <- names(allFits[[bestModel]][[1]]$par)
  parameterCIs <- rep(NA, times = 3 * length(parameters))
  names(parameterCIs) <- unlist(lapply(parameters, function(par)c(par, paste(par, "Lower", sep = ""), paste(par, "Upper", sep = ""))))
  parameterCIs <- rep(list(parameterCIs), length(fitsTLogisticGamma))
  
  cI <- function(par, x, z, opt, dMax){
    pars <- par
    pars[names(opt)] <- opt
    nll <- tLogisticGamma(pars, x, z, dMax)
    return(nll)
  }
  
  for(i in 1:length(fitsTLogisticGamma)){
    parameterVecs <- lapply(fitsTLogisticGamma, function(fit){
      parStart <- fit$par / 100
      parEnd <- fit$par * 20
      parVec <- list(seq(parStart[1], parEnd[1], length.out = 400), seq(parStart[2], parEnd[2], length.out = 400), seq(parStart[3], parEnd[3], length.out = 400), seq(parStart[4], parEnd[4], length.out = 400))
      names(parVec) <- names(fit$par)
      return(parVec)
    })
    optPar <- fitsTLogisticGamma[[i]]$par
    parameterCIs[[i]][names(optPar)] <- optPar
    parameterCIs[[i]][unlist((lapply(names(optPar), function(n)c(paste(n, "Lower", sep = ""), paste(n, "Upper", sep = "")))))] <-
      unlist(lapply(1:length(parameterVecs[[i]]), function(j){
        t(tibble(calculate95CI(parStartVec = optPar[-which(names(optPar) == names(fitsTLogisticGamma[[i]]$par)[[j]])],
                               optVec = parameterVecs[[i]][j],
                               model = cI,
                               optNLL = fitsTLogisticGamma[[i]]$value,
                               x = dataTreatment[[i]][,1][[1]],
                               z = dataTreatment[[i]][,2][[1]],
                               dMax = 1.0)))
      }))
    
  }
  names(parameterCIs) <- names(dataTreatment)
  
  # Write to table
  CIOut <- tibble(Treatment = names(parameterCIs)) %>% 
    bind_cols(bind_rows(parameterCIs))
  write.table(CIOut, paste(dir, "data_analysis/", dataType, "/", dataType, "_", weedSuffix, "_", bestModel, "_parameters.txt", sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  # Create triple plots for intercrops
  dataTreatment <- lapply(1:length(dataTreatment), function(i){
    dat <- dataTreatment[[i]]
    treatment <- names(dataTreatment)[i]
    dat <- dat %>% 
      mutate(Treatment = treatment)
  })
  
  dataTreatment <- lapply(dataTreatment, function(Data){
    Average <- aggregate(Data$Val, list(Data$Time), function(x)mean(x, na.rm = TRUE))
    Data$Average <- Average$x[match(Data$Time, Average$Group.1)]
    return(Data)
  })
  
  geomPoints <- lapply(dataTreatment, function(d){
    return(geom_point(data = d, aes(x = Time, y = Average, col = Treatment[1]), size = 1))
  })
  
  geomFunctions <- lapply(1:length(fitsTLogisticGamma), function(i){
    o <- fitsTLogisticGamma[[i]]
    return(geom_function(fun = function(x)(exp(o$par["r"]*(x - o$par["h"])))/(1 + exp(o$par["r"]*(x - o$par["h"]))),
                         aes(col = dataTreatment[[i]]$Treatment[1]),
                         size = 1.5))
  })
  
  # Create a vector of the most distinctive colours
  qualitativeColors = brewer.pal.info[brewer.pal.info$category == 'qual',]
  colourVector = unlist(mapply(brewer.pal, qualitativeColors$maxcolors, rownames(qualitativeColors)))
  
  values <- colourVector[1:length(fitsTLogisticGamma)]
  names(values) <- unname(sapply(dataTreatment, function(d)d$Treatment[1]))
  labels <- unname(sapply(dataTreatment, function(d)d$Treatment[1]))
  
  ### Get all triple combinations
  treatmentsToCompare <- unique(data$Treatment)
  treatmentsToCompare <- treatmentsToCompare[-which(treatmentsToCompare %in% c("T+", "F+"))] # The increased density sole crops are not compared here
  
  triples <- lapply(treatmentsToCompare, function(treatment){
    if(!treatment %in% c("T", "F")){
      return(c(treatment, "T", "F"))
    }else{
      return(NA)
    }
  })
  triples <- triples[!is.na(triples)]
  names(dataTreatment) <- lapply(dataTreatment, function(l)l$Treatment[1])
  names(geomPoints) <- names(dataTreatment)
  names(geomFunctions) <- names(dataTreatment)
  
  triplePlots <- lapply(triples, function(triple){
    plotMultipleGraphs(geomPoints[triple], geomFunctions[triple], triple, colourVector, XLab = xLab, YLab = yLab, Title = "")
  })
  
  doubles <- list(c("T", "T+"), c("F", "F+"))
  
  doublePlots <- lapply(doubles, function(double){
    plotMultipleGraphs(geomPoints[double], geomFunctions[double], double, c(colourVector[2], colourVector[1]), XLab = xLab, YLab = yLab, Title = "")
  })
  
  combinedPlot <- plotMultipleGraphs(geomPoints[treatmentsToCompare], geomFunctions[treatmentsToCompare], treatmentsToCompare, c(colourVector[1], colourVector[2], colourVector[5], colourVector[4], colourVector[3], colourVector[6]), XLab = xLab, YLab = yLab, Title = "")
  
  pdf(paste(dir, "data_analysis/", dataType, "/", dataType, "_", weedSuffix, "_", bestModel, "_graphs.pdf", sep = ""))
  print(triplePlots)
  print(combinedPlot)
  print(doublePlots)
  dev.off()
  
  return(sort(allAICs))
}


fitCombinations <- function(data, triple, model, startParameters, dMax = 1.0){
  if(is.na(dMax)){
    fitIC <- fitDataToModel(data = data[[triple[1]]], model = model, startParameters = startParameters[[1]])
    fitC <- fitDataToModel(data = data[[triple[2]]], model = model, startParameters = startParameters[[2]])
    fitL <- fitDataToModel(data = data[[triple[3]]], model = model, startParameters = startParameters[[3]])
    
    print(plotFit(data[[triple[1]]], fitIC$par, "tLogistic"))
    print(plotFit(data[[triple[2]]], fitC$par, "tLogistic"))
    print(plotFit(data[[triple[3]]], fitL$par, "tLogistic"))
    
    AICICCL <- ( 2 * length(fitIC$par) + 2 * fitIC$value) + (2 * length(fitC$par) + 2 * fitC$value) + (2 * length(fitL$par) + 2 * fitL$value)
    
    fitCIC <- fitDataToModel(data = bind_rows(data[[triple[1]]], data[[triple[2]]]), model = model, startParameters = startParameters[[4]])
    fitL <- fitDataToModel(data = data[[triple[3]]], model = model, startParameters = startParameters[[3]])
    
    print(plotFit(bind_rows(data[[triple[1]]], data[[triple[2]]]), fitCIC$par, "tLogistic"))
    print(plotFit(data[[triple[3]]], fitL$par, "tLogistic"))
    
    AICCICL <- ( 2 * length(fitCIC$par) + 2 * fitCIC$value) + (2 * length(fitL$par) + 2 * fitL$value)
    
    fitLIC <- fitDataToModel(data = bind_rows(data[[triple[1]]], data[[triple[3]]]), model = model, startParameters = startParameters[[5]])
    fitC <- fitDataToModel(data = data[[triple[2]]], model = model, startParameters = startParameters[[2]])
    
    print(plotFit(bind_rows(data[[triple[1]]], data[[triple[3]]]), fitLIC$par, "tLogistic"))
    print(plotFit(data[[triple[2]]], fitC$par, "tLogistic"))
    
    AICCLIC <- ( 2 * length(fitC$par) + 2 * fitC$value) + (2 * length(fitLIC$par) + 2 * fitLIC$value)
  }else{
    fitIC <- fitDataToModel(data = data[[triple[1]]], model = model, startParameters = startParameters[[1]], dMax = dMax)
    fitC <- fitDataToModel(data = data[[triple[2]]], model = model, startParameters = startParameters[[2]], dMax = dMax)
    fitL <- fitDataToModel(data = data[[triple[3]]], model = model, startParameters = startParameters[[3]], dMax = dMax)
    
    print(plotFit(data[[triple[1]]], fitIC$par, "tLogistic"))
    print(plotFit(data[[triple[2]]], fitC$par, "tLogistic"))
    print(plotFit(data[[triple[3]]], fitL$par, "tLogistic"))
    
    AICICCL <- ( 2 * length(fitIC$par) + 2 * fitIC$value) + (2 * length(fitC$par) + 2 * fitC$value) + (2 * length(fitL$par) + 2 * fitL$value)
    
    fitCIC <- fitDataToModel(data = bind_rows(data[[triple[1]]], data[[triple[2]]]), model = model, startParameters = startParameters[[4]], dMax = dMax)
    fitL <- fitDataToModel(data = data[[triple[3]]], model = model, startParameters = startParameters[[3]], dMax = dMax)
    
    print(plotFit(bind_rows(data[[triple[1]]], data[[triple[2]]]), fitCIC$par, "tLogistic"))
    print(plotFit(data[[triple[3]]], fitL$par, "tLogistic"))
    
    AICCICL <- ( 2 * length(fitCIC$par) + 2 * fitCIC$value) + (2 * length(fitL$par) + 2 * fitL$value)
    
    fitLIC <- fitDataToModel(data = bind_rows(data[[triple[1]]], data[[triple[3]]]), model = model, startParameters = startParameters[[5]], dMax = dMax)
    fitC <- fitDataToModel(data = data[[triple[2]]], model = model, startParameters = startParameters[[2]], dMax = dMax)
    
    print(plotFit(bind_rows(data[[triple[1]]], data[[triple[3]]]), fitLIC$par, "tLogistic"))
    print(plotFit(data[[triple[2]]], fitC$par, "tLogistic"))
  
    AICCLIC <- ( 2 * length(fitC$par) + 2 * fitC$value) + (2 * length(fitLIC$par) + 2 * fitLIC$value)
  }
  return(c(AICICCL, AICCICL, AICCLIC))
}
