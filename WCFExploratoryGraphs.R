## Author: David Kottelenberg
## Institute: Wageningen University & Research
## Last edited: 25-04-2023
## Summary: Definition of functions to generate exploratory graphs of input data of a specific type

library(tidyverse)
library(RColorBrewer)

# Set N treatments
NTreatments <- 19

##
## Plotting functions
##

plotGraph = function(DataPlot, title, i, XLab = "X", YLab = "Y", XMin = NULL, XMax = NULL, YMin = NULL, YMax = NULL){
  ggplot(data = DataPlot) + 
    geom_point(aes(x = Time, y = unname(unlist(DataPlot[,ncol(DataPlot)-2]))), col = "red4", shape = 1) + 
    geom_line(aes(x = Time, y = unname(unlist(DataPlot[,ncol(DataPlot)-1]))), col = "red4", size = 1) + 
    theme_classic(base_size = 25) +
    labs(title = paste(title, i, sep = " "),
         subtitle = gsub("_", "-", DataPlot$Treatment[[1]]),
         y = YLab,
         x = XLab) + 
    aes(xmin = XMin,
        ymin = YMin,
        xmax = XMax,
        ymax = YMax)
}

plotMultipleGraphs <- function(DataBarsList, DataLinesList, DataNames, ColourVector, XLab = "X", YLab = "Y", Title = NULL, TitleSuffix = "", LegendName = "Treatment", XMin = NULL, XMax = NULL, YMin = NULL, YMax = NULL){
  Values <- ColourVector[1:length(DataBarsList)]
  names(Values) <- DataNames
  Labels <- DataNames
  names(Labels) <-  DataNames
  if(is.null(Title)){
    Title <- paste(Labels[length(Labels) - 1], " & ", Labels[length(Labels)], TitleSuffix, sep = "")
  }
  
  Plot <- ggplot() +
    DataLinesList +
    DataBarsList +
    theme_classic(base_size = 25) +
    scale_color_manual(labels = Labels,
                       values = Values) +
    ylim(c(0, 1.0)) +
    labs(title = Title,
         y = YLab,
         x = XLab,
         color = LegendName) + 
    aes(xmin = XMin,
        ymin = YMin,
        xmax = XMax,
        ymax = YMax)
  return(Plot)
}

createExploratoryGraphs <- function(Data, Variable, Dir, OutputDir, ReturnTreatmentAverages = TRUE, ReturnAllTreatments = TRUE, ReturnTriples = TRUE, PDF = TRUE, XLab = "X", YLab = "Y", XMin = NULL, XMax = NULL, YMin = NULL, YMax = NULL){
  Output <- list()
  # Get averages of treatments
  TreatmentAverages <- as_tibble(aggregate(Data[,ncol(Data)], list(Treatment = Data$Treatment, TreatmentN = Data$TreatmentN, Time = Data$Time), function(x)mean(x, na.rm = TRUE)))
  TreatmentSDs <- as_tibble(aggregate(Data[,ncol(Data)], list(Treatment = Data$Treatment, TreatmentN = Data$TreatmentN, Time = Data$Time), function(x)sd(x, na.rm = TRUE)))
  TreatmentCount <- as_tibble(aggregate(Data[,ncol(Data)], list(Treatment = Data$Treatment, TreatmentN = Data$TreatmentN, Time = Data$Time), length))
  # DataTableTime <- unlist(lapply(1:length(unique(Data$Time)), function(i)table(Data$Plot[which(Data$Time == unique(Data$Time)[i])])))
  DataTreatmentAverages <- Data %>% 
    arrange(Time, TreatmentN) %>% 
    mutate(Average = unlist(lapply(1:nrow(TreatmentAverages), function(i)rep(unlist(TreatmentAverages[,ncol(TreatmentAverages)])[i], times = TreatmentCount[i,ncol(TreatmentCount)]))),
           SD = unlist(lapply(1:nrow(TreatmentSDs), function(i)rep(unlist(TreatmentSDs[,ncol(TreatmentSDs)])[i], times = TreatmentCount[i,ncol(TreatmentCount)]))))
  
  if(ReturnTreatmentAverages){
    # Plot treatment average graphs
    PlotGraphsAverage <- lapply(1:NTreatments, function(i){
      DataPlot <- DataTreatmentAverages %>% 
        filter(TreatmentN == i)
      return(plotGraph(DataPlot, "Treatment", i, XLab, YLab, XMin = XMin, XMax = XMax, YMin = YMin, YMax = YMax))
    })
    Output <- append(Output, list(PlotGraphsAverage))
    
    if(PDF){
      # Write to pdf
      pdf(paste(OutputDir, gsub(" ", "_", Variable), "_treatments.pdf", sep = ""))
      print(PlotGraphsAverage)
      dev.off()
    }
  }
  
  if(ReturnAllTreatments | ReturnTriples){
    DataTreatmentAverages <- DataTreatmentAverages %>% 
      mutate(TreatmentN = factor(TreatmentN))
    
    ## Combine all graphs
    DataTreatmentAveragesList <- split(DataTreatmentAverages, DataTreatmentAverages$TreatmentN)
    GeomPoints <- lapply(DataTreatmentAveragesList, function(d){
      return(geom_point(data = d, aes(x = Time, y = unname(unlist(d[,ncol(d)-2])), col = Treatment[1]), size = 1))
    })
    GeomLines <- lapply(DataTreatmentAveragesList, function(d){
      return(geom_line(data = d, aes(x = Time, y = unname(unlist(d[,ncol(d)-1])), col = Treatment[1]), size = 1))
    })
    GeomBars <- lapply(DataTreatmentAveragesList, function(d){
      return(geom_errorbar(data = d, aes(x = Time,
                                         ymin = unname(unlist(d[,ncol(d)-1])) - SD,
                                         ymax = unname(unlist(d[,ncol(d)-1])) + SD, col = Treatment[1]), width = 0.3))
    })
    
    # Create a vector of the most distinctive colours
    QualitativeColors = brewer.pal.info[brewer.pal.info$category == 'qual',]
    ColourVector = unlist(mapply(brewer.pal, QualitativeColors$maxcolors, rownames(QualitativeColors)))
    
    Values <- ColourVector[1:19]
    names(Values) <- unname(sapply(DataTreatmentAveragesList, function(d)d$Treatment[1]))
    Labels <- unname(sapply(DataTreatmentAveragesList, function(d)d$Treatment[1]))
    names(Labels) <- Labels
    
    if(ReturnAllTreatments){
      All <- ggplot() +
        GeomPoints +
        GeomLines +
        theme_classic(base_size = 25) +
        scale_color_manual(labels = Labels,
                           values = Values) +
        labs(title = "All treatments",
             y = YLab,
             x = XLab,
             color = "Treatment") + 
        aes(xmin = XMin,
            ymin = YMin,
            xmax = XMax,
            ymax = YMax)
      
      ## Intercrops crops only
      IC <- ggplot() +
        GeomPoints[1:12] +
        GeomLines[1:12] +
        theme_classic(base_size = 25) +
        scale_color_manual(labels = Labels[1:12],
                           values = Values[1:12]) +
        labs(title = "Intercrop treatments",
             y = YLab,
             x = XLab,
             color = "Treatment") + 
        aes(xmin = XMin,
            ymin = YMin,
            xmax = XMax,
            ymax = YMax)
      
      ## Sole crops only
      SC <- ggplot() +
        GeomPoints[13:19] +
        GeomLines[13:19] +
        theme_classic(base_size = 25) +
        scale_color_manual(labels = Labels[13:19],
                           values = Values[13:19]) +
        labs(title = "Sole crop treatments",
             y = YLab,
             x = XLab,
             color = "Treatment") + 
        aes(xmin = XMin,
            ymin = YMin,
            xmax = XMax,
            ymax = YMax)
      
      Output <- append(Output, list(All, IC, SC))
      
      if(PDF){
        # Write to pdf
        pdf(paste(OutputDir, gsub(" ", "_", Variable), "_all.pdf", sep = ""))
        print(All)
        print(IC)
        print(SC)
        dev.off()
      }
    }
    
    if(ReturnTriples){
      ## Combine IC SC specific crops
      Treatments <- read_delim(paste(Dir, "PlotTreatments.txt", sep = ""), delim = "\t")
      
      ### Get all triple combinations
      Triples <- lapply(unique(DataTreatmentAverages$Treatment), function(treatment){
        first <- strsplit(treatment, "_")[[1]][1]
        second <- strsplit(treatment, "_")[[1]][2]
        if(!is.na(first) & !is.na(second)){
          return(c(treatment, first, second))
        }else{
          return(NA)
        }
      })
      Triples <- Triples[!is.na(Triples)]
      names(DataTreatmentAveragesList) <- lapply(DataTreatmentAveragesList, function(l)l$Treatment[1])
      names(GeomPoints) <- names(DataTreatmentAveragesList)
      names(GeomLines) <- names(DataTreatmentAveragesList)
      names(GeomBars) <- names(DataTreatmentAveragesList)
      
      TriplePlots <- lapply(Triples, function(triple){
        plotMultipleGraphs(GeomBars[triple], GeomLines[triple], triple, ColourVector, XLab = XLab, YLab = YLab, XMin = XMin, XMax = XMax, YMin = YMin, YMax = YMax)
      })
      
      Output <- append(Output, list(TriplePlots))      
      
      if(PDF){
        # Write to pdf
        pdf(paste(OutputDir, gsub(" ", "_", Variable), "_triples.pdf", sep = ""))
        print(TriplePlots)
        dev.off()
      }
    }
  }
  return(Output)
}

##
## End
##