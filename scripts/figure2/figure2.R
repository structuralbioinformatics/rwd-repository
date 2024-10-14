# Load necessary libraries
library(readr)         # For reading CSV files
library(ggplot2)       # For creating plots
library(cowplot)       # For combining multiple ggplot objects into a single plot
library(reshape2)      # For reshaping data, especially for ggplot
library(stringr)       # For string manipulation
library(limma)         # For differential expression analysis of datasets
library(ggrepel)       # For better placement of text labels in ggplot2
library(tidyr)         # For data manipulation, particularly reshaping
library(naniar)        # For handling missing data visualization
library(RColorBrewer)  # For color palettes in plots
library(matrixTests)   # For performing matrix-based statistical tests
library(purrr)         # For functional programming (map, reduce, etc.)

# Set the date for file naming purposes
date = '211126'

# Function to retrieve the gene name from a UniProt ID
getGeneName <- function(uniprotId) {
  geneName <- tryCatch({
    unipObj <- uniprot(uniprotId)  # Fetch UniProt information
    outName <- unipObj$gene[1]     # Extract the first gene name
    return(outName)
  },
  error = function(cond) {
    message("Here's the original warning message:")  # If there's an error, display message
    message(cond)
    outName <- NA
    return(outName)
  })
  return(geneName)
}

# Function to compute summary statistics of the dataset (max, min, mean)
statDescribe <- function(dataset) {
  numericPart <- as.data.frame(dataset[,-1])  # Remove the first column (non-numeric)
  numericPart[!is.finite(as.matrix(numericPart))] <- NA  # Replace non-finite values with NA
  totalMax <- max(numericPart, na.rm = T)  # Maximum value in the dataset
  totalMin <- min(numericPart, na.rm = T)  # Minimum value in the dataset
  totalMean <- mean(rowMeans(numericPart, na.rm = T), na.rm = T)  # Mean of row means
  
  return(c(totalMax, totalMin, totalMean))  # Return the max, min, and mean
}

# Function to identify outlier samples in the dataset using boxplot statistics
getOutlierSamples <- function(dataset) {
  charList <- apply(dataset[, c(2:length(dataset))], 2, mean, na.rm=T)  # Calculate the mean of each column (sample)
  outliersList <- boxplot(charList, plot=FALSE)$out  # Identify outliers using boxplot statistics
  outliersList <- names(outliersList)  # Extract the names of the outliers
  return(unique(outliersList))  # Return unique outlier sample names
}

# Function to generate evaluation graphs and preprocess the dataset
makeEvaluationGraphs <- function(dataName, maxn, date) {
  
  # Load the datasets: raw, processed, spline, and polynomial
  rawDataset <- read_csv(paste0("../../data/figure2/", date, '_', dataName, "_raw.csv"),
                         col_types = cols(UNIPROT = col_character(), .default = col_double()))
  proceDataset <- read_csv(paste0("../../data/figure2/", date, '_', dataName, "_aftproc.csv"),
                           col_types = cols(UNIPROT = col_character(), .default = col_double()))
  splineDataset <- read_csv(paste0("../../data/figure2/", date, '_', dataName, "_spline.csv"),
                            col_types = cols(UNIPROT = col_character(), .default = col_double()))
  polyDataset <- read_csv(paste0("../../data/figure2/", date, '_', dataName, "_poly.csv"),
                          col_types = cols(UNIPROT = col_character(), .default = col_double()))
  
  # Log2 transformation for raw and processed datasets
  rawDataset[, 2:ncol(rawDataset)] <- log2(rawDataset[, 2:ncol(rawDataset)])
  proceDataset[, 2:ncol(proceDataset)] <- log2(proceDataset[, 2:ncol(proceDataset)])
  
  # Identify and remove outlier samples from all datasets
  outliersTotal <- unique(c(
    getOutlierSamples(rawDataset),
    getOutlierSamples(proceDataset),
    getOutlierSamples(splineDataset),
    getOutlierSamples(polyDataset)
  ))
  
  # Remove outliers from datasets
  rawDataset <- rawDataset[, !(names(rawDataset) %in% outliersTotal)]
  proceDataset <- proceDataset[, !(names(proceDataset) %in% outliersTotal)]
  splineDataset <- splineDataset[, !(names(splineDataset) %in% outliersTotal)]
  polyDataset <- polyDataset[, !(names(polyDataset) %in% outliersTotal)]
  
  # Compute and save basic statistics for each dataset (Max, Min, Mean)
  rawDesc <- statDescribe(rawDataset)
  proceDesc <- statDescribe(proceDataset)
  splineDesc <- statDescribe(splineDataset)
  polyDesc <- statDescribe(polyDataset)
  
  # Compile statistics into a data frame and write to a CSV file
  collapsedStats <- data.frame(
    Raw = rawDesc, 
    Processing = proceDesc, 
    Spline = splineDesc, 
    Polynomial = polyDesc,
    row.names = c('Max', 'Min', 'Mean')
  )
  write.csv(collapsedStats, file = paste0(date, '_', dataName, '_datasetdescription.csv'), row.names = TRUE)
  
  # Adjust max number of samples if dataset has fewer than maxn
  maxn <- min(maxn, ncol(splineDataset))
  
  # Subset datasets to the top maxn samples
  rawDataset <- rawDataset[, 1:maxn]
  proceDataset <- proceDataset[, 1:maxn]
  splineDataset <- splineDataset[, 1:maxn]
  polyDataset <- polyDataset[, 1:maxn]
  
  # Create melted datasets for ribbon plot (median, min, max values per sample)
  meltDataset <- function(dataset) {
    data.frame(
      datapoint = 1:(ncol(dataset) - 1),
      median = apply(dataset[, -1], 2, median, na.rm = TRUE),
      min = apply(dataset[, -1], 2, min, na.rm = TRUE),
      max = apply(dataset[, -1], 2, max, na.rm = TRUE)
    )
  }
  
  melt_rawDataset <- meltDataset(rawDataset)
  melt_proceDataset <- meltDataset(proceDataset)
  melt_splineDataset <- meltDataset(splineDataset)
  melt_polyDataset <- meltDataset(polyDataset)
  
  # Function to create ribbon plots for each dataset
  createRibbonPlot <- function(melted_data) {
    ggplot(melted_data, aes(x = 1:nrow(melted_data), y = median)) +
      geom_line() +
      geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.5) +
      theme_cowplot() +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) +
      ylim(-20, 20) +  # Adjust y-axis limits for consistency
      xlab("Sample") + 
      ylab(expression(log2(intensity)))
  }
  
  # Create individual plots for each dataset
  raw_plot <- createRibbonPlot(melt_rawDataset)
  proce_plot <- createRibbonPlot(melt_proceDataset)
  spline_plot <- createRibbonPlot(melt_splineDataset)
  poly_plot <- createRibbonPlot(melt_polyDataset)
  
  # Combine all plots into one grid
  final_plot <- plot_grid(raw_plot, proce_plot, spline_plot, poly_plot, ncol = 4)
  
  # Return the combined plot and the datasets used
  return(list(plot = final_plot, raw = rawDataset, proc = proceDataset, spline = splineDataset, poly = polyDataset))
}

# Function to calculate mean values for each dataset (raw, processed, spline, polynomial)
getMeans <- function(dataset) {
  # Calculate column means for each dataset, excluding the first column (UNIPROT)
  meanValsRaw <- cbind(dataset$raw$UNIPROT, colMeans(dataset$raw[ , 2:ncol(dataset$raw)], na.rm = TRUE))
  meanValsProc <- cbind(dataset$proc$UNIPROT, colMeans(dataset$proc[ , 2:ncol(dataset$proc)], na.rm = TRUE))
  meanValsSpline <- cbind(dataset$spline$UNIPROT, colMeans(dataset$spline[ , 2:ncol(dataset$spline)], na.rm = TRUE))
  meanValsPoly <- cbind(dataset$poly$UNIPROT, colMeans(dataset$poly[ , 2:ncol(dataset$poly)], na.rm = TRUE))
  
  # Return a list with the mean values for each dataset
  return(list(raw = meanValsRaw, proc = meanValsProc, spline = meanValsSpline, poly = meanValsPoly))
}

# Function to create MA plots (M vs A plots) for test and control datasets
makeMAPlots <- function(test, control) {
  # Get mean values for test and control datasets
  meansControl <- getMeans(control)
  meansTest <- getMeans(test)
  
  # Merge control and test datasets by UNIPROT ID (assumed as the first column)
  mergeDatasets <- function(control, test) {
    merge(control, test, by = "V1")
  }
  
  rawMeans <- mergeDatasets(meansControl$raw, meansTest$raw)
  procMeans <- mergeDatasets(meansControl$proc, meansTest$proc)
  splineMeans <- mergeDatasets(meansControl$spline, meansTest$spline)
  polyMeans <- mergeDatasets(meansControl$poly, meansTest$poly)
  
  # Calculate M (difference) and A (average) values for each dataset
  calcMA <- function(means) {
    as.data.frame(cbind(
      M = as.numeric(means$V2.x) - as.numeric(means$V2.y),  # M: Difference between test and control
      A = 1/2 * (as.numeric(means$V2.x) + as.numeric(means$V2.y))  # A: Average of test and control
    ))
  }
  
  rawMA <- calcMA(rawMeans)
  procMA <- calcMA(procMeans)
  splineMA <- calcMA(splineMeans)
  polyMA <- calcMA(polyMeans)
  
  # Function to create a single MA plot
  createMAPlot <- function(maData) {
    ggplot(maData, aes(x = A, y = M)) +
      stat_density_2d(geom = "raster", aes(fill = after_stat(density)), contour = FALSE) +
      scale_fill_gradient(low = "white", high = "deepskyblue4") +
      geom_smooth(color = 'red') + 
      geom_hline(yintercept = 0, color = 'blue') +
      theme_cowplot() +
      theme(legend.position = "none") + 
      ylim(-10, 10) +
      xlab("A") + 
      ylab("M")
  }
  
  # Create MA plots for raw, processed, spline, and polynomial datasets
  rawMAPlot <- createMAPlot(rawMA)
  procMAPlot <- createMAPlot(procMA)
  splineMAPlot <- createMAPlot(splineMA)
  polyMAPlot <- createMAPlot(polyMA)
  
  # Combine all MA plots into a grid
  final_plt <- plot_grid(rawMAPlot, procMAPlot, splineMAPlot, polyMAPlot, ncol = 4)
  
  return(final_plt)
}

# Function to create density plots for each dataset (raw, processed, spline, polynomial)
densityPlotsTest <- function(dataset) {
  # Define color palette for the density plots
  getDensPalette <- colorRampPalette(brewer.pal(9, "Spectral"))
  
  # Function to create a density plot for a given dataset
  createDensityPlot <- function(dataset, paletteLength) {
    datasetMelted <- melt(dataset, id.vars = 'UNIPROT')  # Reshape data for ggplot
    ggplot(datasetMelted, aes(x = value, color = variable)) +
      geom_density() +
      scale_color_manual(values = getDensPalette(paletteLength)) +
      theme_cowplot() +
      theme(legend.position = "none") +
      xlab("Expression") +
      ylab("Density")
  }
  
  # Create density plots for raw, processed, spline, and polynomial datasets
  densPlotRaw <- createDensityPlot(dataset$raw, length(dataset$spline) - 1)
  densPlotProc <- createDensityPlot(dataset$proc, length(dataset$spline) - 1)
  densPlotSpline <- createDensityPlot(dataset$spline, length(dataset$spline) - 1)
  densPlotPoly <- createDensityPlot(dataset$poly, length(dataset$poly) - 1)
  
  # Combine all density plots into a grid
  final_plt <- plot_grid(densPlotRaw, densPlotProc, densPlotSpline, densPlotPoly, ncol = 4)
  
  return(final_plt)
}

# Number of samples to generate
nsamples <- 10000

# Generate evaluation graphs for 'control' and 'PROSTATIC_NEOPLASMS' datasets
control <- makeEvaluationGraphs('control', nsamples, date)
prostcan <- makeEvaluationGraphs('PROSTATIC_NEOPLASMS', nsamples, date)

# Create MA Plot by comparing prostcan and control datasets
maPlot <- makeMAPlots(prostcan, control)

# Create density plots for control and prostcan datasets
densityPlotsControl <- densityPlotsTest(control)
densityPlotsProstCan <- densityPlotsTest(prostcan)

# Combine all plots (control, density plots, and MA plot) into a single plot grid
final_plt <- plot_grid(
  control$plot,           # Control plot
  densityPlotsControl,    # Density plot for control dataset
  prostcan$plot,          # Prostcan plot
  densityPlotsProstCan,   # Density plot for prostcan dataset
  maPlot,                 # MA plot
  ncol = 1                # Arrange plots in one column
)

# Save the final combined plot to a PDF file
ggsave2('figure3.pdf', final_plt, width = 13, height = 16, units = 'in', dpi = 600)

# Save the final combined plot to a PNG file
ggsave2('figure3.png', final_plt, width = 13, height = 16, units = 'in', dpi = 600)
