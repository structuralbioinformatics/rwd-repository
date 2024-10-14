library(ggplot2)       # For creating graphics and visualizations
library(tidyverse)     # For data manipulation and wrangling
library(readr)         # For reading CSV files
library(cowplot)       # For combining multiple ggplot objects into a single plot
library(reshape2)      # For reshaping data, especially for ggplot
library(stringr)       # For string manipulation
library(ggrepel)       # For better placement of text labels in ggplot2
library(tidyr)         # For data manipulation, particularly reshaping
library(naniar)        # For handling missing data visualization
library(RColorBrewer)  # For color palettes in plots
library(matrixTests)   # For performing matrix-based statistical tests
library(purrr)         # For functional programming (map, reduce, etc.)
library(biomaRt)

data_location <- '../../data/figure4/'

# Load necessary libraries


# Set the date for file naming purposes
date = '211126'

# Function to retrieve the gene name from a UniProt ID
convert_uniprot_to_gcode <- function(uniprot_ids) {
  # Use the Ensembl BioMart service to query the human gene dataset
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Query BioMart for mapping Entrez gene IDs to UniProt IDs
  result <- getBM(
    attributes = c("uniprotswissprot", "hgnc_symbol"),  # Attributes to retrieve
    filters = "uniprotswissprot",                            # Filter based on Entrez gene ID
    values = uniprot_ids,                                  # Values to filter (input Entrez IDs)
    mart = ensembl                                        # BioMart object
  )
  
  # Replace empty UniProt IDs with "UNK" (Unknown)
  hgcn_symbols <- result$hgnc_symbol
  hgcn_symbols[hgcn_symbols == ""] <- "UNK"
  
  # Update the result with modified UniProt IDs
  result$hgnc_symbol <- hgcn_symbols
  
  # Create a named vector mapping Entrez IDs to UniProt IDs
  uniprot_map <- setNames(result$hgnc_symbol, result$uniprotswissprot)
  
  # Return the UniProt mapping
  return(uniprot_map)
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
  rawDataset <- read_csv(paste0('../../data/figure4/', date, '_', dataName, "_raw.csv"),
                         col_types = cols(UNIPROT = col_character(), .default = col_double()))
  proceDataset <- read_csv(paste0('../../data/figure4/', date, '_', dataName, "_aftproc.csv"),
                           col_types = cols(UNIPROT = col_character(), .default = col_double()))
  splineDataset <- read_csv(paste0('../../data/figure4/', date, '_', dataName, "_spline.csv"),
                            col_types = cols(UNIPROT = col_character(), .default = col_double()))
  polyDataset <- read_csv(paste0('../../data/figure4/', date, '_', dataName, "_poly.csv"),
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

doDifferentialExpressionAnalysis <- function(dataset, studydesign) {
  geneNames <- dataset[,1]
  dataset <- dataset[,-1]
  controlPop <- dataset[, studydesign == 1]
  testPop <- dataset[, studydesign == 0]
  
  meanExprControl <- rowMeans(controlPop, na.rm = T)
  meanExprTest <- rowMeans(testPop, na.rm = T)
  
  tTestMeans <- matrixTests::row_t_welch(testPop, controlPop)
  tTestMeans <- cbind(geneNames, tTestMeans)
  return(tTestMeans)
}

make_vulcano_plots <- function(control_dataset, test_dataset, dataset) {
  combined <- merge(control_dataset,test_dataset,by="UNIPROT")
  
  isControl <- c(1:(length(control_dataset)-1))
  isControl[] <- 1
  isTest <- c(1:(length(test_dataset)-1))
  isTest[] <- 0
  controlColumn <- c(isControl, isTest)
  
  ttestResult <- doDifferentialExpressionAnalysis(combined, controlColumn)
  geneNamesUprot <- ttestResult$geneNames
  HUGOMAPPER <- convert_uniprot_to_gcode(geneNamesUprot)
  
  ttestResult$geneNames <- HUGOMAPPER[as.character(ttestResult$geneNames)]
  ttestResult$geneNames[is.na(ttestResult$geneNames)] <- ""
  
  ttestResult$HUGOCODE <- ttestResult$geneNames
  
  ####IN-HOUSE VULCANO PLOT
  FC_Cutoff <- 0.5
  ttestResult$diffexpressed <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
  ttestResult$diffexpressed[ttestResult$mean.diff > FC_Cutoff & ttestResult$pvalue < 0.05] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  ttestResult$diffexpressed[ttestResult$mean.diff < -FC_Cutoff & ttestResult$pvalue < 0.05] <- "DOWN"
  
  ttestResult$delabel <- NA
  ttestResult$delabel[ttestResult$diffexpressed != "NO"] <- ttestResult$HUGOCODE[ttestResult$diffexpressed != "NO"]
  
  write.table(ttestResult, paste0(dataset,'_prostCanvsControl_HouseTTestDifferentialExpression.csv'), row.names = T,
              col.names = T, sep = ',')
  vulcplot_house <- ggplot(ttestResult) +
    geom_point(aes(x=mean.diff, y=-log10(pvalue), colour=diffexpressed)) +
    geom_text_repel(aes(x=mean.diff, y=-log10(pvalue), label=delabel, colour=diffexpressed)) +
    scale_color_manual(values=c("#F22222", "#BFBFBF", "#2947F2")) +
    geom_vline(xintercept=FC_Cutoff, linetype = "longdash") + geom_vline(xintercept=-FC_Cutoff, linetype = "longdash") +
    geom_hline(yintercept=-log10(0.05), linetype = "longdash") +
    theme_cowplot() + theme(legend.position = "none") +
    scale_y_continuous(name=expression("-log"[10]*"(p-value)"), limits = c(0, 15)) +
    scale_x_continuous(name=expression("log(FC)"), limits = c(-3, 3))
  
  return(vulcplot_house)
}

# Number of samples to generate
nsamples <- 10000

# Generate evaluation graphs for 'control' and 'PROSTATIC_NEOPLASMS' datasets
control <- makeEvaluationGraphs('control', nsamples, date)
prostcan <- makeEvaluationGraphs('PROSTATIC_NEOPLASMS', nsamples, date)

#Differential Expression
vulc_spline <- make_vulcano_plots(control$spline,prostcan$spline, 'spline')
vulc_poly <- make_vulcano_plots(control$poly,prostcan$poly, 'polynomial')

#Enrichment plot
panther_enrichment <- read_excel(paste0(data_location, "/panther_enrichment.xlsx"))
panther_enrichment = panther_enrichment %>% arrange(`Category name - L0`, `Percent of gene hit against total # genes - L1`)

empty_bar <- 4
to_add <- data.frame( matrix(NA, empty_bar*nlevels(as.factor(panther_enrichment$`Category name - L0`)), ncol(panther_enrichment)) )
colnames(to_add) <- colnames(panther_enrichment)
to_add$`Category name - L0` <- rep(levels(as.factor(panther_enrichment$`Category name - L0`)), each=empty_bar)
panther_enrichment <- rbind(panther_enrichment, to_add)
panther_enrichment <- panther_enrichment %>% arrange(`Category name - L0`)
panther_enrichment$Index <- seq(1, nrow(panther_enrichment))

# Get the name and the y position of each label
label_data <- panther_enrichment
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$Index-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- panther_enrichment %>% 
  group_by(`Category name - L0`) %>% 
  summarize(start=min(Index), end=max(Index) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

enrichment_plot <- ggplot(panther_enrichment, aes(x=as.factor(Index), y=`Percent of gene hit against total # genes - L1`, fill=`Category name - L0`)) +
  geom_bar(aes(x=as.factor(Index), y=`Percent of gene hit against total # genes - L1`, fill=`Category name - L0`), stat="identity", alpha=0.5) +
  
  geom_segment(data=grid_data, aes(x = end, y = 0.8, xend = start, yend = 0.8), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.6, xend = start, yend = 0.6), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.4, xend = start, yend = 0.4), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.2, xend = start, yend = 0.2), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  
  annotate("text", x = rep(max(panther_enrichment$Index),4), y = c(0.2, 0.4, 0.6, 0.8), label = c("20%", "40%", "60%", "80%") , color="grey", size=1.5 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(Index), y=`Percent of gene hit against total # genes - L1`, fill=`Category name - L0`), stat="identity", color='black') +
  
  ylim(-0.3,1.3) +
  theme_minimal() +
  scale_fill_brewer(palette = "Spectral") +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = 'bottom',
    legend.title = element_blank(),
  ) +
  coord_polar() +
  geom_text(data=label_data, aes(x=as.factor(Index), y=`Percent of gene hit against total # genes - L1`+0.1, label=`Category name - L1`, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=1.5, angle= label_data$angle, inherit.aes = FALSE ) +

  geom_segment(data=base_data, aes(x = start, y = -0.05, xend = end, yend = -0.05), colour = "black", size=0.6 , inherit.aes = FALSE )

final_plot <- plot_grid(plot_grid(vulc_spline, vulc_poly, ncol = 2, labels = c("A", "B")),
                        enrichment_plot, ncol = 1, labels = c("", "C"), rel_heights = c(0.4, 1))
ggsave('figure4.png', final_plot, width = 10, height = 12, dpi = 600, units = c('in'))
ggsave('figure4.pdf', final_plot, width = 10, height = 12, dpi = 600, units = c('in'))
