# Load necessary libraries for data manipulation and visualization
library(readr)      # For reading CSV files
library(readxl)     # For reading Excel files
library(reshape2)   # For reshaping data
library(ggplot2)    # For creating visualizations
library(cowplot)    # For enhancing ggplot2 plots with a cohesive layout

# Define the location of the data
data_location <- '../../data/figure1'

# Read in the pre-processed results from an Excel file, specifying the sheet to be used
preprocessingresuls <- read_excel(paste0(data_location, "/PAPER_processingPredictionModels_NoLogAdjustedPValue.xlsx"), 
                                  sheet = "PlotDataset_Manual")

# Rename columns in the dataframe for easier reference
colnames(preprocessingresuls) <- c("Metric", "True Positives (TP)", "True Negatives (TN)", "False Positives (FP)", 
                                   "False Negatives (FN)", "ACC", "PRE", "SNS", "SPC", "Negative Predictive value", 
                                   "F1-Score", "Dataset")  

# Reshape the dataframe from wide to long format using melt function
molten <- melt(preprocessingresuls, id.vars = c('Metric', 'Dataset'))

# Rename the columns of the molten dataframe for clarity
colnames(molten) <- c('Pre-processing', 'Dataset', 'Metric', 'Score')

# Filter out specific metrics that are not needed for the analysis
molten <- subset(molten, !(Metric %in% c('F1-Score',
                                         'True Positives (TP)',
                                         'True Negatives (TN)',
                                         'False Positives (FP)',
                                         'False Negatives (FN)',
                                         'Negative Predictive value')))

# Further filter out pre-processing methods that are not required
molten <- subset(molten, !(`Pre-processing` %in% c('Norm')))

# Create a subset of the molten data for the 'Log 2' pre-processing method
pltA <- molten[molten$`Pre-processing` == 'Log 2',]

# Generate a bar plot for the 'Log 2' pre-processing results
plotA <- ggplot(pltA, aes(x=Dataset, y=Score, fill=Dataset)) + 
  geom_col(aes(fill=Dataset), colour='black') + 
  scale_fill_manual("legend", values = c("Test" = "#69AAD8", "Train" = "#0872BC", 'All' = '#F36A6E')) +
  facet_grid(~`Metric`) + 
  theme_cowplot() +  # Use cowplot theme for improved aesthetics
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")  # Customize axis text and legend

# Create a subset of the molten data for the 'Log 10' pre-processing method
pltB <- molten[molten$`Pre-processing` == 'Log 10',]

# Generate a bar plot for the 'Log 10' pre-processing results
plotB <- ggplot(pltB, aes(x=Dataset, y=Score, fill=Dataset)) + 
  geom_col(aes(fill=Dataset), colour='black') + 
  scale_fill_manual("legend", values = c("Test" = "#69AAD8", "Train" = "#0872BC", 'All' = '#F36A6E')) +
  facet_grid(~`Metric`) + 
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")

# Create a subset of the molten data for the 'Ratio' pre-processing method
pltC <- molten[molten$`Pre-processing` == 'Ratio',]

# Generate a bar plot for the 'Ratio' pre-processing results
plotC <- ggplot(pltC, aes(x=Dataset, y=Score, fill=Dataset)) + 
  geom_col(aes(fill=Dataset), colour='black') + 
  scale_fill_manual("legend", values = c("Test" = "#69AAD8", "Train" = "#0872BC", 'All' = '#F36A6E')) +
  facet_grid(~`Metric`) + 
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")

# Combine the three plots into a single plot layout
plot <- plot_grid(plotA, plotB, plotC, labels = c('A', 'B', 'C'), nrow = 1)

# Save the combined plot as both PNG and PDF files with specified dimensions
save_plot('figure1.png', plot, base_height = 4, base_width = 15)
save_plot('figure1.pdf', plot, base_height = 4, base_width = 15)
