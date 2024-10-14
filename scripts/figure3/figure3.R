# Load necessary libraries
library(readr)          # For reading CSV files
library(dplyr)          # For data manipulation
library(tidyr)          # For data tidying
library(ggplot2)        # For plotting
library(ggVennDiagram)  # For generating Venn diagrams
library(cowplot)        # For enhanced plotting themes and layouts
library(reshape2)       # For reshaping data
library(biomaRt)        # For biomart queries

# Function to convert Entrez gene IDs to UniProt IDs using Ensembl BioMart
convert_entrez_to_uniprot <- function(entrez_ids) {
  # Use the Ensembl BioMart service to query the human gene dataset
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Query BioMart for mapping Entrez gene IDs to UniProt IDs
  result <- getBM(
    attributes = c("entrezgene_id", "uniprotswissprot"),  # Attributes to retrieve
    filters = "entrezgene_id",                            # Filter based on Entrez gene ID
    values = entrez_ids,                                  # Values to filter (input Entrez IDs)
    mart = ensembl                                        # BioMart object
  )
  
  # Replace empty UniProt IDs with "UNK" (Unknown)
  uniprots <- result$uniprotswissprot
  uniprots[uniprots == ""] <- "UNK"
  
  # Update the result with modified UniProt IDs
  result$uniprotswissprot <- uniprots
  
  # Create a named vector mapping Entrez IDs to UniProt IDs
  uniprot_map <- setNames(result$uniprotswissprot, result$entrezgene_id)
  
  # Return the UniProt mapping
  return(uniprot_map)
}

# Function to generate a boxplot for a dataset, coloring based on dataset type
makeBoxplot <- function(dataset, dataset_name) {
  # Reshape the dataset for plotting (convert wide format to long format)
  plottingDataset <- reshape2::melt(dataset, id.vars = c('UNIPROT'))
  colnames(plottingDataset) <- c('Uniprot', 'Sample', 'Intensity')
  
  # Create initial boxplot
  plt <- ggplot(plottingDataset) +
    geom_boxplot(aes(x = Sample, y = Intensity)) +
    theme_cowplot() +
    theme(axis.text.x = element_blank())  # Remove x-axis text
  
  # Get data used in the plot for additional customization
  dat <- ggplot_build(plt)$data[[1]]
  
  # Set boxplot segment color based on dataset type
  if (dataset_name == "raw") {
    colour <- "grey"
  } else if (dataset_name == "macaroon") {
    colour <- "#06d6a0"
  } else if (dataset_name == "yugene") {
    colour <- "#118ab2"
  }
  
  # Add colored segment for median values in boxplot
  plt <- plt + geom_segment(data = dat, aes(x = xmin, xend = xmax,
                                            y = middle, yend = middle),
                            colour = colour, size = 2)
  
  # Return the final plot
  return(plt)
}

# Function to perform differential expression analysis using Limma
differential_expression_analysis_limma <- function(df, group_info) {
  # Ensure the group_info is treated as a factor
  group_info <- factor(group_info)
  
  # Create a design matrix for the linear model
  design <- model.matrix(~ 0 + group_info)
  colnames(design) <- levels(group_info)
  
  # Fit the linear model
  fit <- lmFit(df, design)
  
  # Create contrast matrix to compare the groups
  contrast_matrix <- makeContrasts(
    contrasts = paste(levels(group_info)[1], "-", levels(group_info)[2], sep = ""),
    levels = design
  )
  
  # Apply contrasts and compute statistics using eBayes
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  # Return the differential expression results
  results <- topTable(fit2, sort.by = "logFC", adjust.method = "BH", number = Inf)
  
  return(results)
}

# Set the data location path
data_location <- ("../../data/figure3")

# Read datasets from CSV files
MA_MOD <- read_csv(paste0(data_location, "/GSE102124/20240619_pcaexpression.csv"))
MA_MACAROON <- read_csv(paste0(data_location, "/GSE102124/20240619_pcaexpression_genobablewadjust.csv"))
MA_YUGENE <- read_csv(paste0(data_location, "/GSE102124/20240619_pcaexpression_yugene.csv"))

# Convert Entrez IDs to UniProt for the MOD dataset
MA_MOD_MAPPER <- convert_entrez_to_uniprot(MA_MOD$UNIPROT)
MA_MOD$UNIPROT <- MA_MOD_MAPPER[as.character(MA_MOD$UNIPROT)]
MA_MOD$UNIPROT[is.na(MA_MOD$UNIPROT)] <- "UNK"

# Convert Entrez IDs to UniProt for the YuGene dataset
MA_YUGENE_MAPPER <- convert_entrez_to_uniprot(MA_YUGENE$UNIPROT)
MA_YUGENE$UNIPROT <- MA_YUGENE_MAPPER[as.character(MA_YUGENE$UNIPROT)]
MA_YUGENE$UNIPROT[is.na(MA_YUGENE$UNIPROT)] <- "UNK"

# Generate a density plot for the MACAROON dataset (excluding the first column)
density_dataframe <- stack(MA_MACAROON[,-1])
ggplot(density_dataframe) + geom_density(aes(x = values))

# Create boxplots for the datasets
MA_MOD_boxplt <- makeBoxplot(MA_MOD, "raw")
MA_MACAROON_boxplt <- makeBoxplot(MA_MACAROON, "macaroon")
MA_YUGENE_boxplt <- makeBoxplot(MA_YUGENE, "yugene")

# Combine the boxplots into one figure panel
pannel1 <- plot_grid(MA_MOD_boxplt, MA_MACAROON_boxplt, MA_YUGENE_boxplt, labels = c('A', 'B', 'C'), ncol = 1)

# Read metadata file
METADATA <- read.csv(paste0(data_location, "/METADATA_GSE102124.csv"), header = 0)
METADATA$V2[METADATA$V2 == "prostate, tumor, treated"] <- "treated"
METADATA$V2[METADATA$V2 == "prostate, tumor, untreated"] <- "untreated"
treated_samples <- METADATA$V2 == 'treated'

# Perform differential expression analysis using Limma
MA_MOD_LIMMA <- differential_expression_analysis_limma(MA_MOD, METADATA$V2)
MA_MACAROON_LIMMA <- differential_expression_analysis_limma(MA_MACAROON, METADATA$V2)
MA_YUGENE_LIMMA <- differential_expression_analysis_limma(MA_YUGENE, METADATA$V2)

# Merge the results from the different datasets based on UniProt IDs
joint_stats <- merge(MA_MOD_LIMMA, MA_MACAROON_LIMMA, suffixes = c("_RAW", "_MACAROON"), by.x = 'UNIPROT', by.y = 'UNIPROT')
joint_stats <- merge(joint_stats, MA_YUGENE_LIMMA, suffixes = c("", "_YUGENE"), by.x = 'UNIPROT', by.y = 'UNIPROT')

# Reshape the data for visualization and filter relevant variables
stats_density <- melt(joint_stats, na.rm = FALSE, id = 'UNIPROT')
stats_density <- stats_density[stats_density$variable %in% c('logFC_RAW', 'logFC_MACAROON', 'logFC'),]

# Create a violin plot for log fold change values across datasets
hist_logf <- ggplot(stats_density, aes(x = variable, y = value, fill = variable)) +
  geom_violin() +
  labs(y = 'log FC', x = '') +
  scale_fill_manual(values = c("grey", "#06d6a0", "#118ab2")) +
  theme_cowplot() +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("Raw", "MACAROON", "YuGene"))

# Create scatter plots to compare log fold changes between datasets
scatters <- ggplot(joint_stats) +
  geom_point(aes(x = logFC_RAW, y = logFC_MACAROON), color = '#06d6a0', alpha = 0.5) +
  geom_smooth(aes(x = logFC_RAW, y = logFC_MACAROON), method = lm, color = "#037758", fill = "grey", se = TRUE) +
  geom_point(aes(x = logFC_RAW, y = logFC), color = '#118ab2', alpha = 0.5) +
  geom_smooth(aes(x = logFC_RAW, y = logFC), method = lm, color = "#0a5770", fill = "grey", se = TRUE) +
  labs(y = 'Normalization log FC', x = 'Raw log FC') +
  theme_cowplot() +
  theme(legend.position = "none")

# Filter significant genes and prepare data for Venn diagram
joint_stats_2 <- joint_stats[joint_stats$UNIPROT != "UNK",]
joint_stats_2 <- joint_stats[(joint_stats_2$adj.P.Val_RAW < 0.01 | joint_stats_2$adj.P.Val_MACAROON < 0.01 | joint_stats_2$adj.P.Val < 0.01),]

is_selected_gene_raw <- (rank(joint_stats_2$logFC_RAW) <= 50) | (rank(-joint_stats_2$logFC_RAW) <= 50)
is_selected_gene_macarroon <- (rank(joint_stats_2$logFC_MACAROON) <= 50) | (rank(-joint_stats_2$logFC_MACAROON) <= 50)
is_selected_gene_yugene <- (rank(joint_stats_2$logFC) <= 50) | (rank(-joint_stats_2$logFC) <= 50)

# Create Euler diagram for selected genes
eulerData <- list(A = joint_stats_2$UNIPROT[is_selected_gene_raw],
                  B = joint_stats_2$UNIPROT[is_selected_gene_macarroon],
                  C = joint_stats_2$UNIPROT[is_selected_gene_yugene])
eulerDiagram <- ggVennDiagram(eulerData, category.names = c("Raw", "MACAROON", "YuGene"), label = 'count') +
  theme(legend.position = 'bottom') +
  scale_fill_gradient(low = '#D9DCE1', high = '#495057')

# Create a second panel with scatter plots and Euler diagram
panel2 <- plot_grid(scatters, eulerDiagram, labels = c('E', 'F'), nrow = 1)

# Create a third panel with the violin plot and the second panel
panel3 <- plot_grid(hist_logf, panel2, labels = c('D', ''), ncol = 1)

# Combine all panels into a grid
pannelgrid <- plot_grid(pannel1, panel3, labels = NA, ncol = 2)

# Save the final plot to file
save_plot('figure3.png', pannelgrid, base_height = 8, base_width = 14)
save_plot('figure3.png', pannelgrid, base_height = 8, base_width = 14)
