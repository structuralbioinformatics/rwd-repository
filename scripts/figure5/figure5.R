# Load necessary libraries
library(readr)    # For data reading
library(dplyr)    # For data manipulation
library(ggplot2)  # For creating plots
library(cowplot)  # For combining multiple plots
library(readxl)   # For reading Excel files
library(xlsx)     # For reading and writing Excel files
library(rstatix)  # For statistical tests
library(ggbreak)  # For adjusting scales in plots

# Define the location of the data
data_location <- "../../data/figure5"

# Read age data
age <- read_delim(paste0(data_location, "/datasets_metadata/20220317_age.txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE)
age <- age[age$regression >= 0,]  # Filter data for non-negative ages

# Read BMI data
bmi <- read_delim(paste0(data_location, "/datasets_metadata/20220317_bmi.txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE)
# Classify BMI data into categories
bmi$classification <- 'NA'
bmi$classification[bmi$regression < 18.5] = 'underweight'  # Underweight
bmi$classification[bmi$regression >= 18.5] = 'normal'  # Normal weight
bmi$classification[bmi$regression >= 25] = 'pre-obesity'  # Pre-obesity
bmi$classification[bmi$regression >= 30] = 'obesity class I'  # Obesity Class I
bmi$classification[bmi$regression >= 35] = 'obesity class II'  # Obesity Class II
bmi$classification[bmi$regression >= 40] = 'obesity class III'  # Obesity Class III

# Read height data
height <- read_delim(paste0(data_location, "/datasets_metadata/20220317_height.txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# Read sex data
sex <- read_delim(paste0(data_location, "/datasets_metadata/20220317_sex.txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE)
# Group by sex and calculate patient counts
sex_plot <- sex %>% group_by(regression) %>% 
  summarise(count_patients = n()) %>%
  mutate(percent = count_patients/sum(count_patients))  # Calculate percentage

# Read ethnicity data
ethnicity <- read_delim(paste0(data_location, "/datasets_metadata/20220317_ethnicity.txt"), delim = "\t", escape_double = FALSE, trim_ws = TRUE)
# Group by ethnicity and calculate patient counts
ethnicity_plot <- ethnicity %>% group_by(regression) %>% 
  summarise(count_patients = n()) %>%
  mutate(percent = count_patients/sum(count_patients))  # Calculate percentage

# Create histogram for age distribution
age_cont_hist <- ggplot(age, aes(x=regression, fill=classification)) +
  geom_histogram(bins=40, color = "#000000") +  # Create histogram with 40 bins and black outline
  geom_vline(xintercept = median(age$regression), col = "#000000") +  # Vertical line for median age
  annotate("text", x = median(age$regression)-18, y = 7500, label = paste("Median =", median(age$regression)), angle=0, col = "#000000") +  # Annotate median value
  geom_vline(xintercept = min(age$regression), col = "#000000") +  # Vertical line for minimum age
  annotate("text", x = min(age$regression)+10, y = 5000, label = paste("Min =", min(age$regression)), angle=0, col = "#000000") +  # Annotate minimum value
  geom_vline(xintercept = max(age$regression), col = "#000000") +  # Vertical line for maximum age
  annotate("text", x = max(age$regression)-15, y = 5000, label = paste("Max =", max(age$regression)), angle=0, col = "#000000") +  # Annotate maximum value
  theme_cowplot() +  # Apply cowplot theme
  xlab('Age') +  # Label x-axis
  theme(axis.title.y = element_blank(), legend.position = 'none',) +  # Remove y-axis title and legend
  scale_fill_brewer(palette="YlOrBr")  # Set fill color palette for age classification

# Create histogram for BMI distribution
bmi_cont_hist <- ggplot(bmi, aes(x=regression, fill=classification)) +
  geom_histogram(bins=40, color = "#000000") +  # Create histogram with 40 bins and black outline
  geom_vline(xintercept = median(bmi$regression), col = "#000000") +  # Vertical line for median BMI
  annotate("text", x = median(bmi$regression)+8, y = 1550, label = paste("Median =", round(median(bmi$regression))), angle=0, col = "#000000") +  # Annotate median value
  geom_vline(xintercept = min(bmi$regression), col = "#000000") +  # Vertical line for minimum BMI
  annotate("text", x = min(bmi$regression)+5, y = 900, label = paste("Min =", round(min(bmi$regression))), angle=0, col = "#000000") +  # Annotate minimum value
  geom_vline(xintercept = max(bmi$regression), col = "#000000") +  # Vertical line for maximum BMI
  annotate("text", x = max(bmi$regression)-6, y = 900, label = paste("Max =", round(max(bmi$regression))), angle=0, col = "#000000") +  # Annotate maximum value
  theme_cowplot() +  # Apply cowplot theme
  xlab('BMI') +  # Label x-axis
  theme(axis.title.y = element_blank(), legend.position = 'none',) +  # Remove y-axis title and legend
  scale_fill_brewer(palette="Reds")  # Set fill color palette for BMI classification

# Create histogram for height distribution
height_cont_hist <- ggplot(height, aes(x=regression)) +
  geom_histogram(bins=40, fill="#B781EA", color = "#000000") +  # Create histogram with 40 bins, purple fill and black outline
  geom_vline(xintercept = median(height$regression), col = "#000000") +  # Vertical line for median height
  annotate("text", x = median(height$regression)+15, y = 850, label = paste("Median =", round(median(height$regression))), angle=0, col = "#000000") +  # Annotate median value
  geom_vline(xintercept = min(height$regression), col = "#000000") +  # Vertical line for minimum height
  annotate("text", x = min(height$regression)+13, y = 400, label = paste("Min =", round(min(height$regression))), angle=0, col = "#000000") +  # Annotate minimum value
  geom_vline(xintercept = max(height$regression), col = "#000000") +  # Vertical line for maximum height
  annotate("text", x = max(height$regression)-13, y = 400, label = paste("Max =", round(max(height$regression))), angle=0, col = "#000000") +  # Annotate maximum value
  theme_cowplot() +  # Apply cowplot theme
  xlab('Height') +  # Label x-axis
  theme(axis.title.y = element_blank(), legend.position = 'none',)  # Remove y-axis title and legend

# Create a bar plot for sex distribution
sex_cont_hist <- ggplot(sex_plot) +
  geom_col(aes(x=regression,y=count_patients, fill=regression), color='#000000', position = "dodge2", show.legend = TRUE, alpha = .9) +  # Bar plot with dodging
  scale_y_continuous(limits = c(0, 110000), expand = c(0, 0), breaks = c(0:6) * 25000) +  # Set y-axis limits and breaks
  theme_cowplot() +  # Apply cowplot theme
  xlab('Biological Sex') +  # Label x-axis
  theme(legend.position = 'none', axis.title.y = element_blank(),)  # Remove legend and y-axis title

# Create a polar bar plot for ethnicity distribution
eth_cont_hist <- ggplot(ethnicity_plot) +
  geom_hline(aes(yintercept = y), data.frame(y = c(0:4) * 5000), color = "lightgrey") +  # Add light grey horizontal lines for y-axis reference
  geom_col(aes(x=regression,y=count_patients, fill=regression), position = "dodge2", show.legend = TRUE, alpha = .9, color='#000000') +  # Bar plot with dodging
  coord_polar() +  # Convert to polar coordinates
  scale_y_continuous(limits = c(-5000, 21000), expand = c(0, 0), breaks = c(0:5) * 5000) +  # Set y-axis limits and breaks
  theme_minimal_grid() +  # Apply minimal grid theme
  xlab('Ethnicity') +  # Label x-axis
  theme(legend.position = 'none', axis.title.y = element_blank()) +  # Remove legend and y-axis title
  scale_fill_brewer(palette="Blues")  # Set fill color palette for ethnicity

# Read medical conditions data for parents
medicalconditions_parents <- read_delim(paste0(data_location, "/20220405_medicalconditions_parents.csv"), 
                                        delim = ";", escape_double = FALSE, trim_ws = TRUE)
medicalconditions_parents <- medicalconditions_parents[!is.na(medicalconditions_parents$`Family Name`),]  # Remove rows with NA in Family Name

# Create a bar plot for medical conditions by family name
disease_plot_hist <- ggplot(medicalconditions_parents) +  # Initialize plot
  geom_bar(aes(x=`Family Name`, fill=`Family Name`), alpha=.9, color='#000000') +  # Bar plot for family names
  theme_cowplot() +  # Apply cowplot theme
  scale_y_cut(breaks = c(30000, 150000), which = c(1, 2, 3), scales = c(.2, 0.1, 3)) +  # Cut y-axis at specific breaks
  theme(legend.position = 'none', axis.title.y = element_blank(), axis.text.x = element_text(angle = 40, vjust = 1, hjust=1, size = 9))  # Customize axis labels and theme
disease_plot_hist <- print(disease_plot_hist)  # Print the disease plot histogram

# Create a specific plot for cancer-related diseases
disease_plot_hist_cancer <- ggplot(medicalconditions_parents[medicalconditions_parents$`Family Name` == 'Neoplasms',]) + 
  scale_fill_brewer(palette="Spectral") +  # Set fill color palette for cancer
  geom_bar(aes(x=Names, fill=Names), alpha=.9, color='#000000') +  # Bar plot for names related to Neoplasms
  theme_cowplot() +  # Apply cowplot theme
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1, size = 9))  # Customize axis labels for x-axis

# Read global accuracy report for tag classification
tagclassification_globalaccuracyreport <- read_excel(paste0(data_location, "/20221004_tagclassification_globalaccuracyreport.xlsx"),
                                                     sheet = "Sheet2")
colnames(tagclassification_globalaccuracyreport) <- c("Family Name", "BED Condition", "Number of Models", "B. Accuracy", "CI" )  # Set column names

# Create a boxplot for classification accuracy by family name
conditionplot <- ggplot(data=tagclassification_globalaccuracyreport, aes(x=`Family Name`, y=`B. Accuracy`)) +
  geom_boxplot(aes(fill=`Family Name`), alpha=0.4) +  # Boxplot with transparency
  geom_jitter(color="black", size=0.4, alpha=0.9) +  # Jittered points for better visibility
  theme_cowplot() +  # Apply cowplot theme
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1, size = 9),  # Customize x-axis labels
        legend.position = 'none')  # Remove legend

# Combine all individual plots into a panel layout
pannel_plots_1 <- plot_grid(age_cont_hist, bmi_cont_hist, height_cont_hist,  # Combine age, BMI, and height plots
                            labels=c('A', 'B', 'C'), nrow=1)  # Set labels and row number
pannel_plots_2 <- plot_grid(eth_cont_hist, sex_cont_hist,  # Combine ethnicity and sex plots
                            labels=c('D', 'E'), nrow=1)  # Set labels and row number
pannel_plots_final <- plot_grid(pannel_plots_1, pannel_plots_2, disease_plot_hist, conditionplot,  # Combine previous plots with disease plot and condition plot
                                labels=c('', '', 'F', 'G'), ncol=1, rel_heights = c(1, 1, 2, 2))  # Set final layout and relative heights

# Save the final plot to file
save_plot('figure5.png', pannel_plots_final, base_height = 16, base_width = 11)
save_plot('figure5.pdf', pannel_plots_final, base_height = 16, base_width = 11)
