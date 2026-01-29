#' Course   BIO3010
#' Year     2026
#' 
#' Study    Metabolic changes in Alzheimer patient-derived induced neurons versus non-demented controls
#' ID       ST002213
#' https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST002213&StudyType=MS&ResultType=1
#' 
#' 

# (0) Packages -----------------------------------------------------------------

#install packages 
install.packages(c("readr","ggplot2","pheatmap","reshape2"))

#import packages
library(readr)
library(ggplot2)
library(pheatmap)
library(reshape2)

# (1) Read data negative ion mode data: AN003618_N.txt
data_negative_mode <- read_delim("AN003618_N.txt", delim = "\t", show_col_types = FALSE)

#visualize data in the console - see the first 6 rows of the table
head(data_negative_mode) 


# (2) Volcano plots (AD vs Control) --------------------------------------------

#elongate data
negative_long_data <- melt(data_negative_mode,
                           id.vars = "metabolite_name",
                           variable.name = "sample",
                           value.name = "value")

#visualize elongated data - note that the data is now redesigned
head(negative_long_data) 


#add a new column showing whether it is Control or AD sample
negative_long_data$group <- ifelse(grepl("^Control", negative_long_data$sample), "Control", "AD")

#visualize the new data - node the added column
head(negative_long_data) 


#Perform statistical analysis on the data

##1 create a list of metabolites
metabolite_list <- unique(negative_long_data$metabolite)

##2 import the compute_volcano_stats function from volcano_stats.R file in the same directory
source("volcano_stats.R")

# a custom-made function is added to your environment

##3 run the compute_volcano_stats function which performs statistical analysis
analyzed_neg_mode <- compute_volcano_stats(negative_long_data, metabolite_list)

#visualize new, analyzed dataframe with new columns
head(analyzed_neg_mode) 
#' in this new data set,  "mean_control" is the mean of control group and "mean_ad" is
#' the mean of Alzheimer Disease group for the given metabolite. "log2FC" is the fold-change 
#' of metabolites (log2), e.g. when log2FC = 2, the metabolite is found 4 times more in
#' AD samples. p_value shows the p-value after unpaired t-test.
       



# Now it is time to plot our data.

# Plot1 
ggplot(data = analyzed_neg_mode, aes(x = log2FC, y = -log10(p_value))) + 
   geom_point() 

# Plot2 — add axis names
ggplot(data = analyzed_neg_mode, aes(x = log2FC, y = -log10(p_value))) + 
   geom_point() +
   labs(title = "Volcano plot: AD vs Control (negative mode)",
           x = "log2 Fold Change (AD / Control)",
           y = "-log10 p-value") 

# Plot3 - label data points
ggplot(analyzed_neg_mode, aes(x = log2FC, y = -log10(p_value))) + 
   geom_point() +
   labs(title = "Volcano plot: AD vs Control (negative mode)",
        x = "log2 Fold Change (AD / Control)",
        y = "-log10 p-value") + 
   geom_text(aes(label = metabolite))

# Plot3 - label wanted datapoint   
ggplot(analyzed_neg_mode, aes(x = log2FC, y = -log10(p_value))) + 
   geom_point() +
   labs(title = "Volcano plot: AD vs Control (negative mode)",
        subtitle = "Today is monday",
        x = "log2 Fold Change (AD / Control)",
        y = expression(-log[10]("p-value"))
   ) + 
   geom_text(data=subset(analyzed_neg_mode, log2FC>1), aes(label = metabolite))

# Take a screenshot



#' EXERCISE 1
#' Label three most down-regulated metabolites in AD
#' Add your name and student number to subtitle
#' take a screenshot and upload to Canvas

#' EXERCISE 2
#' Analyze the positive ion-mode data and label most-upregulated metabolites.
#' Add your name and student number to subtitle
#' take a screenshot and upload to Canvas

# END END HERE
