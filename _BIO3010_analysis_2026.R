#' Course   BIO3010
#' Year     2026
#' 
#' Study    Metabolic changes in Alzheimer patient-derived induced neurons versus non-demented controls
#' ID       ST002213
#' https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST002213&StudyType=MS&ResultType=1
#' 
#' 

# (0) Packages -----------------------------------------------------------------

install.packages(c("readr","ggplot2","pheatmap","reshape2"))

library(readr)
library(ggplot2)
library(pheatmap)
library(reshape2)

# (1) Read data (metabolites x samples) ----------------------------------------
data_negative_mode <- read_delim("AN003618_N.txt", delim = "\t", show_col_types = FALSE)

#visualize data
head(data_negative_mode) 


# (2) Volcano plots (AD vs Control) --------------------------------------------

#elongate data
negative_long_data <- melt(data_negative_mode,
                           id.vars = "metabolite_name",
                           variable.name = "sample",
                           value.name = "value")

#visualize elongated data
head(negative_long_data) 


#group data accordingly
negative_long_data$group <- ifelse(grepl("^Control", negative_long_data$sample), "Control", "AD")
head(negative_long_data) #visualize the new added showing sample column


#Perform statistical analysis on the data

##1 create a list of metabolites
metabolite_list <- unique(negative_long_data$metabolite)

##2 import the compute_volcano_stats function from volcano_stats.R file 
source("volcano_stats.R")

##3 compute_volcano_stats performs statistical analysis
analyzed_neg_mode <- compute_volcano_stats(negative_long_data, metabolite_list)
head(analyzed_neg_mode) #visualize new, analyzed dataframe

#Start plotting

# Plot1 
ggplot(analyzed_neg_mode, aes(x = log2FC, y = -log10(p_value))) + 
   geom_point() 

# Plot2 — add axis names
ggplot(analyzed_neg_mode, aes(x = log2FC, y = -log10(p_value))) + 
   geom_point() +
   labs(title = "Volcano plot: AD vs Control (negative mode)",
           x = "log2 Fold Change (AD / Control)",
           y = expression(-log[10]("BH-adjusted p-value"))) 

# Plot3 - label data points
ggplot(analyzed_neg_mode, aes(x = log2FC, y = -log10(p_value))) + 
   geom_point() +
   labs(title = "Volcano plot: AD vs Control (negative mode)",
        x = "log2 Fold Change (AD / Control)",
        y = expression(-log[10]("p-value"))) + 
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

#' EXERCISE
#' Label three most downregulated metabolites in AD
#' Add your name and student number to subtitle
#' take a screenshot and upload to Canvas
# 

# POSSIBLE END HERE


# (3) Heatmaps + correlation heatmaps ------------------------------------------

#Convert to matrices with metabolite names as rownames
neg_matrix <- as.data.frame(data_negative_mode)
rownames(neg_matrix) <- neg_matrix$metabolite_name
neg_matrix$metabolite_name <- NULL
neg_matrix <- as.matrix(neg_matrix)

##visualize the matrix
head(data_negative_mode)
head(neg_matrix)

# Group labels from sample names
neg_samples <- colnames(neg_matrix)
neg_group <- ifelse(grepl("^Control", neg_samples), "Control", "AD")
neg_group <- factor(neg_group, levels = c("Control", "AD"))
neg_annotation <- data.frame(group = neg_group)
rownames(neg_annotation) <- neg_samples
anno_colors <- list(group = c(Control = "black", AD = "red"))

# For metabolite heatmap we want metabolites as rows, samples as columns
neg_log <- log10(neg_matrix + 1)  # metabolites x samples
neg_log_transposed <- t(neg_log)        # samples x metabolites - transposed
neg_log_scaled <- scale(neg_log_transposed)  # metabolites x samples (scaled)

pheatmap(neg_log_scaled,
         annotation_col = neg_annotation,
         annotation_colors = anno_colors,
         clustering_distance_cols = "correlation",
         clustering_distance_rows = "correlation",
         clustering_method = "complete",
         show_colnames = FALSE,
         main = "Metabolomics heatmap (negative mode)")

# Sample–sample correlation (correlate rows of samples x metabolites => use t())
neg_sample_cor <- cor(t(neg_log_scaled), use = "pairwise.complete.obs")
pheatmap(neg_sample_cor,
         annotation_col = neg_annotation,
         annotation_row = neg_annotation,
         annotation_colors = anno_colors,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "complete",
         main = "Sample–sample correlation (negative mode)")

# Metabolite–metabolite correlation (correlate metabolites across samples)
neg_metabolite_cor <- cor(neg_log_scaled, use = "pairwise.complete.obs")  # metabolites x metabolites
pheatmap(neg_metabolite_cor,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "complete",
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Metabolite–metabolite correlation (negative mode)")

###### optional: top N metabolites ----------------------------

N <- 15
score <- rowMeans(abs(neg_metabolite_cor), na.rm = TRUE)
top_metabolites <- names(sort(score, decreasing = TRUE))[1:min(N, length(score))]
cor_top_metabolites <- neg_metabolite_cor[top_metabolites, top_metabolites]

pheatmap(cor_top_metabolites,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,     # adjust if labels overlap
         fontsize_col = 8,
         main = paste0("Top ", N, " metabolite–metabolite"))
