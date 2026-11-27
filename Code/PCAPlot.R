library(DESeq2)
library(ggplot2)
library(stringr)

source("./R/Utils.R")
source("./R/group_analysis.R")

# - - - Specifying the custom colors for the treatments - - - - -
treatment_colors <- c("untreated" = "#999999", "mock"      = "#159D91", "flg22"     = "#D17A35", "nlp24"     = "#B05F91", "Pep1"      = "#7251A3", "SCOOP10"   = "#526A7A")

# - - - Specifying the plotting directory - - - - -
plot_dir <- "../Data/analysis_pca" # set this to the directory you want the plots saved to 

# - - - Loading the file/data - - - - - - -
file <- str_glue("./counts_genes_x_samples.tsv")
my_colData_2 <- build_coldata(file, ref_treatment ="mock", drop_treatments = c("Room1", "Room2")) # select reference and what to drop 
counts_2  <- load_counts(file, samples = rownames(my_colData_2)) # load counts

all(rownames(my_colData_2) %in% colnames(counts_2)) # checks if sample name in colData is in counts
counts_2 <- counts_2[, rownames(my_colData_2), drop = FALSE] # reorder counts to match colData

my_colData_2$temperature <- factor(my_colData_2$temperature)
my_colData_2$treatment <- factor(my_colData_2$treatment)
my_colData_2$time <- factor(my_colData_2$time)

dds_2 <- DESeqDataSetFromMatrix(countData = counts_2, colData = my_colData_2, design = ~ 1) # does not use full 3 way interaction model (specified by "design")


# - - - Run varianceStabilizingTransformation on the data - - - - -
vsd <- vst(dds_2, blind = TRUE) # blind = TRUE for unbiased comparison


# - - - Plot PCA once on all the data - - - - -
pca_all_data <- suppressMessages(plotPCA(vsd, intgroup = c("treatment", "temperature", "time"), returnData = TRUE, ntop = nrow(vsd))) # ntop = number of genes to use for PCA, default is 500

# plotPCA() returns samples as row names
pca_all_data$sample <- rownames(pca_all_data) # assigns the sample names to a new column in the PCA data frame
percentVar <- round(100 * attr(pca_all_data, "percentVar")) # makes a vector of the percentage of variance explained by each principal component



# - - - - Function to make PCA plots for subsets of the data - - - - -
make_pca_plot <- function(vsd_subset, plot_title) {
  # Keep the full-data PCA coordinates for only the requested samples
  pca_data <- pca_all_data[pca_all_data$sample %in% colnames(vsd_subset), drop = FALSE]

  pca_data$time <- factor(pca_data$time, levels = c("t.24", "t.18", "t0", "t0.5", "t2", "t6", "t24"))

  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = treatment, shape = temperature, size = time)) +
    geom_point(alpha = 0.85) +
    scale_color_manual(values = treatment_colors, drop = FALSE) +
    scale_shape_manual(values = c("T21" = 16, "T28" = 17), drop = FALSE) +
    scale_size_manual(
      values = c("t.24" = 1, "t.18" = 2, "t0" = 3, "t0.5" = 4, "t2" = 5, "t6" = 6, "t24" = 7),
      labels = c("−24 h", "−18 h", "0 h", "30 min", "2 h", "6 h", "24 h"),
      name = "Timepoint",
      drop = FALSE) +
    labs(
      title = plot_title,
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance"),
      color = "Treatment",
      shape = "Temperature") +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right")

  filename <- file.path(plot_dir, paste0("PCA_", gsub("[^A-Za-z0-9_-]+", "_", plot_title), ".png"))

  ggsave(filename = filename, plot = p, width = 10, height = 8, dpi = 300)
  p # return the plot
}


# - - - - Making PCA plots for different subsets of the data - - - - -

### All samples - - - - - - - - - -
plot <- make_pca_plot(vsd, "PCA of all RNA-seq samples")

# plot + facet_wrap(~ time) # able to display a separate plot for each coefficient


### All samples except untreated - - - - - - - - -
vsd_treated <- vsd[, vsd$treatment != "untreated"]
pca_treated <- make_pca_plot(vsd_treated, "All mock- and PDAMP-treated samples")
