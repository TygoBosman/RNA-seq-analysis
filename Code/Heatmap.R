library(DESeq2)
library(pheatmap)
library(NbClust)
library(mclust)

dds <- readRDS("../Data/Processed/dds_fitted.rds") # directories have to be changed to your own directories
# res <- results(dds)
# print(res)


# - - Starting settings - - - - - - - - - - - -
p_cutoff <- 0.05
treatment_of_interest <- "Pep1" # not used in this script
time_numeric <- c(0.5, 2, 6, 24)

# - - Treatment vs mock effect at each timepoint at 21C - - - - - - -
coef_by_time <- list(
  "t0.5" = "treatment_Pep1_vs_mock", # Change the treatment tested here to any of the treatments in the model (e.g., flg22, nlp24, SCOOP10)
  "t2"   = c("treatment_Pep1_vs_mock", "treatmentPep1.timet2"),
  "t6"   = c("treatment_Pep1_vs_mock", "treatmentPep1.timet6"),
  "t24"  = c("treatment_Pep1_vs_mock", "treatmentPep1.timet24")
)



# - - Something to do with contrast - - - - - - -
bp <- BiocParallel::SnowParam(
  workers = 6, # The amount of cores used in the calculation
  type = "SOCK"
)

DESeq_time <- system.time({
  res_by_time <- lapply(coef_by_time, function(coefs) {
    DESeq2::results(
      dds,
      contrast = list(coefs),
      parallel = TRUE,
      BPPARAM = bp
    )
  })
}) # Takes 11 minutes to run on 6 cores for 24,003 genes and 4 timepoints



# Assemble log2 fold changes and adjusted p-values into genes x time matrices
lfc_mat <- do.call(
  cbind,
  lapply(res_by_time, function(res) res$log2FoldChange)
)
padj_mat <- do.call(
  cbind,
  lapply(res_by_time, function(res) res$padj)
)

rownames(lfc_mat) <- rownames(res_by_time[[1]])
rownames(padj_mat) <- rownames(res_by_time[[1]])
colnames(lfc_mat) <- names(coef_by_time)
colnames(padj_mat) <- names(coef_by_time)

# Keep genes significant at one or more timepoints and with finite LFCs
significant_at_any_time <- rowSums(
  padj_mat < p_cutoff,
  na.rm = TRUE
) > 0

complete_lfc <- rowSums(is.finite(lfc_mat)) == ncol(lfc_mat)
eligible <- significant_at_any_time & complete_lfc
# sum(eligible) # 14689

if (!any(eligible)) {
  stop("No genes pass the adjusted p-value cutoff at any timepoint.")
}

lfc_eligible <- lfc_mat[eligible, , drop = FALSE]

# Limit the input size for faster EM fitting
gene_means <- rowMeans(lfc_eligible)

gene_var <- rowSums(sweep(lfc_eligible, 1, gene_means, FUN = "-")^2) / (ncol(lfc_eligible) - 1)
names(gene_var) <- rownames(lfc_eligible)

selected_genes <- names(sort(gene_var, decreasing = TRUE))[
  seq_len(min(24003, length(gene_var))) # can edit the amount of genes tested (total = 24,003)
  ]

# Keep the actual log2 fold changes: positive = higher in treatment than mock, negative = reverse
cluster_input <- as.data.frame(lfc_eligible[selected_genes, , drop = FALSE], check.names = FALSE)




# - - Check the number of clusters that is optimal - - - - - - -
x_all <- as.matrix(lfc_eligible)
finite_rows <- apply(x_all, 1, function(x) all(is.finite(x)))
x_all <- x_all[finite_rows, , drop = FALSE]

# Use the top 2,000 variable genes as a pool
gene_var_all <- apply(x_all, 1, var)
names(gene_var_all) <- rownames(x_all)

pool_size <- min(2000L, nrow(x_all))
gene_pool <- names(sort(gene_var_all, decreasing = TRUE))[seq_len(pool_size)]

# Draw a different subset of 500 genes
set.seed(101)  # controls which genes are sampled
genes_subset <- sample(gene_pool, min(500L, length(gene_pool))) # able to change test pool size by editing the 500L number here
x_test <- as.data.frame(x_all[genes_subset, , drop = FALSE])

# - - Run NbClust - - - - - - - - - - - - -
set.seed(1)
nb_subset <- NbClust(
  data = x_test,
  distance = "euclidean",
  min.nc = 2,
  max.nc = min(75, nrow(x_test) - 1), # able to change the max amount of clusters checked here 
  method = "kmeans"
)

print(nb_subset$Best.nc)




# - - hclust - - - - - - - - - - - - - - - - - - - - - -
optimal_k <- 7 # Cluster number to test (tested with NbClust, can be changed to any number of clusters)

# cluster_input: genes in rows, timepoints in columns
x <- as.matrix(cluster_input)
stopifnot(all(is.finite(x)))

# Hierarchical clustering
gene_dist <- dist(x, method = "euclidean")
hc <- hclust(gene_dist, method = "ward.D2")

# Cut the tree into the number of clusters suggested by NbClust
gene_clusters <- cutree(hc, k = optimal_k)

# Cluster sizes and gene assignments
# print(table(gene_clusters))
cluster_assignments <- data.frame(
  gene = rownames(x),
  cluster = gene_clusters
)

# - - Center colour scale at zero; use quantiles to reduce influence of outliers - - - - - - - - - -
x <- as.matrix(cluster_input)
lims <- quantile(as.vector(x), probs = c(0.02, 0.98), na.rm = TRUE)
lim <- max(abs(lims))
if (!is.finite(lim) || lim == 0) lim <- 1



# - - Plot and save the heatmaps - - - - - - - - - -
heatmap_dir <- "../Data/Heatmaps"

hclust_heatmap <- pheatmap::pheatmap(
  x,
  cluster_rows = hc,       # displays the hierarchical dendrogram
  cluster_cols = FALSE,    # keeps timepoints in their existing order
  annotation_row = annotation_row,
  show_rownames = FALSE,
  scale = "none",          # displays the original log2 fold changes
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-lim, lim, length.out = 101),
  treeheight_row = 60,
  main = paste("Pep1 vs mock at T21; k =", optimal_k),
  filename = file.path(
    heatmap_dir,
    "hclust_Pep1_vs_mock_T21_hierarchical_heatmap.png"
  )
)
