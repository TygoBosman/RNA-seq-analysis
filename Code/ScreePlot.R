# - - - This is a continuation of PCAPlot code - - - - -
# # # Rely on the same @vsd as input, as well as the same @plot_dir for the generation of the Scree plot

library(DESeq2)
library(ggplot2)

plot_scree <- function(vsd, standardize = FALSE, n_components = NULL) {
  x <- as.matrix(t(assay(vsd)))
  storage.mode(x) <- "double"

  if (anyNA(x) || any(!is.finite(x))) {
    stop("`Vsd` must contain only finite numeric values.")
  }

  max_components <- min(dim(x))
  if (max_components < 1) {
    stop("`Vsd` must have at least one row and one column.")
  }

  if (is.null(n_components)) {
    n_components <- max_components
  } else {
    if (length(n_components) != 1 || is.na(n_components) || n_components < 1) {
      stop("`n_components` must be a positive number.")
    }
    n_components <- min(as.integer(n_components), max_components)
  }

  # Center columns; optionally scale them to unit variance.
  pca <- prcomp(x, center = TRUE, scale. = standardize) # calculates the PCA 

  # Calculate proportions using the variance from the full PCA.
  variance <- pca$sdev^2 / sum(pca$sdev^2)
  variance <- variance[seq_len(n_components)]

  pcs <- seq_along(variance)
  ticks <- unique(c(
    1,
    if (length(pcs) >= 5) seq(5, length(pcs), by = 5),
    length(pcs)
  ))

  scree_data <- data.frame(
    PC = pcs,
    variance_explained = variance * 100
  )

  p_scree <- ggplot(scree_data, aes(x = PC, y = variance_explained)) +
    geom_col(aes(fill = "Individual"), width = 0.65, colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = c("Individual" = "#4C78A8"), name = NULL) +
    scale_x_continuous(breaks = ticks) +
    labs(
      x = "Principal component number",
      y = "Variance explained (%)",
      title = "PCA scree plot"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )

  # Save the scree plot to the specified directory
  ggsave(
  filename = file.path(plot_dir, "Scree_plot_all_samples.png"), 
  plot = p_scree,
  width = 10,
  height = 8,
  dpi = 300
  )

  p_scree # return the plot
}

scree_plot <- plot_scree(vsd, standardize = FALSE)
