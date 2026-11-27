library(ggplot2)
library(scales)

file <- "../Data/counts_genes_x_samples.tsv"
colData <- build_coldata(file)
counts  <- load_counts(file, samples = rownames(colData))

plot_dir <- "../Data/"

plot_raw_counts <- function(counts) {
  count_df <- as.data.frame(counts, check.names = FALSE)

  # All columns are sample counts; gene IDs are already row names
  count_df[] <- lapply(count_df, function(z) {
    as.numeric(as.character(z))
  })

  count_mat <- as.matrix(count_df)

  if (anyNA(count_mat)) {
    stop("Count columns contain missing or non-numeric values.")
  }

  sample_totals <- data.frame(
    sample = colnames(count_mat),
    total_counts = colSums(count_mat)
  )

  # Arrange samples by total count
  sample_totals <- sample_totals[
    order(sample_totals$total_counts),
    ,
    drop = FALSE
  ]

  # Space sample positions 1.5 units apart
  sample_totals$position <- seq_len(nrow(sample_totals)) * 1.5

  p_raw_counts <- ggplot(
    sample_totals,
    aes(x = position, y = total_counts)
  ) +
    geom_col(fill = "#0072B2", width = 0.9) +
    scale_x_continuous(
      breaks = sample_totals$position,
      labels = sample_totals$sample
    ) +
    coord_flip() +
    labs(
      title = "Total raw counts per sample",
      x = "Sample",
      y = "Total counts"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 7)
    )

  # Add plot height as the number of samples increases
  plot_height <- max(8, 0.22 * nrow(sample_totals))

  ggsave(
    filename = file.path(plot_dir, "Raw_counts.png"),
    plot = p_raw_counts,
    width = 10,
    height = plot_height,
    dpi = 300
  )

  p_raw_counts # returns the plot
}

p_raw_counts <- plot_raw_counts(counts)
p_raw_counts