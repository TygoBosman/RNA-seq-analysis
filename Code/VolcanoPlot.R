library(DESeq2)
library("apeglm")
library(EnhancedVolcano)
library(dplyr)

dds <- readRDS("../Data/Processed/dds_fitted.rds")
# res <- results(dds)
# print(res)

outdir <- "./lfcShrink"

# - - All coefficients - - - - - - - - - - - - - -
coef_names_full <- setdiff(resultsNames(dds), "Intercept")

# Check that these are exact coefficient names in the fitted model
missing_names <- setdiff(coef_names_full, resultsNames(dds))
if (length(missing_names) > 0) {
  stop("These names are not in resultsNames(dds): ",
       paste(missing_names, collapse = ", "))
}

# - - - - Shrink the LFC - - - - - - - - - - - - - - - - - - - - -
bp7 <- BiocParallel::SnowParam(workers = 7, type = "SOCK") # using 7 different cores for parallel processing

lfc_trimmed_results <- setNames(
  vector("list", length(coef_names_full)),
  coef_names_full
)

for (i in seq_along(coef_names_full)) {
  nm <- coef_names_full[i]
  message("Shrinking ", i, "/", length(coef_names_full), ": ", nm)

  timing <- system.time({
    lfc_trimmed_results[[nm]] <- lfcShrink(
      dds,
      coef = nm,
      type = "apeglm", # method for shrinkage, can also use "ashr" or "normal"
      parallel = TRUE,
      BPPARAM = bp7,
      quiet = TRUE
    )
  })

  print(timing)

  safe_name <- gsub("[^A-Za-z0-9_-]", "_", nm)
  saveRDS(
    lfc_trimmed_results[[nm]],
    file.path(outdir, paste0(sprintf("%02d", i), "_", safe_name, ".rds"))
  )
}

saveRDS(
  lfc_trimmed_results,
  file.path(outdir, "all_lfc_results.rds") # Save the shrunk file
)


# - - - - - Plot the things - - - - - - - - - - - - - - - -
p_cutoff <- 1e-6
fc_cutoff <- 2

plot_dir <- "../Data/Volcano_plots" # output directory for the plots

# Store plot objects in a list so they can also be reused in R
volcano_plots <- list()

# Make an empty dictionary to append the sig hit under the name of compared treatment 
sig_hit_all <- list()

for (nm in names(lfc_trimmed_results)) {

  message("Creating volcano plot for: ", nm)

  # Extract the current lfcShrink result
  resLFC <- lfc_trimmed_results[[nm]]

  if (is.null(resLFC)) {
    warning("Skipping ", nm, ": result is NULL")
    next
  }
  # Keep genes passing both cutoffs
  sig_hit <- as.data.frame(resLFC)

  sig_hit <- sig_hit[
    !is.na(sig_hit$padj) &
      !is.na(sig_hit$log2FoldChange) &
      sig_hit$padj < p_cutoff &
      abs(sig_hit$log2FoldChange) >= fc_cutoff,
    ,
    drop = FALSE
  ]

  sig_hit_all[[nm]] <- sig_hit

  # Count significant genes
  n_sig <- nrow(sig_hit)

  # Select up to 10 significant genes with the smallest adjusted p-values
  top_genes <- rownames(sig_hit)[
    head(order(sig_hit$padj), min(10, n_sig))
  ]

  # Create the volcano plot
  volcano_plot <- EnhancedVolcano(
    toptable = resLFC,
    x = "log2FoldChange",
    y = "padj",
    xlim = c(-7.5, 7.5),
    lab = rownames(resLFC),
    selectLab = top_genes,
    pCutoff = p_cutoff,
    FCcutoff = fc_cutoff,
    pointSize = 2.0,
    colAlpha = 0.5,
    drawConnectors = TRUE,
    widthConnectors = 0.4,
    title = paste0(
      nm,
      "\n(fold-change cutoff = ", fc_cutoff,
      ", adjusted p-value cutoff = ", p_cutoff, ")"
    )
  ) +
    ggplot2::labs(
      caption = paste(
        "Number of genes with padj <", p_cutoff,
        "and |log2 FC| >=", fc_cutoff, ":", n_sig
      )
    ) +
    ggplot2::theme(
      plot.caption.position = "plot",
      plot.caption = ggplot2::element_text(hjust = 0)
    )

  # Store the plot in the list
  volcano_plots[[nm]] <- volcano_plot

  # Make the coefficient name safe for use as a filename
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", nm)

  # Save the plot
  ggplot2::ggsave(
    filename = file.path(plot_dir, paste0(safe_name, "_volcano.png")),
    plot = volcano_plot,
    width = 8,
    height = 6,
    units = "in",
    dpi = 300
  )
}

# - - print the TAIR gene names - - - - - - - - - - - - - - 
print(sig_hit_all[["treatment_flg22_vs_mock"]]) # for each of the compared coefficients 


# - - Save the significant genes to a CSV file - - - - - - - - 
keep <- vapply(
  sig_hit_all,
  function(x) !is.null(x) && nrow(x) > 0,
  logical(1)
)

nonempty_sig <- sig_hit_all[keep]

if (length(nonempty_sig) == 0) {
  stop("No significant genes to write.")
}

combined_sig <- do.call(rbind, lapply(names(nonempty_sig), function(nm) {
  df <- as.data.frame(nonempty_sig[[nm]])

  data.frame(
    coefficient = nm,
    gene = rownames(df),
    df,
    row.names = NULL,
    check.names = FALSE
  )
}))

write.csv(
  combined_sig,
  file = file.path(plot_dir, "significant_genes_all_coefficients.csv"),
  row.names = FALSE
)