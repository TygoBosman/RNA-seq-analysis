library(DESeq2)
library(stringr)
library(ggplot2)
library(VennDiagram)

# all directories have to be changed to your own directories
dds       <- readRDS(str_glue("../Data/PDAMP/dds_object.rds"))
contrasts <- readRDS(str_glue("../Data/PDAMP/contrasts_against_mock.rds"))
plot_dir <- "../Data/Venn_diagram"

# Select four contrasts
selected_keys <- c(
  "T21__flg22__t0.5",
  "T21__flg22__t2",
  "T21__flg22__t6",
  "T21__flg22__t24"
)

# Stop early if one of the requested contrasts is missing
missing_keys <- setdiff(selected_keys, names(contrasts$results_list))
if (length(missing_keys) > 0) {
  stop("Missing contrast(s): ", paste(missing_keys, collapse = ", "))
}

# Extract the result tables and give each set a short time-point label
res_by_time <- contrasts$results_list[selected_keys]
names(res_by_time) <- paste0(
  sub("^.*__t", "", selected_keys),
  " h"
)

# Create one gene list per time point, separately for up and down DEGs
up_sets <- lapply(
  res_by_time,
  get_deg_direction,
  direction = "up",
  padj_cutoff = 1e-6,
  lfc_cutoff = 2
)

down_sets <- lapply(
  res_by_time,
  get_deg_direction,
  direction = "down",
  padj_cutoff = 1e-6,
  lfc_cutoff = 2
)

# - - - Change the level of logs to stop INFO messages - - - - -
futile.logger::flog.threshold(
  futile.logger::WARN,
  name = "VennDiagramLogger"
)

# Helper to draw a four-set Venn diagram
make_venn <- function(gene_sets, direction) {
  
  direction_title <- if (direction == "up") {
    "Upregulated"
  } else {
    "Downregulated"
  }
  
  old_level <- futile.logger::flog.threshold()
  futile.logger::flog.threshold(futile.logger::WARN)
  
  tryCatch(
    VennDiagram::venn.diagram(
      x = gene_sets,
      filename = file.path(
        plot_dir,
        paste0("DEG_Venn_T21_flg22_", direction, ".png")
      ),

      main = paste0(
        direction_title,
        " DEGs: flg22 vs mock at T21\n",
        "padj < 1e-6, abs(log2FC) >= 2"
      ),
      main.cex = 1.6,
      main.fontface = "bold",
      main.col = "#253247",
      main.fontfamily = "sans",
      main.pos = c(0.5, 0.93),

      # colours of the circles and their borders
      fill = c("#228FF5", "#F4360B", "#1CE330", "#BB1CE3"),
      col = c("#0B83F4", "#C82C09", "#17BA27", "#9917BA"),
      alpha = 0.4,
      lwd = 2,
      lty = "solid",

      cex = 1.2,
      fontface = "bold",
      fontfamily = "sans",
      label.col = "#253247",

      cat.col = c("#0B83F4", "#C82C09", "#17BA27", "#9917BA"), # colour of the labels 
      cat.cex = 1.1,
      cat.fontface = "bold",
      cat.fontfamily = "sans",
      cat.default.pos = "outer",

      margin = 0.12,
      imagetype = "png",
      height = 2400,
      width = 2800,
      resolution = 300
    ),
    finally = futile.logger::flog.threshold(old_level)
  )
}

# Draw separate diagrams for upregulated and downregulated DEGs
make_venn(up_sets, "up")
make_venn(down_sets, "down")
