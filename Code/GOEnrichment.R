library(DESeq2)
library(clusterProfiler)
library(org.At.tair.db)
library(enrichplot)
library(ggplot2)
library(stringr)

# load the data (directories have to be changed to your own directories)
dds       <- readRDS(str_glue("./dds_object.rds"))
contrasts <- readRDS(str_glue("./contrasts_agains_mock.rds"))

plot_dir <- str_glue("C:/Users/tjcbo/OneDrive/Documenten/uni/MCLS/Year 1/Major Research Project (GSLS-MAJRP)/Code/Data/GO enrichment")


# group analysis
deg_obj <- get_deg_sets(dds, contrasts)

go_results <- setNames(
  lapply(names(deg_obj$sets), function(x) list(up = NULL, down = NULL)),
  names(deg_obj$sets)
)

for (contrast_name in names(deg_obj$sets)) {
  for (direction in c("up", "down")) {

    genes <- deg_obj$sets[[contrast_name]][[direction]]
    if (length(genes) == 0) next

    go <- run_arabidopsis_go(gene_list = genes, ont = "BP")
    go_results[[contrast_name]][[direction]] <- go

    if (!is.null(go$plot)) {
      safe_name <- gsub("[^A-Za-z0-9_-]+", "_", contrast_name)
      file_path <- file.path(
        plot_dir,
        paste0(safe_name, "_", direction, "_BP.png")
      )
    
      ggplot2::ggsave(
        filename = file_path,
        plot = go$plot,
        width = 9,
        height = 6,
        dpi = 300
      )
    }
  }
}