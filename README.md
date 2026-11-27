# RNA-seq analysis
This repository provides various files to analyse an RNA-seq dataset.


### List of files supplied in this repository 
1. GO Enrichment 
2. Heatmap
3. PCA plots
    - Scree plots 
4. Venn diagram
5. Volcano plots


# Required Packages
The code for the initial data analysis can be found in https://github.com/JordiAlonsoEsteve/RNAseqAnalysisWorkshop/tree/main/code/Statistical_Analysis. 
- For all files it is necessary to have [utils.R](https://github.com/JordiAlonsoEsteve/RNAseqAnalysisWorkshop/blob/main/code/Statistical_Analysis/DEA.R) and [DEA.R](https://github.com/JordiAlonsoEsteve/RNAseqAnalysisWorkshop/blob/main/code/Statistical_Analysis/utils.R).
    - It is assumed that the starting file is a count table with TAIR12 gene names for running the DEA.R file. 

- For GOEnrichment [group_analysis.R](https://github.com/JordiAlonsoEsteve/RNAseqAnalysisWorkshop/blob/main/code/Statistical_Analysis/group_analysis.R) is also required.


# Specific functions of each file 
## GO Enrichment 
The GO enrichment file uses 2 functions specified in [group_analysis.R](https://github.com/JordiAlonsoEsteve/RNAseqAnalysisWorkshop/blob/main/code/Statistical_Analysis/group_analysis.R) to filter and assign GO terms based on TAIR identifiers for biological process.

## Heatmap 
The heatmap file supplies code that takes the initial fitted DESeq2 object and calculates the contrasts for the specified coefficients. 
- The current code generates a heatmap for all (significant) genes over time relative to the mock 
Some additional filtering is applied. After which the optimal number of clusters is calculated based on NbClust. 
- Important to note: NbClust provides a list of optimal cluster numbers, with the optimal number being the one selected for by the most tests. However, I have observed that "2" clusters seemed selected for the most in my dataset. I have therefore chosen to use the second most chosen optimal cluster number more often (which was 7 or 8).

The optimal number of clusters can then be specified after which classical hierarchical clustering is applied to the data and a heatmap is generated.

## PCA plot 
The PCA plot code provides code for generating PCA plots based on count data. 
- Important to note: The PCA plot uses the variance stabilized transformed counts and as such requires a different input than the "dds" supplied by [DEA.R](https://github.com/JordiAlonsoEsteve/RNAseqAnalysisWorkshop/blob/main/code/Statistical_Analysis/DEA.R).

The PCA is plotted only once, if any of the coefficients want to be removed for easier visualization of the data, this can be specified using the way shown at the bottom. 

## Scree Plot 
The Scree plot plots the variance explained by each principle component. Each principle component is plotted by variance explained as a barplot. 
- Important to note: This code uses the same @vsd as well as @plot_dir as calculated in PCA plot. This is specified again at the top of the file.

## Venn Diagram 
The Venn diagram code draws a venn diagram displaying DEGs for four different conditions. 
The DEGs are split between (significantly) up- and down-regulated genes. After which, the plots are generated. 
- Important to note: the "futile.logger" lines are only there so your terminal does not get clogged by messages following plot generation.

## Volcano Plot
The volcano plot file shrinks the Log2FC data reducing the impact of extreme fold changes, improving the reliability of the results. 
- "apeglm" is used as the shrinkage method, but "ashr" gives similar results

Next you loop over the entire dataset and create the volcano plots for each condition separately.
Some extra options include:
- The ability to print a list of TAIR gene names for the significant gene hits in one condition.
- Writing a .csv file for all significant genes for each condition 