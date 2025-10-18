# RNA-Seq-RShiny-app-BU

## Project Title: RNA-Seq Analysis RShiny App
### Course: R for Biological Sciences, BU (CDSBF 591)

### Brief Project Description: 
Developed an RShiny app to investigate gene expression in fruit flies exposed to temperature changes. This application walks the users through the data, classic RNA-Seq analysis steps, and interactive visualizations. RNA-Seq count normalization and differential expression analysis were done with DESeq2 from Bioconductor. Gene set enrichment analysis (GSEA) was done with the fgsea library. GSEA provides insight into the biological pathways that are negatively and positively enriched in the experimental condition. The goal of this project was to perform a full RNA-Seq analysis from raw count data and create interactive visualizations to walk the user through the steps of the analysis and data preparation.

### Software and Packages
IDE: RStudio

Libraries: tidyverse, DT, gplots, SummarizedExperiment (Bioconductor), DESeq2 (Bioconductor), shiny, bslib, ggplot2, colourpicker, fgsea, biomaRt

### Data Source and Description
Data Link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150450

Article Reference: Salachan, P. V., & Sørensen, J. G. (2022). Molecular mechanisms underlying plasticity in a thermally varying environment. Molecular ecology, 31(11), 3174–3191. https://doi.org/10.1111/mec.16463 
Description: This dataset consists of the raw bulk RNA-Seq counts for 80 Drosophila melanogaster individuals. These samples consist of a combination of various experimental and demographic conditions: adult, larvae, male, female, replicates, and experimental treatments.

### Data Files
- data/sample_info.csv is a breakdown of the population and experimental characteristics for each of the 80 samples
- data/gene_count_matrix.csv consists of the raw RNA-Seq counts
- data/norm_counts.csv consists of the DESeq2 normalized RNA-Seq counts
- data/de_results.csv shows the results of the differential expression analysis done with DESeq2
- data/fgbn.to.cg.txt is a conversion of the researchers gene annotation format (fly based gene names) to the gene annotation format used for our fgsea (cg)
- data/pone.0259201.s002.gmt is the reference gene matrix for fruit flies to be used in GSEA
- data/fgsea_res.csv shows the results of the functional gene set enrichment analysis done with fgsea

### Code Structure
- app.R is the RShiny application, it consists of all the interactive components, the page setups, and overall organization of the application
- helpers.R consists of the helper functions we use to create any figure, table, or result that is displayed in the app
- metadata_creation.R takes the sample information and extracts the important information that is then stored in a table
- DESeq_norm_res_creation.R is where the count normalization and differential expression analysis is done, which results in the normalized counts matrix and the DE results
- fgsea_creation.R is where the FGSEA is performed and returns the fgsea results
