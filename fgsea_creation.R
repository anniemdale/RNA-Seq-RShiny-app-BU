library('tidyverse')
library('fgsea')
library('biomaRt')

# label the DESeq results based on whether or not the genes are up or down regulated
label_res <- function(deseq2_res) {
  # read in deseq results as a tibble
  deseq_tibble <- as_tibble(read.csv(deseq2_res))
  # create new column and assign values based  on fitting the conditions
  deseq_tibble <- deseq_tibble %>% mutate(volc_plot_status = case_when(
    padj < 0.01 & log2FoldChange > 0.58 ~ "UP",
    padj < 0.01 & log2FoldChange < -0.58 ~ "DOWN",
    padj >= 0.01 | abs(log2FoldChange) < 0.58 ~ "NS"
  ))
  # reorder the columns
  deseq_tibble <- deseq_tibble %>% dplyr::select(gene_id, volc_plot_status, log2FoldChange, padj, baseMean, lfcSE,stat, pvalue)
  # sort rows by padj in ascending order
  deseq_tibble <- deseq_tibble %>% arrange(padj)
  return(deseq_tibble)
}

# rank the results based on log2FC and merge with the gene symbol used in gmt file
make_ranked_log2fc <- function(labeled_results, id2gene_path) {
  # read in txt file without a header
  id <- read.delim(id2gene_path, header=FALSE)
  # rename the columns
  id <- id %>% dplyr::rename(gene_id = V1, symbol = V2)
  # combine the two files by their genes
  combined <- left_join(labeled_results, id, by = "gene_id")
  # arrange the log2FC in descending order
  ranked <- combined %>% arrange(desc(log2FoldChange))
  # create vector
  ranked_vector <- ranked$log2FoldChange
  # set vector names
  names(ranked_vector) <- ranked$symbol
  return(ranked_vector)
}

# run the fgsea using the ranked list
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  # read in gmt file
  gmt <- gmtPathways(gmt_file_path)
  # remove NAs and infinite values
  rnk_list <- rnk_list[is.finite(rnk_list)]
  # call fgsea
  fgseaRes <- fgsea(pathways = gmt,
                    stats = rnk_list,
                    minSize = min_size,
                    maxSize = max_size)
  # return results as a tibble
  fgsea_tibble <- fgseaRes %>% as_tibble()
  fgsea_df <- fgsea_tibble %>% as.data.frame()
  return(fgsea_df)
}

# create labeled results
label <- label_res('data/de_results.csv')
# create ranked vector
ranked <- make_ranked_log2fc(label, 'data/fgbn.to.cg.txt')
# create fgsea object
fgsea <- run_fgsea('data/pone.0259201.s002.gmt', ranked, 15, 500)
# mutate fgsea object for formating reasons
fgsea <- fgsea %>% mutate(leadingEdge=map_chr(leadingEdge, toString))
# write fgsea object to csv
write.csv(fgsea, "data/fgsea_res.csv", row.names=FALSE)
