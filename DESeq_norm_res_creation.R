library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')

#' Here we are running DESeq2 in order to normalize our counts and analyze them to see which 
#' genes are being differentially expressed (up or down regulated) based on the experimental condition

# make summarized experiment object
make_se <- function(counts_csv, metafile_csv, selected_treatment) {
  # read in counts and sample info
  counts <- as_tibble(read.csv(counts_csv))
  meta <- as_tibble(read.csv(metafile_csv))
  # select the columns of interest
  meta <- meta %>% dplyr::select(sample, treatment)
  treatment <- selected_treatment
  # create count matrix
  count_mat <- as.matrix(counts[-1])
  row.names(count_mat) <- counts$gene_id
  # specify the selected timepoints in the meta data
  selected_meta <- meta %>% filter(treatment %in% selected_treatment)
  # set timepoint as a factor
  selected_meta$treatment <- factor(selected_meta$treatment)
  # set the reference level
  selected_meta$treatment <- relevel(selected_meta$treatment, ref = "control")
  # have the count matrix columns only be the ones that are in the subsetted meta data
  count_mat <- count_mat[, selected_meta$sample]
  # create summarized experiment
  sum_exp <- SummarizedExperiment(assays=SimpleList(counts=count_mat),
                                  colData=selected_meta)
  return(sum_exp)
}

# initialize DESeq analysis
create_dds <- function(se, design){
  dds <- DESeqDataSetFromMatrix(countData = assay(se),
                                colData = colData(se),
                                design = design)
  dds <- DESeq(dds)
  return(dds)
}

# normalize the counts using DESeq
norm_dds <- function(se, design) {
  dds <- create_dds(se, design)
  norm <- counts(dds, normalized=TRUE)
  norm <- norm %>% as.data.frame()
  norm <- rownames_to_column(norm, "gene_id")
  return(norm)
}

# extract the DESeq results
res_dds <- function(se, design) {
  dds <- create_dds(se, design)
  res <- as.data.frame(results(dds))
  res <- rownames_to_column(res, "gene_id")
  return(res)
}

# create summarized experiment object
se <- make_se('data/gene_count_matrix.csv', 'data/sample_info.csv', c("control","fluctuating"))
# create normalized counts
norm <- norm_dds(se, design=~treatment)
# create DESeq results object
res <- res_dds(se, design=~treatment)
# write results to csv
write.csv(res, "data/de_results.csv", row.names=FALSE)
# write normalized counts to csv
write.csv(norm, "data/norm_counts.csv", row.names=FALSE)

