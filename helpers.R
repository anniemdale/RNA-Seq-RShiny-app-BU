library('tidyverse')
library('DT')
library('gplots')

# SAMPLE DATA EXPLORATION
# create tibble that summarizes the sample information
sample_summary <- function(df) {
  # extract column names
  columns <- colnames(df)
  # determine class type of each column
  type <- sapply(df, class)
  # extract unique values for each column
  values <- sapply(df, unique)
  # put results into a tibble
  tibble <- tibble(
    'Column Name' = columns,
    Type = type,
    value = values
  )
  # mutate tibble for formatting
  tibble <- tibble %>% mutate(value=map_chr(value, toString))
  return(tibble)
}

# turn sample info into a filterable table
sample_datatable <- function(df) {
  table <- datatable(df)
  return(table)
}

# create histograms/bar plots for the counts of unique values in each column
create_sample_histogram <- function(df, column) {
  histogram <- df %>% 
    ggplot(aes(x=!!sym(column))) +
      geom_bar()
  return(histogram)
}

# NORMALIZED COUNTS DATA EXPLORATION
# filter counts based on variance percentine and number of non-zero genes
filtered_counts <- function(df, var_slider, zero_slider) {
  # convert gene_id column to rownames
  df <- df %>% column_to_rownames(var="gene_id")
  # calculate variance for df
  variance <- sapply(df,var)
  # calculate variance percentiles
  variance_threshold <- quantile(variance, probs = var_slider/100)
  # filter based on variance percentile
  filtered_genes <- df[variance >= variance_threshold, ]
  # get non-zero genes
  non_zero_counts <- apply(df, 1, function(row) sum(row != 0))
  # filter to get zero_slider of non-zero counts
  filter <- head(filtered_genes, n=zero_slider)
  return(filter)
}

# create a summary saying how many genes are being filtered
count_summary_table <- function(df, var_slider, zero_slider) {
  # create filtered object
  filter <- filtered_counts(df, var_slider, zero_slider)
  # count number of samples
  sample_count <- ncol(df) - 1
  # count number of genes
  gene_count <- nrow(df)
  # count number of genes passing filter
  number_pass <- nrow(filter)
  # calculate % genes passing filter
  percent_pass <- (number_pass/gene_count)*100
  # count number of genes failing filter
  number_fail <- gene_count - number_pass
  # calculate % genes failing filter
  percent_fail <- 100 - percent_pass
  # create tibble holding this information
  tibble <- tibble(
    'Number of Samples' = sample_count,
    'Total Number of Genes' = gene_count,
    'Number Genes Passing Filter' = number_pass,
    '% Genes Passing Filter' = percent_pass,
    'Number Genes Failing Filter' = number_fail,
    '% Genes Failing Filter' = percent_fail
  )
  return(tibble)
}

# function to calculate how many zeros there are in a vector
sum_zero <- function(x, na.rm=TRUE) {
  return(sum(x==0, na.rm=TRUE))
}

# create tibble containing stats for all the genes
counts_stats_tibble <- function(df) {
  # create gene_id vector
  gene_id <- df$gene_id
  # convert df gene_id column to rownames
  df <- df %>% column_to_rownames(var="gene_id")
  # calculate median for each gene
  median <- apply(df, 1,median, na.rm=TRUE)
  # calculate variance for each gene
  variance <- apply(df, 1,var)
  # calculate number of zeros for each gene
  zero_counts <- apply(df, 1,sum_zero)
  # scale the variance for plotting
  log_variance <- log10(variance)
  # create tibble based on above calculations
  tibble <- tibble(
    gene_id = gene_id,
    median = median,
    variance = variance,
    log10variance = log_variance,
    zero_counts = zero_counts
  )
  return(tibble)
}

# extract stats for genes passing the filters
filter_stats <- function(df, var_slider, zero_slider) {
  # create count stats object
  counts_stats <- counts_stats_tibble(df)
  # caclculate quantiles for variance threshold
  variance_threshold <- quantile(counts_stats$variance, probs = var_slider/100)
  # only keep the stats the pass the variance threshold
  filter_stats <- counts_stats[counts_stats$variance >= variance_threshold, ]
  # only return the filtered stats and number of genes decided by zero filter
  filtered_stats <- head(filter_stats, n=zero_slider)
  # convert gene_id to rownames
  filtered_stats <- filtered_stats %>% column_to_rownames(var="gene_id")
  return(filtered_stats)
}

# extract stats for genes failing the filters
fail_filter_stats <- function(df, var_slider, zero_slider) {
  count_stats <- counts_stats_tibble(df)
  filtered_stats <- filter_stats(df, var_slider, zero_slider)
  # get count stats for genes that aren't passing the filter
  fail_filter <- count_stats %>% filter(!(count_stats$gene_id %in% filtered_stats$gene_id))
  fail_filter <- fail_filter %>% as.data.frame()
  fail_filter <- fail_filter %>% column_to_rownames(var="gene_id")
  return(fail_filter)
}

# create scatter plot for variance and median for the genes
count_variance_scatter <- function(df, var_slider, zero_slider) {
  # stats that pass the filter
  filtered_stats <- filter_stats(df, var_slider, zero_slider)
  # specify that they are passing
  filtered_stats <- filtered_stats %>% mutate(filter_status='pass')
  # stats that fail the filter
  fail_filter <- fail_filter_stats(df, var_slider, zero_slider)
  # specify that they are failing
  fail_filter <- fail_filter %>% mutate(filter_status='fail')
  # bind the failed and passing stats
  filter_status_df <- rbind(fail_filter,filtered_stats)
  # create scatter plot, color based on passing/failing filter
  variance_scatter <- ggplot(filter_status_df, aes(x=log10variance, y=median, color=filter_status)) +
    geom_point(alpha=0.5) +
    scale_color_manual(values=c("lightblue","mediumvioletred")) +
    labs(title = "Log10(Variance) v. Median")
  return(variance_scatter)
}

# create scatter plot for number of zeros and median for the genes
count_zero_scatter <- function(df, var_slider, zero_slider) {
  filtered_stats <- filter_stats(df, var_slider, zero_slider)
  filtered_stats <- filtered_stats %>% mutate(filter_status='pass')
  fail_filter <- fail_filter_stats(df, var_slider, zero_slider)
  fail_filter <- fail_filter %>% mutate(filter_status='fail')
  filter_status_df <- rbind(fail_filter, filtered_stats)
  zero_scatter <- ggplot(filter_status_df, aes(x=zero_counts, y=median, color=filter_status)) +
    geom_point(alpha=0.5) +
    scale_color_manual(values=c("lightblue","mediumvioletred")) +
    labs(title = "Zero Counts v. Median", x="# zero counts")
  return(zero_scatter)
}

# create clustered heatmap 
count_heatmap <- function(df, var_slider, zero_slider) {
  filter_counts <- filtered_counts(df, var_slider, zero_slider)
  # log scale the counts so there is more contrast in the colors
  filter_counts <- filter_counts %>% mutate(across(everything(), log10))
  # only keep rows that don't contain infinite values
  filter_counts <- filter_counts[is.finite(rowSums(filter_counts)), ]
  # heatmap.2 from gplots requires matrix input
  filter_matrix <- filter_counts %>% as.matrix()
  heatmap <- heatmap.2(filter_matrix, trace="none", xlab="Samples")
  return(heatmap)
}

# create PCA to show clustering in the samples based on filtered counts
plot_pca <- function(count_df, sample_df, first_pc, second_pc, var_slider, zero_slider) {
  # get filtered counts
  counts <- filtered_counts(count_df, var_slider, zero_slider)
  # need to discard gene ids
  counts <- counts %>% rownames_to_column("gene_id")
  counts <- counts %>% dplyr::select(-gene_id)
  # perform PCA on the counts
  pca_results <- prcomp(t(counts))
  # grab just the pcs of choice
  pc_first <- pca_results$x[, first_pc]
  pc_second <- pca_results$x[, second_pc]
  # create a tibble of the PCs in order to plot
  pcs <- tibble(pc1=pc_first, pc2=pc_second)
  # combine with the sample info
  combined <- bind_cols(pcs, sample_df)
  # calculate explained_variance
  explained_variance <- pca_results$sdev^2 / sum(pca_results$sdev^2)
  # get the variance explained by each of the selected PCs
  pc1_var <- round(explained_variance[first_pc] * 100, 1)
  pc2_var <- round(explained_variance[second_pc] * 100, 1)
  # plot PC1 against PC2 colored by the treatment of the samples
  plot <- combined %>% ggplot(aes(x=pc1, y=pc2, color=treatment)) + 
    geom_point(aes(shape=sex)) + 
    labs(title = paste0("PCA Plot of PC", first_pc, " v. PC", second_pc),
         x = paste0("PC",first_pc, ": ", pc1_var, "% variance"),
         y = paste0("PC",second_pc, ": ", pc2_var, "% variance")
    )
  return(plot)
}

# DESEQ RESULTS DATA EXPLORATION
# create sortable table of the deseq results
de_datatable <- function(res_df) {
  table <- datatable(res_df)
  return(table)
}

# create volcano plot of the results based on the researchers' parameters
plot_volcano <- function(res_df) {
  # add column for significance levels
  volcano_data <- res_df %>%
    mutate(significance = ifelse((abs(log2FoldChange) > 0.58 & pvalue < 0.01), "high", "low"))
  # make the plot
  plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color=significance)) +
    geom_point() +
    scale_color_manual(values = c("high" = "red", "low" = "black")) +
    labs(x = "log2FoldChange", y = "-log10(padj)", title = "Volcano Plot DEG Based on Treatment") +
    theme_minimal()
  return(plot)
}

# FGSEA RESULTS DATA EXPLORATION
# get the number of pathways to look at based on the padj value slider
num_pathways <- function(fgsea_res, padj_slider) {
  filtered <- fgsea_res[fgsea_res$padj < padj_slider, ]
  num <- nrow(filtered)
  return(num)
}

# top positive enriched pathways
up_pathways <- function(fgsea_res, padj_slider) {
  # get number of pathways
  num_paths <- num_pathways(fgsea_res, padj_slider)
  # filter fgsea for positive
  top_up <- fgsea_res %>%
    filter(NES > 0) %>%
    arrange(pval) %>%
    # extract number of rows
    slice_head(n=num_paths) %>%
    mutate(direction = "Positive NES")
  return(top_up)
}

# top negative enriched pathways
down_pathways <- function(fgsea_res, padj_slider) {
  num_paths <- num_pathways(fgsea_res, padj_slider)
  # filter fgsea for negative
  top_down <- fgsea_res %>% 
    filter(NES < 0) %>% 
    arrange(pval) %>% 
    slice_head(n = num_paths) %>% 
    mutate(direction = "Negative NES")
  return(top_down)
}

# create horizontal bar plot for the top positive and negative pathways
top_pathways_bar <- function(fgsea_res, padj_slider) {
  # filter out the top up pathways
  topPathwaysUp <- up_pathways(fgsea_res, padj_slider)
  # filter out the top down pathways
  topPathwaysDown <- down_pathways(fgsea_res, padj_slider)
  # bind the two dfs together
  topPathways <- bind_rows(topPathwaysUp, topPathwaysDown) %>%
    arrange(desc(NES))
  #plot pathways with color for positive/negative directions
  plot <- topPathways %>% ggplot(aes(x=NES, y=reorder(pathway,NES), fill=direction)) +
    geom_col() +
    labs (x="Normalized Enrichment Score (NES)", y="", title="FGSEA Results Top Pathways") +
    theme(axis.text.y = element_text(size=6), legend.position="none")
  return(plot)
}

# create table for the the top positive and negative and both pathways
fgsea_table <- function(fgsea_res, padj_slider, value_button) {
  # get top positive pathways
  topPathwaysUp <- up_pathways(fgsea_res, padj_slider)
  # get top negative pathways
  topPathwaysDown <- down_pathways(fgsea_res, padj_slider)
  # get all top pathways 
  topPathways <- bind_rows(topPathwaysUp, topPathwaysDown) %>%
    arrange(desc(NES))
  # return either positive, negative or all based on user selection
  if(value_button=="All") {
    return(topPathways)
  } else if(value_button=="Negative") {
    return(topPathwaysDown)
  } else if(value_button=="Positive") {
    return(topPathwaysUp)
  }
}

# create scatter plot of the pathways passing the filter
fgsea_scatter <- function(fgsea_res, padj_slider) {
  # mutate the fgsea results to denote which pathways are passing/failing the filter
  fgsea <- fgsea_res %>%
    mutate(filter = ifelse(padj < padj_slider, "pass", "fail"))
  # color points based on whether they pass/fail the filter
  plot <- ggplot(fgsea, aes(x = NES, y = -log10(pval), color=filter)) +
    geom_point() +
    scale_color_manual(values = c("pass" = "darkgreen", "fail" = "grey")) +
    labs(x = "Normalized Enrichment Score (NES)", y = "-log10(pval)", title = "FGSEA Scatter Based on Filter") +
    theme_minimal()
  return(plot)
}