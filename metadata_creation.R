library('tidyverse')

#' Here we are creating the metadata or sample information based on the column names which represent Samples.
#' This will be used in the DESeq analysis and normalization.
#' 
#' Sample Names are formatted for ex: CAH1F or CLH1
#' First position specifies treatment
#' Second position specifies lifestage
#' Third position specifies timepoint
#' Fourth position specifies period and replicate
#' Last position specifies sex

# Extract Treatment from Sample Name (C or F)
treatment_from_sample <- function(x) {
  x_sub <- substring(x, 1, 1)
  return(x_sub)
}
# Extract Lifestage from Sample Name (A or L)
lifestage_from_sample <- function(x) {
  x_sub <- substring(x,2,2)
  return(x_sub)
}
# Extract Timepoint from Sample Name (H or L)
timepoint_from_sample <- function(x) {
  x_sub <- substring(x,3,3)
  return(x_sub)
}
# Extract Period from Sample Name (1 through 8)
period_from_sample <- function(x) {
  x_sub <- substring(x,4,4)
  return(x_sub)
}
# Extract Sex from Sample Name (F or M or blank)
sex_from_sample <- function(x) {
  x_sub <- substring(x,5,5)
  return(x_sub)
}
# function to create the sample information based on sample names using the above helper functions
meta_info_from_labels <- function(sample_names) {
  treatment <- sapply(sample_names, FUN=treatment_from_sample)
  lifestage <- sapply(sample_names, FUN=lifestage_from_sample)
  timepoint <- sapply(sample_names, FUN=timepoint_from_sample)
  period <- sapply(sample_names, FUN=period_from_sample)
  sex <- sapply(sample_names, FUN=sex_from_sample)
  # create tibble with columns and information that was extracted
  tibble <- tibble(sample = sample_names,
                   treatment = treatment,
                   lifestage = lifestage,
                   timepoint = timepoint,
                   period = period,
                   sex = sex,
                   replicate = period)
  # redefine the period based on sex and what period was previously defined as
  tibble <- tibble %>% mutate(period = case_when(
    period < 5 & sex=='M' | sex=='F' ~ 2,
    period > 4 & sex=='M' | sex=='F' ~ 3,
    .default = 1
  ))
  # expand treatment names
  tibble <- tibble %>% mutate(treatment = case_when(
    treatment=='C' ~ "control",
    treatment=='F' ~ "fluctuating",
  ))
  # expand lifestage names
  tibble <- tibble %>% mutate(lifestage = case_when (
    lifestage=='A' ~ "adult",
    lifestage=='L' ~ "larvae"
  ))
  # expand timepoint names
  tibble <- tibble %>% mutate(timepoint = case_when(
    timepoint=='H' ~ "high",
    timepoint=='L' ~ "low"
  ))
  # expand sexes and include unknown sex (aka larvae)
  tibble <- tibble %>% mutate(sex = case_when(
    sex=='F' ~ "female",
    sex=='M' ~ "male",
    sex=='' ~ "unknown"
  ))
  return(tibble)
}

# Read in the raw counts
counts <- read.csv('data/gene_count_matrix.csv')
# Keep all columns besides the gene names
counts <- counts[-1]
# Extract the column names - these are the sample names
columns <- colnames(counts)
# create the meta data
meta <- meta_info_from_labels(columns)
# convert to dataframe
meta_df <- as.data.frame(meta)
# write it to a csv file
write.csv(meta_df, "data/sample_info.csv", row.names=FALSE)