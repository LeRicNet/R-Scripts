# Simple Script to Combine FPKM Files into one DataFrame
# Eric Prince
# 2018-04-11
#
# Input:
#   root_dir = top directory for all folders containing FPKM files with the suffix ".sorted".
# Output:
#   FPKM = DataFrame with the first three columns being Ensemble GeneID, the Gene Short Name, 
#           and the locus. With the subsequent columns for each sample.  Sample column names
#           are appended with the sample ID.
#
# Usage:
#   source(/this/file)
#   FPKM <- build.FPKM.file(root_dir)
#

library(tidyverse)
library(magrittr)

build.FPKM.file <- function(root_dir) {
  file_list <- list.files(root_dir, recursive = T, full.names = T, include.dirs = T)
  fpkm_files <- unique(file_list[grepl("sorted$", file_list)])
  
  for (f in 1:length(fpkm_files)) {
    split_name <- unlist(strsplit(fpkm_files[f], "[/]"))
    sample_name <- split_name[length(split_name)]
    sample_name <- sapply(strsplit(sample_name, split = "_"), '[', 1)
    # Remove extra characters after id
    sample_name <- sapply(strsplit(sample_name, split = "+[a-z]"), '[', 1)
    # Exchange hyphen for underscore
    sample_name <- gsub(pattern = "[-]", replacement = ".", sample_name)
    print(paste0("Gathering: ", sample_name))
    
    # Import DataFrame
    fpkm_data <- read.table(fpkm_files[f], header = T)
    fpkm_data <- fpkm_data[,c('gene_id', 'gene_short_name', 'locus', 'FPKM')]
    names(fpkm_data)[length(fpkm_data)] <- paste0("FPKM_", sample_name)
    
    if (!exists("FPKM")) {
      FPKM <- fpkm_data
    } else {
      FPKM <- full_join(FPKM, fpkm_data, by = c('gene_id', 'gene_short_name', 'locus'))
    }
    
    
  }
  return(FPKM)
}


  

