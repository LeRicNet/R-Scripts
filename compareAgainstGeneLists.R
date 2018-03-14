# Compare Gene Significance Against MSigDB Lists
# Author: Eric Prince
# Date: 2018-03-14
#
# DESCRIPTION:
# This script takes in a .csv file exported from Loupe Cell Browser
# containing gene expression by cluster id.  Genes are filtered out to be
# p < pValThresh, and then counts the number of genes in each cluster found
# in a given gene list.

library(tidyverse)
library(magrittr)
library(GSA)

compareAgainstGeneLists <- function(geneset.gmt, dataset.path, 
                                    pValThresh = 0.1,
                                    value) {
  # Inputs:
  #   1.  geneset.gmt
  #         (str) - path to .gmt file [obtained from MSigDB]
  #
  #   2.  dataset.path
  #         (str) - path to .csv exported from Loupe Cell Browser with
  #         the format:
  #           EnsemblID   GeneName `Cluster 1 Avera… `Cluster 1 Log2 Fol… `Cluster 1 P-Val…
  #            <chr>       <chr>                <dbl>                <dbl>             <dbl>
  #          1 ENSG000002… TRAC                  2.16                 3.60     0.00000000602
  #          2 ENSG000001… CD3D                  1.23                 3.86     0.0000000217
  #          3 ENSG000000… IL32                  4.77                 3.45     0.0000000287
  #
  #   3.  pValThresh
  #         (num) - numeric value to filter for p-value
  #   4.  value
  #         (char) - 'raw' or 'pct'.
  #             'raw' returns pure count values of genes in a geneset for each cluster.
  #             'pct' returns the count divided by the number of genes in each respective geneset.
  
  print(paste(">> Importing:", dataset.path))
  dataset <- read_csv(dataset.path)
  df <- df %>%
    select(-contains('Average'))
  
  print(paste(">> Dataset contains:", nrow(df), "genes."))
  
  print(paste(">> GeneSet Location:", geneset.gmt))
  geneset <- GSA.read.gmt(geneset.gmt)
  num_rows <- length(geneset$geneset.names)
  num_cols <- (ncol(df) - 2) / 2
  print(paste(">> K =", num_cols))
  mat <- matrix(0, nrow = num_rows, ncol = num_cols+1)
  
  count_genes <- function(dataset, geneset = geneset) {
    count <- nrow(dataset[dataset$GeneName %in% geneset,])
    return(count)
  }
  
  cluster_data <- list()
  idx_a <- grep('Fold', names(df))
  idx_b <- idx_a + 1
  nms <- grep('Fold', names(df), value = TRUE)
  for (i in c(1:length(nms))) {
    nm <- paste0(sapply(strsplit(nms[i], split = " "), '[', c(1:2)), collapse = "")
    target_seq <- c(1:2, idx_a[i]:idx_b[i])
    tmp_df <- df[,target_seq]
    names(tmp_df) <- c("EnsemblID", "GeneName","Log2FoldChange", "pVal")
    tmp_df <- tmp_df %>%
      filter(pVal < pValThresh)
    print(paste(nm, "contains:", nrow(tmp_df), "genes with p-value <", pValThresh))
    cluster_data[[paste0(nm)]] <- tmp_df
  }
  
  if (value == 'raw') {
    dataset <- cluster_data
    for (i in c(1:num_rows)) {
      gene_set_name <- unlist(geneset$geneset.names[i])
      gene_set_genes <- unlist(geneset$genesets[i])
      mat[i,1] <- gene_set_name
      rng <- c(1:num_cols)
      for (j in rng) {
        k = j + 1
        mat[i,k] <- count_genes(dataset = dataset[[j]], geneset = gene_set_genes)
      }
    }
  } else if (value == 'pct') {
    dataset <- cluster_data
    for (i in c(1:num_rows)) {
      gene_set_name <- unlist(geneset$geneset.names[i])
      gene_set_genes <- unlist(geneset$genesets[i])
      total <- length(gene_set_genes)
      mat[i,1] <- gene_set_name
      rng <- c(1:num_cols)
      for (j in rng) {
        k = j + 1
        mat[i,k] <- count_genes(dataset = dataset[[j]], geneset = gene_set_genes) / total
      }
    }
  }
  
  
  cluster_ids <- c("Gene List", paste0("Cluster", rng))
  out <- as_tibble(mat)
  names(out)[(ncol(out) - length(cluster_ids)):ncol(out)] <- cluster_ids
  
  # Convert values to numeric
  for (i in c(2:ncol(out))) {
    out[,i] <- as.numeric(unlist(out[,i]))
  }
  
  return(out)
}

