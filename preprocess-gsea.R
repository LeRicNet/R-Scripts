# Preprocessing Script For GSEA ================================================
# Author: Eric Prince
# Date: 2018-01-27
#
# Description -
# Gene Set Enrichment Analysis (GSEA) examines differences between defined gene
# sets of two classes (phenotypes) of biological samples.  The input for GSEA is
# four files:
#       1. Expression Data
#             Data acquired from NGS platforms, typically microarray.  Must be
#             of the format where genes and cells are represented as numeric
#             values in a matrix format
#       2. Phenotype Data
#             A .cls file describing the classes and identifying samples to
#             classes.
#       3. Gene Set(s) -
#             This is a simple list of the gene set to be compared for enrich-
#             ment.  This file can contain one gene set, or any number of gene
#             sets.  There are many gene sets available at MSigDB, but they may
#             also be user defined.
#       4. Microarray Chip Annotation
#             A description .chip file for the genes probed in the expression
#             data.
#
# Initialization ---------------------------------------------------------------
library(tidyverse)


# Define Functions -------------------------------------------------------------
GeneratePhenotypeFile <- function(expression_data, class_titles, class_n,
                                  out_file_cls) {
  
  num_samples <- length(names(expression_data)[3:ncol(expression_data)])
  num_classes <- length(class_titles)
  
  write(paste0(num_samples, " ", num_classes, " 1"),
        file = out_file_cls)
  
  write(c(paste0("#"), paste0(class_titles)),
        file = out_file_cls,
        ncolumns = length(class_titles) + 1,
        sep = " ",
        append = TRUE)
  
  pheno_seq <- c()
  for (i in 1:length(class_n)) {
    k = i-1
    y = rep(k, as.numeric(class_n[i]))
    pheno_seq <- c(pheno_seq, y)
  }
  write(paste0(pheno_seq),
        file = out_file_cls,
        ncolumns = length(pheno_seq) + 1,
        sep = " ",
        append = TRUE)
  
}

GenerateChipFile <- function(gene_probes, out_file_chip, 
                             input_gene_format = "standard") {
  if (input_gene_format == "standard") {
    gene_symbols <- gene_probes
  } else if (input_gene_format == "ensembl") {
    require('biomaRt')
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    gene_symbols <- getBM(filters= "ensembl_gene_id",
                          attributes= c("external_gene_name", "description"),
                          values=gene_probes,
                          mart= mart)
  }
  chip <- data_frame('Probe Set ID' = paste0('A', seq(1, length(gene_probes))),
                     'Gene Symbol' = gene_symbols$external_gene_name,
                     'Gene Title' = gene_symbols$description)
  write_tsv(chip, out_file_chip)
}



