library(optparse)
library(Seurat)

# Convert processed Seurat object to csv for import into Inception-v3 classifier.
#
# Eric Prince
# 2018-04-17
#
# This is a command line argument based in R for subsetting a processed Seurat object by a given label.
# For each label in the dataset, the data is subsetted and appended with an ID column with the respective
# label.  The result is a concatanated labeled dataset in csv format. 



###############################################################################################
#                                         M E T H O D S
###############################################################################################

BuildTFTrainInput <- function(file_path, label_id) {  
  # EP: Consider adding transformation.. (i.e. log(x+1))
  object <- readRDS(file_path)
  labels <- levels(as.factor(object@meta.data[,label_id]))
  
  
  for (label in labels) {
    print(paste("Gathering: ", label))  # Pick out each label, and subse the RAW exprs data
    mat <- as.matrix(x = object@raw.data[, WhichCells(object = object,
                                                      ident = label)])
    ids <- rep(label, ncol(mat))
    print(paste("  ", "Number of Cells:", ncol(mat)))
    print(paste("  ", "Number of Genes:", nrow(mat)))
    print("")
    
    mat <- t(rbind(mat, ids))
    
    if (!exists("out")) {
      out <- mat
    } else {
      out <- rbind(out, mat)
    }
  }
  return(out)
}

# Combines the function above, and the write out component
DataPipeline <- function(file_path, label_id, out_path) {
  out <- BuildTFTrainInput(file_path = file_path,
                           label_id = label_id)
  write.table(out, out_path, sep = ",", col.names = FALSE, 
            row.names = FALSE)
}

# Argument Parsing
option_list <- list(
  make_option(c("-d", "--data_dir"), type = "character",
              default = NULL, help = "Dataset path to processed input data",
              metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character",
              default = NULL, help = "Output path for csv",
              metavar = "character"),
  make_option(c("-l", "--label_id"), type = "character",
              default = "ClusterNames_0.6", help = "Column Name for Labels in Processed Data",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if (is.null(opt$data_dir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

###############################################################################################
#                                         M A I N
###############################################################################################
DataPipeline(file_path = opt$data_dir, label_id = opt$label_id,
             out_path = opt$output_dir)

# End (Not Run)