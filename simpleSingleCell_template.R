#==============================================================================================
#
#                      Analyzing single-cell RNA-seq data containing UMI counts
#
#==============================================================================================
setwd('/Users/mijiarui/Nature_Biotechnology_Paper/simpleSingleCell')
library(simpleSingleCell)

readFormat <- function(infile) { 
  # First column is empty.
  metadata <- read.delim(infile, stringsAsFactors=FALSE, header=FALSE, nrow=10)[,-1] 
  rownames(metadata) <- metadata[,1]
  metadata <- metadata[,-1]
  metadata <- as.data.frame(t(metadata))
  
  # First column after row names is some useless filler.
  counts <- read.delim(infile, stringsAsFactors=FALSE, 
                       header=FALSE, row.names=1, skip=11)[,-1] 
  counts <- as.matrix(counts)
  return(list(metadata=metadata, counts=counts))
}

## Load in data, read in counts for the endogenous genes, ERCC spike-in transcripts and mitochondrial genes
endo.data <- readFormat("expression_mRNA_17-Aug-2014.txt")
spike.data <- readFormat("expression_spikes_17-Aug-2014.txt")
mito.data <- readFormat("expression_mito_17-Aug-2014.txt")


## We also need to rearrange the columns for the mitochondrial data, as the order is not consistent with the 
## other files.
m <- match(endo.data$metadata$cell_id, mito.data$metadata$cell_id)
mito.data$metadata <- mito.data$metadata[m,]
mito.data$counts <- mito.data$counts[,m]

## In this particular dataset, some genes are represented by multiple rows corresponding to alternative genomic 
## locations. We sum the counts for all rows corresponding to a single gene for ease of interpretation.
raw.names <- sub("_loc[0-9]+$", "", rownames(endo.data$counts))
new.counts <- rowsum(endo.data$counts, group=raw.names, reorder=FALSE)
endo.data$counts <- new.counts

## The counts are then combined into a single matrix for constructing a  SingleCellExperiment object. For 
## convenience, metadata for all cells are stored in the same object for later access.
library(SingleCellExperiment)
all.counts <- rbind(endo.data$counts, mito.data$counts, spike.data$counts)
sce <- SingleCellExperiment(list(counts=all.counts), colData=endo.data$metadata)
dim(sce)
