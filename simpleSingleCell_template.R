#==============================================================================================
#
#                      Analyzing single-cell RNA-seq data containing UMI counts
#
#==============================================================================================
setwd('/Users/mijiarui/Nature_Biotechnology_Paper/simpleSingleCell')
library(simpleSingleCell)

#################################### Setting up the data ####################################

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

## We add gene-based annotation identifying rows that correspond to each class of features. We also determine the Ensembl identifier for each row.
## Specifying the nature of each row.
nrows <- c(nrow(endo.data$counts), nrow(mito.data$counts), nrow(spike.data$counts))
is.spike <- rep(c(FALSE, FALSE, TRUE), nrows); length(is.spike)
is.mito <- rep(c(FALSE, TRUE, FALSE), nrows); length(is.mito)
isSpike(sce, "Spike") <- is.spike

## Adding Ensembl IDs.
library(org.Mm.eg.db)
ensembl <- mapIds(org.Mm.eg.db, keys=rownames(sce), keytype="SYMBOL", column="ENSEMBL")
rowData(sce)$ENSEMBL <- ensembl
sce 



#################################### Quality control on the cells ####################################

# The original authors of the study have already removed low-quality cells prior to data publication. 
# Nonetheless, we compute some quality control metrics with scater (McCarthy et al. 2017) to check whether the 
# remaining cells are satisfactory.
library(scater)
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=is.mito)) 
sce # 查看sce，发现最大的区别在于colData从原来的10列变成现在的65列


## We examine the distribution of the QC metrics across all cells (Figure 1). The library sizes here are at 
## least one order of magnitude lower than observed in the 416B dataset. This is consistent with the use of UMI 
## counts rather than read counts, as each transcript molecule can only produce one UMI count but can yield many 
## reads after fragmentation. In addition, the spike-in proportions are more variable than observed in the 416B 
## dataset. This may reflect a greater variability in the total amount of endogenous RNA per cell when many cell 
## types are present.
names(colData(sce))
par(mfrow=c(2,2), mar=c(5.1, 4.1, 0.1, 0.1))
hist(sce$total_counts/1e3, xlab="Library sizes (thousands)", main="",  # 使用$符号可以从SingleCellExperiment里面提取所有colData的信息，相当于对细胞进行质控
     breaks=20, col="grey80", ylab="Number of cells")
hist(sce$total_features, xlab="Number of expressed genes", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
hist(sce$pct_counts_Spike, xlab="ERCC proportion (%)",
     ylab="Number of cells", breaks=20, main="", col="grey80")
hist(sce$pct_counts_Mt, xlab="Mitochondrial proportion (%)", 
     ylab="Number of cells", breaks=20, main="", col="grey80")


## We remove small outliers for the library size and the number of expressed features, and large outliers for 
## the spike-in proportions. Again, the presence of spike-in transcripts means that we do not have to use the 
## mitochondrial proportions.
## 这边nmads的参数含义是“A numeric scalar, specifying the minimum number of MADs away from median required for 
## a value to be called an outlier.” MAD的意思是中位绝对偏差；我们滤除reads counts数，检测的gene数偏少，和那些
## spike-in比例偏高的细胞
libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce$total_features, nmads=3, type="lower", log=TRUE)
spike.drop <- isOutlier(sce$pct_counts_Spike, nmads=3, type="higher")


## Removal of low-quality cells is then performed by combining the filters for all of the metrics. The majority 
## of cells are retained, which suggests that the original quality control procedures were generally adequate.
sce <- sce[,!(libsize.drop | feature.drop | spike.drop)]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), 
           BySpike=sum(spike.drop), Remaining=ncol(sce))
### We could improve our cell filtering procedure further by setting batch in isOutlier to one or more known 
### factors, e.g., mouse/plate of origin. As previously mentioned, this would avoid inflation of the MAD and 
### improve power to remove low-quality cells. However, for simplicity, we will not do this as sufficient 
### quality control has already been performed.

#################################### Cell cycle classifications ####################################
library(scran)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran")) # 目前有人和小鼠的数据，暂无zebrafish的数据
assignments <- cyclone(sce, mm.pairs, gene.names=rowData(sce)$ENSEMBL) # 这一步很耗时
table(assignments$phase)
plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16)


