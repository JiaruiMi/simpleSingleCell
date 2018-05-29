#==============================================================================================
#
#             0.Analyzing single-cell RNA sequencing data from droplet-based protocols
#
#==============================================================================================
######################################### Overview ##########################################
# 可以在10X的官网“https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k”
# 下载数据

####################################### Setting up the data #######################################
# Reading in a sparse matrix
## We load in the raw count matrix using the read10xCounts() function from the DropletUtils package. 
## This will create a SingleCellExperiment object where each column corresponds to a cell barcode.
untar("pbmc4k_raw_gene_bc_matrices.tar.gz", exdir="pbmc4k")

## 我们看到10X的数据解压以后是一个文件夹，进入文件夹有三个文件分别是“barcodes.tsv", "genes.tsc", "matrix.mtx"
## 一般10X的数据可以使用Seurat(适合大量细胞)，当然在细胞数有限的情况下考虑使用SingleCellExperiment对象
## 和Scater包也是可以一样处理的，在这里我们使用scater
library(DropletUtils)
fname <- "pbmc4k/raw_gene_bc_matrices/GRCh38"
sce <- read10xCounts(fname, col.names=TRUE)
sce # assays里面只有counts一个slots

# Here, each count represents the number of unique molecular identifiers (UMIs) assigned to a gene 
# for a cell barcode. Note that the counts are loaded as a sparse matrix object - specifically, a 
# dgCMatrix instance from the Matrix package. This avoids allocating memory to hold zero counts, 
# which is highly memory-efficient for low-coverage scRNA-seq data.
# 一般10X的方法是和UMI序列一起配合使用的。
class(counts(sce))





#==============================================================================================
#
#                    1.Analyzing single-cell RNA-seq data containing UMI counts
#
#==============================================================================================
## updated in 2018-05-25, R version: 3.5.0
setwd('/Users/mijiarui/Nature_Biotechnology_Paper/simpleSingleCell')
library(simpleSingleCell)

######################################### Overview ##########################################

## In this workflow, we examine a heterogeneous dataset from a study of cell types in the mouse 
## brain (Zeisel et al. 2015). This contains approximately 3000 cells of varying types such as 
## oligodendrocytes, microglia and neurons. Individual cells were isolated using the Fluidigm C1 
## microfluidics system (Pollen et al. 2014) and library preparation was performed on each cell 
## using a UMI-based protocol. After sequencing, expression was quantified by counting the number 
## of UMIs mapped to each gene.

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
par(mfrow = c(1,1), mar = c(5,5,3,2))
plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16)

#################################### Examining gene-level metrics ####################################
# Figure 3 shows the most highly expressed genes across the cell population in the brain dataset. This is 
# mostly occupied by spike-in transcripts, reflecting the use of spike-in concentrations that span the entire 
# range of expression. There are also a number of constitutively expressed genes, as expected.
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotQC(sce, type = "highest-expression", n=50) + fontsize # 这步非常耗时，导出图片也非常耗时

# Gene abundance is quantified by computing the average count across all cells (Figure 4). As previously 
# mentioned, the UMI count is generally lower than the read count.
ave.counts <- calcAverage(sce, use_size_factors=FALSE)
hist(log10(ave.counts), breaks=100, main="Histogram of log-average counts for all genes in the brain dataset", 
     col="grey", xlab=expression(Log[10]~"average count"))

# We save the average counts into the SingleCellExperiment object for later use. We also remove genes that have 
# average counts of zero, as this means that they are not expressed in any cell.
rowData(sce)$ave.count <- ave.counts
to.keep <- ave.counts > 0
sce <- sce[to.keep,]
summary(to.keep)


############################## Normalization of cell-specific biases ##############################
# For endogenous genes, normalization is performed using the computeSumFactors function as previously described. 
# Here, we cluster similar cells together and normalize the cells in each cluster using the deconvolution method 
# (Lun, Bach, and Marioni 2016). This improves normalization accuracy by reducing the number of DE genes between 
# cells in the same cluster. Scaling is then performed to ensure that size factors of cells in different clusters 
# are comparable.

# We use a average count threshold of 0.1 to define high-abundance genes to use during normalization. 
# This is lower than the default threshold of min.mean=1 in  computeSumFactors, reflecting the fact that UMI 
# counts are generally smaller than read counts.
clusters <- quickCluster(sce, min.mean=0.1, method="igraph")
## quickCluster uses distances based on Spearman’s rank correlation for clustering. This ensures that scaling 
## biases in the counts do not affect clustering, but yields very coarse clusters and is not recommended for 
## biological interpretation.
## For large datasets, using method="igraph" in quickCluster will speed up clustering. This uses a graph-based 
## clustering algorithm - see ?buildSNNGraph for more details.

sce <- computeSumFactors(sce, cluster=clusters, min.mean=0.1)
## Only a rough clustering is required to avoid pooling together very different cell types in  computeSumFactors. 
## The function is robust to a moderate level of differential expression between cells in the same cluster.

summary(sizeFactors(sce))


# Compared to the 416B analysis, more scatter is observed around the trend between the total count and size 
# factor for each cell (Figure 5). This is consistent with an increased amount of DE between cells of different 
# types, which compromises the accuracy of library size normalization (Robinson and Oshlack 2010). In contrast, 
# the size factors are estimated based on median ratios and are more robust to the presence of DE between cells.
plot(sizeFactors(sce), sce$total_counts/1e3, log="xy",
     ylab="Library size (thousands)", xlab="Size factor")

# We also compute size factors specific to the spike-in set, as previously described.
sce <- computeSpikeFactors(sce, type="Spike", general.use=FALSE)

# Finally, normalized log-expression values are computed for each endogenous gene or spike-in transcript using 
# the appropriate size factors.
sce <- normalize(sce) # assays增加了logcounts这个slot


############################## Modelling and removing technical noise ##############################
var.fit <- trendVar(sce, parametric=TRUE, loess.args=list(span=0.4))
var.out <- decomposeVar(sce, var.fit)
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
     ylab="Variance of log-expression")
points(var.out$mean[isSpike(sce)], var.out$total[isSpike(sce)], col="red", pch=16)
curve(var.fit$trend(x), col="dodgerblue", add=TRUE, lwd=2)
chosen.genes <- order(var.out$bio, decreasing=TRUE)[1:10]
plotExpression(sce, rownames(var.out)[chosen.genes], 
               alpha=0.05, jitter="jitter") + fontsize
sce <- denoisePCA(sce, technical=var.fit$trend, approximate=TRUE)
ncol(reducedDim(sce, "PCA"))
sce # 增加了reducedDimNames

########################## Data exploration with dimensionality reduction ##########################
# We perform dimensionality reduction on the denoised PCs to check if there is any substructure. 
# Cells separate into clear clusters in the t-SNE plot (Van der Maaten and Hinton 2008) in Figure 8, 
# corresponding to distinct subpopulations. This is consistent with the presence of multiple cell types in the 
# diverse brain population. We increase the perplexity to favour visualization of the overall structure at the 
# expense of local scale.
sce <- runTSNE(sce, use_dimred="PCA", perplexity=50, rand_seed=1000)
tsne1 <- plotTSNE(sce, colour_by="Neurod6") + fontsize
tsne2 <- plotTSNE(sce, colour_by="Mog") + fontsize
multiplot(tsne1, tsne2, cols=2)

# The PCA plot is less effective at separating cells into many different clusters (Figure 9). 
# This is because the first two PCs are driven by strong differences between specific subpopulations, 
# which reduces the resolution of more subtle differences between some of the other subpopulations. 
# Nonetheless, some substructure is still visible.
pca1 <- plotReducedDim(sce, use_dimred="PCA", colour_by="Neurod6") + fontsize
pca2 <- plotReducedDim(sce, use_dimred="PCA", colour_by="Mog") + fontsize
multiplot(pca1, pca2, cols=2)

## For both methods, we colour each cell based on the expression of a particular gene. This is a 
## useful strategy for visualizing changes in expression across the lower-dimensional space. 
## It can also be used to characterise each cluster if the selected genes are known markers for 
## particular cell types. For example, Mog can be used to identify clusters corresponding to 
## oligodendrocytes.

########################## Clustering cells into putative subpopulations ##########################
# The reduced dimension coordinates are used to cluster cells into putative subpopulations. 
# We do so by constructing a shared-nearest-neighbour graph (Xu and Su 2015), in which cells 
# are the nodes and edges are formed between cells that share nearest neighbours. Clusters are 
# then defined as highly connected communities of cells within this graph, using methods from 
# the igraph package. This is more efficient than forming a pairwise distance matrix for 
# hierarchical clustering of large numbers of cells.
# 共享最近邻聚类算法: Decreasing the number of neighbours k in buildSNNGraph will reduce the 
# connectivity of the graph. This will generally result in the formation of smaller clusters 
# (Xu and Su 2015), which may be desirable if greater resolution is required.
snn.gr <- buildSNNGraph(sce, use.dimred="PCA")
cluster.out <- igraph::cluster_walktrap(snn.gr)
my.clusters <- cluster.out$membership
table(my.clusters)

# The modularity score provides a global measure of clustering performance for community 
# detection methods. Briefly, it compares the number of within-cluster edges to the expected 
# number under a null model of random edges. A high modularity score (approaching the maximum of 
# 1) indicates that the detected clusters are enriched for internal edges, with relatively few 
# edges between clusters.

# Notice that we do not run library(igraph), but instead use igraph:: to extract methods from 
# the package. This is because igraph contains a normalize method that will override its 
# counterpart from scater, resulting in some unusual bugs.
igraph::modularity(cluster.out)


# We further investigate the clusters by examining the total weight of edges for each pair of 
# clusters. For each pair, the observed total weight is compared to what is expected under a 
# null model, similar to the modularity calculation. Most clusters contain more internal links 
# than expected (Figure 10), while links between clusters are fewer than expected. This indicates 
# that we successfully clustered cells into highly-connected communities.
mod.out <- clusterModularity(snn.gr, my.clusters, get.values=TRUE)
ratio <- log10(mod.out$observed/mod.out$expected + 1)
library(pheatmap)
pheatmap(ratio, cluster_rows=FALSE, cluster_cols=FALSE, 
         color=colorRampPalette(c("white", "blue"))(100))

# We visualize the cluster assignments for all cells on the t-SNE plot in Figure 11. Adjacent 
# cells are generally assigned to the same cluster, indicating that the clustering procedure was 
# applied correctly.
sce$cluster <- factor(my.clusters)
plotTSNE(sce, colour_by="cluster") + fontsize


# An alternative approach is to use graph-based visualizations such as force-directed layouts 
# (Figure 12). These are appealing as they directly represent the relationships used during 
# clustering. However, convergence tends to be slow for large graphs, so some tinkering with 
# niter= may be required to ensure that the results are stable.
set.seed(2000)
reducedDim(sce, "force") <- igraph::layout_with_fr(snn.gr, niter=5000)
plotReducedDim(sce, colour_by="cluster", use_dimred="force")


########################## Detecting marker genes between subpopulations ##########################
# We use the findMarkers function with direction="up" to identify upregulated marker genes for 
# each cluster. As previously mentioned, we focus on upregulated genes as these can quickly 
# provide positive identification of cell type in a heterogeneous population. We examine the 
# table for cluster 1, in which log-fold changes are reported between cluster 1 and every other 
# cluster. The same output is provided for each cluster in order to identify genes that 
# discriminate between clusters.
markers <- findMarkers(sce, my.clusters, direction="up")
marker.set <- markers[["1"]]
head(marker.set[,1:8], 10) # only first 8 columns, for brevity

# We save the list of candidate marker genes for further examination, using compression to 
# reduce the file size.
gzout <- gzfile("brain_marker_1.tsv.gz", open="wb")
write.table(marker.set, file=gzout, sep="\t", quote=FALSE, row.names=FALSE)
close(gzout)

# Figure 13 indicates that most of the top markers are strongly DE in cells of cluster 1 compared 
# to some or all of the other clusters. We can use these markers to identify cells from cluster 1
# in validation studies with an independent population of cells. A quick look at the markers suggest 
# that cluster 1 represents interneurons based on expression of Gad1 and Slc6a1 (Zeng et al. 2012), 
# differing from closely related cells in cluster 10 by virtue of high Synpr expression.
top.markers <- rownames(marker.set)[marker.set$Top <= 10]
plotHeatmap(sce, features=top.markers, columns=order(my.clusters),
            colour_columns_by="cluster", cluster_cols=FALSE, 
            center=TRUE, symmetric=TRUE, zlim=c(-5, 5))


################################## Concluding remarks ##################################
# Having completed the basic analysis, we save the SingleCellExperiment object with its associated 
# data to file. This is especially important here as the brain dataset is quite large. If further 
# analyses are to be performed, it would be inconvenient to have to repeat all of the pre-processing 
# steps described above.
saveRDS(file="brain_data.rds", sce)


#==============================================================================================
#
#                   2. Analyzing single-cell RNA-seq data containing read counts
#
#==============================================================================================

############################################ Overview ############################################

# In this workflow, we use a relatively simple dataset (Lun et al. 2017) to introduce most of the 
# concepts of scRNA-seq data analysis. This dataset contains two plates of 416B cells (an immortalized 
# mouse myeloid progenitor cell line), processed using the Smart-seq2 protocol (Picelli et al. 2014). 
# A constant amount of spike-in RNA from the External RNA Controls Consortium (ERCC) was also added to 
# each cell’s lysate prior to library preparation. High-throughput sequencing was performed and the 
# expression of each gene was quantified by counting the total number of reads mapped to its exonic 
# regions. Similarly, the quantity of each spike-in transcript was measured by counting the number of 
# reads mapped to the spike-in reference sequences. 

##################################### Setting up the data #####################################
# One matrix was generated for each plate of cells used in the study.
# unzip("E-MTAB-5522.processed.1.zip"). Unzip once. 这里我们关注标明“Calero”的数据集
# Reading in the count tables for each of the two plates.
plate1 <- read.delim("counts_Calero_20160113.tsv", 
                     header=TRUE, row.names=1, check.names=FALSE)
plate2 <- read.delim("counts_Calero_20160325.tsv", 
                     header=TRUE, row.names=1, check.names=FALSE)
head(plate1[1:5,1:5])
gene.lengths <- plate1$Length # First column is the gene length.
plate1 <- as.matrix(plate1[,-1]) # Discarding gene length (as it is not a cell).
plate2 <- as.matrix(plate2[,-1])
rbind(Plate1=dim(plate1), Plate2=dim(plate2))

# We combine the two matrices into a single object for further processing. This is done after 
# verifying that the genes are in the same order between the two matrices.
stopifnot(identical(rownames(plate1), rownames(plate2)))
all.counts <- cbind(plate1, plate2)

# For convenience, the count matrix is stored in a SingleCellExperiment object from the SingleCellExperiment 
# package. This allows different types of row- and column-level metadata to be stored alongside the 
# counts for synchronized manipulation throughout the workflow.
library(SingleCellExperiment)
sce <- SingleCellExperiment(list(counts=all.counts))
rowData(sce)$GeneLength <- gene.lengths
sce

# We identify the rows corresponding to ERCC spike-in transcripts from the row names. We store this 
# information in the SingleCellExperiment object for future use. This is necessary as spike-ins 
# require special treatment in downstream steps such as normalization.
isSpike(sce, "ERCC") <- grepl("^ERCC", rownames(sce))
summary(isSpike(sce, "ERCC"))

# This dataset is slightly unusual in that it contains information from another set of spike-in 
# transcripts, the Spike-In RNA Variants (SIRV) set. For simplicity, we will only use the ERCC 
# spike-ins in this analysis. Thus, we must remove the rows corresponding to the SIRV transcripts 
# prior to further analysis, which can be done simply by subsetting the  SingleCellExperiment object.
is.sirv <- grepl("^SIRV", rownames(sce))
sce <- sce[!is.sirv,] 
summary(is.sirv)


metadata <- read.delim("E-MTAB-5522.sdrf.txt", check.names=FALSE, header=TRUE)
m <- match(colnames(sce), metadata[["Source Name"]]) # Enforcing identical order.
stopifnot(all(!is.na(m))) # Checking that nothing's missing.
metadata <- metadata[m,]
head(colnames(metadata))


colData(sce)$Plate <- factor(metadata[["Factor Value[block]"]])
pheno <- metadata[["Factor Value[phenotype]"]]
levels(pheno) <- c("induced", "control")
colData(sce)$Oncogene <- pheno
table(colData(sce)$Oncogene, colData(sce)$Plate)





library(org.Mm.eg.db)
symb <- mapIds(org.Mm.eg.db, keys=rownames(sce), keytype="ENSEMBL", column="SYMBOL")
rowData(sce)$ENSEMBL <- rownames(sce)
rowData(sce)$SYMBOL <- symb
head(rowData(sce))

new.names <- rowData(sce)$SYMBOL
missing.name <- is.na(new.names)
new.names[missing.name] <- rowData(sce)$ENSEMBL[missing.name]
dup.name <- new.names %in% new.names[duplicated(new.names)]
new.names[dup.name] <- paste0(new.names, "_", rowData(sce)$ENSEMBL)[dup.name]
rownames(sce) <- new.names
head(rownames(sce))

library(TxDb.Mmusculus.UCSC.mm10.ensGene)
location <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene, keys=rowData(sce)$ENSEMBL, 
                   column="CDSCHROM", keytype="GENEID")
rowData(sce)$CHR <- location
summary(location=="chrM")


##################################### Quality control on the cells #####################################
library(scater)
mito <- which(rowData(sce)$CHR=="chrM")
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=mito))
head(colnames(colData(sce)), 10)


sce$PlateOnco <- paste0(sce$Oncogene, ".", sce$Plate)
multiplot(
  plotColData(sce, y="total_counts", x="PlateOnco"),
  plotColData(sce, y="total_features", x="PlateOnco"),
  plotColData(sce, y="pct_counts_ERCC", x="PlateOnco"),
  plotColData(sce, y="pct_counts_Mt", x="PlateOnco"),
  cols=2)


par(mfrow=c(1,3))
plot(sce$total_features, sce$total_counts/1e6, xlab="Number of expressed genes",
     ylab="Library size (millions)")
plot(sce$total_features, sce$pct_counts_ERCC, xlab="Number of expressed genes",
     ylab="ERCC proportion (%)")
plot(sce$total_features, sce$pct_counts_Mt, xlab="Number of expressed genes",
     ylab="Mitochondrial proportion (%)")



libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", 
                          log=TRUE, batch=sce$PlateOnco)
feature.drop <- isOutlier(sce$total_features, nmads=3, type="lower", 
                          log=TRUE, batch=sce$PlateOnco)


spike.drop <- isOutlier(sce$pct_counts_ERCC, nmads=3, type="higher",
                        batch=sce$PlateOnco)


keep <- !(libsize.drop | feature.drop | spike.drop)
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
           BySpike=sum(spike.drop), Remaining=sum(keep))


sce$PassQC <- keep
saveRDS(sce, file="416B_preQC.rds")
sce <- sce[,keep]
dim(sce)



################################### Classification of cell cycle phase ###################################
set.seed(100)
library(scran)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", 
                                package="scran"))
assignments <- cyclone(sce, mm.pairs, gene.names=rowData(sce)$ENSEMBL)

par(mfrow = c(1,1))
plot(assignments$score$G1, assignments$score$G2M, 
     xlab="G1 score", ylab="G2/M score", pch=16)

sce$phases <- assignments$phases
table(sce$phases)


################################### Examining gene-level expression metrics ###################################
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotQC(sce, type = "highest-expression", n=50) + fontsize

ave.counts <- calcAverage(sce, use_size_factors=FALSE)
hist(log10(ave.counts), breaks=100, main="", col="grey80", 
     xlab=expression(Log[10]~"average count"))

demo.keep <- ave.counts >= 1
filtered.sce <- sce[demo.keep,]
summary(demo.keep)


num.cells <- nexprs(sce, byrow=TRUE)
smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells", 
              xlab=expression(Log[10]~"average count"))


to.keep <- num.cells > 0
sce <- sce[to.keep,]
summary(to.keep)


############################## Normalization of cell-specific biases ##############################
sce <- computeSumFactors(sce)
summary(sizeFactors(sce))

plot(sce$total_counts/1e6, sizeFactors(sce), log="xy",
     xlab="Library size (millions)", ylab="Size factor",
     col=c("red", "black")[sce$Oncogene], pch=16)
legend("bottomright", col=c("red", "black"), pch=16, cex=1.2,
       legend=levels(sce$Oncogene))

sce <- computeSpikeFactors(sce, type="ERCC", general.use=FALSE)

sce <- normalize(sce)


########################### Modelling the technical noise in gene expression ###########################
var.fit <- trendVar(sce, parametric=TRUE, block=sce$Plate,
                    loess.args=list(span=0.3))
var.out <- decomposeVar(sce, var.fit)
head(var.out)


plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
     ylab="Variance of log-expression")
curve(var.fit$trend(x), col="dodgerblue", lwd=2, add=TRUE)
cur.spike <- isSpike(sce)
points(var.out$mean[cur.spike], var.out$total[cur.spike], col="red", pch=16)

chosen.genes <- order(var.out$bio, decreasing=TRUE)[1:10]
plotExpression(sce, features=rownames(var.out)[chosen.genes]) + fontsize


################################ Removing the batch effect ################################
library(limma)
assay(sce, "corrected") <- removeBatchEffect(logcounts(sce), 
                                             design=model.matrix(~sce$Oncogene), batch=sce$Plate)
assayNames(sce)


############################## Denoising expression values using PCA ##############################
sce <- denoisePCA(sce, technical=var.fit$trend, assay.type="corrected")
dim(reducedDim(sce, "PCA")) 


sce2 <- denoisePCA(sce, technical=var.fit$trend, 
                   assay.type="corrected", value="lowrank") 
assayNames(sce2)



######################## Data exploration with dimensionality reduction #######################
plotReducedDim(sce, use_dimred="PCA", ncomponents=3, 
               colour_by="Oncogene") + fontsize

plotReducedDim(sce, use_dimred="PCA", ncomponents=3, 
               colour_by="Plate") + fontsize



run_args <- list(rand_seed=100, use_dimred="PCA")
out5 <- plotTSNE(sce, run_args=c(run_args, perplexity=5),
                 colour_by="Oncogene") + fontsize + ggtitle("Perplexity = 5")

out10 <- plotTSNE(sce, run_args=c(run_args, perplexity=10),
                  colour_by="Oncogene") + fontsize + ggtitle("Perplexity = 10")

out20 <- plotTSNE(sce, run_args=c(run_args, perplexity=20),
                  colour_by="Oncogene") + fontsize + ggtitle("Perplexity = 20")

multiplot(out5, out10, out20, cols=3)


######################## Clustering cells into putative subpopulations ########################
pcs <- reducedDim(sce, "PCA")
my.dist <- dist(pcs)
my.tree <- hclust(my.dist, method="ward.D2")

library(dynamicTreeCut)
my.clusters <- unname(cutreeDynamic(my.tree, distM=as.matrix(my.dist), 
                                    minClusterSize=10, verbose=0))

table(my.clusters, sce$Plate)


sce$cluster <- factor(my.clusters)
plotTSNE(sce, run_args=list(use_dimred="PCA", perplexity=20, 
                            rand_seed=200), colour_by="cluster") + fontsize


library(cluster)
clust.col <- scater:::.get_palette("tableau10medium") # hidden scater colours
sil <- silhouette(my.clusters, dist = my.dist)
sil.cols <- clust.col[ifelse(sil[,3] > 0, sil[,1], sil[,2])]
sil.cols <- sil.cols[order(-sil[,1], sil[,3])]
plot(sil, main = paste(length(unique(my.clusters)), "clusters"), 
     border=sil.cols, col=sil.cols, do.col.sort=FALSE) 


markers <- findMarkers(sce, my.clusters, block=sce$Plate)


marker.set <- markers[["1"]]
head(marker.set, 10)


write.table(marker.set, file="416B_marker_1.tsv", sep="\t", 
            quote=FALSE, row.names=FALSE)



top.markers <- rownames(marker.set)[marker.set$Top <= 10]
plotHeatmap(sce, features=top.markers, columns=order(sce$cluster), 
            colour_columns_by=c("cluster", "Plate", "Oncogene"),
            cluster_cols=FALSE, center=TRUE, symmetric=TRUE, zlim=c(-5, 5))
            


#################################### Concluding remarks ####################################
saveRDS(file="416B_data.rds", sce)





#==============================================================================================
#
#                3. Further strategies for analyzing single-cell RNA-seq data
#
#==============================================================================================

######################################### Overview #########################################
# The previous workflows focused on analyzing single-cell RNA-seq data with “standard” procedures. 
# However, a number of alternative parameter settings and strategies can be used at some steps of the 
# workflow. This workflow describes a few of these alternative settings as well as the rationale 
# behind choosing them instead of the defaults.


##################################### Quality control on cells #####################################


# Checking for discarded cell types
## We can diagnose loss of distinct cell types during QC by looking for differences in gene expression 
## between the discarded and retained cells (Figure 1). If the discarded pool is enriched for a 
## certain cell type, we should observe increased expression of the corresponding marker genes. No 
## systematic upregulation of genes is apparent in the discarded pool in Figure 1, indicating that 
## the QC step did not inadvertently filter out a cell type in the 416B dataset.
library(SingleCellExperiment)
sce.full.416b <- readRDS("416B_preQC.rds")

library(scater)
suppressWarnings({
  lost <- calcAverage(counts(sce.full.416b)[,!sce.full.416b$PassQC])
  kept <- calcAverage(counts(sce.full.416b)[,sce.full.416b$PassQC])
})
logfc <- log2((lost+1)/(kept+1))
head(sort(logfc, decreasing=TRUE), 20)


# Each point represents a gene, with spike-in and mitochondrial transcripts in red and blue respectively.
plot(lost, kept, xlab="Average count (discarded)", 
     ylab="Average count (retained)", log="xy", pch=16, 
     main = 'Average counts across all discarded and retained cells in the 416B dataset')
is.spike <- isSpike(sce.full.416b)
points(lost[is.spike], kept[is.spike], col="red", pch=16)
is.mito <- rowData(sce.full.416b)$is_feature_control_Mt
points(lost[is.mito], kept[is.mito], col="dodgerblue", pch=16)


# By comparison, a more stringent filter in the PBMC dataset would remove the previously identified 
# platelet population (see the previous workflow). This manifests in Figure 2 as a shift to the 
# bottom-right for a number of genes, including PF4 and PPBP.
sce.pbmc <- readRDS("pbmc_data.rds")
wrong.keep <- sce.pbmc$total_counts >= 1000
suppressWarnings({
  lost <- calcAverage(counts(sce.pbmc)[,!wrong.keep])
  kept <- calcAverage(counts(sce.pbmc)[,wrong.keep])
})
logfc <- log2((lost+1)/(kept+1))
head(sort(logfc, decreasing=TRUE), 20)


plot(lost, kept, xlab="Average count (discarded)", 
     ylab="Average count (retained)", log="xy", pch=16)
platelet <- c("PF4", "PPBP", "SDPR")
points(lost[platelet], kept[platelet], col="orange", pch=16)


# Using PCA-based outliers
# Another strategy is to perform a principal components analysis (PCA) based on the quality metrics 
# for each cell, e.g., the total number of reads, the total number of features and the proportion of 
# mitochondrial or spike-in reads. Outliers on a PCA plot may be indicative of low-quality cells that 
# have aberrant technical properties compared to the (presumed) majority of high-quality cells. This 
# is demonstrated below on a brain cell dataset from Tasic et al. (2016), using functions from the 
# scater package (McCarthy et al. 2017).

# Obtaining the dataset.
library(scRNAseq)
data(allen)

# Setting up the data.
sce.allen <- as(allen, "SingleCellExperiment")
assayNames(sce.allen) <- "counts"
isSpike(sce.allen, "ERCC") <- grep("ERCC", rownames(sce.allen))

# Computing the QC metrics and running PCA.
library(scater)
sce.allen <- calculateQCMetrics(sce.allen)
sce.allen <- runPCA(sce.allen, use_coldata=TRUE, detect_outliers=TRUE)
table(sce.allen$outlier)


################################## Normalizing based on spike-in coverage ##################################


###################################### Detecting highly variable genes #####################################


################################# Advanced modelling of the technical noise ################################

# Loading the saved object.
sce.416B <- readRDS("416B_data.rds") 

# Repeating the trendVar() call.
var.fit <- trendVar(sce.416B, parametric=TRUE, block=sce.416B$Plate,
                    loess.args=list(span=0.3))

matplot(var.fit$means, var.fit$vars, col=c("darkorange", "forestgreen"))

tmp.416B <- sce.416B
tmp.416B$log_size_factor <- log(sizeFactors(sce.416B))
plotColData(tmp.416B, x="Plate", y="log_size_factor")

sce.416B.2 <- normalize(sce.416B, size_factor_grouping=sce.416B$Plate)
comb.out <- multiBlockVar(sce.416B.2, block=sce.416B.2$Plate,
                          trend.args=list(parametric=TRUE, loess.args=list(span=0.4)))

head(comb.out[,1:6])

par(mfrow=c(1,2))
is.spike <- isSpike(sce.416B.2)
for (plate in levels(sce.416B.2$Plate)) {
  cur.out <- comb.out$per.block[[plate]]
  plot(cur.out$mean, cur.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
       ylab="Variance of log-expression", main=plate)
  curve(metadata(cur.out)$trend(x), col="dodgerblue", lwd=2, add=TRUE)
  points(cur.out$mean[is.spike], cur.out$total[is.spike], col="red", pch=16)
}

lfit <- trendVar(sce.416B, design=model.matrix(~sce.416B$Plate))
######################### Identifying correlated gene pairs with Spearman’s rho #########################


######################### Using parallel analysis to choose the number of PCs #########################
set.seed(1000)
npcs <- parallelPCA(sce.416B, assay.type="corrected", 
                    subset.row=comb.out$bio > 0, value="n")
as.integer(npcs)
