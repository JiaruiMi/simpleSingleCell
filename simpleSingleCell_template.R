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

## Here, each count represents the number of unique molecular identifiers (UMIs) assigned to a gene 
## for a cell barcode. Note that the counts are loaded as a sparse matrix object - specifically, a 
## dgCMatrix instance from the Matrix package. This avoids allocating memory to hold zero counts, 
## which is highly memory-efficient for low-coverage scRNA-seq data.
## 一般10X的方法是和UMI序列一起配合使用的。一般来讲droplet-based methods是没有spike-in RNA的
class(counts(sce))

# Annotating the rows
## We relabel the rows with the gene symbols for easier reading. This is done using the  
## uniquifyFeatureNames() function, which ensures uniqueness in the case of duplicated or missing 
## symbols.
library(scater)
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)
head(rownames(sce))

## We also identify the chromosomal location for each gene. The mitochondrial location is particularly 
## useful for later quality control.
library(EnsDb.Hsapiens.v86)
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce)$ID, 
                   column="SEQNAME", keytype="GENEID")
rowData(sce)$CHR <- location
summary(location=="MT")


#################################### Calling cells from empty droplets ####################################
# An interesting aspect of droplet-based data is that we have no prior knowledge about which droplets 
# (i.e., cell barcodes) actually contain cells, and which are empty. Thus, we need to call cells from 
# empty droplets based on the observed expression profiles. This is not entirely straightforward as 
# empty droplets can contain ambient (i.e., extracellular) RNA that can be captured and sequenced. An 
# examination of the distribution of total counts suggests a fairly sharp transition between barcodes 
# with large and small total counts (Figure 1), probably corresponding to cell-containing and empty 
# droplets respectively.
bcrank <- barcodeRanks(counts(sce))

# Only showing unique points for plotting speed.
uniq <- !duplicated(bcrank$rank)
par(mfrow = c(1,1))
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
     xlab="Rank", ylab="Total UMI count", cex.lab=1.2,
     main = 'Total UMI count for each barcode in the PBMC dataset, plotted against its rank (in decreasing order of total counts)')

abline(h=bcrank$inflection, col="darkgreen", lty=2)
abline(h=bcrank$knee, col="dodgerblue", lty=2)

legend("bottomleft", legend=c("Inflection", "Knee"), 
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)


# We use the emptyDrops() function to test whether the expression profile for each cell barcode is 
# significantly different from the ambient pool (Lun et al. 2018). Any significant deviation indicates 
# that the barcode corresponds to a cell-containing droplet. We call cells at a false discovery rate 
# (FDR) of 1%, meaning that no more than 1% of our called barcodes should be empty droplets on average.
set.seed(100)
e.out <- emptyDrops(counts(sce))
sum(e.out$FDR <= 0.01, na.rm=TRUE)

# emptyDrops() computes Monte Carlo p-values, so it is important to set the random seed to obtain 
# reproducible results. The number of Monte Carlo iterations also determines the lower bound for 
# the _p_values. If any non-significant barcodes are TRUE for Limited, we may need to increase the 
# number of iterations to ensure that they can be detected.
table(Sig=e.out$FDR <= 0.01, Limited=e.out$Limited)
# We then subset our SingleCellExperiment object to retain only the detected cells.
# using which() to automatically remove NAs.
sce <- sce[,which(e.out$FDR <= 0.01)]
## emptyDrops() assumes that cell barcodes with total UMI counts below a certain threshold (default 
## of 100) correspond to empty droplets, and uses them to estimate the ambient expression profile. By 
## definition, these barcodes cannot be cell-containing droplets and are excluded from the hypothesis 
## testing, hence the NAs in the output. 


#################################### Quality control on the cells ####################################
# The previous step only distinguishes cells from empty droplets, but makes no statement about the 
# quality of the cells. It is entirely possible for droplets to contain damaged or dying cells, which 
# need to be removed prior to downstream analysis. We compute some QC metrics using 
# calculateQCMetrics() (McCarthy et al. 2017) and examine their distributions in Figure 2.
sce <- calculateQCMetrics(sce, feature_controls=list(Mito=which(location=="MT")))
par(mfrow=c(1,3))
hist(sce$log10_total_counts, breaks=20, col="grey80",
     xlab="Log-total UMI count")
hist(sce$log10_total_features_by_counts, breaks=20, col="grey80",
     xlab="Log-total number of expressed features")
hist(sce$pct_counts_Mito, breaks=20, col="grey80",
     xlab="Proportion of reads in mitochondrial genes")

# Ideally, we would remove cells with low library sizes or total number of expressed features as 
# described previously. However, this would likely remove cell types with low RNA content, 
# especially in a heterogeneous PBMC population with many different cell types. Thus, we use a more 
# relaxed strategy and only remove cells with large mitochondrial proportions, using it as a proxy 
# for cell damage. (Keep in mind that droplet-based datasets usually do not have spike-in RNA.)
high.mito <- isOutlier(sce$pct_counts_Mito, nmads=3, type="higher")
sce <- sce[,!high.mito]
summary(high.mito)

## The above justification for using a more relaxed filter is largely retrospective. In practice, 
## we may not know a priori the degree of population heterogeneity and whether it manifests in the 
## QC metrics. We recommend performing the analysis first with a stringent QC filter, and then 
## relaxing it based on further diagnostics (see here for an example). This is motivated by the fact 
## that low-quality cells can often yield misleading results, clustering separately from other cells; 
## containing genes that appear to be strongly “upregulated” due to the presence of very small size 
## factors; containing genes that are strongly downregulated due to the loss of RNA upon cell damage; 
## and distorting the characterization of population heterogeneity during variance estimation or PCA.


##################################### Examining gene expression #####################################
# The average expression of each gene is much lower here compared to the previous datasets (Figure 3). 
# This is due to the reduced coverage per cell when thousands of cells are multiplexed together for 
# sequencing.
ave <- calcAverage(sce)
rowData(sce)$AveCount <- ave
par(mfrow = c(1,1))
hist(log10(ave), col="grey80")

# The set of most highly expressed genes is dominated by ribosomal protein and mitochondrial genes (Figure 4), as expected.
plotHighestExprs(sce)


################################## Normalizing for cell-specific biases ################################
# We apply the deconvolution method to compute size factors for all cells (Lun, Bach, and Marioni 2016). 
# We perform some pre-clustering to break up obvious clusters and avoid pooling cells that are very 
# different.
library(scran)
clusters <- quickCluster(sce, method="igraph", min.mean=0.1,
                         irlba.args=list(maxit=1000)) # for convergence.
table(clusters)

sce <- computeSumFactors(sce, min.mean=0.1, cluster=clusters)
summary(sizeFactors(sce))

# The size factors are well correlated against the library sizes (Figure 5), indicating that capture 
# efficiency and sequencing depth are the major biases.
plot(sce$total_counts, sizeFactors(sce), log="xy")

# Finally, we compute normalized log-expresion values. There is no need to call  computeSpikeFactors() 
# here, as there are no spike-in transcripts available.
sce <- normalize(sce)


#################################### Modelling the mean-variance trend ####################################
# The lack of spike-in transcripts complicates the modelling of the technical noise. One option is to 
#assume that most genes do not exhibit strong biological variation, and to fit a trend to the 
# variances of endogenous genes. However, this assumption is generally unreasonable for a heterogeneous 
# population. Instead, we assume that the technical noise is Poisson and create a fitted trend on that 
# basis using the makeTechTrend() function.
new.trend <- makeTechTrend(x=sce)

# We estimate the variances for all genes and compare the trend fits in Figure 6. The Poisson-based 
# trend serves as a lower bound for the variances of the endogenous genes, consistent with non-zero 
# biological components.
# The blue line represents the mean-dependent trend fitted to the variances, while the red line 
# represents the Poisson noise.
fit <- trendVar(sce, use.spikes=FALSE, loess.args=list(span=0.05))
plot(fit$mean, fit$var, pch=16)
curve(fit$trend(x), col="dodgerblue", add=TRUE)
curve(new.trend(x), col="red", add=TRUE)

# We decompose the variance for each gene using the Poisson-based trend, and examine the genes with 
# the highest biological components.
fit0 <- fit
fit$trend <- new.trend
dec <- decomposeVar(fit=fit)
top.dec <- dec[order(dec$bio, decreasing=TRUE),] 
head(top.dec)

# We can plot the genes with the largest biological components, to verify that they are indeed highly 
# variable (Figure 7). Distributions of normalized log-expression values for the top 10 genes with the 
# largest biological components in the PBMC dataset. Each point represents the log-expression value in 
# a single cell.
plotExpression(sce, features=rownames(top.dec)[1:10])


#################################### Dimensionality reduction ####################################
# We use the denoisePCA() function with the assumed Poisson technical trend, to choose the number of 
# dimensions to retain after PCA. The red line represents the chosen number of PCs.
sce <- denoisePCA(sce, technical=new.trend, approx=TRUE)
ncol(reducedDim(sce, "PCA"))

plot(attr(reducedDim(sce), "percentVar"), xlab="PC",
     ylab="Proportion of variance explained")
abline(v=ncol(reducedDim(sce, "PCA")), lty=2, col="red")

# Examination of the first few PCs already reveals some strong substructure in the data (Figure 9).
# Pairwise PCA plots of the first three PCs in the PBMC dataset, constructed from normalized 
# log-expression values of genes with positive biological components.
# Each point represents a cell, coloured by the log-number of expressed features.
plotPCA(sce, ncomponents=3, colour_by="log10_total_features_by_counts")

# This is recapitulated with a t-SNE plot (Figure 10). Again, note that we set use_dimred= to perform 
# t-SNE on the denoised expression matrix. Each point represents a cell and is coloured according to 
# the log-number of expressed features.
sce <- runTSNE(sce, use_dimred="PCA", perplexity=30, rand_seed=100)
plotTSNE(sce, colour_by="log10_total_features_by_counts")

#################################### Clustering with graph-based methods ####################################
# We build a shared nearest neighbour graph (Xu and Su 2015) and use the Walktrap algorithm to 
# identify clusters.\
snn.gr <- buildSNNGraph(sce, use.dimred="PCA")
clusters <- igraph::cluster_walktrap(snn.gr)
sce$Cluster <- factor(clusters$membership)
table(sce$Cluster)

# We look at the ratio of the observed and expected edge weights to confirm that the clusters are 
# modular. (We don’t look at the modularity score itself, as that varies by orders of magnitudes 
# across clusters and is difficult to interpret.) Figure 11 indicates that most of the clusters 
# are well seperated, with few strong off-diagonal entries.
cluster.mod <- clusterModularity(snn.gr, sce$Cluster, get.values=TRUE)
log.ratio <- log2(cluster.mod$observed/cluster.mod$expected + 1)

library(pheatmap)
pheatmap(log.ratio, cluster_rows=FALSE, cluster_cols=FALSE, 
         color=colorRampPalette(c("white", "blue"))(100))

# We examine the cluster identities on a t-SNE plot (Figure 12) to confirm that different clusters 
# are indeed separated.
plotTSNE(sce, colour_by="Cluster")



######################################## Marker gene detection ######################################
# We detect marker genes for each cluster using findMarkers(). Again, we only look at upregulated 
# genes in each cluster, as these are more useful for positive identification of cell types in a 
# heterogeneous population.
markers <- findMarkers(sce, clusters=sce$Cluster, direction="up")
# We examine the markers for cluster 1 in more detail. The upregulation of genes such as PF4 and 
# PPBP suggests that cluster 1 contains platelets or their precursors.
marker.set <- markers[["1"]]
head(marker.set[,1:8], 10) # only first 8 columns, for brevity

# This is confirmed in Figure 13, where the transcriptional profile of cluster 1 is clearly distinct 
# from the others.
chosen <- rownames(marker.set)[marker.set$Top <= 10]
plotHeatmap(sce, features=chosen, exprs_values="logcounts", 
            zlim=5, center=TRUE, symmetric=TRUE, cluster_cols=FALSE,
            colour_columns_by="Cluster", columns=order(sce$Cluster))

######################################## Concluding remarks ########################################
saveRDS(sce, file="pbmc_data.rds")


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

# 注意，在上游数据处理的时候：
# Some feature-counting tools will report mapping statistics in the count matrix (e.g., the number of 
# unaligned or unassigned reads). While these values can be useful for quality control, they would be 
# misleading if treated as gene expression values. Thus, they should be removed (or at least moved to 
# the colData) prior to further analyses.


##################################### Setting up the data #####################################
# One matrix was generated for each plate of cells used in the study.
# unzip("E-MTAB-5522.processed.1.zip"). Unzip once. 这里我们关注标明“Calero”的数据集
# Reading in the count tables for each of the two plates.
setwd('/Users/mijiarui/Nature_Biotechnology_Paper/simpleSingleCell')
plate1 <- read.delim("counts_Calero_20160113.tsv", 
                     header=TRUE, row.names=1, check.names=FALSE)
plate2 <- read.delim("counts_Calero_20160325.tsv", 
                     header=TRUE, row.names=1, check.names=FALSE)
head(plate1[1:5,1:5])
gene.lengths <- plate1$Length # First column is the gene length.
plate1 <- as.matrix(plate1[,-1]) # Discarding gene length (as it is not a cell).
plate2 <- as.matrix(plate2[,-1])
rbind(Plate1=dim(plate1), Plate2=dim(plate2)) # 两个文件中均有96个细胞，检测46703个基因

# We combine the two matrices into a single object for further processing. This is done after 
# verifying that the genes are in the same order between the two matrices.
stopifnot(identical(rownames(plate1), rownames(plate2)))
all.counts <- cbind(plate1, plate2)

# For convenience, the count matrix is stored in a SingleCellExperiment object from the SingleCellExperiment 
# package. This allows different types of row- and column-level metadata to be stored alongside the 
# counts for synchronized manipulation throughout the workflow.
library(SingleCellExperiment)
sce <- SingleCellExperiment(list(counts=all.counts)) # 构建SingleCellExperiment对象的时候最少只需要输入矩阵信息
rowData(sce)$GeneLength <- gene.lengths
sce # assays中只有一个slot，是'counts'

# We identify the rows corresponding to ERCC spike-in transcripts from the row names. We store this 
# information in the SingleCellExperiment object for future use. This is necessary as spike-ins 
# require special treatment in downstream steps such as normalization.
# 一个很小的注意事项，在使用正则表达式去除ERCC spike-in的时候：
# Be aware of using the ^ERCC regular expression for human data where the row names of the count 
# matrix are gene symbols. An ERCC gene family actually exists in human annotation, so this would 
# result in incorrect identification of genes as spike-in transcripts. This problem can be avoided by 
# publishing count matrices with standard identifiers (e.g., Ensembl, Entrez).
isSpike(sce, "ERCC") <- grepl("^ERCC", rownames(sce)) # isSpike是sce这个SingleCellExperiment对象中的一个单元
row.names(sce)[grepl("^ERCC", rownames(sce))] # 可以用这句代码来查看ERCC的标签
summary(isSpike(sce, "ERCC")) # 或者用table函数也可以

# This dataset is slightly unusual in that it contains information from another set of spike-in 
# transcripts, the Spike-In RNA Variants (SIRV) set. For simplicity, we will only use the ERCC 
# spike-ins in this analysis. Thus, we must remove the rows corresponding to the SIRV transcripts 
# prior to further analysis, which can be done simply by subsetting the  SingleCellExperiment object.
is.sirv <- grepl("^SIRV", rownames(sce))
rownames(sce)[is.sirv] # 可以用这句代码来查看以‘SIRV’开头的另一套spike-in标签
sce <- sce[!is.sirv,] 
summary(is.sirv)
sce # SIRV不是常规使用的spike-in，所以在SingleCellExperiment对象中并没有对应的存储单元


# Incorporating cell-based annotation
# We load in the metadata for each library/cell from the sdrf.txt file. It is important to check that 
# the rows of the metadata table are in the same order as the columns of the count matrix. Otherwise, 
# incorrect metadata will be assigned to each cell. metadata的列名和matrix的行名必须完全一致，尤其是顺序。
metadata <- read.delim("E-MTAB-5522.sdrf.txt", check.names=FALSE, header=TRUE)
head(metadata) # 在这个文件中，metadata的行名存储在列名为“Source Name"(实际为第一列)
m <- match(colnames(sce), metadata[["Source Name"]]) # Enforcing identical order，检测是否匹配
stopifnot(all(!is.na(m))) # Checking that nothing's missing.
metadata <- metadata[m,]
head(colnames(metadata))
dim(metadata)  # 我们看到metadata中的信息（列数）非常多，注意目前只是存储在一个名为metadata的中间变量


# We only retain relevant metadata fields to avoid storing unnecessary information in the  colData 
# of the SingleCellExperiment object. In particular, we keep the plate of origin (i.e., block) and 
# phenotype of each cell. The second field is relevant as all of the cells contain a CBFB-MYH11 
# oncogene, but the expression of this oncogene is only induced in a subset of the cells.
colData(sce)$Plate <- factor(metadata[["Factor Value[block]"]])
pheno <- metadata[["Factor Value[phenotype]"]]
levels(pheno) <- c("induced", "control") # 根据因子型变量的首字母在字母表中的顺序
colData(sce)$Oncogene <- pheno
table(colData(sce)$Oncogene, colData(sce)$Plate)

# Incorporating gene-based annotation
# Feature-counting tools typically report genes in terms of standard identifiers from Ensembl or 
# Entrez. These identifiers are used as they are unambiguous and highly stable. However, they are 
# difficult to interpret compared to the gene symbols which are more commonly used in the literature. 
# Given the Ensembl identifiers, we obtain the corresponding gene symbols using annotation packages 
# like org.Mm.eg.db. 定量分析软件一般使用Ensembl或者Entrez基因id，但是为了方便阅读，需要转换成gene symbol
# 这一步转化工作可以使用注释包，例如：org.Mm.eg.db; 如果是斑马鱼的，则是org.Dr.eg.db
library(org.Mm.eg.db)
symb <- mapIds(org.Mm.eg.db, keys=rownames(sce), keytype="ENSEMBL", column="SYMBOL") # 类似于键值对的概念；
## 提示信息：有可能多个Ensembl id会对应同一个gene symbol或者有缺失值
rowData(sce)$ENSEMBL <- rownames(sce)
rowData(sce)$SYMBOL <- symb
head(rowData(sce))
table(is.na(rowData(sce)$SYMBOL))  # 有相当一部分Ensembl id不能对应到gene symbol，显示为NA


# It is often desirable to rename the row names of sce to the gene symbols, as these are easier to 
# interpret. However, this requires some work to account for missing and duplicate symbols. The code 
# below will replace missing symbols with the Ensembl identifier and concatenate duplicated symbols 
# with the (unique) Ensembl identifiers.
new.names <- rowData(sce)$SYMBOL
missing.name <- is.na(new.names)
new.names[missing.name] <- rowData(sce)$ENSEMBL[missing.name]
dup.name <- new.names %in% new.names[duplicated(new.names)]
new.names[dup.name] <- paste0(new.names, "_", rowData(sce)$ENSEMBL)[dup.name]
rownames(sce) <- new.names
head(rownames(sce))

# We also determine the chromosomal location for each gene using the TxDb.Mmusculus.UCSC.mm10.ensGene 
# package. This will be useful later as several quality control metrics will be computed from rows 
# corresponding to mitochondrial genes.
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
location <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene, keys=rowData(sce)$ENSEMBL, 
                   keytype="GENEID", column="CDSCHROM")
rowData(sce)$CHR <- location
summary(location=="chrM")


##################################### Quality control on the cells #####################################
# Defining the quality control metrics
# Cells with small library sizes are of low quality as the RNA has not been efficiently captured 
# (i.e., converted into cDNA and amplified) during library preparation.

# The number of expressed features in each cell is defined as the number of features with non-zero 
# counts for that cell. Any cell with very few expressed genes is likely to be of poor quality as the 
# diverse transcript population has not been successfully captured.

# The proportion of reads mapped to spike-in transcripts is calculated relative to the library size 
# for each cell. High proportions are indicative of poor-quality cells, where endogenous RNA has been 
# lost during processing (e.g., due to cell lysis or RNA degradation). The same amount of spike-in 
# RNA to each cell, so an enrichment in spike-in counts is symptomatic of loss of endogenous RNA. 
# spike-in所占的比例过高提示cell quality比较差

# In the absence of spike-in transcripts, the proportion of reads mapped to genes in the mitochondrial 
# genome can also be used. High proportions are indicative of poor-quality cells 
# (Islam et al. 2014; Ilicic et al. 2016), possibly because of loss of cytoplasmic RNA from perforated cells. 
# The reasoning is that mitochondria are larger than individual transcript molecules and less likely to escape 
# through tears in the cell membrane.
# 在没有spike-in的情况下，使用线粒体基因的转录本也可以用于cell quality control；过高比例的线粒体基因转录本
# 提示细胞破损。

# For each cell, we calculate these quality control metrics using the calculateQCMetrics function 
# from the scater package (McCarthy et al. 2017). These are stored in the row- and column-wise 
# metadata of the SingleCellExperiment for future reference.
library(scater)
mito <- which(rowData(sce)$CHR=="chrM") # 给出线粒体基因转录本的行号
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=mito))
head(colnames(colData(sce)), 10)

# The distributions of these metrics are shown in Figure 1, stratified by oncogene induction status 
# and plate of origin. The aim is to remove putative low-quality cells that have low library sizes, 
# low numbers of expressed features, and high spike-in (or mitochondrial) proportions.
sce$PlateOnco <- paste0(sce$Oncogene, ".", sce$Plate)
multiplot(                  # multiplot是ggplot2下的函数
  plotColData(sce, y="total_counts", x="PlateOnco"),
  plotColData(sce, y="total_features", x="PlateOnco"),
  plotColData(sce, y="pct_counts_ERCC", x="PlateOnco"),
  plotColData(sce, y="pct_counts_Mt", x="PlateOnco"),
  cols=2)

# Generally, they will be in rough agreement, i.e., cells with low total counts will also have low 
# numbers of expressed features and high ERCC/mitochondrial proportions. Clear discrepancies may 
# correspond to technical differences between batches of cells (see below) or genuine biological 
# differences in RNA content.
par(mfrow=c(1,3))
plot(sce$total_features, sce$total_counts/1e6, xlab="Number of expressed genes",
     ylab="Library size (millions)")
plot(sce$total_features, sce$pct_counts_ERCC, xlab="Number of expressed genes",
     ylab="ERCC proportion (%)")
plot(sce$total_features, sce$pct_counts_Mt, xlab="Number of expressed genes",
     ylab="Mitochondrial proportion (%)")

# Identifying outliers for each metric
# 有的时候阈值的设定非常的tricky，我们不妨就假设大部分的cell都是high quality的，那么在离群值范围的
# 细胞我们就认为是低质量的细胞，需要filter掉。Outliers are defined based on the median absolute 
# deviation (MADs) from the median value of each metric across all cells. We remove cells with 
# log-library sizes that are more than 3 MADs below the median log-library size. A log-transformation 
# improves resolution at small values, especially when the MAD of the raw values is comparable to or 
# greater than the median. We also remove cells where the log-transformed number of expressed genes 
# is 3 MADs below the median value. 我们使用的是中位绝对偏差，我们讲检测到的library size(counts)和gene
# 数在3个MADs以下的排除掉
libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower",  # isOutlier()是scater包的内置函数
                          log=TRUE, batch=sce$PlateOnco)
feature.drop <- isOutlier(sce$total_features, nmads=3, type="lower", 
                          log=TRUE, batch=sce$PlateOnco)

# 我们也要滤除spike-in比例过高的细胞，为了更加清楚的展示(因为发现离群值中的大值)，我们不对ERCC的percentage
# 进行log转换
spike.drop <- isOutlier(sce$pct_counts_ERCC, nmads=3, type="higher",
                        batch=sce$PlateOnco)

# Subsetting by column will retain only the high-quality cells that pass each filter described above. 
# We examine the number of cells removed by each filter as well as the total number of retained cells. 
# Removal of a substantial proportion of cells (> 10%) may be indicative of an overall issue with data 
# quality.
keep <- !(libsize.drop | feature.drop | spike.drop)
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
           BySpike=sum(spike.drop), Remaining=sum(keep))

# We then subset the SingleCellExperiment object to retain only the putative high-quality cells. 
# We also save the original object to file for later use.
sce$PassQC <- keep
saveRDS(sce, file="416B_preQC.rds")
sce # before cell QC，192个细胞
sce <- sce[,keep]
dim(sce) # after cell QC，183个细胞


################################### Classification of cell cycle phase ###################################
# 英国 Sanger 研究院的 Teichmann 等人就开发了一款单细胞隐藏 可变模型 (single-cell latent variable model, 
# scLVM) 和 Cyclone 软件。Cyclone 软 件可以利用机器学习技术和统计学的方法， 将细胞周期信息与单细胞 RNA 
# 测序数据结合起来，来帮助我们判断哪些基因表达信号与细胞 周期的哪个阶段有关。对于任何单细胞的 RNA 测序数
# 据，使用 Cyclone 软件就能够使其与每一个细 胞在细胞周期中所处的阶段一一对应。两者的 scLVM 就采用了在整个
# 细胞周期中表达程度高 度可变的基因作为研究对象，来明确基因表达 与细胞周期的关系。他们确定了某一种决定细
# 胞周期的因子，也发现了在细胞发育或分化的 整个细胞周期中，推动细胞转化的因子，而且 还发现了一些独特的亚
# 群 (subpopulations)。设计到训练集和测试集：We use the prediction method described by Scialdone et al. 
# (2015) to classify cells into cell cycle phases based on the gene expression data. Using a training dataset, 
# the sign of the difference in expression between two genes was computed for each pair of genes. Pairs with 
# changes in the sign across cell cycle phases were chosen as markers. Cells in a test dataset can then be 
# classified into the appropriate phase, based on whether the observed sign for each marker pair is consistent 
# with one phase or another.
# This approach is implemented in the cyclone function from the scran package. The package contains a 
# pre-trained set of marker pairs for mouse data, which we can load in the the readRDS function. We use the 
# Ensembl identifiers for each gene in our dataset to match up with the names in the pre-trained set of gene pairs.
set.seed(100)
library(scran)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", # 加载训练集
                                package="scran"))
assignments <- cyclone(sce, mm.pairs, gene.names=rowData(sce)$ENSEMBL)  # 判断测试集细胞的细胞周期状态
assignments # 结果是一个列表，包含分类的属性(phases)，score和normalized score
# Each cell is assigned a score for each phase, with a higher score corresponding to a higher probability 
# that the cell is in that phase. We focus on the G1 and G2/M scores as these are the most informative 
# for classification. 
par(mfrow = c(1,1))
plot(assignments$score$G1, assignments$score$G2M, 
     xlab="G1 score", ylab="G2/M score", pch=16)

# Cells are classified as being in G1 phase if the G1 score is above 0.5 and greater than the G2/M 
# score; in G2/M phase if the G2/M score is above 0.5 and greater than the G1 score; and in S phase 
# if neither score is above 0.5. Here, the vast majority of cells are classified as being in G1 phase. 
# We save these assignments into the SingleCellExperiment object for later use.
sce$phases <- assignments$phases
table(sce$phases)

# 可视化的展示
pheatmap::pheatmap(t(assignments$scores))
pheatmap::pheatmap(t(assignments$normalized.scores))

## 有几个注意点：
## The classifier may not be accurate for data that are substantially different from those used in 
## the training set, e.g., due to the use of a different protocol. In such cases, users can construct 
## a custom classifier from their own training data using the sandbag function. This will also be 
## necessary for other model organisms where pre-trained classifiers are not available. 也就是说对于特别
## 的建库方法，使用默认的训练集并不合适；亦或是其它模式动物，可以考虑使用sandbag函数来构建自己的训练集。

## 在进行细胞周期检测之前，不要将低表达gene滤除：Do not filter out low-abundance genes before applying 
## cyclone. Even if a gene is not expressed in any cell, it may still be useful for classification if 
## it is phase-specific. Its lack of expression relative to other genes will still yield informative pairs, 
## and filtering them out would reduce power.

################################### Examining gene-level expression metrics ###################################
# We examine the identities of the most highly expressed genes (Figure 4). This should generally be dominated by 
# constitutively expressed transcripts, such as those for ribosomal or mitochondrial proteins. The presence of 
# other classes of features may be cause for concern if they are not consistent with expected biology. For example, 
# a top set containing many spike-in transcripts suggests that too much spike-in RNA was added during library preparation, 
# while the absence of ribosomal proteins and/or the presence of their pseudogenes are indicative of suboptimal alignment.
# 通常来讲，高表达的基因是和细胞生物学功能相关，或者是线粒体内或者核糖体基因。过好的spike-in不是好事情。
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotQC(sce, type = "highest-expression", n=50) + fontsize  # 这一步可能会卡，最好先清一清内存空间

# Filtering out low-abundance genes：
# Several metrics can be used to define low-abundance genes. The most obvious is the average count for 
# each gene, computed across all cells in the dataset. We calculate this using the calcAverage() function, 
# which also performs some adjustment for library size differences between cells. We typically observe a 
# peak of moderately expressed genes following a plateau of lowly expressed genes (Figure 5).
# A minimum threshold can be applied to this value to filter out genes that are lowly expressed. The example 
# below demonstrates how we could remove genes with average counts less than 1. The number of TRUE values in 
# demo.keep corresponds to the number of retained rows/genes after filtering.
ave.counts <- calcAverage(sce, use_size_factors=FALSE)
hist(log10(ave.counts), breaks=100, main="", col="grey80", 
     xlab=expression(Log[10]~"average count"))

demo.keep <- ave.counts >= 1
filtered.sce <- sce[demo.keep,]
summary(demo.keep)

# We also examine the number of cells that express each gene. This is closely related to the average 
# count for most genes, as expression in many cells will result in a higher average (Figure 6). Genes 
# expressed in very few cells are often uninteresting as they are driven by amplification artifacts 
# (though they may also also arise from rare populations). We could then remove genes that are 
# expressed in fewer than n cells. “Intensity of colour corresponds to the number of genes at any given location.”
num.cells <- nexprs(sce, byrow=TRUE)
smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells", 
              xlab=expression(Log[10]~"average count"))

# We remove genes that are not expressed in any cell to reduce computational work in downstream steps. 
# Such genes provide no information and would be removed by any filtering strategy.
to.keep <- num.cells > 0
sce <- sce[to.keep,]
summary(to.keep)


############################## Normalization of cell-specific biases ##############################
# Using the deconvolution method to deal with zero counts
# Single-cell data can be problematic for these bulk data-based methods (DESeq2 and edgeR normalization) due 
# to the dominance of low and zero counts. To overcome this, we pool counts from many cells to increase the 
# count size for accurate size factor estimation (Lun, Bach, and Marioni 2016). Pool-based size factors are 
# then “deconvolved” into cell-based factors for cell-specific normalization.
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
# Fitting a trend to the spike-in variances
# Variability in the observed expression values across genes can be driven by genuine biological heterogeneity 
# or uninteresting technical noise. To distinguish between these two possibiltiies, we need to model the 
# technical component of the variance of the expression values for each gene. We do so using the set of 
# spike-in transcripts, which were added in the same quantity to each cell. Thus, the spike-in transcripts 
# should exhibit no biological variability, i.e., any variance in their counts should be technical in 
# origin.

# We use the trendVar() function to fit a mean-dependent trend to the variances of the log-expression values for the spike-in transcripts. We set block= to block on the plate of origin for each cell, to ensure that technical differences between plates do not inflate the variances. Given the mean abundance of a gene, the fitted value of the trend is then used as an estimate of the technical component for that gene. The biological component of the variance is finally calculated by subtracting the technical component from the total variance of each gene with the decomposeVar function.
var.fit <- trendVar(sce, parametric=TRUE, block=sce$Plate,
                    loess.args=list(span=0.3))
var.out <- decomposeVar(sce, var.fit)
head(var.out)

# We visually inspect the trend to confirm that it corresponds to the spike-in variances (Figure 8)). 
# The wave-like shape is typical of the mean-variance trend for log-expression values. A linear 
# increase in the variance is observed as the mean increases from zero, as larger variances are 
# possible when the counts increase. At very high abundances, the effect of sampling noise decreases 
# due to the law of large numbers, resulting in a decrease in the variance.
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
     ylab="Variance of log-expression")
curve(var.fit$trend(x), col="dodgerblue", lwd=2, add=TRUE)
cur.spike <- isSpike(sce)
points(var.out$mean[cur.spike], var.out$total[cur.spike], col="red", pch=16)

# We check the distribution of expression values for the genes with the largest biological components. 
# This ensures that the variance estimate is not driven by one or two outlier cells (Figure 9).
chosen.genes <- order(var.out$bio, decreasing=TRUE)[1:10]
plotExpression(sce, features=rownames(var.out)[chosen.genes]) + fontsize


################################ Removing the batch effect ################################
library(limma)
assay(sce, "corrected") <- removeBatchEffect(logcounts(sce), 
                                             design=model.matrix(~sce$Oncogene), batch=sce$Plate)
assayNames(sce)


############################## Denoising expression values using PCA ##############################
# Once the technical noise is modelled, we can use principal components analysis (PCA) to remove random technical noise.
# We assume that biological processes involving co-regulated groups of genes will account for the 
# most variance in the data. If this is the case, this process should be represented by one or more 
# of the earlier PCs. In contrast, random technical noise affects each gene independently and will be 
# represented by later PCs. The denoisePCA() function removes later PCs until the total discarded 
# variance is equal to the sum of technical components for all genes used in the PCA.

# 在减少了背景噪声（某些基因）的基础上，对数据进行进一步降维（一般降维到原来维度的4-7%），以进一步去除背景噪音
# denoisePCA() will only use genes that have positive biological components, i.e., variances greater 
# than the fitted trend. This guarantees that the total technical variance to be discarded will not 
# be greater than the total variance in the data. 所以上面modeling technical noise一部是必须先执行的。

# No filtering is performed on abundance here, which ensures that PCs corresponding to rare subpopulations 
# can still be detected. Discreteness is less of an issue as low-abundance genes also have lower variance, 
# thus reducing their contribution to the PCA.
sce <- denoisePCA(sce, technical=var.fit$trend, assay.type="corrected")
dim(reducedDim(sce, "PCA")) 

# It is also possible to obtain a low-rank approximation of the original expression matrix, capturing 
# the variance equivalent to the retained PCs. This is useful for denoising prior to downstream 
# procedures that require gene-wise expression values.
sce2 <- denoisePCA(sce, technical=var.fit$trend, 
                   assay.type="corrected", value="lowrank") 
assayNames(sce2)

######################## Data exploration with dimensionality reduction #######################
# We visualize the relationships between cells by constructing pairwise PCA plots for the first three 
# components (Figure 10). Cells with similar expression profiles should be located close together in 
# the plot, while dissimilar cells should be far apart. In this case, we observe a clear separation 
# of cells based on the oncogene induction status, consistent with the expected effects on the transcriptome.
# 不同的颜色可以区分(induced or control)
plotReducedDim(sce, use_dimred="PCA", ncomponents=3, 
               colour_by="Oncogene") + fontsize

# By comparison, we observe no clear separation of cells by batch (Figure 11). This indicates that 
# our batch correction step using removeBatchEffect() was successful. 不同颜色(batch)应该搅合在一起
plotReducedDim(sce, use_dimred="PCA", ncomponents=3, 
               colour_by="Plate") + fontsize


# Another widely used approach for dimensionality reduction is the t-stochastic neighbour embedding 
# (t-SNE) method (Van der Maaten and Hinton 2008). t-SNE tends to work better than PCA for separating 
# cells in more diverse populations. This is because the former can directly capture non-linear 
# relationships in high-dimensional space, whereas the latter must represent them on linear axes. 
# However, this improvement comes at the cost of more computational effort and requires the user to
# consider parameters such as the random seed and perplexity (see comments). We demonstrate the 
# generation of t-SNE plots in Figure 12, using the low-rank approximation of the data to take 
# advantage of the denoising step.
run_args <- list(rand_seed=100, use_dimred="PCA")
out5 <- plotTSNE(sce, run_args=c(run_args, perplexity=5),
                 colour_by="Oncogene") + fontsize + ggtitle("Perplexity = 5")

out10 <- plotTSNE(sce, run_args=c(run_args, perplexity=10),
                  colour_by="Oncogene") + fontsize + ggtitle("Perplexity = 10")

out20 <- plotTSNE(sce, run_args=c(run_args, perplexity=20),
                  colour_by="Oncogene") + fontsize + ggtitle("Perplexity = 20")

multiplot(out5, out10, out20, cols=3)

## There are many other dimensionality reduction techniques that we do not consider here but could 
## also be used, e.g., multidimensional scaling, diffusion maps. These have their own advantages and 
## disadvantages – for example, diffusion maps (see plotDiffusionMap) place cells along a continuous 
## trajectory and are suited for visualizing graduated processes like differentiation

## t-SNE is a stochastic method, so users should run the algorithm several times to ensure that the 
## results are representative. Scripts should set a seed (via the  rand_seed argument) to ensure 
## that the chosen results are reproducible. It is also advisable to test different settings of the 
## “perplexity” parameter as this will affect the distribution of points in the low-dimensional 
## space. A good guide on how to interpret t-SNE plots can be found at http://distill.pub/2016/misread-tsne/.
######################## Clustering cells into putative subpopulations ########################
# 大致策略是先用PCA降维去噪，计算距离矩阵，根据距离矩阵来进行clustering。其中使用Ward’s criterion可减少
# 同一个cluster的variance。然后使用动态剪切树的方法来进一步减少cluster的数量。然后在经过PCA降维的数据基础上
# 进行tSNE可视化以进一步区分不同的细胞，cluster是用距离矩阵和剪切数来定义的。
# 这种在去除背景噪音基因的基础上，再进行两次数据降维的方法（PCA降维后再用tSNE）非常适合可视化。
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

# Detecting marker genes
# 1, Once putative subpopulations are identified by clustering, we can identify marker genes for each 
# cluster using the findMarkers function. This fits a linear model to the log-expression values for 
# each gene using limma (Ritchie et al. 2015). The aim is to test for DE in each cluster compared to 
# the others while blocking on uninteresting factors such as the plate of origin in design. The top 
# DE genes are likely to be good candidate markers as they can effectively distinguish between cells 
# in different clusters.
# 2, findMarkers can also be directed to find genes that are DE between the chosen cluster and all 
# other clusters. This should be done by setting pval.type="all", which defines the p-value for each 
# gene as the maximum value across all pairwise comparisons involving the chosen cluster. Combined 
# with direction="up", this can be used to identify unique markers for each cluster. However, this 
# is sensitive to overclustering, as unique marker genes will no longer exist if a cluster is split 
# into two smaller subclusters.
markers <- findMarkers(sce, my.clusters, block=sce$Plate)
markers

markers_strict <- findMarkers(sce, my.clusters, block=sce$Plate, pval.type= "all", direction = 'up')
markers_strict # 看一下fdr
# To construct a marker set for cluster 1 from the top 10 genes of each comparison, one would filter 
# marker.set to retain rows with Top less than or equal to 10. Other statistics are also reported for 
# each gene, including the adjusted p-values (see below) and the log-fold changes relative to every 
# other cluster.
marker.set <- markers[["1"]]
head(marker.set, 10)

# We save the list of candidate marker genes for further examination.
write.table(marker.set, file="416B_marker_1.tsv", sep="\t", 
            quote=FALSE, row.names=FALSE)

# We visualize the expression profiles of the top candidates to verify that the DE signature is robust 
# (Figure 15). Most of the top markers have strong and consistent up- or downregulation in cells of 
# cluster 1 compared to some or all of the other clusters. A cursory examination of the heatmap 
# indicates that cluster 1 contains oncogene-induced cells with strong downregulation of DNA replication 
# and cell cycle genes. This is consistent with the potential induction of senescence as an 
# anti-tumorigenic response (Wajapeyee et al. 2010). A more comprehensive investigation of the 
# function of these markers can be performed with gene set enrichment analyses, e.g., using kegga 
# or goana from limma.
top.markers <- rownames(marker.set)[marker.set$Top <= 10]
plotHeatmap(sce, features=top.markers, columns=order(sce$cluster), 
            colour_columns_by=c("cluster", "Plate", "Oncogene"),
            cluster_cols=FALSE, center=TRUE, symmetric=TRUE, zlim=c(-5, 5))
            
## 找到单一一个细胞类型特异表达的gene有的时候标准过于苛刻，在部分subpopulation存在DE就可以认为marker
## Many of the markers in Figure 15 are not uniquely up- or downregulated in the chosen cluster. 
# Testing for unique DE tends to be too stringent as it overlooks important genes that are expressed 
# in two or more clusters. For example, in a mixed population of CD4+-only, CD8+-only, double-positive 
# and double-negative T cells, neither Cd4 or Cd8 would be detected as subpopulation-specific markers 
# because each gene is expressed in two subpopulations. With our approach, both of these genes will 
# be picked up as candidate markers as they will be DE between at least one pair of subpopulations. 
# A combination of markers can then be chosen to characterize a subpopulation, which is more flexible 
# than trying to find uniquely DE genes.

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
