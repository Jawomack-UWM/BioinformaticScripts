#Project 2 DE Analysis
#Justin Womack

#########Making sure all requirements installed#################
# Indicate package repositories to R
repositories <- c("https://cloud.r-project.org",
                  "https://bioconductor.org/packages/3.7/bioc",
                  "https://bioconductor.org/packages/3.7/data/annotation",
                  "https://bioconductor.org/packages/3.7/data/experiment",
                  "https://www.stats.ox.ac.uk/pub/RWin",
                  "http://www.omegahat.net/R",
                  "https://R-Forge.R-project.org",
                  "https://www.rforge.net",
                  "https://cloud.r-project.org",
                  "http://www.bioconductor.org",
                  "http://www.stats.ox.ac.uk/pub/RWin")

# Package list to download
packages <- c("vsn", "UpSetR", "gplots", "NMF", "org.Hs.eg.db",
              "pheatmap", "tximport", "readr","edgeR", "biomaRt",
              "VennDiagram", "plyr", "dplyr","DESeq2", "AnnotationDbi",
              "Biobase", "ensembldb", "ggpubr", "ggplot2", "limma", "magrittr")

# Install and load missing packages
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  install.packages(new.packages, repos = repositories)
}
lapply(packages, require, character.only = TRUE)

#############################################
#####Project Begin - Data Transformation#####
#############################################

setwd('C:\\Users\\Justi\\Desktop\\Project2')
samples <- read.table('samples.txt')

# The 'conditions' (Itraconazole of DMSO treatment) are in the 10th column
conditions <- samples[[10]]
sampleTable <- data.frame(condition = conditions)

# get the table of read counts
readcounts <- read.table('featureCounts_results.txt', header = TRUE)

# the gene IDs should be stored as row.names
row.names(readcounts) <- readcounts$Geneid

# exclude all columns that do not contain read counts
readcounts <- readcounts[, -c(1:6)]
head(readcounts)

# give meaningful sample names 
names(readcounts) <- samples$V5
head(readcounts)

######################################################
####Preparing DESeqDataSet for use with DESeq2########
######################################################
dds_star <- DESeqDataSetFromMatrix(readcounts, sampleTable, ~condition)
dim(dds_star)

# remove genes without any counts
dds_star <- dds_star[ rowSums(counts(dds_star)) > 0, ]
dim(dds_star)

# Inspect DESeq dataset 
colData(dds_star) %>% head

# investigate different library sizes 
colSums(counts(dds_star))


######################################################
##########Normalization and Transformation############
######################################################

# calculate the size factor and add it to the data set
dds_star <- estimateSizeFactors(dds_star)
sizeFactors(dds_star)

# retrieve the normalized read counts
counts.sf_normalized <- counts(dds_star, normalized = TRUE)

# inspect counts.sf_normalized
head(counts.sf_normalized)

# transform size-factor normalized read counts
# to log2 scale using a pseudocount of 1
log.norm.counts <- log2(counts.sf_normalized + 1)

# first, boxplots of non-transformed read counts (one per sample)
boxplot(counts.sf_normalized,
        notch = TRUE,
        main = "untransformed read counts",
        ylab = "read counts")

# box plots of log2-transformed read counts
boxplot(log.norm.counts,
        notch = TRUE,
        main = "log2 -transformed read counts",
        ylab = "log2(read counts)")

plot(log.norm.counts[, 1:2], cex = 0.1,
     main = "Normalized log2 (read counts)")

#check for heteroskedasticity
msd_plot <- meanSdPlot(log.norm.counts, ranks = FALSE, plot = FALSE)
msd_plot$gg + ggtitle("sequencing depth normalized log2(read counts)") +
  ylab("standard deviation")

###Trying to reduce heterosket####

# obtain regularized log-transformed values
DESeq.rlog <- rlog(dds_star, blind = TRUE)
rlog.norm.counts <- assay(DESeq.rlog)
head(rlog.norm.counts)

# mean-sd plot for rlog-transformed data
msd_plot <- meanSdPlot(rlog.norm.counts, ranks = FALSE, plot = FALSE)
msd_plot$gg + ggtitle("rlog-transformed read counts") +
  ylab("standard deviation")

plot(rlog.norm.counts[, 1:2], cex = 0.1,
     main = "Normalized log2(read counts)")

###########FUN WITH PCA#############

cor(counts.sf_normalized, method = "pearson")
# cor() calculates the correlation between columns of a matrix
distance.m_rlog <- as.dist(1 - cor(rlog.norm.counts, method = "pearson"))
# plot() can directly interpret the output of hclust()
plot(hclust(distance.m_rlog),
     labels = colnames(rlog.norm.counts),
     main = "rlog-transformed read counts\ndistance: Pearson correlation")

#pc with base graphics
pc <- prcomp(t(rlog.norm.counts))
plot(pc$x[, 1], pc$x[, 2],
     col = colData(dds_star)[, 1],
     main = "PCA of seq.depth normalized\n and rlog-transformed read counts")

# PCA with Deseq (looks much prettier)
P <- plotPCA(DESeq.rlog)
# plot cosmetics
P <- P + theme_bw() + ggtitle("rlog-transformed counts")
print(P)

##################################################
########DIFFERENTIAL GENE ANALYSIS################
##################################################

# DESeq uses the levels of the condition
# to determine the order of the comparison
str(colData(dds_star)$condition)
# set 'new_born' as the first-level factor
colData(dds_star)$condition <- relevel(colData(dds_star)$condition, "DMSO")

# Running DGE analysis using the DESeq() function
dds2 <- DESeq(dds_star)
DGE.results <- results(dds2, independentFiltering = TRUE, alpha = 0.05)
summary(DGE.results)
head(DGE.results)
table(DGE.results$padj < 0.05)
rownames(subset(DGE.results, padj < 0.05))

# sequencing depth normalization between the samples
dds3 <- estimateSizeFactors(dds_star)
dds3 <- estimateDispersions(dds3)
dds3 <- nbinomWaldTest(dds3)


############################
#######PLOTTING DE##########
############################

#hist of pvalues
hist(DGE.results$pvalue,
     col = "grey",
     border = "white",
     xlab = "",
     ylab = "",
     main = "frequencies of p-values")
#realtionship between expression change and condition
plotMA(DGE.results,
       alpha = 0.05,
       main = "Imazoldole vs. DMSO",
       ylim = c(-4,4))

# Volcano plot for a threshold of adjusted pval=0.05 and logFC=7
with(DGE.results,
     plot(log2FoldChange, -log10(padj), pch = 20,
          main = "Volcano plot", xlim = c(-10,10)))
with(subset(DGE.results, padj < 0.05),
     points(log2FoldChange,-log10(padj), pch = 20, col = "blue"))

with(subset(DGE.results, padj < 0.05 & abs(log2FoldChange) > 7),
     points(log2FoldChange, -log10(padj), pch = 20, col = "red"))

#heatmap fun
# sort the results according to the log2FoldChange
DGE.results.sorted <- DGE.results[order(DGE.results$log2FoldChange), ]

# identify genes with the desired adjusted p-value or fold change
DGEgenes <- rownames(subset(DGE.results.sorted , padj < 0.0001))
DGEgenes <- rownames(subset(DGE.results.sorted, abs(log2FoldChange) > 5))
length(DGEgenes)

hm.mat_DGEgenes <- log.norm.counts[DGEgenes, ]

# scale the read counts per gene to emphasize
# the sample type-specific differences
aheatmap(hm.mat_DGEgenes,
         Rowv = TRUE,
         Colv = TRUE,
         distfun = "euclidean",
         hclustfun = "average",
         scale = "row")

