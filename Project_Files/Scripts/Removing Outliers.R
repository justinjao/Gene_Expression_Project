
# Introduction to Script --------------------------------------------------

#This script is for cleaning the pre-filtered dataset. It plots a dendrogram
#that visualizes any obvious outliers, which you can use to remove if need be

#Loading Datasets and Packages --------------------------------------------------------
##imports prior data and required packages

library(here)

load(here("Datasets/RData_Files/Data_Preprocessing/data1.Rdata"))

load(here("Datasets/RData_Files/Data_Preprocessing/data2.Rdata"))

library(WGCNA)


# Visualize Hierarchical Clustering to Identify Obvious Outliers ----------

#perform hierarchical clustering
sampleTree <- hclust(dist(filteredLog2rpkm), method = "average")

#create a plot to visualize dendrogram

sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
pdf(here("Output/Data_Preprocessing/Sample_Clustering_Visualizing_Outliers.pdf"), width = 10, height = 10)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

#can add horizontal red line to show cutoff point
#in this case, ACHC35, ACSTR52 & ACSTR60 look like outliers, 
#so we can pick a suitable number to remove them
abline(h = 86, col = "red")

dev.off()
# Removing Potential Outliers ---------------------------------------------

#split the tree into 2 clusters, one group for all the clusters below and one for above the line
clust <- cutreeStatic(sampleTree, cutHeight = 86, minSize = 10)

#displays number of clusters in each group
table(clust)

# We keep clust 1
keepSamples <- (clust==1)
filteredLog2rpkm_outliers_removed <- filteredLog2rpkm[keepSamples, ]
nGenes <- ncol(filteredLog2rpkm)
nSamples <- nrow(filteredLog2rpkm)

#visualize cut dendrogram again, if need be
sampleTree_outliers_removed <- hclust(dist(filteredLog2rpkm_outliers_removed), method = "average")

pdf(here("Output/Data_Preprocessing/Clustering_Outliers_Removed.pdf"), width = 10, height = 10)
plot(sampleTree_outliers_removed, main = "Clustering w/ Outliers Removed", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

dev.off()

save(filteredLog2rpkm, filteredLog2rpkm_outliers_removed, 
     file = here("Datasets/RData_Files/Data_Preprocessing/data.RData"))

