# Script Comments ---------------------------------------------------------
## New script for identifying relevant traits and relating to modules

library(here)

# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
load(here("Datasets/RData_Files/Data_Preprocessing/data.RData"))


# Load network data saved in the second part.

load(here("Datasets/RData_Files/Network_Construction/data.Rdata"))

#read in traits data:
traits <- read.csv(here("Datasets/data.csv"))



# Wrangle Traits Data -----------------------------------------------------

#match sample names and trait data:
Samples <- rownames(filteredLog2rpkm_outliers_removed)
traitRows <- match(Samples, traits$SampleName)

traits <- cbind(traits, group = paste(traits$Var1, traits$Var2, traits$Var3, sep = "_"))

datTraits <- traits[traitRows, ]

datTraits<- as.data.frame(datTraits, )


# Calculate Eigengenes ----------------------------------------------------


# Define numbers of genes and samples
nGenes = ncol(filteredLog2rpkm_outliers_removed);
nSamples = nrow(filteredLog2rpkm_outliers_removed);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(filteredLog2rpkm_outliers_removed, moduleColors)$eigengenes



# Dendrogram for Eigengenes -----------------------------------------------

MEList = moduleEigengenes(filteredLog2rpkm_outliers_removed, colors = moduleColors) 
MEs = MEList$eigengenes

#write out MEs combined with trait data for later use:

write.csv(cbind(datTraits, MEs), here("Output/Bicor/Quantifying_Modules_and_ANOV_Prep/data.csv"))

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)

pdf(file = here("Output/Bicor/Quantifying_Modules_and_ANOV_Prep/data.pdf"), 
    wi = 9, he = 6)

plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# Plot the cut line into the dendrogram
abline(h=0.25, col = "red")

dev.off()

## number of genes in each module
table(moduleColors)

write.csv(table(moduleColors), 
          file = here("Output/Bicor/Quantifying_Modules_and_ANOV_Prep/data.csv"), 
          row.names = F)
############# merge module pink and red

MEDissThres = 0.25

# Call an automatic merging function
merge = mergeCloseModules(filteredLog2rpkm_outliers_removed, labels2colors(moduleLabels), cutHeight = MEDissThres, verbose = 5)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

table(mergedColors)
#replot
sizeGrWindow(12, 9)

par(mar=c(1,1,1,1))
pdf(file = here("Output/Bicor/Quantifying_Modules_and_ANOV_Prep/data.pdf"))

plotDendroAndColors(geneTree, mergedColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "All Samples Cluster Dendrogram")
dev.off()

###

# Relabeling Categorical Data as Numerical --------------------------------


#change groupings to numeric

datTraits$Treatment[datTraits$Treatment == "CTRL"] <- 0
datTraits$Treatment[datTraits$Treatment == "AVar1"] <- 1
datTraits$Treatment[datTraits$Treatment == "BVar1"] <- 2
datTraits$Treatment[datTraits$Treatment == "CVar1"] <- 3
datTraits$Treatment[datTraits$Treatment == "DVar1"] <- 4

datTraits$BrainRegion[datTraits$BrainRegion == "Var2.1"] <- 5
datTraits$BrainRegion[datTraits$BrainRegion == "Var2.2"] <- 6
datTraits$BrainRegion[datTraits$BrainRegion == "Var2.3"] <- 7

datTraits$Hemisphere.[datTraits$Hemisphere. == "Var3.1"] <- 8
datTraits$Hemisphere.[datTraits$Hemisphere. == "Var3.2"] <- 9


#numeric traits only
numerictraits <- datTraits[,4:6]
rownames(numerictraits) <- datTraits$SampleName




#We now have the expression data in the variable datExpr, and the corresponding clinical traits 
#in the variable datTraits. Before we continue with network construction and module detection, 
#we visualize how the clinical traits relate to the sample dendrogram.

# Re-cluster samples

datExpr <- filteredLog2rpkm_outliers_removed

sampleTree2 = hclust(dist(datExpr), method = "average")


##Var2

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors_brain = numbers2colors(as.numeric(numerictraits$Var2), signed = FALSE);

# Plot the sample dendrogram and the colors underneath
pdf(file = here("Output/Bicor/Quantifying_Modules_and_ANOV_Prep/data.pdf"), 
    width = 12, height = 9)

plotDendroAndColors(sampleTree2, traitColors_brain,
                    groupLabels = names(numerictraits),
                    main = "Sample Dendrograms")
dev.off()

##Treatment

traitColors_treatment = numbers2colors(as.numeric(numerictraits$Treatment), signed = FALSE);

# Plot the sample dendrogram and the colors underneath.
pdf(file = here("Output/Bicor/Quantifying_Modules_and_ANOV_Prep/data.pdf"), 
    width = 12, height = 9)

plotDendroAndColors(sampleTree2, traitColors_treatment,
                    groupLabels = names(numerictraits),
                    main = "Sample Dendrogram Across Var1 Treatments")
dev.off()

##Var3

traitColors_var3 = numbers2colors(as.numeric(numerictraits$Var3), 
                                        colors = blueWhiteRed(50), signed = FALSE);

# Plot the sample dendrogram and the colors underneath.
pdf(file = here("Output/Bicor/Quantifying_Modules_and_ANOV_Prep/Sample_Clustering_data.pdf"), 
    width = 12, height = 9)

plotDendroAndColors(sampleTree2, traitColors_var3,
                    groupLabels = names(numerictraits),
                    main = "Sample Dendrogram Across Var3")
dev.off()


#write out MEs together with trait info:


write.csv(MEs, file = here("Output/Bicor/Quantifying_Modules_and_ANOV_Prep/data.csv"))

save(MEs0, MEList, MEs, MEDiss, METree, merge, mergedColors, datTraits, numerictraits, 
     file = here("Datasets/RData_Files/Quantifying_Modules_and_ANOV_Prep/data.RData"))
