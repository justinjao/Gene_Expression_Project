# Script Comments ---------------------------------------------------------
## Justin Revised Dataprep scripts for WGCNA Input. This performs necessary
## preprocessing steps to wrangle the data into the right format

# Loading Packages --------------------------------------------------------

library(ggplot2)
library(dplyr)
library(cowplot)
library(RColorBrewer)
library(gplots)
library(ggrepel)
library(tidyr)
library(readxl)
library(tidyverse)
library(here)

#Import Datasets & Workspaces -------------------------------------------------

##read in count matrix containing expression info
UMIcollapsed_countmatrix <- as.data.frame(read_csv(here("Datasets/data.csv")))


##read in clinical info dataset containing Var1 treatment 
##information for later corrrelation
clinical_info <- as.data.frame(read_xlsx(here("Datasets/data2.xlsx")))

##change column names to be valid (basically remove spaces etc.)
colnames(clinical_info) <- make.names(colnames(clinical_info))


#Organizing Datasets -----------------------------------------------------

##change row names of clinical info to sample IDs
rownames(clinical_info) <- clinical_info$Samlple.ID

#make counts into matrix and convert to dataframe
countmatrix <- UMIcollapsed_countmatrix[,2:ncol(UMIcollapsed_countmatrix)]
countmatrix <- as.data.frame(countmatrix)

##make rownames into GENE IDs
rownames(countmatrix) <- UMIcollapsed_countmatrix$Geneid

##Get order of clinical info
colorder <- clinical_info$Samlple.ID

##convert order into character format
colorder <- as.character(colorder)

colorder

colnames(countmatrix)
#Rearrange expression info to match count matrix order
countmatrix <- countmatrix[,colorder]
colnames(countmatrix)

#convert back to matrix
countmatrix <- as.matrix(countmatrix)


##Making DGE List Object in EdgeR -------------------------------------------
##compile a differential gene expression list object
##this combines the clinical info with the gene expression info together

##load necessary library
library(edgeR)
library(limma)


##concatenate Var2 and Var1 columns to create grouping
clinical_info$group <- paste(clinical_info$Var2,clinical_info$Var1, sep="_")

##create the differential gene expression list object, grouping by "group"
RG <- DGEList(counts = countmatrix, group=clinical_info$group)

# use dimensions to see number of unique transcripts
dim(RG)
#53801    60


# Filtering via CPM -------------------------------------------------------

# Filter for one count per million in at least 6 libraries

#ensure that at least 
keep <- rowSums(cpm(RG)>1)>=4
RGcpm <- RG[keep,]
dim(RGcpm)
# 16928    60

table(rowSums(RGcpm$counts==0)==4)

#create vector of ensemble gene IDs for biomart
geneid <- rownames(RGcpm)


# Reading in Gene Info ----------------------------------------------------

library(biomaRt)

#select mouse gene dataset
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
#head(listAttributes(human),100)
#attributes <- listAttributes(mouse)
#attributes[grep("exon", attributes$name),]

#extract genetic data off of biomart
genes <- getBM(attributes = c("ensembl_gene_id","chromosome_name","start_position","end_position","mgi_symbol","external_gene_name","description"), 
               filter= "ensembl_gene_id",
               values = geneid, 
               mart = mouse)

#calculate gene length
genes$genelength <- abs(genes$end_position - genes$start_position) 

#remove duplicates:keeps only first entry
genes <- genes[!duplicated(genes$ensembl_gene_id),]

#match order: df[match(target, df$name),]
genes <- genes[match(geneid,genes$ensembl_gene_id),]

#save(genes,file="mm10_Ensemble_geneInfo.RData")


# Adding Gene Info to DEGList Object --------------------------------------

#add into DEGlist
RGcpm$genes <- genes


# Filtering Plot ----------------------------------------------------------
##create a density plot to view the effect of filtering
library(RColorBrewer)
nsamples <- ncol(RGcpm)

colourCount = nsamples
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
fill=getPalette(colourCount)

#plot:
pdf(here("Output/Data_Preprocessing/FilteringCPM_plots.pdf"), h=6,w=10)
par(mfrow=c(1,2))

#prefilter:
lcpm <- cpm(RG, log=TRUE, prior.count=2)


plot(density(lcpm[,1]), col=fill[1], lwd=2, ylim=c(0,0.5), las=2, 
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=fill[i], lwd=2)
}


#filtered data
#og-CPM of zero threshold (equivalent to a CPM value of 1) used in the filtering ste
lcpm <- cpm(RGcpm, log=TRUE, prior.count=2)
plot(density(lcpm[,1]), col=fill[1], lwd=2, ylim=c(0,0.5), las=2, 
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=fill[i], lwd=2)
}
#legend("topright", Samples, text.col=fill, bty="n")
dev.off()


#reset library sizes
RGcpm$samples$lib.size <- colSums(RGcpm$counts)


#plot library sizes
pdf(here("Output/Data_Preprocessing/LibrarySizes.pdf"),w=30,h=8)
barplot(RGcpm$samples$lib.size,names=colnames(RGcpm),las=2)
# Add a title to the plot
title("Barplot of library sizes")
dev.off()



# Full Model with Repeated Measure Correction -----------------------------

#main effect of treatment: four groups

#make design matrix
clinical_info$Var1 <- factor(clinical_info$Condition, levels=c("CTRL","AAA", "BBB", "CCC", "DDD"))

clinical_info$Var2 <- factor(clinical_info$Var2)

clinical_info$Var3 <- factor(clinical_info$Var3)

clinical_info$SampleName <- as.character(clinical_info$SampleName)


#set design matrix
design <- model.matrix(~0+Var2+Var1, data=clinical_info) #if using a 0 intercept must set up contrasts

colnames(design)


# Normalization TMM and Voom ----------------------------------------------

#Normalize library
RGnorm=calcNormFactors(RGcpm,method =c("TMM")) #TMM normalization for library size/composition

RGnorm$samples$group


#estimate dispersion
RGnorm <- estimateGLMCommonDisp(RGnorm,design,verbose=TRUE)##RG change
#Disp = 0.01928 , BCV = 0.1388 

RGnorm <- estimateGLMTrendedDisp(RGnorm,design)##RG change
RGnorm <- estimateGLMTagwiseDisp(RGnorm,design)##RG change

pdf(here("Output/Data_Preprocessing/plotBCV.pdf"))

plotBCV(RGnorm)##RG change

dev.off()



# Write Out FPKM Table ----------------------------------------------------


#get normalized count matrix:
nc <- cpm(RGnorm, normalized.lib.sizes=T,log=F, prior.count=1)
nc <-as.data.frame(nc)

#RPKM
rpkm <- rpkm(RGnorm, gene.length=RGnorm$genes$genelength, normalized.lib.sizes=TRUE, log=FALSE, prior.count=0.25)

#log2RPKM
log2rpkm <- rpkm(RGnorm, gene.length=RGnorm$genes$genelength, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.25)

#save data
save(nc,rpkm,log2rpkm, clinical_info,RGnorm,file=here("Datasets/RData_Files/Data_Preprocessing/CountInputData.RData"), version = 2)


# Data Wrangling for WGCNA ------------------------------------------------

#transform rows= samples, col=genes

rpkm_t <- t(rpkm)


#BiocManager::install("WGCNA")
library(WGCNA)

################## process data ########################
#detect genes with missing values or samples with missing values
gsg = goodSamplesGenes(rpkm_t, verbose = 3);
gsg$allOK # FALSE: Excluding 67 genes from the calculation due to too many missing samples or zero variance.

#remove the genes with too many missing values
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(rpkm_t)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(rpkm_t)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  rpkm_t = rpkm_t[gsg$goodSamples, gsg$goodGenes]
}

gsg = goodSamplesGenes(rpkm_t, verbose = 3);
gsg$allOK 

########################################################################################################################
#Filtering
########################################################################################################################
#gene with RPKM value of 0.25 or higher in at least one sample.
#max(col)>=.25
#https://stackoverflow.com/questions/24212739/how-to-find-the-highest-value-of-a-column-in-a-data-frame-in-r/24212879
colMax <- function(data) sapply(data, max, na.rm = TRUE)

rpkm_t <- as.data.frame(rpkm_t)

filteredrpkm_t <- rpkm_t[, (colMax(rpkm_t) >= .25)]
#from 16861genes
#to 12956 genes


#log2 transform:
Log2rpkm <- log2(filteredrpkm_t +.25)

# median absolute deviation 
#remove if = 0 https://support.bioconductor.org/p/65124/
colMad <- function(data) sapply(data, mad, na.rm = TRUE)

filteredLog2rpkm = Log2rpkm[, colMad(Log2rpkm) != 0]
#11626 genes for input to network

save(filteredLog2rpkm, file=here("Datasets/RData_Files/Data_Preprocessing/filteredLog2RPKM.networkinput.Rdata"), version = 2)
