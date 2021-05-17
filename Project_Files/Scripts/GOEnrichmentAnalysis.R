# Load Necessary Packages and Datasets -------------------------------------------------
#run the following scripts as needed if you don't have the necessary packages installed:

#install.packages("here")
#install.packages(here("Datasets/Enrichment_Analysis_Collections/anRichment_1.01-2.tar.gz"), repos = NULL, type = "source")
#install.packages(here("Datasets/Enrichment_Analysis_Collections/anRichmentMethods_0.90-1.tar.gz", repos = NULL, type = "source"))
#install.packages(here("Datasets/Enrichment_Analysis_Collections/SCSBrainCellTypeCollection_1.00.tgz", repos = NULL, type = "source"))
#install.packages(here("Datasets/Enrichment_Analysis_Collections/BrainDiseaseCollection_1.00.tgz", repos = NULL, type = "source"))

library(here)
library(dplyr)
library(anRichment)
library(anRichmentMethods)


options(stringsAsFactors = FALSE)


load(here("sample_dataset.RData"))


# Getting Entrez IDs from BiomaRt -----------------------------------------
##these scripts query the entrez ID from biomaRt.It takes a while to run
##uncomment if you need to reconstruct from the source the entrez IDs. 
##Otherwise, you can use the previously retrieved dataset that is loaded below


# 
# BiocManager::install("biomaRt")
# 
# library(biomaRt)
## get mouse ensemble data
# mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

##query ensemble and entrez ids
# genes <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
#                filter= "ensembl_gene_id",
#                values = datExpr_grey_rm$Ensemble_IDs,
#                mart = mouse)

##remove duplicates
# genes_duplicates_rm <- genes[!duplicated(genes$ensembl_gene_id),]

##reorder dataset to ensure gene order matches
# genes_duplicates_rm_2 <- genes_duplicates_rm[match(datExpr_grey_rm$Ensemble_IDs, 
#                                                    genes_duplicates_rm$ensembl_gene_id),]
# 
##check that the old dataset and queried gene list is equal in order
# all.equal(datExpr_grey_rm$Ensemble_IDs, genes_duplicates_rm_2$ensembl_gene_id)
# 
##combine
# merged_gene_ids <- cbind(genes_duplicates_rm_2, datExpr_grey_rm$Modules)

##set column names
# colnames(merged_gene_ids) <- c("Ensemble_ID",
#                                "Entrez_ID",
#                                "Module")

##save for reference just in case
# write.csv(merged_gene_ids, 
#           file = here("file.csv"),
#           row.names = FALSE)


# Wrangling Entrez IDs ----------------------------------------------------

##import entrez ID dataset
merged_gene_ids <-read.csv(here("dataset.csv"))

datExpr_entrez <- cbind(Entrez_ID = merged_gene_ids$Entrez_ID, 
                        datExpr_grey_rm)

datExpr_entrez_rm_na <- datExpr_entrez[!is.na(datExpr_entrez$Entrez_ID),]

row.names(datExpr_entrez_rm_na) <- datExpr_entrez_rm_na$Entrez_ID

entrez_id <- datExpr_entrez_rm_na$Entrez_ID


# Building GO Collection -----------------------------------------------------------
#build GO collection
GOColl <- buildGOcollection(organism = "mouse")

#run enrichment analysis
GOenrichment <- enrichmentAnalysis(classLabels = datExpr_entrez_rm_na$Modules,
                                   identifiers = entrez_id,
                                   refCollection = GOColl,
                                   useBackground = "intersection",
                                   threshold = 0.05,
                                   thresholdType = "FDR",
                                   geneSeparator = ",",
                                   getOverlapEntrez = FALSE,
                                   getOverlapSymbols = TRUE,
                                   maxReportedOverlapGenes = 500,
                                   ignoreLabels = "grey")

collectGarbage()

#extract enrichment table
GOenrichment_table <- GOenrichment$enrichmentTable 

#save files
write.csv(GOenrichment_table, 
          file = here("dataset.csv"),
          row.names = FALSE)

#subset only xxx data
GOenrichment_xxx <- GOenrichment_table %>%
  filter(class == "xxx")

#save files
write.csv(GOenrichment_xxx,
        file = here("dataset"),
         row.names = FALSE)





# Building Internal Collection -----------------------------------------------------
#build Internal collection
InternalColl <- internalCollection(organism = "mouse")

#run enrichment analysis
Internalenrichment <- enrichmentAnalysis(classLabels = datExpr_entrez_rm_na$Modules,
                                   identifiers = entrez_id,
                                   refCollection = InternalColl,
                                   useBackground = "intersection",
                                   threshold = 0.05,
                                   thresholdType = "FDR",
                                   geneSeparator = ",",
                                   getOverlapEntrez = FALSE,
                                   getOverlapSymbols = TRUE,
                                   maxReportedOverlapGenes = 500,
                                   ignoreLabels = "grey")

collectGarbage()

#extract enrichment table
Internalenrichment_table <- Internalenrichment$enrichmentTable 

#save files
write.csv(Internalenrichment_table, 
          file = here("dataset.csv"),
          row.names = FALSE)

#subset only xxx data
Internalenrichment_xxx <- Internalenrichment_table %>%
  filter(class == "xxx")

#save files
write.csv(Internalenrichment_xxx,
          file = here("Output/Bicor/Enrichment_Analysis/Enrichment_Using_All_Modules/xxx_Only/InternalColl_Enrichment_xxx.csv"),
          row.names = FALSE)

# NCBI BioSystems Collection ----------------------------------------------
#build NCBI collection
NCBIColl <- BioSystemsCollection(organism = "mouse")

#run enrichment analysis
NCBIenrichment <- enrichmentAnalysis(classLabels = datExpr_entrez_rm_na$Modules,
                                   identifiers = entrez_id,
                                   refCollection = NCBIColl,
                                   useBackground = "intersection",
                                   threshold = 0.05,
                                   thresholdType = "FDR",
                                   geneSeparator = ",",
                                   getOverlapEntrez = FALSE,
                                   getOverlapSymbols = TRUE,
                                   maxReportedOverlapGenes = 500,
                                   ignoreLabels = "grey")

collectGarbage()

#extract enrichment table
NCBIenrichment_table <- NCBIenrichment$enrichmentTable 

#save files
write.csv(NCBIenrichment_table, 
          file = here("dataset.csv"),
          row.names = FALSE)

#subset only xxx data
NCBIenrichment_xxx <- NCBIenrichment_table %>%
  filter(class == "xxx")

#save files
write.csv(NCBIenrichment_xxx,
          file = here("Output/Bicor/Enrichment_Analysis/Enrichment_Using_All_Modules/xxx_Only/NCBIColl_Enrichment_xxx.csv"),
          row.names = FALSE)



# Molecular Signatures Database Collection --------------------------------
#build MSDB collection
## this takes annoyingly long to run interactively so I just saved it as an Rdata file
## to load it in and work quickly. If it needs to be built from scratch again, will need the XML from the website. See the 
## anRichment tutorial if this is needed, in the MSDB section. The file is too big to store.


#save(MSDBColl, file = here("dataset.RData"))

load(here("dataset.RData"))
     
#run enrichment analysis
MSDBenrichment <- enrichmentAnalysis(classLabels = datExpr_entrez_rm_na$Modules,
                                   identifiers = entrez_id,
                                   refCollection = MSDBColl,
                                   useBackground = "intersection",
                                   threshold = 0.05,
                                   thresholdType = "FDR",
                                   geneSeparator = ",",
                                   getOverlapEntrez = FALSE,
                                   getOverlapSymbols = TRUE,
                                   maxReportedOverlapGenes = 500,
                                   ignoreLabels = "grey")

collectGarbage()

#extract enrichment table
MSDBenrichment_table <- MSDBenrichment$enrichmentTable 

#save files
write.csv(MSDBenrichment_table, 
          file = here("dataset.csv"),
          row.names = FALSE)

#subset only xxx data
MSDBenrichment_xxx <- MSDBenrichment_table %>%
  filter(class == "xxx")

#save files
write.csv(MSDBenrichment_xxx,
          file = here("dataset.csv"),
          row.names = FALSE)



# SCS Brain Cell Types Collection -----------------------------------------
library(SCSBrainCellTypeCollection)

#build SCSBrainCellType collection
SCSBrainCellTypeColl <- SCSBrainCellTypeCollection(organism = "mouse")

#run enrichment analysis
SCSBrainCellTypeenrichment <- enrichmentAnalysis(classLabels = datExpr_entrez_rm_na$Modules,
                                   identifiers = entrez_id,
                                   refCollection = SCSBrainCellTypeColl,
                                   useBackground = "intersection",
                                   threshold = 0.05,
                                   thresholdType = "FDR",
                                   geneSeparator = ",",
                                   getOverlapEntrez = FALSE,
                                   getOverlapSymbols = TRUE,
                                   maxReportedOverlapGenes = 500,
                                   ignoreLabels = "grey")

collectGarbage()

#extract enrichment table
SCSBrainCellTypeenrichment_table <- SCSBrainCellTypeenrichment$enrichmentTable 

#save files
write.csv(SCSBrainCellTypeenrichment_table, 
          file = here("data.csv"),
          row.names = FALSE)

#subset only xxx data
SCSBrainCellTypeenrichment_xxx <- SCSBrainCellTypeenrichment_table %>%
  filter(class == "xxx")

#save files
write.csv(SCSBrainCellTypeenrichment_xxx,
          file = here("data.csv"),
          row.names = FALSE)




# Brain Disease Collection ------------------------------------------------
library(BrainDiseaseCollection)


#build BrainDisease collection
BrainDiseaseColl <- BrainDiseaseCollection(organism = "mouse")

#run enrichment analysis
BrainDiseaseenrichment <- enrichmentAnalysis(classLabels = datExpr_entrez_rm_na$Modules,
                                   identifiers = entrez_id,
                                   refCollection = BrainDiseaseColl,
                                   useBackground = "intersection",
                                   threshold = 0.05,
                                   thresholdType = "FDR",
                                   geneSeparator = ",",
                                   getOverlapEntrez = FALSE,
                                   getOverlapSymbols = TRUE,
                                   maxReportedOverlapGenes = 500,
                                   ignoreLabels = "grey")

collectGarbage()

#extract enrichment table
BrainDiseaseenrichment_table <- BrainDiseaseenrichment$enrichmentTable 

#save files
write.csv(BrainDiseaseenrichment_table, 
          file = here("data.csv"),
          row.names = FALSE)

#subset only xxx data
BrainDiseaseenrichment_xxx <- BrainDiseaseenrichment_table %>%
  filter(class == "xxx")

#save files
write.csv(BrainDiseaseenrichment_xxx,
          file = here("data.csv"),
          row.names = FALSE)


#Custom Microglial Collection ----------------------------------------
#run analysis for custom made microglial collection

#import microglia collection
load(here("data.RData"))

#run enrichment analysis
CustomMGenrichment <- enrichmentAnalysis(classLabels = datExpr_entrez_rm_na$Modules,
                                   identifiers = entrez_id,
                                   refCollection = MGCollection,
                                   useBackground = "intersection",
                                   threshold = 0.05,
                                   thresholdType = "FDR",
                                   geneSeparator = ",",
                                   getOverlapEntrez = FALSE,
                                   getOverlapSymbols = TRUE,
                                   maxReportedOverlapGenes = 500,
                                   ignoreLabels = "grey")

collectGarbage()

#extract enrichment table
CustomMGenrichment_table <- CustomMGenrichment$enrichmentTable 

#save files
write.csv(CustomMGenrichment_table, 
          file = here("data.csv"),
          row.names = FALSE)

#subset only xxx data
CustomMGenrichment_xxx <- CustomMGenrichment_table %>%
  filter(class == "xxx")

#save files
write.csv(CustomMGenrichment_xxx,
          file = here("data.csv"),
          row.names = FALSE)


