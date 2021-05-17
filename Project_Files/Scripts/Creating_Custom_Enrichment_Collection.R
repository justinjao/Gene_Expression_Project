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
library(biomaRt)


options(stringsAsFactors = FALSE)

##read in relevant data
custom_list <- read.csv(here("Datasets/Enrichment_Analysis_Collections/data.csv"))


# Getting Entrez IDs for Custom xxx Enrichment Collection from BiomaRt -----------------------------------------
##these scripts query the entrez ID of the custom xxx dataset from biomaRt.

## character vector of relevant groups to keep
relevant_groups <- c("xxx",
                     "xxx yyy",
                     "bbb ddd",
                     "vvv")


##subset the custom list for only the groups we want
custom_list <- custom_list %>% 
    filter(custom_list$groups %in% relevant_groups)

## generate a vector of all the unique ensemble IDs we have, to match the entrez IDs to.
## once we retrieve all the entrez IDs, we can use this to match to our original list
## we create a different vector here instead of subsetting the original dataset because
## we want to keep the duplicates, but if we submit a query with duplicates to BioMart, this makes
## later matching of IDs a little more difficult

unique_ensemble_ids <- custom_list$ensembl_gene_id[!duplicated(custom_list$ensembl_gene_id)]

## get mouse ensemble data
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

##query ensemble and entrez ids of custom dataset
genes <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
               filter= "ensembl_gene_id",
               values = unique_ensemble_ids,
               mart = mouse)


##remove duplicates from our entry
genes_duplicates_rm <- genes[!duplicated(genes$ensembl_gene_id),]

## after merging, strangely, there are a few genes in our custom list that aren't on Ensemble.
## the problem IDs are listed here:

#get list of IDs that are in our original custom list but not in our queried list from ensemble
problem_ids <- custom_list$ensembl_gene_id[!(custom_list$ensembl_gene_id %in% genes_duplicates_rm$ensembl_gene_id)]


# trial <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
#                filter= "ensembl_gene_id",
#                values = problem_ids,
#                mart = mouse)

## Thus, we remove them from our original list to ensure that they match.
custom_list_fixed <- custom_list %>% 
  filter(!custom_list$ensembl_gene_id %in% problem_ids)


## create vector containing matching indices for the gene IDs 
## (i.e. the position of the unique gene IDs in the original list)
MatchList <- match(custom_list_fixed$ensembl_gene_id, genes_duplicates_rm$ensembl_gene_id,)

## use the matching indices to modify the gene IDs taken from biomart to generate a new dataframe containing the ensemble IDs
## in the same order as the original custom list we have
matched_gene_ids <- genes_duplicates_rm[match(custom_list_fixed$ensembl_gene_id, genes_duplicates_rm$ensembl_gene_id,),]


#check that the ensemble IDs between the original custom list and the new generated IDs are equal. If true, this would mean
## that the entrez IDs would be in the correct order as well
all.equal(matched_gene_ids$ensembl_gene_id, custom_list_fixed$ensembl_gene_id)

## add a new column containing the entrez IDs to the original custom list
custom_list_entrez <- custom_list_fixed %>% 
    mutate(Entrez_ID = matched_gene_ids$entrezgene_id)



##save gene IDs for reference just in case
write.csv(matched_gene_ids,
          file = here("Datasets/Enrichment_Analysis_Collections/data.csv"),
          row.names = FALSE)

#import data again if needed
matched_gene_ids <- read.csv(here("Datasets/Enrichment_Analysis_Collections/data.csv"))

###now we remove the NAs
custom_list_entrez_na_rm <- custom_list_entrez[!is.na(custom_list_entrez$Entrez_ID),]

## generate unique list of genesets for later use
unique_list <- custom_list_entrez_na_rm[!duplicated(custom_list_entrez_na_rm$listname),]


# Generating geneSetInfo --------------------------------------------------

##create dataframe for gene sets
geneSetInfo <- data.frame(ID = unique_list$listname,
                          Name = unique_list$listname,
                          ShortName = unique_list$listname,
                          Description = unique_list$description,
                          Source = unique_list$source,
                          Groups = unique_list$groups,
                          InternalClassification = unique_list$groups,
                          organism = "mouse",
                          evidence = "ND")


# Generating geneContent --------------------------------------------------

## create dataframe of gene content
geneContent <- data.frame(Entrez = custom_list_entrez_na_rm$Entrez_ID,
                          ID = custom_list_entrez_na_rm$listname,
                          GeneSetName = custom_list_entrez_na_rm$listname,
                          evidence = "ND",
                          source = custom_list_entrez_na_rm$source)




# Generating GroupInfo ----------------------------------------------------

## create dataframe of group information 
groupDF <- data.frame(name = unique(unique_list$groups),
                      description = c("xxx genes",
                                      "xxx dd genes",
                                      "yyy genes",
                                      "DD genes"),
                      source = "genelist",
                      Parents = "",
                      AlternateNames = unique(unique_list$groups))
                      

# Generating Collection from Dataframes -----------------------------------------------------------
##generate collection
XXCollection <- collectionFromDataFrames(geneSetInfoDF = geneSetInfo,
                         geneSetContentDF = geneContent,
                         groupDF = groupDF)

##save data for use in GO enrichment
save(XXCollection, file = here("Datasets/Enrichment_Analysis_Collections/xxx.RData"))

