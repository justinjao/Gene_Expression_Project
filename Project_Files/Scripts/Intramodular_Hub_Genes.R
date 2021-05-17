# Script Comments -------------------------------------------------
## this is the script for analyzing:
## - intramodular connecrtivity vs gene significance
## - module membership vs gene significance
## - intramodular connectivity vs module membership
## - identification of hub genes within the color module
## - modifying the module membership vs gene significance plot to include hub genes

## for loops are used to perform the same analysis for each treatment condition


# Load Packages & Necessary Datasets --------------------------------------
library(dplyr)
library(WGCNA)
library(here)

#read in  data
datExpr_orig <- read.csv("data.csv")

load(here("data.RData"))

#read in treatment information
order_df <- read.csv(here("data.csv"))

# Data Wrangling ----------------------------------------------------------
#change rownames to ensemble IDs
rownames(datExpr_orig) <- datExpr_orig$Ensemble_IDs

#get module color names
modNames = substring(names(MEs), 3)

#Get Module Assignments
colorh1 <- as.character(datExpr_orig$Modules)

#remove unecessary columns
datExpr <- datExpr_orig[,4:ncol(datExpr_orig)]

#transpose matrix to get it in the right format (genes on columns, samples on rows)
datExpr_t <- data.frame(t(datExpr))

#make sure gene names transfer to columns
colnames(datExpr_t) <- rownames(datExpr)


# Calculating Intramodular Connectivity -----------------------------------------------

#create adjacency matrix and analyze connectivity of different modules
ADJ1=abs(cor(datExpr_t,use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ1, colorh1)
head(Alldegrees1)

# Calculating Gene Module Membership ------------------------------
#calculate gene module membership
geneModuleMembership = as.data.frame(cor(t(datExpr), MEs, use = "p"));

#Calculating Gene Significance ------------------------------


#binarizes treatment conditions to allow for correlation analysis
#(sets value to 1 if it is the treatment condition, otherwise sets to 0)
#see here for more info(https://peterlangfelder.com/2018/11/25/working-with-categorical-variables/)
bin <- binarizeCategoricalColumns(order_df$Treatment_Condition, dropFirstLevelVsAll = F)


#calculate gene significance for each treatment condition
geneTraitSignificance   = as.data.frame(cor(datExpr_t, bin$data.Ctrl.vs.all, use = "p"));
geneTraitSignificancexxxxA = as.data.frame(cor(datExpr_t, bin$data.Axxxx.vs.all, use = "p"));
geneTraitSignificancexxxxB = as.data.frame(cor(datExpr_t, bin$data.Bxxxx.vs.all, use = "p"));
geneTraitSignificancexxxxC = as.data.frame(cor(datExpr_t, bin$data.Cxxxx.vs.all, use = "p"));
geneTraitSignificancexxxxD = as.data.frame(cor(datExpr_t, bin$data.Dxxxx.vs.all, use = "p"));

#get unique color labels
colorlevels=unique(colorh1)

# Module Membership vs Gene Significance ----------------------------------
## create a scatter plot showing module membership vs gene significance for genes
## in the color module, for each treatment condition

#set module to color
module = "color"
column = match(module, modNames);
#subset genes to only be color module genes
moduleGenes = colorh1==module;

#set graphical parameters
sizeGrWindow(9, 9);
par(mfrow = c(1,1))

#for loop to run the same analysis for each treatment condition
for (Treatment in c("Ctrl", "xxxxA", "xxxxB", "xxxxC", "xxxxD")) {
  
  #saves unique file for each treatment condition
  pdf(here(paste0("Output/Bicor/Intramodular_Gene_Significance_Analysis/",
                 Treatment,
                 "/",
                 Treatment,
                 "_ModuleMembership_vs_GeneSignificance.pdf")))
  
  #create call for gene significance values that varies with treatment condition
  geneTraitSignificance <- as.name((paste0("geneTraitSignificance", Treatment)))
  
  #create scatterplot
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(eval(geneTraitSignificance)[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = paste("Gene Significance For", Treatment),
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  
  dev.off()
  
}

# Intramodular Connectivity vs Gene Significance --------------------------

## create plot of intramodular connectivity vs gene significance

#for loop to run the same analysis for each treatment condition
for (Treatment in c("Ctrl", "xxxxA", "xxxxB", "xxxxC", "xxxxD")) {
  
  #create unique file for each treatment condition  
  pdf(here(paste0("Output/Bicor/Intramodular_Gene_Significance_Analysis/",
                 Treatment,
                 "/",
                 Treatment,
                 "_IntramodularConnectivity_vs_GeneSignificance.pdf")))
  
  #set graphical parameters
  par(mfrow=c(4,2))
  par(mar = c(4,5,3,1))
  
  #create gene significance call that varies with treatment condition
  geneTraitSignificance <- as.name((paste0("geneTraitSignificance", Treatment)))
  
  #creates a scatter plot for each module
  for (i in c(1:length(colorlevels)))
  {
    #subsets the module color for each loop
    whichmodule=colorlevels[[i]];
    
    #subsets the genes to be plotted based on the color
    restrict1 = (colorh1==whichmodule);
    
    #makes the plot
    verboseScatterplot(Alldegrees1$kWithin[restrict1],
                       abs(unlist(eval(geneTraitSignificance))[restrict1]), col=colorh1[restrict1],
                       main=whichmodule,
                       xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
  }
  
  dev.off()
  
}


# Writing GS vs Intramodular Connectivity into Table ----------------------

## basically calculates the correlation and p values for gene signficance
## vs intramodular connectivity. Written to have identical values to the plots
## that were generated in the previous section. Also automatically writes it out
## to a table for each treatment condition 

#create empty dataframe with headers to be filled in via for loop. Has headers
#for correlation and p values for each treatment condition

cor_table <- data.frame(Module      = colorlevels,
                        Ctrl_Cor     = numeric(8),
                        Ctrl_Cor_P   = numeric(8),
                        xxxxA_Cor   = numeric(8),
                        xxxxA_Cor_P = numeric(8),
                        xxxxB_Cor   = numeric(8),
                        xxxxB_Cor_P = numeric(8),
                        xxxxC_Cor   = numeric(8),
                        xxxxC_Cor_P = numeric(8),
                        xxxxD_Cor   = numeric(8),
                        xxxxD_Cor_P = numeric(8))

#for loop to run analysis for each treatment condition
for (Treatment in c("Ctrl", "xxxxA", "xxxxB", "xxxxC", "xxxxD")) {
  
  #calculates gene significance for each iteration
  geneTraitSignificance <- as.name((paste0("geneTraitSignificance", Treatment)))
  
  #subsets the relevant genes for each module colour
  for (i in c(1:length(colorlevels))) {
    whichmodule=colorlevels[[i]];
    
    restrict1 = (colorh1==whichmodule);
  
    #calculate correlation value
  IMMvsGS_cor <- signif(cor(unlist(abs(eval(geneTraitSignificance)))[restrict1], 
                            Alldegrees1$kWithin[restrict1]),2)
   #calculate p value
  IMMvsGS_corp <- signif(corPvalueStudent(IMMvsGS_cor, 
                                          sum(eval(geneTraitSignificance)[colorh1==colorlevels[[i]],]
                                              & Alldegrees1$kWithin[colorh1==colorlevels[[i]]])),2)
  
  #write out value for each module to relevant treatment condition columns
  cor_table[,paste0(Treatment, "_Cor")][i] <- IMMvsGS_cor
  
  cor_table[,paste0(Treatment, "_Cor_P")][i] <- IMMvsGS_corp
  }
  
}

#write out to csv
write.csv(cor_table, row.names = F,
          here("data.csv"))
# Calculate SignedKME (Module Membership) ---------------------------------
#recalculate MM

#get Module Eigengene
datME <- moduleEigengenes(datExpr_t, colorh1)$eigengenes

#calculate connectivity of each gene to all module eigengenes (defining connectivity
#to modular genes across all genes)
datKME=signedKME(datExpr_t, datME, outputColumnName="MM.")

# Identifying Hub Genes through Filtering by GS & MM ----------------------
## this identifies hub genes by filtering via gene significance and module membership
## and retrieves gene information from biomart

library(biomaRt)


#get mouse biomart
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

#retrieve gene names of color genes
color_genes <- as.character(filter(datExpr_orig, Modules == "color")$Ensemble_IDs)

#retrieve gene information from biomart of all color genes
color_genes_description <- getBM(attributes = c("description", "mgi_symbol", "mgi_id", "ensembl_gene_id"),
                                filters = "ensembl_gene_id",
                                values = color_genes,
                                mart = mouse)

#write out to csv
write.csv(color_genes_description, 
          file = here("Output/Bicor/Intramodular_Gene_Significance_Analysis/data.csv"))


#filter genes based on minimum value of gene significance and 
#module membership for each treatment condition

#set cutoff values
GS.cutoff <- 0.7
MM.cutoff <- 0.8

#create dataframe indicating whether a gene was a hub gene in each treatment condition
HubGeneStatus <- data.frame("Hub_Gene_Status_Ctrl"= as.vector(abs(geneTraitSignificanceCtrl)> GS.cutoff   & abs(datKME$MM.color)>MM.cutoff),
                            "Hub_Gene_Status_xxxxA"= as.vector(abs(geneTraitSignificancexxxxA)> GS.cutoff & abs(datKME$MM.color)>MM.cutoff),
                            "Hub_Gene_Status_xxxxB"= as.vector(abs(geneTraitSignificancexxxxB)> GS.cutoff & abs(datKME$MM.color)>MM.cutoff),
                            "Hub_Gene_Status_xxxxC"= as.vector(abs(geneTraitSignificancexxxxC)> GS.cutoff & abs(datKME$MM.color)>MM.cutoff),
                            "Hub_Gene_Status_xxxxD"= as.vector(abs(geneTraitSignificancexxxxD)> GS.cutoff & abs(datKME$MM.color)>MM.cutoff))




#check to see how many hub genes were identified in each treatment condition
#aiming for 5-10 genes ideally
sum(HubGeneStatus$Hub_Gene_Status_Ctrl)   #3
sum(HubGeneStatus$Hub_Gene_Status_xxxxA) #0
sum(HubGeneStatus$Hub_Gene_Status_xxxxB) #9
sum(HubGeneStatus$Hub_Gene_Status_xxxxC) #0
sum(HubGeneStatus$Hub_Gene_Status_xxxxD) #0


# Module Membership vs Intramodular Connectivity --------------------------
## create a scatter plot to visualize module membership vs intramodular connectivity
## ideally, in conditions of interest, we would see a positive correlation

#for loop to run the same analysis for each treatment condition
for (Treatment in c("Ctrl", "xxxxA", "xxxxB", "xxxxC", "xxxxD")) {
  
  #create unique file for each treatment condition
  pdf(here(paste0("Output/Bicor/Intramodular_Gene_Significance_Analysis/",
                 Treatment,
                 "/",
                 Treatment,
                 "_ModuleMembership_vs_IntramodularConnectivity.pdf")))
  
  par(mfrow=c(2,2))
  # We choose 4 modules to plot: turquoise, blue, brown, green.
  # For simplicity we write the code out explicitly for each module.
  which.color="color";
  restrictGenes=colorh1==which.color
  verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                     (datKME[restrictGenes, paste0("MM.", which.color)])^6,
                     col=which.color,
                     xlab="Intramodular Connectivity",
                     ylab="(Module Membership)^6")
  which.color="blue";
  restrictGenes=colorh1==which.color
  verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                     (datKME[restrictGenes, paste0("MM.", which.color)])^6,
                     col=which.color,
                     xlab="Intramodular Connectivity",
                     ylab="(Module Membership)^6")
  which.color="brown";
  restrictGenes=colorh1==which.color
  verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                     (datKME[restrictGenes, paste0("MM.", which.color)])^6,
                     col=which.color,
                     xlab="Intramodular Connectivity",
                     ylab="(Module Membership)^6")
  
  which.color="turquoise";
  restrictGenes=colorh1==which.color
  verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                     (datKME[restrictGenes, paste0("MM.", which.color)])^6,
                     col=which.color,
                     xlab="Intramodular Connectivity",
                     ylab="(Module Membership)^6")
  
  dev.off()
  
}

# Assembling New Spreadsheet ------------------------------------------------
## here we assemble the new spreadsheet containing information on the color genes,
## and their corresponding gene significance values, module membership, and hub gene status
geneTraitSignificance_df <- data.frame("Gene_SignificanceCtrl" = abs(geneTraitSignificanceCtrl),
                                       "Gene_SignificancexxxxA"= abs(geneTraitSignificancexxxxA),
                                       "Gene_SignificancexxxxB" = abs(geneTraitSignificancexxxxB),
                                       "Gene_SignificancexxxxC"= abs(geneTraitSignificancexxxxC),
                                       "Gene_SignificancexxxxD" = abs(geneTraitSignificancexxxxD))

#set column names
colnames(geneTraitSignificance_df) <- c("Gene_SignificanceCtrl", 
                                        "Gene_SignificancexxxxA",
                                        "Gene_SignificancexxxxB",
                                        "Gene_SignificancexxxxC",
                                        "Gene_SignificancexxxxD")



#create filters for color genes to only have the color genes in the final spreadsheet
color_filter <- rownames(geneTraitSignificance_df) %in% color_genes_description$ensembl_gene_id

color_filter_imc <- rownames(Alldegrees1) %in% color_genes_description$ensembl_gene_id

color_filter_mm <- rownames(datKME) %in% color_genes_description$ensembl_gene_id

#assemble final spreadsheet
color_genes_complete_df <- cbind(color_genes_description, 
                                HubGeneStatus[color_filter,],
                                geneTraitSignificance_df[color_filter,],
                                "Intramodular_Connectivity" = Alldegrees1$kWithin[color_filter_imc],
                                datKME[color_filter_mm,])


#write out to csv
write.csv(color_genes_complete_df, here("Output/Bicor/Intramodular_Gene_Significance_Analysis/color_Genes_with_All_Conditions.csv"))       

# Gene Significance vs Module Membership ggplot --------------------------

## more detailed scatter plot of gene significance vs module membership 
## to include hub gene information and filtering cutoffs

library(ggplot2)
library(ggrepel)
library(ggpubr)


# Brief Aside on ggplot List ----------------------------------------------


## ggplot reconstructs the plot everytime it is called. Because of this,
## we can't simply just assign the new plot as a variable using assign() and a for loop
## (with the x or y values of the plot varying based on treatment condition), 
## because when you later call it again, it will reconstruct the plot, calling the x and y values
## and it will only run based on the last treatment condition (i.e. the previous set up where
## the 'Treatment' variable changed in for loops and we used it in the x/y variables won't work).
## Technically we could bypass this by using aes_string() instead of aes() in ggplot, since a string won't cause
## gpplot to reconstruct the variable (I think), but best practice would be to store
## the plots in a list so they remain unchanged and to call them from plots later

## see here: https://stackoverflow.com/questions/45430325/assign-loop-not-working-for-plot

# Actual ggplot -----------------------------------------------------------

plot_list <- list()
#for loop to run the same analysis for each treatment condition
for (Treatment in c("Ctrl", "xxxxA", "xxxxB", "xxxxC", "xxxxD")) {
  
  #create gene significance call that varies with treatment condition
  Gene_Significance <- paste0("Gene_Significance", Treatment)
  
  #calculate correlation between gene significance and module membership
  GSvsMM_cor <- signif(cor(color_genes_complete_df[,Gene_Significance],
                        color_genes_complete_df$MM.color, use = "p"), 2)
  
  #calculate p value of correlation
  GSvsMM_cor_p <- signif(corPvalueStudent(GSvsMM_cor, sum(color_genes_complete_df[,Gene_Significance] & 
                                                     color_genes_complete_df$MM.color)), 2)
  
  #create detailed plot
  
       plotx <- ggplot(color_genes_complete_df, aes_string(x = "MM.color", y = paste0("Gene_Significance",Treatment))) + 
           geom_point(aes_string(color = paste("Hub_Gene_Status_", Treatment, sep = ""))) +
           scale_colour_manual(values=c("#A9A9A9","#d42a1e"))+
           geom_hline(aes(yintercept = GS.cutoff), linetype = "dashed") +
           geom_vline(aes(xintercept= MM.cutoff), linetype = "dashed") +
           geom_label_repel(data=filter(color_genes_complete_df, eval(parse(text=paste0("Hub_Gene_Status_", Treatment)))),
                            aes(label = mgi_symbol), size = 2, label.padding = 0.2)+
           theme_classic()+
           labs(x = "Module Membership in color Module",
                y = "Gene Significance",
                col = "Hub Gene Status",
                title = paste("Gene Significance vs Module Membership",Treatment),
                subtitle = paste0("cor = ", GSvsMM_cor, ",", " p = ", GSvsMM_cor_p))+
           theme(plot.title = element_text(face = "bold"),
                 plot.subtitle = element_text(face = "italic"))
         
       plot_name <- paste0(Treatment, "_ggplot")
       plot_list[[plot_name]] <- plotx
         
       #we could use assign here if we used aes_string(), but best to use the list instead
       #assign(x= paste0(Treatment, "_ggplot"), value = plotx, envir = globalenv())
  
  #save plot
  ggsave(here(paste0("Output/Bicor/Intramodular_Gene_Significance_Analysis/",
                    Treatment,"/",Treatment,
                    "_GeneSignificance_vs_ModuleMembership_HubGenes.pdf")),
         plot = plot_list[[paste0(Treatment, "_ggplot")]],
         width = 9,height = 5)
}

#collate figures into main figure
main_hub_figure <- ggarrange(plot_list$Ctrl_ggplot, 
                             plot_list$xxxxB_ggplot,
                             ncol = 2, nrow = 1,
                             common.legend = T,
                             align = "h",
                             labels = "AUTO")

#save to pdf
ggsave(plot = main_hub_figure, 
       filename = here("Output/Bicor/Intramodular_Gene_Significance_Analysis/Figure.pdf"), 
       width = 11, height = 5)


#run the same thing but just comment out the labels in the ggplot script for supplemental
supp_hub_figure <- ggarrange(plot_list$Ctrl_ggplot, 
                plot_list$xxxxA_ggplot, 
                plot_list$xxxxB_ggplot,
                plot_list$xxxxC_ggplot,
                plot_list$xxxxD_ggplot,
                ncol = 2, nrow = 3,
                common.legend = T,
                align = "h",
                labels = "AUTO")


#save to pdf
ggsave(plot = supp_hub_figure, 
       filename = here("Output/Bicor/Intramodular_Gene_Significance_Analysis/Figure.pdf"), 
       width = 11, height = 11)

dev.off()

