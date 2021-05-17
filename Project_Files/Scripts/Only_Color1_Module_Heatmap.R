# Load Required Packages & Datasets --------------------------------------------------

library(here)
library(tidyverse)
library(WGCNA)
library(pheatmap)
library(RColorBrewer)


#load in expression dataset
load(here("Datasets/RData_Files/Data_Preprocessing/data1.RData"))
load(here("Datasets/RData_Files/Quantifying_Modules_and_ANOV_Prep/data2.RData"))

#rename for readability
datExpr <- filteredLog2rpkm_outliers_removed

#get ensemble data
ensemble <- colnames(datExpr)

#transpose matrix
datExpr <- t(datExpr)

#convert to dataframe for ease of working
datExpr <-as.data.frame(datExpr)
#bind module colours to dataframe
datExpr <- cbind(Ensemble_IDs = ensemble,
                 Modules = merge$colors,
                 datExpr)


#remove grey module genes and order genes by modules
datExpr_grey_rm <- datExpr %>% 
  filter(Modules != "grey") %>% 
  arrange(Modules)

datExpr_grey_rm$Modules <- factor(datExpr_grey_rm$Modules)

rownames(datExpr_grey_rm) <- datExpr_grey_rm$Ensemble_IDs

datExpr_only_color1 <- datExpr %>% 
  filter(Modules == "color1")

# Only color1 Module --------------------------------------------------------

order <- data.frame(colnames(datExpr_only_color1[,3:ncol(datExpr_only_color1)]))

#add info from clinical dataset
order <- cbind(order, 
               Var1 = datTraits$Var1.,
               Var2 = datTraits$Var2)

colnames(order) <- c("Name", "Var1", "Var2")

#remove AC and digits from name, and treatment info
order_df <- order %>% 
  mutate(Name = str_sub(Name, start = 3, end = str_length(Name)),
         Var1 = str_replace(Var1, "8", "left"),
         Var1 = str_replace(Var1, "9", "right"),
         Var3 = gsub("[[:digit:]]", "", Name),
         Var2 = substr(Var2, start = 1, 
                                      stop = nchar(as.character(Var2)) - 2))


#convert treatment conditions to factor
order_df$Var2 <- as.factor(order_df$Var2)





# Heatmap Generation ------------------------------------------------------

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)

annotation_col <- order_df

rownames(annotation_col) <- order_df$Name

annotation_row <- data.frame(datExpr_only_color1$Modules)

rownames(annotation_row) <- rownames(datExpr_only_color1)

colnames(annotation_row) <- "Module"



qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

colfunc_blue <- colorRampPalette(c("lightblue", "darkblue", "cyan"))

colfunc_yellow <-colorRampPalette(c("yellow2","goldenrod"))

Var1 <- c("blue","turquoise2","darkgoldenrod2", "color1","purple")

names(Var1) <- unique(annotation_col$Var2)



# Var2
Var2 <- colfunc_blue(length(unique(annotation_col$Var3)))

names(Var2) <- unique(annotation_col$Var3)

# Modules
Var3 <- as.character(unique(datExpr_only_color1$Modules))
names(Var3) <-  unique(datExpr_only_color1$Modules)


# Var1
Var4 <- colfunc_yellow(length(unique(annotation_col$Var1)))
names(Var4) <- unique(annotation_col$Var1)

anno_colors <- list(Var1 = Var4,
                    Var2 = Var1, 
                    Var3 = Var2, 
                    Module = Var3)

datExpr_heatmap_df <- subset(datExpr_only_color1, 
                             select = -c(Ensemble_IDs, Modules))
annotation_col <- subset(annotation_col,
                         select= -c(Name))


colnames(datExpr_heatmap_df) <- gsub("AC", "", colnames(datExpr_heatmap_df))



##only color1
pdf(file = here("Output/Bicor/Log2RPKM_Module_Heatmap/data.pdf"))
pheatmap(datExpr_heatmap_df,
         cluster_row = F,
         cluster_cols = F,
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         show_rownames = F,
         show_colnames = T,
         annotation_names_row=T,
         annotation_names_col = T,
         color = my_palette, 
         fontsize = 10,
         fontsize_row=6, 
         fontsize_col = 6,
         annotation_colors = anno_colors, 
         scale = c("row"),
         main = "Gene Expression Modules in Brain")

dev.off()

