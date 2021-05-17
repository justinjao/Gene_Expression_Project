# Load Required Packages & Datasets --------------------------------------------------

library(here)
library(tidyverse)
library(WGCNA)
library(pheatmap)
library(RColorBrewer)


#load in expression dataset
load(here("Datasets/RData_Files/Data_Preprocessing/data.RData"))

load(here("Datasets/RData_Files/Quantifying_Modules_and_ANOV_Prep/data.RData"))

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

datExpr_all <- arrange(datExpr, Modules)

write.csv(datExpr_all, here("Output/Bicor/Log2RPKM_Module_Heatmap/data.csv"))



#remove grey module genes and order genes by modules
datExpr_grey_rm <- datExpr %>% 
    filter(Modules != "grey") %>% 
    arrange(Modules)

datExpr_grey_rm$Modules <- factor(datExpr_grey_rm$Modules)

rownames(datExpr_grey_rm) <- datExpr_grey_rm$Ensemble_IDs

write.csv(datExpr_grey_rm, 
          file = here("Output/Bicor/Log2RPKM_Module_Heatmap/data.csv"))

datExpr_only_color <- datExpr %>% 
    filter(Modules == "color")

# Creating Levels for the Heatmap -----------------------------------------


#These are the levels, for reference
# treatment_levels <- c("XXXA1",
#                       "XXXB",
#                       "XXXC",
#                       "XXXD",
#                       "XXXE")
# 
# SSS_levels <- c("AAA",
#                 "BBB",
#                 "CCC")
# 
# GGG_levels <- c("DDD",
#                 "EEE")



#trial <- separate(datTraits, group, c("Treatment", "No", "SSS", "GGG"))

order <- data.frame(colnames(datExpr_grey_rm[,3:ncol(datExpr_grey_rm)]))

#add info from clinical dataset
order <- cbind(order, 
               GGG = datTraits$GGG,
               Treatment_Condition = datTraits$Treatment_Condition)

colnames(order) <- c("Name", "GGG", "Treatment_Condition")

#remove AC and digits from name, and treatment info
order_df <- order %>% 
  mutate(Name = str_sub(Name, start = 3, end = str_length(Name)),
         SSS = str_replace(SSS, "8", "left"),
         SSS = str_replace(SSS, "9", "right"),
         GGG = gsub("[[:digit:]]", "", Name),
         GGG = substr(Treatment_Condition, start = 1, 
                                     stop = nchar(as.character(Treatment_Condition)) - 2))


#convert treatment conditions to factor
order_df$Treatment_Condition <- as.factor(order_df$Treatment_Condition)





# Heatmap Generation ------------------------------------------------------

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)

annotation_col <- order_df

rownames(annotation_col) <- order_df$Name

annotation_row <- data.frame(datExpr_grey_rm$Modules)

rownames(annotation_row) <- rownames(datExpr_grey_rm)

colnames(annotation_row) <- "Module"



qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

colfunc_blue <- colorRampPalette(c("lightblue", "darkblue", "cyan"))

colfunc_yellow <-colorRampPalette(c("yellow2","goldenrod"))

Var1 <- c("blue","turquoise2","darkgoldenrod2", "deepcolor2","purple")

names(Var1) <- unique(annotation_col$Treatment_Condition)



#GGG section
Var2 <- colfunc_blue(length(unique(annotation_col$GGG)))

names(Var2) <- unique(annotation_col$GGG)

#YYY
Var3 <- as.character(unique(datExpr_grey_rm$YYY))
names(Var3) <-  unique(datExpr_grey_rm$YYY)


# SSS 
Var4 <- colfunc_yellow(length(unique(annotation_col$SSS)))
names(Var4) <- unique(annotation_col$SSS)

anno_colors <- list(SSS = Var4,
                    Treatment_Condition = Var1, 
                    GGG = Var2, 
                    YYY = Var3)

datExpr_heatmap_df <- subset(datExpr_grey_rm, 
                             select = -c(Ensemble_IDs, Modules))
annotation_col <- subset(annotation_col,
                         select= -c(Name))


colnames(datExpr_heatmap_df) <- gsub("AC", "", colnames(datExpr_heatmap_df))

pdf(file = here("Output/Bicor/Log2RPKM_Module_Heatmap/data.pdf"), 
    wi = 8, he = 10)

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
         fontsize_col = 8,
         annotation_colors = anno_colors, 
         scale = c("row"),
         main = "Gene Expression Modules in Brain")
dev.off()



save(datExpr_grey_rm, datExpr_only_color, 
     file = here("Datasets/RData_Files/Log2RPKM_Module_Heatmap/data.RData"))



