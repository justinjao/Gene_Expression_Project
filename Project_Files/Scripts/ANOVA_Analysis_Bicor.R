# ANOVA Analysis Script Comments ------------------------------------------

#This script contains code split up into the following sections:
# 1. ANOVA on the module eigengenes + posthoc
# 2. Heatmap for the ANOVA of module eigengenes
# 3. ANOVA on the log2RPKM Values
# 4. Posthoc for the log2RPKM


# Loading Required Packages & Datasets ----------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

library(nlme)
library(lsmeans)
library(here)


#load expression dataset
load(here("Datasets/RData_Files/Quantifying_Modules_and_ANOV_Prep/data.RData"))
load(here("Datasets/RData_Files/Data_Preprocessing/data.RData"))

#load new summarized data
summarized_samples_df <- read.csv(here("Datasets/data.csv"))

# load clinical trait data
traits <- read.csv(here("Datasets/data.csv"))


# Wrangling Dataset for ANOVA on Module Eigengenes -------------------------------------------------------


# match sample names and trait data
Samples <- rownames(filteredLog2rpkm_outliers_removed)
traitRows <- match(Samples, traits$SampleName)

#create new column specifying unique group name
traits <- cbind(traits, group = paste(traits$Var1_Condition, traits$Var2, traits$Var3., sep = "_"))

#subset data to only keep relevant columns
datTraits <- traits[traitRows, ]

#set it to a dataframe
datTraits <- as.data.frame(datTraits)

#set rownames
rownames(datTraits) <- datTraits$SampleName

datTraits$subject <- rownames(datTraits)
MEs$subject <- rownames(MEs)

#combine trait info with the module eigengene info, by sample
combined_trait_ME_df <- merge(datTraits, MEs, by = "subject")

#tidy data by gathering it
combined_trait_ME_df_2 <- combined_trait_ME_df %>% 
    gather(module, ME, 9:ncol(combined_trait_ME_df))

#get list of modules
mod <- unique(combined_trait_ME_df_2$module)



#write out for future use
write.csv(combined_trait_ME_df_2,here("Output/Bicor/ANOVA_Results/data.csv"))


# ANOVA + Posthoc on Module Eigengenes-------------------------------------------------------------------

#create empty output locations for anova and posthoc
ME_anova_out <- NULL
ME_posthoc <- NULL

#for loop to run ANOVA across all modules
for (i in mod) {
  print(i)
  
  #select data: for each loop, it runs across only 1 module
  ME_data <- combined_trait_ME_df_2 %>% dplyr::filter(module == i)
  
  #setup the linear model
  ME_lme <- lme(ME ~ Var1*Var2+Var3., ~1|subject, data = ME_data)
  #run ANOVA on linear model
  ME_anova <- anova(ME_lme, type = "marginal") #marginal gives Type 3 SS for ANOVA
  
  #set as dataframe
  ME_anova <- as.data.frame(ME_anova)
  #set condition column to be the rownames
  ME_anova$condition <- rownames(ME_anova)
  #remove the first column
  ME_anova <- ME_anova[-1,]
  
  #wrangle dataset to final format
  ME_anova2 <- ME_anova %>% gather(stats,values,1:4) %>%
    mutate(col = paste(condition,stats, sep="_")) %>%
    dplyr::select(col,values) %>%
    spread(col,values)
  
  #paste current module into "module column"
  ME_anova2$module <- paste(i)
  
  #add new case to the output location
  ME_anova_out <- rbind(ME_anova_out,ME_anova2)
  

## POSTHOC starts here ##
  #create a reference grid model object
  ME_refgrid <- ref.grid(ME_lme)
  
  #creates a fitted model using the reference grid, based on Var1 and Var2 interactions
  ME_lsmeans <- lsmeans(ME_refgrid, ~Var1|Var2)
 
  #summarize paired comparisons
  ME_lsmeans_summary <- summary(pairs(ME_lsmeans, adjust = "none"))
  ME_lsmeans_summary <- as.data.frame(ME_lsmeans_summary)
  
  #posthocHSD = contrast(refgrid, method = "pairwise",adjust = "none")
  #outsum <- as.data.frame(summary(posthocHSD))
  
  
  #get only control vs Var999
  # outsum <- outsum[grepl("CONTROL", outsum$contrast),]
  
  #correct for multiple comparisons
  ME_lsmeans_summary$padjust <- p.adjust(ME_lsmeans_summary$p.value, method = "BH")
  
  #paste module information
  ME_lsmeans_summary$module <- paste(i)
  
  #save to output
  ME_posthoc <- rbind(ME_posthoc,ME_lsmeans_summary)
  
  
}

#write to file
write.csv(ME_anova_out, here("Output/Bicor/ANOVA_Results/data.csv"))
write.csv(ME_posthoc, here("Output/Bicor/ANOVA_Results/data.csv"))





# ANOVA with Mean Log2RPKM Values ----------------------------------------------
#this is the ANOVA for the log2RPKM values. Since we use the summarized data,
#which only uses info from the pink module with the average log2RPKM values across all genes,
#we don't need to loop through any modules, so this is just a single ANOVA

# run lme model with Var2 and condition as interaction effets
MeanLog2RPKM_lme <- lme(Mean_Log2RPKM ~ Var2*Var1, ~1|SampleName, data = summarized_samples_df)

#run ANOVA
MeanLog2RPKM_anova <- anova(MeanLog2RPKM_lme, type = "marginal") #marginal gives Type 3 SS for ANOVA

#reformat data
MeanLog2RPKM_anova <- as.data.frame(MeanLog2RPKM_anova)



#write anova to output file
write.csv(MeanLog2RPKM_anova, here("Output/Bicor/ANOVA_Results/Log2RPKM_Averages_ANOVA/data.csv"))


# Posthoc for Mean Log2RPKM Values ----------------------------------------------------
#runs the posthoc analysis. Similar to previous posthoc for the ME data

#create reference grid
MeanLog2RPKM_refgrid <- ref.grid(MeanLog2RPKM_lme)

#create fitted model from reference grid
MeanLog2RPKM_lsmeans <- lsmeans(MeanLog2RPKM_refgrid, ~Var2|Var3)

#summarize paired comparisons
MeanLog2RPKM_lsmeans_summary <- summary(pairs(MeanLog2RPKM_lsmeans, adjust = "none"))
MeanLog2RPKM_lsmeans_summary <- as.data.frame(MeanLog2RPKM_lsmeans_summary)

#posthocHSD = contrast(refgrid, method = "pairwise",adjust = "none")
#outsum <- as.data.frame(summary(posthocHSD))


#get only control vs Var999
# outsum <- outsum[grepl("CONTROL", outsum$contrast),]

#correct for multiple comparisons
MeanLog2RPKM_lsmeans_summary$padjust <- p.adjust(MeanLog2RPKM_lsmeans_summary$p.value, method = "BH")

write.csv(MeanLog2RPKM_lsmeans_summary, here("Output/Bicor/ANOVA_Results/Log2RPKM_Averages_ANOVA/data.csv"))


# PHEATMAP GENERATION -----------------------------------------------------

library(pheatmap)

#select relevant columns
anova_heatmap_df <- ME_anova_out %>% 
    select(`Var3._p-value`,
           `Var2_p-value`,
           `Var1_p-value`,
           `Var1:Var2_p-value`)

#set rownames
rownames(anova_heatmap_df) <- ME_anova_out$module

#rename columns for readability
colnames(anova_heatmap_df) <- c("Var3", "Var2", "Var1", "Var1 X Var2")

#create tiny constant to add to ensure values are nonzero
constant <- 0.00000000000000001

#add constant
non_zero_anova_heatmap_df <- anova_heatmap_df + constant

#generate heatmap, plotting the -log10 values for colour, and plotting the actual p values
#the effect is that the lower p values have higher colour expression (to convey greater significance)
pdf(file = here("Output/Bicor/ANOVA_Results/Heatmap_of_ANOVA_P_Values_(non-scaled).pdf"), width = 6)
pheatmap(-log10(non_zero_anova_heatmap_df),
         cluster_cols = F,
         scale = "none", display_numbers = signif(anova_heatmap_df, digits = 3),
         main = "ANOVA of Gene Modules Across Samples (P Values)",
        annotation_names_col = T,
        number_color = "black",
        fontsize_number = 10,
        angle_col = 45)
dev.off()

