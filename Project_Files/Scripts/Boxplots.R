# Boxplot Script Comments ----------------------------------------------------------
# this is the script for making boxplots


# Load Required Packages and Datasets -------------------------------------

library(here)
library(cowplot)
library(tidyverse)
datTrait <- readxl::read_xlsx(here("Datasets/VarX_DNARNAinfo.xlsx"))
datExpr <- read.csv(here("Output/Bicor/Log2RPKM_Module_Heatmap/Log2RPKM_With_Modules.csv"))


# Wrangling Dataset -------------------------------------------------------
# we first have to wrangle the dataset into a tidier format that is amenable
# to working with ggplot, and making boxplots

#select only color modules
datExpr_color <- datExpr %>% 
    filter(Modules == "color")

#set rownames to ensemble ids
rownames(datExpr_color) <- datExpr_color$Ensemble_IDs

#remove the unnecessary columns
datExpr_color <- select(datExpr_color, -1, -Modules, -Ensemble_IDs)

##transpose dataframe for ease of access
datExpr_color <- as.data.frame(t(datExpr_color))

#get sample IDs
sample_ids <- rownames(datExpr_color)

#add sample IDs to dataset
datExpr_color <- cbind(sample_ids, datExpr_color)

#create matching vector of sample IDs
matching_vector <- datTrait$`Samlple ID` %in% sample_ids

#match traits data to the Log2RPKM dataset
datTrait_ordered <- datTrait[matching_vector, c("SampleName", "Var1", "Var2", "Var3")]

#add trait information to original dataset
datExpr_color <- cbind(datTrait_ordered, datExpr_color)


# Tidying Dataset -----------------------------------------------------------

# here, we use gather to tidy the dataset. This makes it longer instead of wider,
# ensuring each specific case is in one row
tidy_df <- datExpr_color %>% 
  gather(key = "Ensemble_ID", value = "Log2RPKM", -c(1:5))

# here, we specify that the column Var1 should be treated as a factor, and we
# list the levels in the desired order. This is necessary as without it CTRL
# goes last, whereas we want it to be first (since it is our control)
# https://community.rstudio.com/t/how-to-manually-order-x-axis-on-bar-chart/9601/2
tidy_df$Var1 <- factor(tidy_df$Var1, 
                            levels = c("CTRL", "AVarX", "BVarX", "CVarX", "DVarX"))

#this was the old grouping used
#grouped_samples_df <- group_by(tidy_df,  SampleName)

#we next create a grouped dataframe via sample, Var2, and Var1
grouped_samples_df <- group_by(tidy_df,  SampleName, `Var2`, Var1)

#summarize the values to calculate the mean log2RPKM values across all genes for each sample
summarized_samples_df <- summarise(grouped_samples_df, "Mean_Log2RPKM" = mean(`Log2RPKM`))

#write this out to a spreadsheet for use in other work
write.csv(summarized_samples_df, here("Datasets/data.csv"), )

#this is old code for the old grouping. Delete if not used.
# summarized_samples_df <- cbind(Var1 = c("AVarX","AVarX","AVarX","AVarX",
#                               "BVarX","BVarX","BVarX","BVarX",
#                               "CVarX","CVarX","CVarX","CVarX",
#                               "DVarX","DVarX","DVarX","DVarX",
#                               "CTRL","CTRL","CTRL","CTRL"), summarized_samples_df)


#again, we specify the order of the factor here to make sure CTRL goes first
summarized_samples_df$Var1 <- factor(summarized_samples_df$Var1, 
                            levels = c("CTRL", "AVarX", "BVarX", "CVarX", "DVarX"))

# Boxplot for All Genes in All Var2s ------------------------------------------------------

pdf(here("/Output/Bicor/Boxplots/data.pdf"))

#create a color blind friendly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#creates the boxplot for all the genes in all Var2s
ggplot(tidy_df, aes(x = Var1, y = `Log2RPKM`, group = Var1)) +
  stat_summary(geom = "boxplot", 
               fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
               position = "dodge", aes(fill=Var1)) +
  scale_fill_manual(values = cbPalette)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(name = "VarX Treatment") +
  scale_y_continuous(name = "Log2RPKM")+
 theme_cowplot(font_size = 12) +
  ggtitle("Log2RPKM of All Genes in All Areas Across VarX Treatment")



# Boxplot for Averaged Genes ---------------------------------------------------
#this is the boxplot script for all the genes averaged across samples

ggplot(summarized_samples_df, aes(x = Var1, y = `Mean_Log2RPKM`, group = Var1)) +
  stat_summary(geom = "boxplot", 
               fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
               position = "dodge", aes(fill=Var1)) +
  geom_point()+
  scale_fill_manual(values = cbPalette)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(name = "VarX Treatment") +
  scale_y_continuous(name = "Log2RPKM")+
  theme_cowplot(font_size = 12) +
  ggtitle("Log2RPKM with Genes Averaged Across VarX Treatment")

dev.off()


# Boxplot Across Var2s and VarX Treatment --------------------------

#this is with all the genes but split by Var2
pdf(here("/Output/Bicor/Boxplots/data.pdf"), width = 15)

ggplot(tidy_df, aes(x = Var1, y = `Log2RPKM`, group = Var1)) +
  facet_wrap(~`Var2`, scales = "free", ncol=3)+
  stat_summary(geom = "boxplot", 
               fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
               position = "dodge", aes(fill=Var1)) +
  scale_fill_manual(values = cbPalette)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(name = "VarX Treatment") +
  scale_y_continuous(name = "Log2RPKM")+
  theme_cowplot(font_size = 12) +
  ggtitle("Log2RPKM of All Genes By Region & Across VarX Treatment")




# Boxplot for Averaged Genes Across Var2s -------------------------

#this is with the averaged genes but split by Var2
ggplot(summarized_samples_df, aes(x = Var1, y = `Mean_Log2RPKM`, group = Var1)) +
  facet_wrap(~`Var2`, scales = "free", ncol=3)+
  stat_summary(geom = "boxplot", 
               fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
               position = "dodge", aes(fill=Var1)) +
  geom_point()+
  scale_fill_manual(values = cbPalette)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(name = "VarX Treatment") +
  scale_y_continuous(name = "Log2RPKM")+
  theme_cowplot(font_size = 12) +
  ggtitle("Log2RPKM of Averaged Genes By Region & Across VarX Treatment")

dev.off()

