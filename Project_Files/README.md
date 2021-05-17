# Gene Expression Project

Justin Jao's contributions to the Gene Expression Project.

## General Notes

- With the exception of the Network_Construction.R scripts and bash scripts (which were written to be be executed on a cluster), all code can be run on a local machine.

- The repository is set up with the **R Project** and **here** package workflow in mind. As such, if the entire repository is cloned, one need only to open the R Project titled **"Gene_Expression_Project.Rproj"**, and any script can then be opened and executed without needing to change any directories or import any datasets. Therefore, when executing any code, one should always start by opening the project file.

- **Do not move/delete any files or folders unless the appropriate scripts have been changed**. To enable the analysis to run without having to edit any code, each script calls on the required datasets and saves its output in specific folders. If the directories are modified (by moving or deleting any files or folders) without editing the affected scripts, you could break a script and it might not run.


- Raw datasets and RData files for use in the analysis can be found in *Datasets*.
  - If you're trying to figure out how a dataset/RData file was made, follow along the script workflow until you find the dataset you're looking for. The folders the datasets are saved in usually correspond to a script name as well.


- Results from the analysis can be found in *Output*.
  - Note: the Pearson folder isn't fully updated with the results, as we eventually stopped using the pearson data in our analysis.
  - Results are saved in individual folders, partitioned according to each step of the analysis.
  - If you rerun any scripts which save files, they'll automatically overwrite the files (as they are all set to save in specific locations) unless you change the location (e.g. by changing the saved filename) in which the file is written to in the script.


- All the scripts used can be found in *Scripts*.
  - Each step of the analysis has been split up into different scripts. To run a script, open the **"Gene_Expression_Project.Rproj"** file, and then open the script desired. This will ensure that the scripts can call the necessary files with the right file path.


- Any writings about the project are in *Documentation*

## Current Script Workflow
In order of execution, the analysis is currently run as follows:

- **Revised_Dataprep.R**
  - Initial data preprocessing, such as filtering genes, converting to log2RPKM, etc.
- **Removing Outliers.R**
  - Removes sample outliers in the expression dataset.
- **Network_Construction_bicor.R**
  - (run on cluster) constructs the network using bicorrelation settings.
- **bicor_v3.sh**
  - (run on cluster) is the bash script to be executed together with the Network_Construction_bicor.R
- **Quantifying_Modules_and_ANOVA_Prep_Bicor.R**
  - Ensures modules are appropriate size (merging if needed), and prepping the data for ANOVA analysis, as well as producing sample clustering graphs by the 3 different variable groupings in the experiment (Var1, Var2, Var3).
- **ANOVA_Analysis_Bicor.R**
  - Runs ANOVA analysis and posthoc comparisons.
- **Log2RPKM_Module_Heatmap.R**
  - Generates the Log2RPKM heatmap for the different modules, without scaling.
- **Only_Color1_Module_Heatmap.R**
  - Generates the Log2RPKM heatmap but only for the Color1 modules, to increase readability.
- **Creating_Custom_Enrichment_Collection.R**
  - Converts csv of manually curated microglia relevant genes (compiled via previous work) and converts it into a gene set collection via the anRichment R package, for enrichment analysis.
- **GOEnrichmentAnalysis.R**
  - Performs enrichment analysis for different collections, including the custom curated gene set.
- **Boxplots.R**
  - Creates boxplots visualizing both log2RPKM(+025) values and averaged log2RPKM(+0.25) values across Var1, split by Var2.
- **Intramodular_Hub_Genes.R**
  - Runs the hub gene analysis, calculating gene significance, module membership, intramodular connectivity and identifying hub genes for each Var1. Also makes plots showcasing these stats, and a table.


If you have any further questions about the code or project, you may contact Justin Jao (jao.c.justin@gmail.com)
