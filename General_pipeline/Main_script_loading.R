## Main script loading is the main script for analysis of LiP-MS thermal denaturation data
## The script uploads the Spectronaut data and transforms it for python fitting and exports it.
## Afterwards, the data has to be fit using python. After the fitting, the data is again loaded to the same script 
## The script then determines the effect of small molecules on binding, aggregation and thermal stability

# Load in the required packages
library(plotly)
library(ggplot2)
library(magrittr)
library(dplyr)
library(data.table)
library(ggpubr)
library(spatstat)
library("rstudioapi")
library(seqinr)
library(reldist)

# Additional package requirements:
# Important: while functions from reshape2 and pyr are needed, you should not load the two packages, as the scripts do not work well otherwise!!
# fclust, reshape2, plyr, kit, tidyr, stringr

## Set the working directory to the location of the script
setwd(dirname(getActiveDocumentContext()$path))

## Load in the supporting functions
source("Spectronaut_to_python_function.R")
source("find_significant_peptides.R")
source("python_to_list_function.R")
source("fuzzy_clustering.R")
source("peptide_fasta_matching.R")
source("aggregate_scores.R")
source("area_calculation.R")

### Import Spectronaut search data and transform them to python-input-format
path_file <- read.csv("~/Documents/Phd/Collaborations/110/Path_directions.csv")

fasta_ecoli <- read.fasta("~/Documents/Phd/Collaborations/110/200320_human_PK_iRT.fasta", "AA",
                          as.string = TRUE, set.attributes = FALSE)

# Define small molecules that you wish to analyse
all_compounds <- c("ATP")

# Define the common path, where you would want to save your data for fitting
# You should then keep the same structure of files such that the data from python can be loaded directly
common_path <- "~/Documents/Phd/Collaborations/110/Fitting/"

### Import Spectronaut search data and min-max scale them to python-input-format
# Path file should include the experiment number and the small molecule
# The files will be saved in the folder: Common_path/Experiment_number/Small_molecule/python_import/Small_molecule_trypticity.csv
# Example for TMAO: Common_path/023/TMAO/python_import/TMAO_FT.csv
for(i in all_compounds){
  Spectronaut_to_python(i, path_file,common_path = common_path,fasta = "human")}

### Transform files from python #########
# Define the empty list so it can be used 
# Define path_file to connect the small molecules and experiment names
all_results <- list()
for(i in all_compounds){
  all_results <- python_to_list(defined_list = all_results, path_file = path_file, sm = i,common_path = common_path)
}

# Data can be saved as RDS - this is the file frequently used in the post-analysis 
#saveRDS(all_results, "/Volumes/PMONIKA/data/Combined_results_other.rds")

### Clustering analysis #########
# The function uses fuzzy clustering to cluster the data
# k = number of clusters, show_clusters = should you plot the clusters
# cutoff - cutoff for determining whether the cluster shows aggregation behaviour or not. The 1.75 was chosen by manual validation to best descriminate between unfolding profiles and aggregation profiles
clusters <- list()
for(i in all_compounds) {
  cluster <- fuzzy_clusters_two_conditions(all_results, i, control = "Control", k = 15, show_clusters = TRUE, cutoff = 1.75)
  clusters[[i]] <- cluster
}

#saveRDS(clusters,paste(common_path, "/clusters.rds", sep=""))

### Significant analysis ########
# Identify significantly stabilised peptides
significant_all <- list()

for(i in all_compounds) {

  significant <- find_significant_peptides(all_results, clusters = clusters,sm = i, control = "Control")
  significant_all[[i]] <- significant

}

#### Matching peptides to proteins + peptide positions #########
all_peptides <- lapply(significant_all, function(x) x$scores$Peptide) %>% unlist() %>% unique()

matched_df <- peptide_fasta_matching(all_peptides, fasta_ecoli)

#### Get AA-level score + summarise protein score ########

for(i in all_compounds) {
  significant_all <- score_aggregation(significant_all, sm = i, matched_df = matched_df, pep_level_quantile = 0.75, prot_level_quantile = 0.75)
}

#saveRDS(significant_all,file = paste(common_path, "/significant_all.rds", sep = ""))
