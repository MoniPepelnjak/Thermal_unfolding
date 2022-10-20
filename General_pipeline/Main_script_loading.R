## Spectronaut to python main script 
## Short script that can 

library(plotly)
library(ggplot2)
library(magrittr)
library(dplyr)
library(data.table)
library(ggpubr)
library(spatstat)

source("Spectronaut_to_python_function.R")
source("post_analysis_functions.R")
source("Testing_differential_plots_functions.R")
source("Clustering_functions_two_conditions_2.R")
source("Significance_analysis_functions.R")
source("python_to_list_function.R")

### Import Spectronaut search data and transform them to python-input-format
path_file <- read.csv("~/polybox/SemPro_Picotti/001_ATP/Path_directions.csv")
all_compounds <- c("ATP")
common_path <- "~/polybox/SemPro_Picotti/001_ATP/"

### Import Spectronaut search data and min-max scale them to python-input-format
for(i in all_compounds){
  Spectronaut_to_python(i, path_file,common_path = common_path)}

### Transform files from python #########

all_results <- list()
for(i in all_compounds){
  all_results <- python_to_list(defined_list = all_results, path_file = path_file, sm = i,common_path = common_path)
}

saveRDS(all_results, "/Volumes/PMONIKA/data/Combined_results_other.rds")

### Clustering analysis #########

all_results <- readRDS("~/polybox/SemPro_Picotti/001_ATP/018/ATP/Combined_results_ATP.rds")[c("ATP")]

clusters <- list()
for(i in all_compounds) {
  cluster <- fuzzy_clusters_two_conditions(all_results, i, control = "Control", k = 15, show_clusters = TRUE, cutoff = 1.75)
  clusters[[i]] <- cluster
}

saveRDS(clusters,"~/polybox/SemPro_Picotti/001_ATP/018/ATP/211012clusters_ATP_k8.rds")

### Significant analysis ########

significant_all <- list()

for(i in all_compounds) {
  significant <- find_significant_peptides(all_results, clusters = clusters,sm = i, control = "Control")
  significant_all[[i]] <- significant
}

saveRDS(significant_all,file = "~/polybox/SemPro_Picotti/001_ATP/018/ATP/211210significant_allATP_old.rds")


#### Matching peptides to proteins + peptide positions #########
all_peptides <- lapply(significant_all, function(x) x$scores$Peptide) %>% unlist() %>% unique()

fasta_ecoli <- read.fasta("211129_Ecoli_proteome_fasta.fasta", "AA",
                          as.string = TRUE, set.attributes = FALSE)

matched_df <- peptide_fasta_matching(all_peptides, fasta_ecoli)

#### Get AA-level score + summarise protein score ########

for(i in all_compounds) {
  significant_all <- score_aggregation_both_quantiles(significant_all, sm = i, matched_df = matched_df, pep_level_quantile = 0.75, prot_level_quantile = 0.75)
}




