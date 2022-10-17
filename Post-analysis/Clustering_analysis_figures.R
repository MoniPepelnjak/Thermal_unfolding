# Analysis of cluster shapes
# Goal of the script: How does cluster shape change within the conditions
# Author: Monika Pepelnjak

# Load required scripts
source("/Users/moni/Documents/Phd/Scripts/Clustering_functions_two_conditions_2.R")
source("/Users/moni/Documents/Phd/Scripts/Clustering_functions_two_experiments.R")
source(("/Users/moni/Documents/Phd/Scripts/post_analysis_functions.R"))
source(("/Users/moni/Documents/Phd/Scripts/Non_precipitator_calculation_functions.R"))

all_results <- list()
all_results[["Trehalose"]] <- readRDS("/Volumes/My Passport for Mac/data/All_data_1M.rds")[["Trehalose"]]
TPP_all <- readRDS("/Volumes/My Passport for Mac/data/Python_TPP_fits_all.rds")
Savitski <- read.csv("/Users/moni/Documents/Phd/databases/Savitski_melting_temp.csv")

Precipitators <- define_precipitators_control(TPP_all,Savitski = Savitski, both = TRUE, cutoff=0.5)

# load libraries
library(magrittr)
library(dplyr)
library(ggplot2)

# Set the seed for reproducible clustering
set.seed(123)
### 
sm_all <- c("TMAO",  "Betaine", "Glycerol", "Proline", "Glucose", "Trehalose")

# Clustering of all the conditions takes a long time and is not recommended to do each time
# Therefore, this part can be run only once and the result saved. For subsequent plotting and analysis pre-calculated data can be uploaded

# output_cluster_all <- NULL
# for(sm in sm_all){
#   all_results <- list()
#   all_results[[sm]] <-  readRDS("/Volumes/My Passport for Mac/data/All_data_1M.rds")[[sm]]
#   
#   output_cluster <- cluster_analysis_all(all_results, c(sm), k=15)
#   
#   output_cluster_all <- rbind(output_cluster, output_cluster_all)
#   write.csv(output_cluster, paste("/Users/moni/Documents/Phd/Experiments/TD_analysis/", sm, "_05_cutoff_Changes_clusters.csv", sep=""))
# }

# Alternatively, if the files have been produced before, use this part of the code to just import the result of clustering
output_cluster_all <- NULL
for(sm in sm_all){
  output_cluster <- read.csv( paste("/Users/moni/Documents/Phd/Experiments/TD_analysis/", sm, "_Changes_clusters.csv", sep=""))
  output_cluster_all <- rbind(output_cluster, output_cluster_all)

}

# Count the changes in different conditions 
cluster_combined_relative <- output_cluster_all %>%
  unique() %>%
  group_by(Change, Experiment) %>% # Count the number of changing peptides
  dplyr::summarise(n=n()) %>%
  group_by(Experiment) %>%
  mutate(Percentage = n/sum(n)) %>% # Calculate the percentage
  ungroup() %>%
  mutate(Experiment = factor(Experiment, levels = sm_all)) %>% # Reorder osmolytes
  mutate(Change = factor(Change, levels = c("No change", "Two step start","Flip", "Two step stop"))) %>%
  na.omit() %>%
  mutate(Experiment = factor(Experiment))

# Plot the changes for figure S2A
Figure_S2A <- ggplot(cluster_combined_relative, aes(x=Experiment, y=Percentage, fill=Change)) +
  geom_bar(stat = "identity", col="black") +
  #  facet_wrap(~Experiment, scales="free") +
  scale_fill_manual(values = c("#F3EEE9", "#C4B8AF", "#7F7063", "#332922")) +
  theme_classic() +
  ylab("Percentage of peptides") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12, color="black"))

### Clustering 
# Perform the clustering for a single condition
# In the paper, trehalose condition was used as we have the highest number of good peptides

# Cluster the Trehalose dataset, into 20 clusters, only the Control condition
fuzzy_output <- fuzzy_clusters_one(all_results, "Trehalose", "Control", "Control", 20)
fuzzy_whole_plot <- fuzzy_output$whole_plot 

# Group the peptides into the ones with increasing / decreasing or two step profile
fuzzy_whole_plot$description <- ifelse(fuzzy_whole_plot$manual_order %in% c(1, 2, 3), "down", 
         ifelse(fuzzy_whole_plot$manual_order %in% c(4, 5, 6), "up", "two-step"))

# Count the peptides in each cluster
fuzzy_whole_count <- fuzzy_whole_plot %>%
  dplyr::select(Peptide_con, cluster) %>%
  unique() %>%
  group_by(cluster) %>%
  summarise(n = n())

# Order the clusters based on the number of peptides per cluster
fuzzy_whole_plot %<>%
  plyr::join(., fuzzy_whole_count, by="cluster") %>% 
  dplyr::arrange(-n) 

# This is hard-coded for exactly this experiment such that the clusters are nicely ordered for the plot
# This would have to be automated to fit other datasets!
reorder <- c("682", "884", "666", "730", "788", "819", "828", "627",  "710", "671", "504", "559", 
             "939", "650", "668", "866", "716", "864", "817", "740")
fuzzy_whole_plot$order <- as.factor(as.character(fuzzy_whole_plot$n))
fuzzy_whole_plot$order <-factor(as.character(fuzzy_whole_plot$n), levels = reorder)

# Plot the plot for Figure 1C
clusters_plot <- ggplot(fuzzy_whole_plot, aes(x=temp, y=y, col=description, alpha=`Membership degree`, group = Peptide_con)) +
    geom_line() +
    scale_alpha_continuous(range=c(0.00000001, 0.1)) +
  scale_color_manual(values = c("#9b9b7a", "#f4d2b9", "#B48F8F")) +
    facet_wrap(~order, strip.position = "top") +
  theme_minimal() +
  theme(text = element_text(size=15), 
        legend.position = "none", 
        strip.background = element_blank(), 
        panel.background = element_blank(), 
        panel.grid = element_blank(),
        axis.ticks = element_line(), 
        axis.title = element_text(size=18),
        #panel.background = element_rect(fill=NA)
        ) +
  # theme(axis.line = element_blank(), panel.grid.minor = element_blank(),
  #       panel.grid.major = element_blank(), legend.position = "none", 
  #       strip.background = element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) + 
  xlab(expression('Temperature ['*degree*C*']')) +
  ylab("Scaled intensity") +
  ylim(c(-0.2,1.2))

ggsave("/Users/moni/Documents/Phd/Osmolyte_paper/Figure1C_ordered.pdf", clusters_plot, height = 5, width = 6)

# Count the number of peptides per main cluster
main_plot <- fuzzy_output$data %>%
  mutate(description = ifelse(Cluster %in% c(1, 2, 3), "down",
                 ifelse(Cluster %in% c(4, 5, 6), "up", "two-step"))) %>%
  plyr::join(., Precipitators, by="Protein") %>% # Join with the data from the Savitski paper to analyse the precipitators
  group_by(trypticity, description,Precipitator_Savitski) %>%
  na.omit() %>% # Remove the proteins that are not in the Savitski dataset
  summarise(n = n()) %>% # Count the numbers per trypticity and precipitators
  ungroup() %>%
  group_by(trypticity, Precipitator_Savitski) %>%
  mutate(percentage = n/sum(n)) %>%
  na.omit()

main_plot$Precipitator <- ifelse(main_plot$Precipitator_Savitski == "Precipitator", "P", "NP")

# Plot the cluster distribution for precipitators and non-precipitators

Figure_1F <- ggplot(main_plot[main_plot$trypticity == "FT",], aes(x=percentage, y=Precipitator, fill=description)) +
  geom_bar(stat = "identity", col="black") +
  scale_fill_manual(values = c("#9b9b7a", "#f4d2b9", "#B48F8F")) +
  scale_x_continuous(labels = scales::percent)+
  theme_classic() +
  ylab("") +
  xlab("Percentage of peptides") +
  theme(text = element_text(size=20),
        axis.text = element_text(size=20),
        #axis.title.y=element_blank(), 
        legend.position = "none")
ggsave("/Users/moni/Documents/Phd/Osmolyte_paper/Figure1F.pdf", Figure_1F, height = 2, width = 8)

Figure_S1G <- ggplot(main_plot[main_plot$trypticity == "HT",], aes(x=percentage, y=Precipitator, fill=description)) +
  geom_bar(stat = "identity", col="black") +
  scale_fill_manual(values = c("#9b9b7a", "#f4d2b9", "#B48F8F")) +
  scale_x_continuous(labels = scales::percent)+
  theme_classic() +
  ylab("") +
  xlab("Percentage of peptides") +
  theme(text = element_text(size=20),
        axis.text = element_text(size=20),
        #axis.title.y=element_blank(), 
        legend.position = "none")
ggsave("/Users/moni/Documents/Phd/Osmolyte_paper/FigureS1G.pdf", Figure_1F, height = 2, width = 8)

# Calculate the enrichment for unfolding profile in the NP group for FT peptide
# Again, this is hardcoded for this exact dataset
matrix_fisher <- matrix(c(133,2593,(2+22), (2148 +1084)), ncol=2, byrow=TRUE)
fish <- fisher.test(matrix_fisher)
fish$p.value

# Combined distribution (all proteins) for HT and FT peptides
main_plot_2 <- fuzzy_output$data %>%
  mutate(description = ifelse(Cluster %in% c(1, 2, 3), "down",
                              ifelse(Cluster %in% c(4, 5, 6), "up", "two-step"))) %>%
  plyr::join(., Precipitators, by="Protein") %>%
  na.omit() %>%
  group_by(trypticity, description) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(trypticity) %>%
  mutate(percentage = n/sum(n)) %>%
  na.omit()

main_plot_2 %<>%
  ungroup() %>%
  mutate(sum_all = sum(n)) %>%
  group_by(description) %>%
  mutate(n_percentage = sum(n)/sum_all)

# Plot the figure 1D
Figure_1D <- ggplot(main_plot_2, aes(y=trypticity, x=percentage, fill=description)) +
  geom_bar(stat = "identity", col="black") +
  #facet_wrap(~Precipitator_Savitski, scales="free", ncol=1) +
  scale_fill_manual(values = c("#9b9b7a", "#f4d2b9", "#B48F8F")) +
  theme_minimal() +
  ylab("Trypticity") +
  xlab("") +
  scale_x_continuous(labels = scales::percent)+
  theme_classic() +
  ylab("") +
  xlab("Percentage of peptides") +
  theme(text = element_text(size=20), 
        axis.text = element_text(size=20),
        #axis.title.y=element_blank(), 
        legend.position = "none")
ggsave("/Users/moni/Documents/Phd/Osmolyte_paper/Figure1D.pdf", Figure_1D, height = 2, width = 8)

