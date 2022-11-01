# Script for analysis of aggregation data and comparison of TPP and LiP results
# Author: Monika Pepelnjak
# Goal of the script: Goal of the script is to analyse the TPP and LiP data and compare the results
# This script was used to produce the figure 5 in the Osmolyte paper

# Data needed for the script:
  # * LiP-MS and TPP data - the output of the fitting as well as the post-processing analysis
  #     Look at the "Main script" to see how the data structures are generated from Spectronaut searches
  # * Combined calculated features for all proteins of the proteome
  # * Data from Savitski publication doi: 10.15252/msb.20188242., supplementary table 3
  # Aggrescan predictions
  # * Dataset of ribosomal proteins

library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path))

# Define the common path
path <- "/Volumes/My Passport for Mac/Final_data/"

# Load in the functions
# Calculate whether the protein precipitates based on our TPP and Savitski published dataset
source("Supportive_functions/Non_precipitator_calculation_functions.R")
# Mostly plotting functions for thermal denaturation data
source("Supportive_functions/post_analysis_functions.R")
# Script to compare the features between two protein groups
source("Supportive_functions/Structural_predictions_diff_test_functions.R")


# Load the libraries
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(ComplexHeatmap)

# Load in the necessary data:
# Load in the LiP-MS data
significant_all <- readRDS(paste(path, "LiP/Significant_both_quantiles.rds", sep=""))
all_results <- readRDS(paste(path, "LiP/All_data_1M.rds", sep=""))

# Load in ribosomal proteins to exclude them from the analysis
ribosomal_proteins <- read.csv(paste(path,  "Databases/ribosomal_proteins.csv", sep="")) %>%
  filter(grepl("ribosomal protein", Protein.names), 
         Status == "reviewed")

# Load in the TPP data
TPP_all <- readRDS(paste(path, "TPP/Python_TPP_fits_all.rds", sep=""))
TPP_significant <- readRDS(paste(path, "TPP/Significant_TPP.rds", sep=""))

# Load in the protein-feature predictions
All_combined_predictions <- read.csv(paste(path, "Databases/Ecoli_combined_predictions.csv", sep=""))
GO_codes <- read.csv(paste(path, "Databases/GO_terms_ecoli.csv", sep=""), stringsAsFactors=F)

# Load in the data from Savitski publication doi: 10.15252/msb.20188242., supplementary table 3
Savitski <- read.csv(paste(path, "Databases/Savitski_melting_temp.csv", sep=""))
Agrescan_list <-readRDS(paste(path, "Databases/Aggrescan_annotations.RDS", sep=""))

# Define colors used
# Color palette includes the names of the small molecule
my_palette = c(Control = "#8D8D92",
               control = "#8D8D92",
               Betaine = "#8FAF7E", 
               Trehalose = "#C18D73", 
               TMAO = "#678EAF", 
               Glycerol = "#997A8D", 
               Proline = "#DFBA82",
               proline = "#DFBA82", 
               glucose = "#BB7676",
               Glucose = "#BB7676",
               Proline0 = "#F9DEB7",
               Low_PK = "#A38AAB")

# Define the precipitators based on the Savitski dataset
Precipitators <- define_precipitators_control(TPP_all,Savitski = Savitski, both = TRUE, cutoff=0.5)
 
  
# Define the proteins that are changing their aggregation behavior
df_all <- NULL
comparisons <- NULL

# List of osmolytes that you want to test
sm_list <- c("TMAO", "glucose", "proline")
for(sm in sm_list){
  # Isolate only the measured points at the end of the gradient, such that aggregation profile can be compared
  comparisons_all <- TPP_all[[sm]]$python_fit %>%
    .[.$t %in% c(76, 72.5, 68.6) & .$type == "measured" & .$condition != "full model", ]
  
  comparisons_all %<>%
    dplyr::select(peptide, t, y, condition) %>%
    na.omit() %>%
    group_by(peptide, condition) %>%
    mutate(n = n()) %>%
    # Limit the proteins to the ones that have at least 3 values at the selected temperatures
    filter(n >= 3) %>%
    ungroup() %>%
    group_by(condition) %>%
    mutate(n_condition = n()) %>%
    # Limit the proteins to the ones that have at least 3 values at the selected temperatures
    filter(n_condition >= 3) %>%
    ungroup() %>%
    group_by(peptide) %>%
    mutate(n_conditions = length(unique(condition))) %>%
    # Limit the proteins to the ones that have two unique conditions (such that proteins only present in control or osmolyte are removed)
    filter(n_conditions == 2) %>%
    # Perform t-test to identify which proteins are significantly changing
    summarise(pvalue = t.test(y[condition== "control"], y[condition == sm])$p.value, 
              FC = mean(y[condition == sm])/mean(y[condition== "control"])) %>%
    ungroup() %>%
    # Correct for multiple hypothesis testing
    mutate(qvalue = p.adjust(pvalue, method="BH"))
  
  # Calculate log2 FC
  comparisons_all$log_FC <- log2(comparisons_all$FC)
  
  # Define whether the change is significant or not
  comparisons_all$is_significant <- comparisons_all$qvalue < 0.05 & comparisons_all$log_FC < (-1)
  # Define whether the aggregation is promoted or prevented
  comparisons_all$Aggregation <- ifelse(comparisons_all$qvalue < 0.05 & comparisons_all$log_FC < (-1), "Promoted", 
                                        ifelse(comparisons_all$qvalue < 0.05 & comparisons_all$log_FC > (1), "Delayed", 
                                               "Not significant"))
  
  # Check in the "Precipitator" data set to see whether the protein precipitates or not
  comparisons_all$Precipitator <- comparisons_all$peptide %in% Precipitators$Protein[Precipitators$Precipitator_TPP == "Precipitator"]

  # Define proteins that have promoted / delayed aggregation
  promoted <- length(comparisons_all$peptide[comparisons_all$log_FC < -1 & comparisons_all$qvalue < 0.05])
  prevented <- length(comparisons_all$peptide[comparisons_all$log_FC > 1 & comparisons_all$qvalue < 0.05])
  
  # Combine the data into the data sets
  df <- data.frame(condition = c(sm, sm),
                   change = c("Promoted", "Prevented"), 
                   n = c(promoted, prevented))
  df_all <- rbind(df,df_all)
  comparisons_all$Condition <- sm
  comparisons <- rbind(comparisons_all, comparisons)
}

# Rename some of the variables for a nice plot
df_all$Change <- ifelse(df_all$change == "Prevented", "- Aggreg.", "+ Aggreg.")
df_all$Change <- factor(df_all$Change, levels = c( "- Aggreg.", "+ Aggreg."))
df_all$condition <- dplyr::recode(df_all$condition , "TMAO" = "TMAO",  "glucose" =  "Glucose", "proline" = "Proline")
df_all$condition <- factor(df_all$condition , levels =(c("TMAO", "Proline", "Glucose")))

# Plot figure 5C
Figure_5C <- ggplot(df_all, aes(x=condition, y=n, fill=Change)) +
  geom_bar(stat="identity", position="dodge", col="black") +
  scale_fill_manual(values=c("- Aggreg." ="#B9B2A2", "+ Aggreg." = "#494439")) +
  ylab("Number of proteins") + 
  theme_classic() +
  theme(legend.text=element_text(size=15),
        legend.position = c(0.85, 0.8),
        axis.text.x = element_text(size=15, angle=45, vjust=1, hjust=1, color="black"), 
        axis.text.y = element_text(size=15,  color="black"), 
        axis.title = element_text(size=15),
        axis.title.x = element_blank(),
        legend.title = element_blank()) +
  xlab("")
#ggsave(paste(path, "Figures/Figure_5C.pdf", sep=""), height=4, width=4)

# Isolate the TMAO part, as TMAO shows the most interesting results
comparisons_all <- comparisons[comparisons$Condition == "TMAO",] %>%
  dplyr::rename("Protein" = "peptide") 

comparisons_all %<>%
  plyr::join(., All_combined_predictions, by="Protein")

comparisons_all %<>%
  plyr::join(., Precipitators, by="Protein")

comparisons_TMAO <- comparisons_all %>%
  #dplyr::rename("Protein" = "peptide") %>%
  #plyr::join(., Precipitators, by="Protein") %>% 
  mutate(decreased = qvalue < 0.05 & log_FC > 1 )

comparisons_TMAO$increased <- comparisons_TMAO$qvalue < 0.05 & comparisons_TMAO$log_FC < (-1)

Figure_S5C <- ggplot(comparisons_TMAO, aes(x=increased, y = log2(FC_TPP), fill= increased))+
  geom_violin(alpha=0.5) +
  scale_fill_manual(values=c("#B9B2A2", "#6E6654")) +
  geom_boxplot(width=0.2) +
  stat_compare_means() +
  theme_classic() +
  xlab("Promoted aggregation") +
  ylab("Log2(FC)") +
  theme(legend.position = "none")

## volcano plot for TMAO
Figure_S5B <- ggplot(comparisons_all, aes(x=log_FC, y = -log10(qvalue), col = Aggregation)) +
  geom_point(size=1) +
  xlab("Log2(FC)") +
  ylab("-log10(p-value)") + 
  geom_hline(yintercept = -log10(0.05), lty="dashed") +
  geom_vline(xintercept = c(-1,1), lty="dashed") +
  scale_color_manual(values = c("Delayed" = "darkblue", "Promoted" = "indianred3", "Not significant" = "azure3")) +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size=3)))+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=15), 
        legend.text = element_text(size=12), 
        legend.title = element_text(size=13),
        legend.position = c(0.8, 0.8))

df_all_LIP <- NULL
sm_LIP <- c("TMAO", "Betaine", "Proline", "Glycerol", "Glucose", "Trehalose")
comparisons_LIP_all <- NULL
for(sm in sm_LIP){
  #sm <- "TMAO"
  good_peptide <- all_results[[sm]]$python_score$peptide[all_results[[sm]]$python_score$Control > 0.5]
  
  comparisons_LIP <- all_results[[sm]]$python_fit %>%
    .[.$t %in% c(76, 72.5) & .$type == "measured" & .$condition != "full model" & .$trypticity == "FT" & .$peptide %in% good_peptide, ]
  
  comparisons_LIP %<>%
    dplyr::select(peptide, t, y, condition) %>%
    unique() %>%
    group_by(peptide, condition) %>%
    mutate(n = n()) %>%
    filter(n >= 3) %>%
    ungroup() %>%
    group_by(condition) %>%
    mutate(n_condition = n()) %>%
    filter(n_condition >= 3) %>%
    ungroup() %>%
    group_by(peptide) %>%
    mutate(n_conditions = length(unique(condition))) %>%
    filter(n_conditions == 2) %>%
    summarise(pvalue = t.test(y[condition== "Control"], y[condition == sm])$p.value, 
              FC = mean(y[condition == sm])/mean(y[condition== "Control"])) %>%
    ungroup() %>%
    mutate(qvalue = p.adjust(pvalue, method="BH"))
  
  comparisons_LIP %<>%
    dplyr::rename("Peptide" = "peptide") 
  
  comparisons_LIP %<>%
    plyr::join(., significant_all[[sm]]$scores[,c("Peptide", "Protein")], by="Peptide")
  
  comparisons_LIP$log_FC <- log2(comparisons_LIP$FC)
  
  comparisons_LIP$is_significant <- comparisons_LIP$qvalue < 0.05 & comparisons_LIP$log_FC < (-1)
  
  prevented <- length(comparisons_LIP$Peptide[comparisons_LIP$log_FC < -1 & comparisons_LIP$qvalue < 0.05])
  prevented_proteins <- length(unique(comparisons_LIP$Protein[comparisons_LIP$log_FC < -1 & comparisons_LIP$qvalue < 0.05]))
  
  promoted <- length(comparisons_LIP$Peptide[comparisons_LIP$log_FC > 1 & comparisons_LIP$qvalue < 0.05])
  promoted_proteins <- length(unique(comparisons_LIP$Protein[comparisons_LIP$log_FC > 1 & comparisons_LIP$qvalue < 0.05]))
  
  #out_LIP[[sm]] <- comparisons_LIP
  all <- length(comparisons_LIP$Peptide)
  proteins_all <- length(unique(comparisons_LIP$Protein))
  
  df <- data.frame(condition = c(sm, sm),
                   Proteolysis = c("Decreased", "Increased"), 
                   n = c(promoted, prevented),
                   n_all = c(all,all), 
                   n_protein = c(promoted_proteins, prevented_proteins), 
                   n_all_protein = c(proteins_all,proteins_all)
  )
  df_all_LIP <- rbind(df,df_all_LIP)
  comparisons_LIP_all[[sm]] <- comparisons_LIP}

df_all_LIP$n_percentage <- df_all_LIP$n/df_all_LIP$n_all
df_all_LIP$n_percentage_protein <- df_all_LIP$n_protein/df_all_LIP$n_all_protein

ggplot(df_all_LIP, aes(x=condition, y=n_percentage, fill=Proteolysis)) +
  geom_bar(stat="identity", position="dodge", col="black") +
  scale_fill_manual(values=c("#B9B2A2", "#494439")) +
  scale_y_continuous(labels = scales::percent) +
  ylab("Percentage of peptides") +
  theme_classic()

ggsave(paste(path,  "Figures/Barplot_LiP_proteolysis_peptide_level.pdf", sep=""))

df_all_LIP$condition <- factor(df_all_LIP$condition, levels = sm_list)
Figure_5H <- ggplot(df_all_LIP, aes(x=condition, y=n_percentage_protein, fill=Proteolysis)) +
  geom_bar(stat="identity", position="dodge", col="black") +
  scale_fill_manual(values=c("Increased" = "#B9B2A2", "Decreased" = "#494439")) +
  scale_y_continuous(labels = scales::percent) +
  ylab("Percentage of proteins") +
  theme_classic() +
  theme(legend.text=element_text(size=15),
        legend.position = c(0.2, 0.8),
        legend.title = element_text(size=15),
        axis.text.x = element_text(size=15, angle=45, vjust=1, hjust=1, color="black"), 
        axis.text.y = element_text(size=15,  color="black"), 
        axis.title = element_text(size=15),
        axis.title.x = element_blank()) +
  labs(fill = "Digestion")
ggsave(paste(path, "Figures/Figure_5H.pdf", sep=""),Figure_5H, height=4, width=5)

### Frr examples:
plot_peptide_adjusted <- function(all_results, sm_all, control = "control", peptide, factor_levels, pos="top_right"){
  #sm_all <- sm_list
  test <- data.frame(matrix(nrow=0, ncol = (ncol(all_results[[sm_all[1]]]$python_fit)+1)))
  colnames(test) <- c(colnames(all_results[[sm_all[1]]]$python_fit), "exp")
  
  for(sm in sm_all){

      test1 <-  all_results[[sm]][["python_fit"]][all_results[[sm]][["python_fit"]]$peptide %in% peptide, ] %>%
      filter(condition != "full model") %>%
      mutate(exp=sm) %>% 
      unique()
    
    test <- rbind(test, test1)
  }
  
  my_palette = c(Control = "#8D8D92",
                 control = "#8D8D92",
                 Betaine = "#8FAF7E", 
                 Trehalose = "#C18D73", 
                 TMAO = "#678EAF", 
                 Glycerol = "#997A8D", 
                 Proline = "#DFBA82",
                 proline = "#DFBA82", 
                 glucose = "#BB7676",
                 Glucose = "#BB7676",
                 Proline0 = "#F9DEB7",
                 Low_PK = "#A38AAB")
  test$exp <- recode(test$exp, "control" = "Control", "glucose"  ="Glucose", "proline" = "Proline")
  
  Labels <- test %>%
    dplyr::select(exp) %>%
    unique()
  
  position <- list()
  position[["left_bottom"]] <- c(-Inf, -Inf, -0.1, -0.5)
  position[["top_right"]] <- c(Inf, Inf, 1, 1.5)
  
  
  test$condition <- recode(test$condition, "control" = "Control", "glucose"  ="Glucose", "proline" = "Proline")
  test$exp <- factor(test$exp, levels = factor_levels)
  Labels$exp <- factor(Labels$exp, levels = factor_levels)
  plot3 <- ggplot(test) +  
    geom_ribbon(data = test[test$type == "fitted",], aes(ymin=conflik_lower, ymax=conflik_upper, x= t, fill=condition, group=interaction(exp, condition)), alpha=0.2) +
    geom_point(data = test[test$type == "measured",], aes(y=y, x=t, fill=condition), size=1, pch=21, col="black") + 
    geom_line(data = test[test$type == "fitted",], aes(y=y, x=t, col=condition, group=interaction(exp,condition)), lwd=0.8) +
    scale_color_manual(values = my_palette) +
    scale_fill_manual(values = my_palette )+
    #  geom_smooth(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), method="loess") +
    facet_wrap(~exp) +
    geom_text(data =Labels,  aes(label = exp), x = position[[pos]][1], y = position[[pos]][2], hjust = position[[pos]][3], vjust = position[[pos]][4], size=4.5) +
    #theme_minimal() +
    theme(axis.title = element_text(size = 12, color="black"),
          text=element_text(color="black"),
          axis.text = element_text(size=10, color="black"),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.background = element_blank(), axis.line = element_line(), panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(), legend.position = "none") +
    
    ylab("Scaled abundance") +
    xlab(paste("Temperature ",  "[", "\U00B0", "C]"))
    #theme_bw() +
    #theme(legend.position="none") +
  return(plot3)
  
}


Figure_5E <- plot_peptide_adjusted(TPP_all, sm_all = c("TMAO", "glucose", "proline"), control = "control", peptide =  "P0A805", factor_levels = c("TMAO", "Proline", "Glucose"), pos="left_bottom") +
  facet_wrap(~exp, scales="free", ncol=1) + 
  ylim(c(-0.15, 1.29)) +
  ylab("Abundance")
ggsave(paste(path, "Figures/Figure_5E.pdf", sep=""), Figure_5E, height=5, width=2)


sm_list <-  c("TMAO", "Betaine", "Glycerol", "Proline", "Trehalose", "Glucose")

Figure_5F1 <- plot_peptide_adjusted(all_results, sm_all = sm_list, control = "Control", peptide =  "ASPSLLDGIVVEYYGTPTPLR", factor_levels = sm_list) +
  facet_wrap(~exp, scales="free", ncol=2) + 
  ylim(c(-0.15, 1.3)) +
  ylab("Abundance")

ggsave(paste(path, "Figures/Figure_5F1.pdf", sep=""), Figure_5F1, height=3.2, width=3)

Figure_5F2 <- plot_peptide_adjusted(all_results, sm_all = sm_list, control = "Control", peptide =  "ASDLGLNPNSAGSDIR", factor_levels = sm_list) +
  facet_wrap(~exp, scales="free", ncol=2) + 
  ylim(c(-0.255, 1.25)) +
  ylab("Abundance")

ggsave(paste(path, "Figures/Figure_5F2.pdf", sep=""), Figure_5F2, height=3.2, width=3)

## Agrescan figure 
prot <- "P0A805"

ag_test <- Agrescan_list[[prot]] %>%
  dplyr::select(Position,a4vAHS ) %>%
  `colnames<-`(c("Position", "Value")) %>%
  mutate(b = Value > 0) %>%
  mutate(Value = with(rle(.$b), rep(lengths >= 5 & values, lengths))) %>%
  dplyr::select(-b) %>%
  mutate(  Protein =prot, 
           Measure = "Aggrescan")

LIP_test <- significant_all$TMAO$AA_level %>%
  filter(Protein %in% prot) %>%
  ungroup() %>%
  dplyr::select(Protein, positions, AA_aggregation) %>%
  `colnames<-`(c("Protein", "Position", "Value")) %>%
  mutate(Measure = "LIP") 

LIP_test <- rbind(LIP_test, ag_test)%>%
  mutate(Position = as.numeric(Position))

LIP_test$Color <- ifelse(LIP_test$Measure == "LIP" & LIP_test$Value > 0, "LIP", 
                         ifelse(LIP_test$Measure == "LIP" & LIP_test$Value <= 0, "Negative", 
                                ifelse(LIP_test$Measure == "Aggrescan" & LIP_test$Value > 0, "Aggrescan", "Negative")))


LIP_test$Measure <- recode(LIP_test$Measure, "LIP" = "LiP")
LIP_test$Color <- recode(LIP_test$Color, "LIP" = "LiP")

Figure_5F3 <- ggplot(LIP_test, aes(x=Position, y=Measure, fill=Color, col=Color)) +
  geom_tile(height = 0.8, lwd=1.5) +
  theme_classic() +
  scale_fill_manual(values=c("LiP" = as.character(my_palette["TMAO"]), "Aggrescan" = "#494439", "Negative" ="gray75")) +
  scale_color_manual(values=c("LiP" = as.character(my_palette["TMAO"]), "Aggrescan" = "#494439", "Negative" ="gray75"), guide="none") +
  theme(axis.title = element_text(size = 20, color="black"),
        legend.text = element_text(size=18), 
        legend.title = element_text(size=20),
        text=element_text(color="black"),
        axis.title.y=element_blank(),
        axis.text = element_text(size=20, color="black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.background = element_blank(), axis.line = element_line(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  labs(fill = "Agg.")
  
ggsave(paste(path, "Figures/Figure_5F3.pdf", sep=""), Figure_5F3, height=1.75, width=12)

##### Barplot for stabilisation ####
sm_TPP <- c("TMAO", "glucose", "proline")
percentage <- list()
all_stabilised <- NULL
for(sm in sm_TPP){
  #sm <- "TMAO"
  good_proteins <- TPP_all[[sm]]$python_score$peptide[TPP_all[[sm]]$python_score$control > 0.5 & TPP_all[[sm]]$python_score[sm] > 0.5]
  good_proteins <- good_proteins[good_proteins %in% Precipitators$Protein[Precipitators$Precipitator_TPP == "Precipitator"]]
  percentage[[sm]] <- length(TPP_significant[[sm]]$Stabilisation$peptide[TPP_significant[[sm]]$Stabilisation$peptide %in% good_proteins & TPP_significant[[sm]]$Stabilisation$Stabilisation > 0])/length(good_proteins)
  
  subset <- TPP_significant[[sm]]$Stabilisation[TPP_significant[[sm]]$Stabilisation$peptide %in% good_proteins,] %>%
    mutate(Condition = sm)
  
  subset$Condition <- ifelse(sm == "TMAO", "TMAO", stringr::str_to_title(sm) )
  
  all_stabilised <- rbind(all_stabilised, subset)
  
}
all_stabilised$Condition <- factor(all_stabilised$Condition, levels=c("TMAO", "Glucose", "Proline"))
Figure_S5A <-  ggplot(all_stabilised[all_stabilised$Stabilisation > 0,], aes(x=Condition, y=Stabilisation, fill=Condition)) +
  geom_boxplot(col="gray20") +
  scale_fill_manual(values=my_palette[c("TMAO", "Glucose", "Proline")]) +
  stat_compare_means(comparisons = list(c("TMAO", "Glucose"), c("Glucose", "Proline")), label="p.signif") +
  theme_classic()+
  theme(legend.position = "none") + 
  ylab("Stabilisation score") +
  xlab("")
ggsave(paste(path, "Figures/Boxplot_stabilisation_TPP.pdf", sep=""))

barplot_df <- data.frame(Condition = c("TMAO", "Glucose", "Proline"),
                         Percentage = unlist(percentage))
barplot_df$Condition <- factor(barplot_df$Condition, levels=c("TMAO", "Glucose", "Proline"))

Figure_S5B <- ggplot(barplot_df, aes(x=Condition, y=Percentage, fill=Condition)) +
  geom_bar(stat="identity", col="gray20") +
  scale_fill_manual(values=my_palette[c("TMAO", "Glucose", "Proline")]) +
  theme_classic() +
  theme(legend.position = "none") +
  ylab("Percentage stabilised proteins") +
  xlab("")
ggsave(paste(path, "Figures/Barplot_stabilisation_TPP.pdf", sep=""))

sm_LIP <- c("TMAO", "Glucose", "Proline")
df_all <- NULL
for(i in 1:length(sm_TPP)){
  
  good_proteins <- TPP_all[[sm_TPP[i]]]$python_score$peptide[TPP_all[[sm_TPP[i]]]$python_score$control > 0.5 & TPP_all[[sm_TPP[i]]]$python_score[sm_TPP[i]] > 0.5]
  good_proteins <- good_proteins[good_proteins %in% Precipitators$Protein[Precipitators$Precipitator_TPP == "Precipitator"]]
  stabilised <- TPP_significant[[sm_TPP[i]]]$Stabilisation$peptide[TPP_significant[[sm_TPP[i]]]$Stabilisation$peptide %in% good_proteins & TPP_significant[[sm_TPP[i]]]$Stabilisation$Stabilisation > 0]
  
  LIP_good <- significant_all[[sm_LIP[i]]]$Protein_level$Protein
  LIP_stabilised <- significant_all[[sm_LIP[i]]]$Protein_level$Protein[significant_all[[sm_LIP[i]]]$Protein_level$Protein_stabilisation > 0]
  
  all_intersect <- intersect(good_proteins, LIP_good)
  Both <- intersect(all_intersect, intersect(LIP_stabilised, stabilised))
  None <- all_intersect[!all_intersect %in% LIP_stabilised & !all_intersect %in% stabilised]
  Only_LIP <- intersect(LIP_stabilised, all_intersect)
  Only_LIP <- Only_LIP[!Only_LIP %in% stabilised]
  Only_TPP <- intersect(stabilised, all_intersect)
  Only_TPP <- Only_TPP[!Only_TPP %in% LIP_stabilised]
  
  df <-  data.frame(All = length(all_intersect), 
                    Both = length(Both), 
                    None = length(None), 
                    Only_TPP = length(Only_TPP), 
                    Only_LIP = length(Only_LIP), 
                    Condition = sm_LIP[i])
  
  df_all <- rbind(df_all, df)
  
}
df_all %<>%
  reshape2::melt(., id.var=c("Condition", "All")) %>%
  mutate(Percentage = value / All)

df_all %<>%
  mutate(Agreement = variable %in% c("Both", "None")) %>%
  group_by(Condition, Agreement) %>%
  mutate(Sum = sum(Percentage))
df_all$variable <- recode(df_all$variable, "Only_TPP" = "Only TPP", "Only_LIP" = "Only LiP")
df_all$variable <- factor(df_all$variable, levels = c("Only TPP", "Only LiP", "Both", "None"))
c("#9b9b7a", "#f4d2b9", "#B48F8F")
Figure_5B <- ggplot(df_all, aes(x=Condition)) +
  geom_bar(stat="identity",  aes(fill = variable, y=Percentage)) +
  geom_bar(stat="summary",aes(color=Agreement, y=Sum), alpha=0, lwd=0.5, show.legend = F) +
  theme_classic() +
  scale_fill_manual(values = c(None = "#9b9b7a", Both = "#5D5D38","Only TPP" = "#B48F8F", "Only LiP" = "#915A5A")) +
  scale_color_manual(values = c("black", "black")) +
  xlab("") +
  theme(legend.text = element_text(size=12),
        legend.title = element_blank(),
        axis.text.x = element_text(size=15, angle=45, vjust=1, hjust=1, color="black"), 
        axis.text.y = element_text(size=15,  color="black"), 
        axis.title = element_text(size=15),
        axis.title.x = element_blank())
ggsave(paste(path, "Figures/Figure_5B.pdf", sep=""),Figure_5B, height=4, width=4 )

names <- c("cysE", "guaC", "trpB", "clpP", "hisD")
pairs_TPP <- c("P0A9D4", "P60560", "P0A879", "P0A6G7", "P06988")
pairs_LIP <- c("EVVEEAYAADPEMIASAACDIQAVR", "VGIGPGSVCTTR", "TNQVLGQALLAK", "SLEQIER", "SFNTIIDWNSCTAEQQR")

for(i in 1:length(names)){
  pdf(paste(path, "Figures/TPP_LiP_agg_example_", names[i], ".pdf", sep="" ), height = 3, width = 6)
  new_labels <- c("TMAO" = "TPP profile")
  pT <- plot_peptide_grid(TPP_all, "TMAO", control = "control", pairs_TPP[i]) +
    theme(panel.grid = element_blank(), title = element_blank()) +
    #labs(title = "") +
    facet_wrap(~exp, labeller = labeller(exp = new_labels))
  new_labels <- c("TMAO" = "LiP profile")
  pL <- plot_peptide_grid(all_results, c("TMAO"), control = "control", pairs_LIP[i]) +
    theme(panel.grid = element_blank(), title = element_blank()) +
    #labs(title = "") +
    facet_wrap(~exp, labeller = labeller(exp = new_labels))
  
  grid.arrange(pT, pL, nrow = 1, top=names[i])
  dev.off()
}

TMAO_box <- significant_all$TMAO$Protein_level %>%
  plyr::join(., comparisons_all, by="Protein")
TMAO_box$Decreased_Aggregation <- ifelse(TMAO_box$Aggregation == "Delayed", TRUE, FALSE)

TMAO_box <- TMAO_box[!is.na(TMAO_box$Decreased_Aggregation),]
Figure_S5C <- ggplot(TMAO_box[TMAO_box$n_peptide > 1,], aes(x=Decreased_Aggregation, y=Protein_stabilisation, fill=Decreased_Aggregation)) +
  geom_violin(alpha=0.5) +
  geom_boxplot(width=0.2) +
  stat_compare_means() +
  scale_fill_manual(values = c("#B9B2A2", "#494439")) + 
  theme_classic() +
  theme(legend.position = "none") +
  ylab("LiP Protein Stabilisation") +
  xlab("TPP decreased aggregation")
ggsave(paste(path, "Figures/Decreased_aggregation_LiP_stabilisation.pdf", sep=""))


### Enrichment analysis ###
library(topGO)
GO_list <- GO_codes$Gene.ontology.IDs %>% strsplit(., "; ") %>% 
  `names<-`( GO_codes$X...Entry)

signif_table_all <- NULL
ontology <- c("BP", "CC", "MF")
for(i in 1:length(ontology)) {

  geneNames <- comparisons_all$Protein
  
  myInterestingGenes <-  comparisons_all$Protein[comparisons_all$Aggregation == "Promoted"]
  
  geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList) <- geneNames
  
  GOdata <- new("topGOdata", ontology = ontology[i], allGenes = geneList,
                annot = annFUN.gene2GO, gene2GO = GO_list)
  
  resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  #?runTest
  allRes <- GenTable(GOdata, classicFisher = resultFisher,
                     orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 100)
  allRes$classicFisher <- as.numeric(allRes$classicFisher)
  allRes$corrected <- p.adjust(allRes$classicFisher, "BH")
  signif_table <- allRes %>%
    filter( .$classicFisher < 0.01 )
  signif_table$Analysis <-  ontology[i]
  signif_table_all <- rbind(signif_table_all, signif_table)
  #table_output <- rbind(table_output, signif_table)
  }

signif_table_all <- signif_table_all[order(signif_table_all$classicFisher, decreasing = TRUE),]
signif_table_all$Term <- factor(signif_table_all$Term, levels=signif_table_all$Term)
ggplot(signif_table_all, aes(y=Term, x=-log10(classicFisher), fill=Analysis)) +
  geom_bar(stat="identity", col="gray20") +
  scale_fill_manual(values = c("#B9B2A2", "#968A71"  , "#494439")) + 
  theme_classic() +
  ylab("GO term") +
  xlab("-log10(pvalue)")
ggsave(paste(path, "Figures/GO_enrichment_aggregatiors.pdf", sep=""))

## Characteristics analysis for figure 5D:
group1 <- comparisons_all$Protein[comparisons_all$Aggregation == "Promoted"]
group2 <- comparisons_all$Protein[comparisons_all$Aggregation != "Promoted"]

group1 <- group1[!group1 %in% ribosomal_proteins$X...Entry]
group2 <- group2[!group2 %in% ribosomal_proteins$X...Entry]

output_test <- Test_two_groups(All_combined_predictions, group1, group2, "Target", "Not target")

ggplot_frame <- output_test$All_values %>%
  reshape2::melt(., id.vars = c("Protein", "group")) %>%
  #filter(variable == "") %>%
  mutate(Condition = "Condition 1")

output_significant <- output_test$Significance_test %>%
  mutate(Condition = "Condition 1")

# Remove redundant features
redundant <- c()
redundant <- c(unique( as.character(output_significant$variable[grepl("RASA_Average", output_significant$variable)])),
               "X.1", "ASAquick_rawscore", "RG_len", "disordered_percentage")

output_significant %<>%
  filter(!variable  %in% redundant) 

# Correct for multiple hypothesis testing
output_significant$qvalue <- p.adjust(output_significant$p_value, "BH")

# Filter for most significant differences
output_plot <- output_significant %>%
  filter(!variable  %in% redundant) %>%
  group_by(variable) %>%
  mutate(any_signif = any(qvalue < 0.01) ) %>%
  filter(any_signif == TRUE)

# Calculate whether the difference is positive/negative
ggplot_test <- ggplot_frame %>%
  filter(variable %in% output_plot$variable) %>%
  group_by(Condition, variable) %>%
  summarise(FC = median(value[group=="Target"], na.rm=TRUE) - median(value[group=="Not target"], na.rm=TRUE))

# Make two matrices for heatmap
matrix_values <- output_plot %>%
  reshape2::dcast(variable~Condition, value.var = "qvalue") %>%
  `rownames<-`(.$variable) %>%
  dplyr::select(-variable) %>%
  as.matrix()

matrix_values2 <- ggplot_test %>%
  reshape2::dcast(variable~Condition, value.var = "FC") %>%
  `rownames<-`(.$variable) %>%
  dplyr::select(-variable) %>%
  as.matrix()

matrix_values2[matrix_values > 0.01] <- 0

matrix_values <- -log10(matrix_values)
matrix_values[(10^(-matrix_values)) > 0.01] <- 0

matrix_values2_sign <- sign(matrix_values2)
matrix_colored <- matrix_values * matrix_values2_sign
f2 = circlize::colorRamp2(c(-max(abs(matrix_colored)), 0,  max(abs(matrix_colored))), c("#3F4DA3",  "white", "#A93F3F"))
rownames(matrix_colored) <- stringr::str_replace(rownames(matrix_colored), pattern = "_", replacement = " ")
rownames(matrix_colored)[rownames(matrix_colored) == "DRNApredDNAscore"] <- "DRNApredDNA"
rownames(matrix_colored)[rownames(matrix_colored) == "DRNApredRNAscore"] <- "DRNApredRNA"
rownames(matrix_colored)[rownames(matrix_colored) == "SCRIBERscore"] <- "SCRIBER"
rownames(matrix_colored)[rownames(matrix_colored) == "DisoRNAscore"] <- "DisoRNA"
rownames(matrix_colored)[rownames(matrix_colored) == "hydrophobicity"] <- "Hydrophobicity"
rownames(matrix_colored)[rownames(matrix_colored) == "protein length"] <- "Protein length"

library(ComplexHeatmap)
#pdf(paste(path, "Figures/Figure_3G_colored.pdf", sep=""), height=8, width = 4.5)
HM <- Heatmap(matrix_colored, name = "-log10(q-value)", col = f2,
              cluster_rows = TRUE, cluster_columns = FALSE,
              rect_gp = grid::gpar(col = "gray40", lwd = 1),
              row_names_gp = grid::gpar(size = 10, position = "right"),
              row_names_side = "right",
              
              show_column_dend = FALSE, show_row_dend = FALSE,
              
              heatmap_legend_param = list(
                legend_direction = "horizontal", 
                legend_width = unit(3, "cm")),
              
)
draw(HM, heatmap_legend_side="bottom")

