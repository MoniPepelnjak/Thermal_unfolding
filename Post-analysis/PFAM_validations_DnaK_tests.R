### PFAM domain mapping ####
### The script is used to perform the domain-level analysis as shown in Figure 6
  # Differential analysis for stabilisation of different domains
  # Example plot for stabilisation of DnaK
  # Analysis of melting temperatures of DnaK 

# Load in the information about the domains - from Pfam
PFAM_domains <- read.csv("/Volumes/My Passport for Mac/SemesterProject/Databases/PFAM/Proteins_Regions_pfam.csv")

# Load in the significance analysis 
significant_all <- readRDS("/Volumes/My Passport for Mac/data/Significant_both_quantiles.rds")

# Load in the gene to uniprot annotation 
Gene_names <- read.csv("/Users/moni/Documents/Phd/databases/Gene_name_ids.csv")

# Define the color palette
my_palette = c(Betaine = "#B5C3AD", 
               Trehalose = "#d4a993", 
               TMAO = "#87A9C6", 
               Glycerol = "#a997a2", 
               Proline = "#edd9bb", 
               Glucose = "#C78C8C")

### Libraries ###
library(magrittr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ComplexHeatmap)

# Select proteins with at least two domains 
PFAM_domains %<>%
  unique() %>%
  filter(Type.1 == "Domain") %>%
  group_by(Protein) %>%
  mutate(n = n()) %>%
  filter(n >= 2) %>%
  mutate(Start = as.numeric(Start), 
         End = as.numeric(End))

sig_test_all <- NULL

# Define all the osmolytes to test
sm_list <- c("TMAO", "Proline", "Trehalose", "Glucose", "Betaine", "Glycerol")
TMAO_test_all <- NULL
for(sm in sm_list){
  TMAO_test <- significant_all[[sm]]$AA_level %>%
    .[.$Protein %in% PFAM_domains$Protein,] # Only select the proteins with the two domains or more
  
  TMAO_test %<>%
    plyr::join(.,PFAM_domains, by="Protein" ) %>%
    group_by(positions, Description) %>%
    mutate(in_domain = ifelse((positions >= Start &  positions <= End), Description, NA)) # Check whether the position is in the domain
  
  TMAO_test %<>%
    na.omit()
  
  TMAO_test %<>%
    ungroup() %>%
    dplyr::select(Peptide, Protein, Description, AA_score) %>%
    group_by(Peptide) %>%
    mutate(AA_score = quantile(AA_score, 0.75)) %>% # Calculate peptide level score
    unique() %>%
    group_by(Protein, Description) %>%
    mutate(n_per_domain = n()) %>%
    filter(n_per_domain >=3 ) %>% # At least 3 peptides per domain
    group_by(Protein) %>%
    mutate(n_unique = length(unique(Description))) %>%
    filter(n_unique >= 2) %>% # At least two conditions
    mutate(Condition = sm)
  
  TMAO_test_all <- rbind(TMAO_test_all, TMAO_test)
  
  # Test for proteins with two domains
  TMAO_test_two_domains <- TMAO_test %>%
    ungroup() %>%
    filter(n_unique == 2) %>%
    group_by(Protein) %>%
    mutate(p = wilcox.test(AA_score[Description == unique(Description)[1]], AA_score[Description == unique(Description)[2]])$p.value) # test for significant difference
  
  # Test for the proteins with more than 2 domains
  TMAO_test_more_domains <- TMAO_test %>%
    ungroup() %>%
    filter(n_unique > 2) %>%
    group_by(Protein) %>%
    mutate(p = kruskal.test(AA_score ~ Description)$p.value) 
  
  # Clean up the df, correct for multiple hypothesis testing and combine the to data frames
  TMAO_test_two_domains_unique <- TMAO_test_two_domains %>%
    dplyr::select(Protein, p) %>%
    ungroup() %>%
    unique() %>%
    na.omit() 
  
  TMAO_test_more_domains_unique <- TMAO_test_more_domains %>%
    dplyr::select(Protein, p) %>%
    ungroup() %>%
    unique() %>%
    na.omit() 
  
  
  sig_test <- rbind(TMAO_test_two_domains_unique, TMAO_test_more_domains_unique) %>%
    mutate(Condition = sm) %>%
    mutate(padj = p.adjust(p, method="BH"))
  
  # export the significance data frame
  sig_test_all <- rbind(sig_test_all,sig_test)
  
}

# Check the number of proteins at which we could test for differential behaviour
length(unique(sig_test_all$Protein))

# Select for proteins with significant cutoff at least for one protein
sig_test_all_sub <- sig_test_all %>%
  mutate(significant = padj < 0.05) %>%
  group_by(Protein) %>%
  mutate(any_sig = any(significant)) %>%
  filter(any_sig == TRUE)

# Assign the colors for heatmap
sig_test_all_sub$p_col <- ifelse(sig_test_all_sub$padj > 0.1, 1, sig_test_all_sub$padj)

# Transform the data to plot with a complex heatmap function
sig_test_all_sub_wide <- sig_test_all_sub %>%
  plyr::join(.,Gene_names, by="Protein" ) %>%
  dplyr::select(-significant) %>%
  dplyr::ungroup() %>%
  base::unique() %>%
  reshape2::dcast(Gene_names~Condition, value.var = "p_col") %>%
  na.omit() %>%
  `rownames<-`(.$Gene_names) %>%
  .[,-1] %>%
  as.matrix()

sig_test_all_sub_wide <- -log10(sig_test_all_sub_wide)
f2 = circlize::colorRamp2(seq(0, max(sig_test_all_sub_wide), length = 3), c("white",  "#B9B2A2", "#494439"))

pdf("/Users/moni/Documents/Phd/Osmolyte_paper/Figure_6A.pdf", height=7.5, width=2.8)
HM <- ComplexHeatmap::Heatmap(sig_test_all_sub_wide, cluster_rows = TRUE, cluster_columns = TRUE, col = f2,
                        row_names_gp = grid::gpar(fontsize = 10), show_row_names = TRUE, raster_device = "CairoPNG", 
                        #left_annotation = row_ha, 
                        show_row_dend = FALSE,
                        show_column_dend = FALSE,
                        rect_gp = grid::gpar(col = "gray20", lwd = 1),
                        heatmap_legend_param = list(
                          title = "-log10 (q-value)", direction = "horizontal", 
                          position = "bottom",legend_width = unit(3, "cm")

                        ))
draw(HM, heatmap_legend_side="bottom")
dev.off()

# Plot examples for DnaK
TMAO_test_all_check <- TMAO_test_all %>%
  group_by(Peptide, Condition) %>%
  ungroup() %>%
  dplyr::select(Peptide, Protein, Description, Condition, AA_score) %>%
  unique()

sub_test <- TMAO_test_all_check[TMAO_test_all_check$Protein %in% "P0A6Y8",] %>%
  unique() %>%
  filter(Condition %in% c("Proline", "TMAO"))

set.seed(10)

# Data frame for plot titles
Labels_examples <- sub_test %>%
  dplyr::select(Condition) %>%
  unique()

# Plot figure 6C
Figure_6C <- ggplot(data = sub_test, aes(x=Description, y=AA_score)) +
  geom_boxplot(data = sub_test, aes(x=Description, y=AA_score, fill=Description), outlier.alpha = 0,  col="black") +
  geom_jitter(data = sub_test, aes(x=Description, y=AA_score, fill=Description), width = 0.1, pch=21, col="black", size=1)+
    facet_wrap(~Condition, scales="free_y") +
  ylim(c(min(sub_test$AA_score), max(sub_test$AA_score)+5)) +
  scale_fill_manual(values=c("#E2DEA1", "#5A7032")) +
  geom_text(data =Labels_examples,  aes(label = Condition), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size=5.5) +
  stat_compare_means(comparisons = list(c("ATP-binding", "Substrate-binding")), label = "p.signif", size=8) +
  #theme_minimal() +
  theme(axis.title = element_text(size = 15, color="black"),
        text=element_text(color="black"),
        axis.text = element_text(size=15, color="black"),
        axis.text.x = element_text(size=15, color="black", angle=45, vjust=1, hjust=1),
        
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.background = element_blank(), axis.line = element_line(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), legend.position = "none", 
        axis.title.x = element_blank()) +
  ylab("Stabilisation score") +
  xlab("")

ggsave("/Users/moni/Documents/Phd/Osmolyte_paper/Figure_6C.pdf", Figure_6C, height=4.5, width=5)

#### DnaK extract peptide info ###
## This part was fit with the python script
## Load in the data from previously fit, manually selected peptides

DnaK_parameters_betaine <- read.csv("/Users/moni/Documents/Phd/Collaborations/C6/LIP_Fit_Proline_output_parameters.csv")[,-1] %>%
  mutate(Condition = "Proline")

DnaK_parameters_TMAO <- read.csv("/Users/moni/Documents/Phd/Collaborations/C6/LIP_Fit_TMAO_output_parameters.csv")[,-1] %>%
  mutate(Condition = "TMAO")

DnaK_parameters_trehalose_control <- read.csv("/Users/moni/Documents/Phd/Collaborations/C6/LIP_Fit_Trehalose_control_output_parameters.csv")[,-1] %>%
  mutate(Condition = "Control")

DnaK_parameters_control <- read.csv("/Users/moni/Documents/Phd/Collaborations/C6/LIP_Fit_DnaK_output_parameters.csv")[,-1] %>%
  mutate(Condition = "Control")

DnaK_parameters <-  rbind(DnaK_parameters_betaine, DnaK_parameters_control) %>%
  rbind(., DnaK_parameters_TMAO) 

selected_peptides_to_map <- read.csv("/Users/moni/Documents/Phd/Collaborations/C6/Selected_peptides_DnaK.csv") %>%
  rbind(., read.csv("/Users/moni/Documents/Phd/Collaborations/C6/Selected_peptides_DnaK_Trehalose.csv"))

# Filter and re-organise the calculated melting temperatures
DnaK_parameters %<>%
  filter(Tm > 310) %>% # filter too low melting temperatures (the melting temperatures should be in our experimental temperature range)
  `colnames<-`(c("Peptides", "Tm", "H", "Condition")) %>% 
  plyr::join(., selected_peptides_to_map, by="Peptides") # Join with peptides for domain annotation

DnaK_parameters$Domain <- recode(DnaK_parameters$Domain, "ATP" = "ATP-binding", "Substrate" = "Substrate-binding")

# DF for plot naames
Labels_examples <- DnaK_parameters %>%
  dplyr::select(Domain) %>%
  unique()

DnaK_parameters%<>%
  unique()

# Plot 6D
Figure_6D <- ggplot(data=DnaK_parameters, aes(x=Condition, y=Tm-273.15)) +
  geom_boxplot(data=DnaK_parameters, aes(x=Condition, y=Tm-273.15, fill=Condition), outlier.alpha = 0, col="black") +
  geom_jitter(data=DnaK_parameters, aes(x=Condition, y=Tm-273.15, fill=Condition), width=0.1, pch=21, col="black", size=1) +
  geom_text(data =Labels_examples,  aes(label = Domain), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size=5.5) +
  facet_wrap(~Domain, scales="free") +
  ylim(c(min(DnaK_parameters$Tm-273.15), max(DnaK_parameters$Tm-273.15)+5))+
  scale_fill_manual(values=my_palette) +
  stat_compare_means(comparisons = list(c("Control", "Proline"), c("Control", "TMAO")), label = "p.signif", size=8, vjust=0.5) +
  theme_classic() +
  theme(axis.title = element_text(size = 15, color="black"),
        text=element_text(color="black"),
        axis.text = element_text(size=15, color="black"),
        axis.text.x = element_text(size=15, color="black", angle=45, vjust=1, hjust=1),
        
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.background = element_blank(), axis.line = element_line(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), legend.position = "none", 
        axis.title.x = element_blank()) +
  ylab("Melting temperature") +
  xlab("Condition")

ggsave("/Users/moni/Documents/Phd/Osmolyte_paper/Figure_6D.pdf", Figure_6D, height=4.5, width=5)
