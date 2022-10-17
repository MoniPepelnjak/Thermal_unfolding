## Correlation FT/HT peptides
# Author: Monika Pepelnjak

# Goal of the script: 
  # Match HT peptides with the "parent" FT peptide
  # Calculate correlation between FT and HT peptides
  # Plot the matched FT/HT examples

# Data needed
  # Python fits
  # Significant analysis list

# Load in the supporting script with the functions
source("FT_peptide_function.R")

# Load in the data
LIP <- readRDS("/Volumes/My Passport for Mac/data/All_data_1M.rds")
significant_all <- readRDS("/Volumes/My Passport for Mac/data/Significant_1M.rds")

# Libraries
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggpubr)

# Select one condition where you would want to perform the matching
# For osmolyte paper, trehalose condition was chosen, as we had the higherst number of high quality peptides
# The script could be optimised in a way to combine the best fits accross different experiments
sm <- "Trehalose"

# Match FT to HT peptides, use cutoff at least 0.5 otherwise the data is too noisy to obtain meaningful results
matched_peptides <- matching_FT_HT(LIP, sm=sm,control = "Control", 0.5) 

# Plot the example peptides, in this case a "weird" HT peptide was removed just for visual representation
matched_peptides_sub <- matched_peptides %>%
  filter(HT_peptide != "IGLNLPSGEMGR")

# Plot the matches with a very specific selected peptides
# This plot needs adjustment if different peptides/data is selected. This includes the size of export, as otherwise the plot does not align nicely
pep_matches_sequence <- plot_matches(matched_peptides_sub, FT_peptide_match = "ITIGLNLPSGEMGR") +
  scale_y_discrete(labels=c("ITIGLNLPSGEMGR" = "FT", "NLPSGEMGR" = "HT",
                            "LPSGEMGR" = "HT", "LNLPSGEMGR" = "HT", "ITIGLNLPSGEM" = "HT"))

ggsave("/Users/moni/Documents/Phd/Thesis/Chapter1/Mapped_sequence_HT_FT_matches.pdf",pep_matches_sequence, width=2.30, height=3.36)

# Plot curve examples
pairs_fits <- plot_HT_pairs(LIP = LIP,sm="Trehalose", plot_which = "Control", matched_peptides = matched_peptides_sub, FT_peptide_sel = "ITIGLNLPSGEMGR" )
ggsave("/Users/moni/Documents/Phd/Thesis/Chapter1/Mapped_profiles_HT_FT_matches.pdf",pairs_fits, width=4, height=3.36)

# restructure the df a bit
matched_peptides_2 <- matched_peptides %>%
  unique() %>%
  filter(matched == "HT_match") %>%
  dplyr::select(FT_peptide, HT_peptide) %>%
  `colnames<-`(c("FT_peptide", "peptide"))

# Select a couple of points from the fit
temp <- LIP[[sm]]$python_fit$t[ LIP[[sm]]$python_fit$type == "fitted"] %>% sort() %>% unique() %>% .[seq(from = 1, to = length(.), by=5 )]

# Selected the fitted points from the python fits - for HT peptides
HT_peptides <- LIP[[sm]]$python_fit %>%
  .[.$peptide %in% c(matched_peptides_2$peptide),] %>% # select only matched peptides
  .[.$type == "fitted" & .$condition == "Control" & .$t %in% temp, ] %>% # select the fitted values for control at previously specified temperatures
  plyr::join(matched_peptides_2, by="peptide") %>% # join with the data set to get the matched ft peptide
  mutate(FT_temp = paste(FT_peptide, as.character(round(t, 2)), sep="%%%")) %>%
  dplyr::select(peptide, FT_temp, t, y)

# Selected the fitted points from the python fits - for FT peptides
FT_peptides <- LIP[[sm]]$python_fit %>%
  .[.$peptide %in% c(matched_peptides_2$FT_peptide),] %>%
  .[.$type == "fitted" & .$condition == "Control" & .$t %in% temp, ] %>%
  mutate(FT_temp = paste(peptide, as.character(round(t, 2)), sep="%%%")) %>%
  dplyr::select(peptide, FT_temp, y) %>%
  `colnames<-`(c("FT_peptide", "FT_temp", "y_FT"))

# Combine the FT and HT peptide data frames
FT_HT <- plyr::join(HT_peptides, FT_peptides, by="FT_temp") 

FT_HT_sum <- FT_HT %>%
  group_by(peptide) %>%
  summarise(cor = cor(y, y_FT)) %>% # calculate correlation for each peptide
  mutate(abs_cor = abs(cor)) %>%
  mutate(data = "data")%>%
  mutate(is_cor = abs(cor) > 0.8)%>%
  mutate(peptide_length = nchar(peptide)) %>%
  mutate(long = peptide_length < 8 ) 

FT_HT_sum %<>%
  group_by(peptide) %>%
  mutate(Last_AA = substr(peptide, peptide_length, peptide_length)) %>%
  mutate(C_tryptic = Last_AA %in% c("K", "R"))

# Plot the histogram of correlations
ggplot(FT_HT_sum, aes(x=cor, fill="#B6AFA8")) +
  geom_histogram(col="black", binwidth = 0.1) +
  scale_fill_manual(values=c("#A1948D")) +
  theme_classic() + 
  #facet_wrap(~C_tryptic, scales="free_y") +
  ylab("Peptide count") +
  xlab("FT-HT correlation") +
  theme(legend.position = "none", 
        text = element_text(size=20), 
        axis.title = element_text(size=18)) 
