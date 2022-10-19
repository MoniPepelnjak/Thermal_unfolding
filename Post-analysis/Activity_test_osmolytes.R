### Activity test all compounds ###
### The goal of this script is to analyse the effect of temperature on PK activity

# Load in the libraries
library(gridExtra)
library(magrittr)
library(ggplot2)
library(dplyr)

# Load in the data 
Activity_data_full <- read.delim("/Users/moni/Documents/Phd/Experiments/024/20210818_Full_report.tsv")
Annotation <- read.csv("/Users/moni/Documents/Phd/Experiments/024/024_annotation_DIA.csv")

# Select for E. coli, FT peptides that are short, no missed cleavages, and abundant
Activity_data <- Activity_data_full %>% 
  unique() %>%
  mutate(pep_len = nchar(PEP.StrippedSequence)) %>%
  filter(pep_len <= 10, 
         PG.FastaFiles == "200331_ecoli", 
         PEP.NrOfMissedCleavages == 0, 
         PEP.DigestType....Trypsin.P. == "Specific") %>%
  filter(PEP.Quantity > 1000) #%>%

# Combine with annotation file
Activity_data %<>%
  plyr::join(., Annotation, by="R.FileName")

# Remove missing values
Activity_data %<>%
  na.omit()

# Select for osmolytes you want to test
sm_list <- c("TMAO", "Betaine", "Glycerol", "Proline", "Trehalose", "Glucose")

# Remove oxidised peptides
Activity_data <- Activity_data[!grepl(pattern = "Oxidation", Activity_data$EG.ModifiedSequence),]

# Calculate differential abundance
Activity_data_low_all <- NULL
for(sm in sm_list){
  Activity_data_compound <- Activity_data %>%
    filter(Condition %in% c(sm, "Control"), # select control and osmolyte
           replicate %in% c(1, 3,4)) # select three replicates

  Activity_data_compound %<>%
    group_by(PEP.StrippedSequence, Condition) %>%
    mutate(mean_fc = mean(log2(PEP.Quantity))) %>%
    group_by(EG.ModifiedSequence) %>%
    mutate( n = n()) %>%
    filter(n >= 6) %>% # select only for full profiles 
    mutate(p_value = t.test(log2(PEP.Quantity)~factor(Condition))$p.value) %>%
    mutate(FC = mean(log2(PEP.Quantity[Condition == sm]))-mean(log2(PEP.Quantity[Condition == "Control"])))
  
  # Correct for multiple hypothesis testing
  Activity_data_low <- Activity_data_compound %>%
    dplyr::select(EG.ModifiedSequence, FC, p_value, pep_len) %>%
    unique() %>%
    ungroup() %>%
    mutate(qvalue = p.adjust(p_value, "BH"),
           log2_FC = FC) %>%
    mutate(Condition = sm)
  
  Activity_data_low_all <- rbind(Activity_data_low, Activity_data_low_all)

}

# Calculate which are significantly changing peptides
Activity_data_summary <- Activity_data_low_all %>%
  ungroup() %>%
  group_by(Condition) %>%
  summarise(All_peptides = n(), 
            All_significant = sum(abs(log2_FC) > 1 & qvalue < 0.05)) %>%
  mutate(percentage = All_significant/All_peptides)

# plot the volcano plot
ggplot(Activity_data_low_all, aes(x=log2_FC, y=-log10(qvalue), col = abs(log2_FC))) +
  geom_point() + 
  geom_hline(yintercept = -log10(0.05), lty="dashed") +
  geom_vline(xintercept = c(-1,1), lty="dashed") +
  scale_color_gradient2(low="azure3", mid="#CA8888", high="#CF1414", midpoint = 4) + 
  xlim(c(-max(abs(Activity_data_low$log2_FC)), max(abs(Activity_data_low$log2_FC)))) +
  #xlim(c(-8, 8)) +
  #theme_minimal() +
  facet_wrap(~Condition) +
  theme(axis.line = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), 
        legend.position = "", 
        panel.background = element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  ylab("-log10(qvalue)") +
  xlab("log2(FC)")

# Make a long format of the summary
Activity_data_summary_long <- Activity_data_summary %>%
  reshape2::melt(id.var=c("Condition", "percentage"))

# Re-factorise the osmolytes
Activity_data_summary_long$Change <- ifelse(Activity_data_summary_long$variable == "All_peptides", "Not changing", "Significant")
Activity_data_summary_long$Condition <- factor(Activity_data_summary_long$Condition, levels=sm_list)

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

# Plot figure S1D
ggplot(Activity_data_summary_long, aes(x=Condition, y=value, fill=Condition)) +
  geom_bar(stat="identity", aes(alpha = Change), col="black") +
  scale_fill_manual(values=my_palette[sm_list], guide="none") +
  scale_alpha_manual(values = c(0.5,1)) +
  geom_label(aes(label = paste(round(percentage*100, 1), "%", sep="") , y = 2000)) +
  theme_classic() +
  #ylim(c(0,45))+
  #geom_hline(yintercept = mean(all_boxplot$Protein_stabilisation), linetype = 2, alpha=0.5, lwd=1)+ # Add horizontal line at base mean
  theme(#legend.position = "none",
        axis.text.x = element_text(size=15, angle=45, vjust=1, hjust=1, color="black"), 
        axis.text.y = element_text(size=15,  color="black"), 
        axis.title = element_text(size=15),
        axis.title.x = element_blank(),
        legend.position = "bottom", 
        legend.text = element_text(size=15), 
        legend.title = element_blank()) +  
  xlab("") +
  ylab("Number of peptides")

ggsave("/Users/moni/Documents/Phd/Osmolyte_paper/Sub_sigures/Figure_S1D.pdf", height = 4, width=4.5)
