### PK Activity test temperature ###
### The goal of this script is to analyse the effect of temperature on PK activity

# Load in the libraries
library(gridExtra)
library(magrittr)
library(ggplot2)
library(dplyr)

# Load in the data
Activity_data_full <- read.delim("/Users/moni/Documents/Phd/Experiments/014/Peptides_normalized_tryptic.tsv")

# Perform filtering to the data 
Activity_data <- Activity_data_full %>% 
  unique() %>%
  mutate(pep_len = nchar(EG.ModifiedSequence)-2) %>% 
  filter(pep_len <= 10) %>% # only select short peptides such that we can assume they have no structure
  group_by(EG.ModifiedSequence, variable1) %>%
  mutate(mean_fc = mean(normalized_value)) %>%
  group_by(EG.ModifiedSequence) %>%
  mutate( n = n()) %>% 
  filter(n == 20) %>%  # only select the "full profiles"
  filter(variable1 < 76) %>% # remove the 76 - we have seen that the activity is reduced by then - we aim to catch the highest difference
  # take the first two temperatures as low and last two as high temperatures --> we only have duplicates and need more points for statistics
  mutate(groups = ifelse(variable1 < 41, "LOW",
                        ifelse(variable1 > 64, "HIGH", NA)) ) %>%
  filter(groups %in% c("LOW", "HIGH")) %>%  # only filter for the temperatures in those two groups
  # calculate the differential abundance and FC
  mutate(p_value = t.test(normalized_value~factor(groups))$p.value) %>% 
  mutate(FC = mean(normalized_value[groups == "HIGH"]) - mean(normalized_value[groups == "LOW"]))

# Correct for multiple hypothesis testing
Activity_data_low <- Activity_data %>%
  dplyr::select(EG.ModifiedSequence, FC, p_value) %>%
  unique() %>%
  ungroup() %>%
  mutate(qvalue = p.adjust(p_value, "BH"))

# Define whether the change is singnificant
Activity_data_low$Change <- ifelse(Activity_data_low$qvalue < 0.05 & abs(Activity_data_low$FC) > 1, "Significant", "Not significant")

# plot volcano plot
p1 <- ggplot(Activity_data_low, aes(x=FC, y=-log10(qvalue), col=Change)) +
  geom_point() + 
  geom_hline(yintercept = -log10(0.05), lty="dashed") +
  geom_vline(xintercept = c(-1,1), lty="dashed") +
  scale_color_manual(values = c("Significant" = "indianred3", "Not significant" = "gray85") ) +
  xlim(c(-max(abs(Activity_data_low$FC)), max(abs(Activity_data_low$FC)))) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), legend.position = "none", 
        axis.text = element_text(size=15), 
        axis.title=element_text(size=15)) +
  xlab("log2(FC)") +
  ylab("-log10(adj. p-value)")

# Summarise the percentage of significantly changing proteins
Activity_data_summary <- Activity_data_low %>%
  mutate(is_significant = abs(FC)> 1 & qvalue < 0.05) %>%
  summarise(Significant = sum(is_significant), 
            Not_significant = sum(!is_significant)) %>%
  mutate(percentage = Significant/(Significant + Not_significant )) %>%
  reshape2::melt(., id.var = "percentage")

# Define as significant or non-significant
Activity_data_summary$Change <- ifelse(Activity_data_summary$variable == "Significant", "Significant", "Not significant")

# Plot the summary plot
p2 <- ggplot(Activity_data_summary, aes(x="Temperature", fill=Change, y=value)) +
  geom_bar(stat="identity", col="black") +
  scale_fill_manual(values = c("Significant" = "indianred3", "Not significant" = "gray85") ) +
  geom_label(aes(label = paste(round(percentage*100, 1), "%", sep="") , y = 2000), fill="indianred3") +
  theme_classic() +
  ylab("Number of peptides") +
  xlab("") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), legend.position = "none", 
        axis.text = element_text(size=15), 
        axis.title=element_text(size=15), 
        legend.text = element_text(size=12), 
        legend.title = element_text(size=12)) +
  theme(legend.position = "")

# Combine the plots in the grid
pdf("/Users/moni/Documents/Phd/Osmolyte_paper/Sub_sigures/Figure_S1C.pdf", height = 4, width=6.5)
grid.arrange(p1, p2,
  widths = c(2, 0.8), nrow=1)
dev.off()

