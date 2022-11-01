### Frr aggregation validations ###
## Script to analyse the precipitation of Frr and Gnd
## Script produces the figures used in the Figure 5 and Figure S5 of the Osmolyte paper

library(magrittr)
library(ggplot2)
library(dplyr)

####
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

# Temperature assays
# Monitor the aggregation of Frr in control condition or in presence of 1 M TMAO, Glycerol or Betaine at different temperatures
Temp_assay <- read.csv("/Users/moni/Documents/Phd/Collaborations/C6/Frr_aggregation_assays/Control_TMAO_Betaine_Glycerol_Temp_assay.csv")

Temp_assay_2 <- Temp_assay %>%
  group_by(Condition) %>%
  mutate(Y_scaled = abs(Absorbance - mean(Absorbance[Sample == "Negative"]) ))  %>%
  filter(Sample != "Negative") %>% 
  #group_by(Replicate) %>%
  mutate(Y_percentage = Y_scaled/Y_scaled[Y_scaled == max(Y_scaled)])

p1 <- ggplot(Temp_assay_2, aes(x=Temperature, y=Y_percentage)) +
  stat_summary(aes(col = Condition), fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", width=0.2) +
  stat_summary(aes(fill = Condition), fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="ribbon", alpha=0.2) +
  stat_summary(aes(col = Condition), fun=mean, geom="point", size=3) +
  stat_summary(aes(col = Condition), fun=mean, geom="line") +
  theme_classic() +
  scale_color_manual(values = c(Control = "azure3", my_palette[c("Betaine", "TMAO", "Glycerol")])) +
  scale_fill_manual(values = c(Control = "azure3", my_palette[c("Betaine", "TMAO", "Glycerol")])) +
  ylab("Relative absorbance") + 
  xlab("Temperature")

ggsave(paste("/Users/moni/Documents/Phd/Collaborations/C6/Frr_aggregation/Frr_aggregation_temperature_dependance.pdf", sep=""), plot = p1, height = 4, width=5)

# Aggregation assay at 76?
# Protein concentration: 0.3 mg/ml
# Time of incubation: 5 minutes
Agg_assay <- read.csv("/Users/moni/Documents/Phd/Collaborations/C6/Frr_aggregation_assays/All_compounds_Bradford_assay.csv")

Background = mean(Temp_assay$Absorbance[Temp_assay$Sample == "Negative"])

Agg_assay %<>%
  mutate(Bradford = Bradford - Background) %>%
  group_by(Condition) %>%
  mutate(Y_scaled = abs(Bradford/mean(Bradford[Sample == "All"]) ))  %>%
  filter(Sample == "Soluble")

median_line <- Agg_assay %>%
  group_by(Condition) %>%
  summarise(Median = median(Y_scaled))
sm_list <- c("Control", "TMAO", "Betaine", "Glycerol", "Proline", "Glucose", "Trehalose")
Agg_assay$Condition <- base::factor(Agg_assay$Condition, levels=sm_list)
median_line$Condition <- factor(median_line$Condition, levels = sm_list)
p2 <- ggplot() +
  #geom_boxplot(outlier.alpha = 0) + 
  geom_crossbar(data=median_line, aes(ymin = Median, ymax= Median, y = Median, x = Condition, col=Condition),
                size=0.5,col="black", width = .5) + 
  geom_jitter(data = Agg_assay, aes(x=Condition, fill = Condition, y= Y_scaled), width = 0.2, pch=21, col="black") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=15, angle=45, vjust=1, hjust=1, color="black"), 
        axis.text.y = element_text(size=15,  color="black"), 
        axis.title = element_text(size=15),
        axis.title.x = element_blank()) +
  scale_color_manual(values = c(Control = "gray75", my_palette[sm_list])) +
  scale_fill_manual(values = c(Control = "gray75", my_palette[sm_list])) +

  ylab("Relative protein concentration") + 
  xlab("Compound")

ggsave(paste("/Users/moni/Documents/Phd/Osmolyte_paper/Figure5G.pdf", sep=""), plot = p2, height = 4, width=3.5)

# Gnd resolubilisation assay for control and betaine
Gnd_assay <- read.csv("/Users/moni/Documents/Phd/Collaborations/C6/Gnd_aggregation/Other_detergent_resulibilisation_2.csv")

Gnd_assay %<>%
  group_by(Condition) %>%
  mutate(scaled = (Value-min(Value))/(max(Value)-min(Value)))

ggplot(Gnd_assay, aes(x=Det_concentration, y=scaled)) +
  stat_summary(aes(col = Condition), fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", width=0.01) +
  stat_summary(aes(fill = Condition), fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="ribbon", alpha=0.2) +
  stat_summary(aes(col = Condition), fun=mean, geom="point", size=3) +
  stat_summary(aes(col = Condition), fun=mean, geom="line") +
  theme_classic() +
  scale_color_manual(values = c(Control = "azure3", my_palette[c("Betaine")])) +
  scale_fill_manual(values = c(Control = "azure3", my_palette[c("Betaine")])) +
  ylab("Resolubilised protein concentration") + 
  xlab("Detergent concentration [%]")

ggplot(Gnd_assay, aes(x=Det_concentration, y=Value)) +
  stat_summary(aes(col = Condition), fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", width=0.01) +
  stat_summary(aes(fill = Condition), fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="ribbon", alpha=0.2) +
  stat_summary(aes(col = Condition), fun=mean, geom="point", size=3) +
  stat_summary(aes(col = Condition), fun=mean, geom="line") +
  theme_classic() +
  scale_color_manual(values = c(Control = "azure3", my_palette["Betaine"])) +
  scale_fill_manual(values = c(Control = "azure3", my_palette["Betaine"])) +
  ylab("Resolubilised protein concentration") + 
  xlab("Detergent concentration [%]")

### proteolysis_test for Gnd 
proteolysis_data <- read.csv("/Users/moni/Documents/Phd/Collaborations/C6/Aggregation_proteolysis/220127_proteolysis_test_aggregate.csv")

p3 <- ggplot(proteolysis_data, aes(x=Time, y=y_scaled)) +
  stat_summary(aes(group = Condition), fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", width=0.5, col="black") +
  stat_summary(aes(fill = Condition), fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="ribbon", alpha=0.25) +
  stat_summary(aes(col = Condition), fun=mean, geom="line", lwd=1.2) +
  stat_summary(aes(group = Condition), fun=mean, geom="line", col="black", lwd=0.2) +
  stat_summary(aes(fill = Condition), fun=mean, geom="point", size=3, col="black", pch=21) +
  theme_classic() +
  theme(legend.text = element_text(size=15),
        legend.position = c(0.8, 0.8),
        legend.title = element_blank(),
        axis.text.x = element_text(size=15, color="black"), 
        axis.text.y = element_text(size=15,  color="black"), 
        axis.title = element_text(size=15)) +
  scale_color_manual(values = c(Control = "grey75", my_palette[c("Betaine", "Proline")])) +
  scale_fill_manual(values = c(Control = "grey75", my_palette[c("Betaine", "Proline")])) +
  ylab("Relative protein concentration") + 
  xlab("PK digestion time [min]")
ggsave(paste("/Users/moni/Documents/Phd/Osmolyte_paper/Figure5I.pdf", sep=""), plot = p3, height = 4, width=4)




