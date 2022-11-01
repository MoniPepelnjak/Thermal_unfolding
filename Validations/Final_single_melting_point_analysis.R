### Final single melting point analysis
## The goal of the script is to extract the data from the fitted values 

# Load libraries
library(ggpubr)
library(ggplot2)
library(dplyr)
library(magrittr)

my_palette = c(Control = "#8D8D92",
               Betaine = "#B5C3AD", 
               Trehalose = "#d4a993", 
               TMAO = "#87A9C6", 
               Glycerol = "#a997a2", 
               Proline = "#edd9bb", 
               Glucose = "#C78C8C")

# Load in the path-directions for the experiments
path = "/Volumes/My Passport for Mac/Final_data"

#path = "/Users/moni/Documents/Phd/Collaborations/C6/Thermal_shift/Fitted_data/"
Path_directions <- read.csv(paste(path, "/Validations/C6/Thermal_shift/Path_directions_validations.csv", sep=""))
significant_all <- readRDS(paste(path, "/LiP/Significant_both_quantiles.rds", sep=""))

Path_directions %<>%
  filter(Rbsk_validation == FALSE)

Proteins_of_interest <- c("Eno", "Gnd", "RbsK")
Uniprot_code <- c("P0A6P9", "P00350", "P0A9J6")
sm_list <- c("TMAO", "Betaine", "Glycerol", "Proline", "Trehalose", "Glucose")
Output_Parameters <- NULL
for(i in 1:length(Proteins_of_interest)){
  experiment <- Path_directions$Experiment[Path_directions$Sample == Proteins_of_interest[i] & Path_directions$Condition_is_TMAO == FALSE ]
  annotation <- read.csv(paste(path, "/Validations/C6/Thermal_shift/", experiment, "_Annotation.csv", sep="")) %>%
    filter(Sample %in% Proteins_of_interest[i])
  
  Parameters <- read.csv(paste(path, "/Validations/C6/Thermal_shift/Fitted_data/", experiment, "_", Proteins_of_interest[i], "_output_parameters.csv" , sep=""))[,-1]
  
  Parameters %<>%
    plyr::join(.,annotation, by="Position" )%>%
    group_by(Position) %>%
    mutate(Conc = ifelse("Concentration" %in% names(.), Concentration, 0) ) %>%
    ungroup() %>%
    filter(Condition != "TMAO", 
           Conc %in% c(0,1))
  
  experiment_TMAO <- Path_directions$Experiment[Path_directions$Sample == Proteins_of_interest[i] & Path_directions$Condition_is_TMAO == TRUE ]
  annotation_TMAO <- read.csv(paste(path, "/Validations/C6/Thermal_shift/", experiment_TMAO, "_Annotation.csv", sep="")) %>%
    filter(Sample %in% Proteins_of_interest[i])
  
  Parameters_TMAO <- read.csv(paste(path, "/Validations/C6/Thermal_shift/Fitted_data/", experiment_TMAO, "_", Proteins_of_interest[i], "_output_parameters.csv" , sep=""))[,-1]

  Parameters_TMAO %<>%
    plyr::join(.,annotation_TMAO, by="Position" )%>%
    group_by(Position) %>%
    mutate(Conc = ifelse("Concentration" %in% names(.), Concentration, 0) ) %>%
    ungroup() %>%
    filter(Condition %in% c("TMAO", "Control", "Glucose", "Trehalose"), 
           Conc %in% c(0,1))
  
  Parameters  %<>%
    mutate(Tm_C = Tm - 273.15) %>%
    mutate(Tm_control = mean(Tm_C[Condition == "Control"]), 
           sd_control = sd(Tm_C[Condition == "Control"])) %>%
    group_by(Condition) %>%
    mutate(Tm_mean = mean(Tm_C)) %>% 
    mutate(sd_sample = sd(Tm_C)) %>%
    ungroup() %>%
    mutate(error_propagated = sqrt(sd_control^2 + sd_sample^2)) %>%
    mutate(dT = Tm_mean - Tm_control) %>%
    mutate(Protein = Uniprot_code[i]) %>%
    dplyr::select(Protein, Sample, Condition, Tm_C, Tm_control, sd_control, Tm_mean, sd_sample, dT, error_propagated)
  
  Parameters_TMAO  %<>%
    mutate(Tm_C = Tm - 273.15) %>%
    mutate(Tm_control = mean(Tm_C[Condition == "Control"]), 
           sd_control = sd(Tm_C[Condition == "Control"])) %>%
    group_by(Condition) %>%
    mutate(Tm_mean = mean(Tm_C)) %>% 
    mutate(sd_sample = sd(Tm_C)) %>%
    ungroup() %>%
    mutate(error_propagated = sqrt(sd_control^2 + sd_sample^2)) %>%
    mutate(dT = Tm_mean - Tm_control)%>%
    mutate(Protein = Uniprot_code[i]) %>%
    dplyr::select(Protein, Sample, Condition, Tm_C, Tm_control, sd_control, Tm_mean, sd_sample, dT, error_propagated)
  
  
  Parameters <- rbind(Parameters, Parameters_TMAO) 
  
  Output_Parameters <- rbind(Output_Parameters, Parameters)
  
}

Stability_scores <- NULL
significant_all$Glycerol$Protein_level <- significant_all$Glycerol$Protein_level %>%
  dplyr::select(-mean_p_score)

for(j in names(significant_all)) {
  sub <- significant_all[[j]]$Protein_level[significant_all[[j]]$Protein_level$Protein %in% Uniprot_code,] %>%
    mutate(Condition = j)
  
  Stability_scores <- rbind(Stability_scores, sub)
}

Output_Parameters %<>%
  plyr::join(., Stability_scores, by=c("Protein", "Condition")) %>%
  unique()

Output_unique <- Output_Parameters %>%
  dplyr::select(Protein, Sample, Condition, Protein_stabilisation, dT, error_propagated, Protein_aggregation, Protein_binding, Coverage) %>%
  unique()

R2 <- round(summary(lm(Output_unique$dT~Output_unique$Protein_stabilisation))$adj.r.squared, 2)

Output_unique_R <- Output_unique %>%
group_by(Sample) %>%
  summarise(R2 = round(summary(lm(Protein_stabilisation~dT))$adj.r.squared, 2))
R2 <- data.frame(R=R2)
#tiff("/Users/moni/Documents/Phd/Experiments/Final_data/All_compounds_correlation_combined.tiff", units="in", width=4, height=4, res=300)
cairo_pdf(paste(path, "/Validations/Figures/Figure_3F_bottom.pdf", sep=""), height=4, width=4.5)
Figure_3F <- ggplot(Output_unique[Output_unique$Condition != "Control",], aes(x=Protein_stabilisation, y=dT)) +
  geom_smooth(method="lm", col="gray40", lty=20) +
  geom_errorbar(aes(ymin = (dT - error_propagated), ymax=(dT + error_propagated )), col="black") +
  geom_errorbar(aes(ymin = (dT - error_propagated), ymax=(dT + error_propagated )), col="black") +
  geom_point(aes(fill=Condition, shape=Sample, colour = Condition), colour="black", size=4.5) + 
  scale_shape_manual(values=c(21, 22, 24))+
  scale_fill_manual(values = my_palette[sm_list]) +  
  #geom_text(data =Labels_rbsk,  aes(label = Experiment), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size=5.5) +
  #theme_minimal() +
  geom_text(x=Inf, y=-Inf, hjust=1, vjust=-0.5, data = R2, aes(label=paste("R2 = ", R, sep="")), size=7) +
  theme(axis.title = element_text(size = 22, color="black"),
        text=element_text(color="black"),
        axis.text = element_text(size=22, color="black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.background = element_blank(), axis.line = element_line(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), legend.position = "none") +
  #scale_color_manual(values = my_palette[1:6]) +  
  theme(panel.grid = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), panel.background = element_blank()) +
  ylab(paste("\U2206Tm",  "[", "\U00B0", "C]")) +
  xlab("Stability score") 
print(Figure_3F)
dev.off()

cor.test(Output_unique[Output_unique$Condition != "Control",]$Protein_stabilisation, Output_unique[Output_unique$Condition != "Control",]$dT)
cor.test(Output_unique[Output_unique$Condition != "Control" & Output_unique$Protein == "P0A9J6",]$Protein_stabilisation, Output_unique[Output_unique$Condition != "Control" & Output_unique$Protein == "P0A9J6",]$dT)
Labels_prot <- Output_unique %>%
  dplyr::select(Sample) %>%
  unique()
cairo_pdf(paste(path, "/Validations/Figures/Figure_3E_bottom.pdf", sep=""), height=4, width=10)
all_shapes <- ggplot(Output_unique[Output_unique$Condition != "Control",], aes(x=Protein_stabilisation, y=dT)) +
  geom_smooth(method="lm", col="gray40", lty=20) +
  geom_errorbar(aes(ymin = (dT - error_propagated), ymax=(dT + error_propagated )), col="black") +
  geom_errorbar(aes(ymin = (dT - error_propagated), ymax=(dT + error_propagated )), col="black") +
  geom_point(aes(fill=Condition, shape=Sample, colour = Condition), colour="black", size=4.5) + 
  scale_fill_manual(values = my_palette[sm_list], guide="none") +  
  scale_shape_manual(values=c(21, 22, 24))+
  geom_text(data =Labels_prot,  aes(label = Sample), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size=8) +
  facet_wrap(~Sample, scales="free_y") +
  ylim(c(0,14)) +
  theme(axis.title = element_text(size = 22, color="black"),
        text=element_text(color="black"),
        axis.text = element_text(size=22, color="black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.background = element_blank(), axis.line = element_line(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), legend.position = "none") +
  geom_text(x=Inf, y=-Inf, hjust=1, vjust=-0.5, data = Output_unique_R, aes(label=paste("R2 = ", R2, sep="")), size=7) +
  theme(panel.grid = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), panel.background = element_blank()) +
  ylab(paste("\U2206Tm",  "[", "\U00B0", "C]")) +
  xlab("Stability score")
all_shapes
dev.off()

circ <- ggplot(Output_unique[Output_unique$Condition != "Control",], aes(x=Protein_stabilisation, y=dT)) +
  geom_smooth(method="lm", col="gray40", lty=20) +
  geom_errorbar(aes(ymin = (dT - error_propagated), ymax=(dT + error_propagated )), col="black") +
  geom_errorbar(aes(ymin = (dT - error_propagated), ymax=(dT + error_propagated )), col="black") +
  geom_point(aes(fill=Condition, colour = Condition), pch = 21,colour="black", size=4) + 
  scale_fill_manual(values = my_palette[sm_list]) +  
  #scale_shape_manual(values=c(21, 22, 24))+
  facet_wrap(~Sample) +
  geom_label(x=20, y=2.5, data = Output_unique_R, aes(label=paste("R2 = ", R2, sep=""))) +
  theme(panel.grid = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), panel.background = element_blank()) +
  ylab(paste("\U2206Tm",  "[", "\U00B0", "C]")) +
  xlab("Stability score")

# Extract the legend. Returns a gtable
leg <- get_legend(circ, position = "right")
leg_shapes <- get_legend(all_shapes, position = "right")
# Convert to a ggplot and print
tiff(paste(path, "/Validations/Figures/Color_palette_right.tiff", sep=""), units="in", width=1, height=2, res=300)
as_ggplot(leg)
dev.off()

tiff(paste(path, "/Validations/Figures/Color_palette_shape_right.tiff",sep=""), units="in", width=1, height=1.5, res=300)
as_ggplot(leg_shapes)
dev.off()

ggplot(Output_unique[Output_unique$Condition != "Control",], aes(x=Protein_aggregation, y=dT)) +
  geom_smooth(method="lm", col="gray40", lty=20) +
  geom_errorbar(aes(ymin = (dT - error_propagated), ymax=(dT + error_propagated )), col="black") +
  geom_errorbar(aes(ymin = (dT - error_propagated), ymax=(dT + error_propagated )), col="black") +
  geom_point(aes(fill=Condition, shape=Sample, colour = Condition), colour="black", size=5) + 
  scale_fill_manual(values = my_palette[1:6]) +  
  scale_shape_manual(values=c(21, 22, 24))+
  scale_color_manual(values=my_palette[1:6]) +
  geom_label(x=20, y=2.5,aes(label=paste("R2 = ", R2, sep=""))) +
  theme(panel.grid = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), panel.background = element_blank()) +
  ylab(paste("\U2206Tm",  "[", "\U00B0", "C]")) +
  xlab("Stability score")

ggplot(Output_unique[Output_unique$Condition != "Control",], aes(x=Coverage, y=dT)) +
  geom_smooth(method="lm", col="gray40", lty=20) +
  geom_point(aes(fill=Condition, shape=Sample), colour="black",pch=21, size=5) + 
  scale_fill_manual(values = my_palette[1:6]) +  geom_errorbar(aes(ymin = (dT - error_propagated), ymax=(dT + error_propagated ), col = Condition)) +
  scale_color_manual(values=my_palette[1:6]) +
  geom_label(x=20, y=2.5,aes(label=paste("R2 = ", R2, sep=""))) +
  theme(panel.grid = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), panel.background = element_blank()) +
  ylab(paste("\U2206Tm",  "[", "\U00B0", "C]")) +
  xlab("Stability score")

export_table <- Output_unique %>%
  dplyr::select(Protein, Sample, Condition, dT, error_propagated) %>%
  unique()

write.csv(export_table, paste(path, "/Validations/C6/Thermal_shift/dT_export_proteins.csv", sep=""))

  