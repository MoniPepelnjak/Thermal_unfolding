### Melting point analysis for RbsK
# The script takes the melting point data from python fitting script and plots the curves as a function of ribose concentration
# Produces the Figure 4C from the Osmolyte paper

# Path directions - a file for all the validation data, to see where to take the fits from
Path_directions <- read.csv("/Users/moni/Documents/Phd/Collaborations/C6/Path_directions_validations.csv")
path = "/Users/moni/Documents/Phd/Collaborations/C6/Thermal_shift/Fitted_data/"

# LiP-MS data for comparison of the data
significant_TMAO <- readRDS("/Volumes/My Passport for Mac/data/Significant_1M.rds")

library(magrittr)
library(dplyr)
library(ggplot2)
library(drc)

Proteins_of_interest <- c("RbsK")
Uniprot_code <- c("P0A9J6")
Output_Parameters <- NULL
rbsk_validation <- TRUE

if(rbsk_validation) {
  Path_directions %<>%
    filter(Rbsk_validation == "All_conc")
}

for(i in 1:length(Proteins_of_interest)){
  experiment <- Path_directions$Experiment[Path_directions$Sample == Proteins_of_interest[i] & Path_directions$Condition_is_TMAO == FALSE ]
  
  annotation <- read.csv(paste("/Users/moni/Documents/Phd/Collaborations/C6/Thermal_shift/", experiment, "_Annotation.csv", sep="")) %>%
    filter(Sample %in% Proteins_of_interest[i])
  
  Parameters <- read.csv(paste(path, experiment, "_", Proteins_of_interest[i], "_output_parameters.csv" , sep=""))[,-1]
  Parameters_2 <- read.csv(paste(path, "13", "_", Proteins_of_interest[i], "_output_parameters.csv" , sep=""))[,-1]
  
  annotation_2 <- read.csv(paste("/Users/moni/Documents/Phd/Collaborations/C6/Thermal_shift/", "13", "_Annotation.csv", sep="")) %>%
    filter(Sample %in% Proteins_of_interest[i])
  
  Parameters %<>%
    plyr::join(.,annotation, by="Position" )
  
  Parameters_2 %<>%
    plyr::join(.,annotation_2, by="Position" )

  Parameters  %<>%
    mutate(Tm_C = Tm - 273.15) %>%
    group_by(Concentration) %>%
    mutate(Tm_control = mean(Tm_C[Condition == "Control"]), 
           sd_control = sd(Tm_C[Condition == "Control"])) %>%
    group_by(Condition, Concentration) %>%
    mutate(Tm_mean = mean(Tm_C)) %>% 
    mutate(sd_sample = sd(Tm_C)) %>%
    ungroup() %>%
    mutate(error_propagated = sqrt(sd_control^2 + sd_sample^2)) %>%
    mutate(dT = Tm_mean - Tm_control) %>%
    mutate(Protein = Uniprot_code[i]) %>%
    dplyr::select(Protein, Sample, Concentration, Condition, Tm_C, Tm_control, sd_control, Tm_mean, sd_sample, dT, error_propagated)
  
  Parameters_2  %<>%
    mutate(Tm_C = Tm - 273.15) %>%
    group_by(Concentration) %>%
    mutate(Tm_control = mean(Tm_C[Condition == "Control"]), 
           sd_control = sd(Tm_C[Condition == "Control"])) %>%
    group_by(Condition, Concentration) %>%
    mutate(Tm_mean = mean(Tm_C)) %>% 
    mutate(sd_sample = sd(Tm_C)) %>%
    ungroup() %>%
    mutate(error_propagated = sqrt(sd_control^2 + sd_sample^2)) %>%
    mutate(dT = Tm_mean - Tm_control) %>%
    mutate(Protein = Uniprot_code[i]) %>%
    dplyr::select(Protein, Sample, Concentration, Condition, Tm_C, Tm_control, sd_control, Tm_mean, sd_sample, dT, error_propagated)
  
   
  Output_Parameters <- rbind(Parameters_2, Parameters)
  
}

Stability_scores <- NULL

for(j in names(significant_TMAO)) {
  sub <- significant_TMAO[[j]]$Protein_level[significant_TMAO[[j]]$Protein_level$Protein %in% Uniprot_code,] %>%
    mutate(Condition = j)
  
  Stability_scores <- rbind(Stability_scores, sub)
}

ggplot(Output_Parameters[Output_Parameters$Condition %in% c("TMAO", "Glycerol", "Glucose"),], aes(x=Concentration, y=dT, col=(Concentration), shape=Condition)) +
  geom_point(size=3) +
  facet_wrap(~Condition)

out_data_raw <- NULL
out_predictions <- NULL
sm_list <- c("TMAO", "Glycerol", "Glucose")
for(sm in sm_list){
  sm
  compound_y <- Output_Parameters$Tm_C[Output_Parameters$Condition == sm]
  compound_x <- Output_Parameters$Concentration[Output_Parameters$Condition == sm]
  
  control_y <- Output_Parameters$Tm_C[Output_Parameters$Condition == "Control"]
  control_x <- Output_Parameters$Concentration[Output_Parameters$Condition == "Control"]
  
  m1 <- drm(compound_y ~ compound_x, fct = MM.3())
  m2 <- drm(control_y ~ control_x, fct = MM.3())
  
  # new dose levels as support for the line
  newdata_sm <- expand.grid(conc=seq(0, 6, length=600))
  data_raw_sm <- m1$data[,1:2] %>%
    `colnames<-`(c("x", "y"))
  data_raw_sm$Condition <- sm
  
  newdata_control <- expand.grid(conc=seq(0, 6, length=600))
  data_raw_control <- m2$data[,1:2] %>%
    `colnames<-`(c("x", "y"))
  data_raw_control$Condition <- "Control"
  
  # predictions and confidence intervals
  pm <- predict(m1, newdata=newdata_sm, interval="confidence")
  pm_control <- predict(m2, newdata=newdata_control, interval="confidence")
  
  # new data with predictions
  newdata_sm$p <- pm[,1]
  newdata_sm$pmin <- pm[,2]
  newdata_sm$pmax <- pm[,3]
  newdata_sm$Condition <- sm
  
  newdata_control$p <- pm_control[,1]
  newdata_control$pmin <- pm_control[,2]
  newdata_control$pmax <- pm_control[,3]
  newdata_control$Condition <- "Control"
  
  data_raw <- rbind(data_raw_sm, data_raw_control)
  data_raw$Experiment <- sm
  newdata <- rbind(newdata_control, newdata_sm)
  newdata$Experiment <- sm
  
  out_data_raw <- rbind(out_data_raw, data_raw)
  out_predictions <- rbind(out_predictions, newdata)

}

# my_palette = c(Betaine = "#B5C3AD", 
#                Trehalose = "#d4a993", 
#                TMAO = "#87A9C6", 
#                Glycerol = "#a997a2", 
#                Proline = "#edd9bb", 
#                Glucose = "#C78C8C",
#                Betaine_2 = "#8E9C86",
#                TMAO_2 = "#718CA2",
#                Glycerol_2 = "#86697B",
#                Proline_2 = "#D8BF9A", 
#                Control = "gray65")

my_palette = c(Control = "gray65",
               control = "#gra65",
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

# plot curve
library(ggplot2)
# need to shift conc == 0 a bit up, otherwise there are problems with coord_trans
#ryegrass$conc0 <- ryegrass$conc
#ryegrass$conc0[ryegrass$conc0 == 0] <- 0.5
# plotting the curve
Labels_rbsk <- out_data_raw %>%
  dplyr::select(Experiment) %>%
  unique()

Figure_4C <- ggplot(out_data_raw, aes(x = x, y = y)) +
  #geom_point(aes(col=Condition)) + 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", width=0.2, aes(color=Condition, group=Condition), col="black") +
  stat_summary(fun=mean, geom="point",aes(fill=Condition), size=2, col="black", pch=21)  +
  geom_ribbon(data=out_predictions, aes(x=conc, y=p, ymin=pmin, ymax=pmax, fill=Condition), alpha=0.5) +
  geom_line(data=out_predictions, aes(x=conc, y=p, col=Condition), lwd=0.8) +
  facet_wrap(~Experiment, scales="free_y")+
  #ylim(c(0, 75)) +
  xlab("Ribose concentration (mM)") + ylab("Tm of RbsK") +
  scale_fill_manual(values = my_palette[c("Glycerol", "Control", "TMAO", "Glucose")]) +
  scale_color_manual(values = my_palette[c("Glycerol", "Control", "TMAO", "Glucose")]) +
  geom_text(data =Labels_rbsk,  aes(label = Experiment), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size=5.5) +
  #theme_minimal() +
  theme(axis.title = element_text(size = 15, color="black"),
        text=element_text(color="black"),
        axis.text = element_text(size=15, color="black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.background = element_blank(), axis.line = element_line(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), legend.position = "none") +
  ylim(c(min(Output_Parameters$Tm_C), max(Output_Parameters$Tm_C+1))) +
   theme(panel.background = element_blank(), 
         axis.line = element_line())

Labels_rbsk$one <- "One"
Figure_4C_short <- ggplot(out_data_raw[out_data_raw$Experiment %in% c("Glucose", "TMAO"),], aes(x = x, y = y)) +
  #geom_point(aes(col=Condition)) + 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", width=0.2, aes(color=Condition, group=Condition), col="black") +
  stat_summary(fun=mean, geom="point",aes(fill=Condition), size=2, col="black", pch=21)  +
  geom_ribbon(data=out_predictions[out_predictions$Experiment %in% c("Glucose", "TMAO"), ], aes(x=conc, y=p, ymin=pmin, ymax=pmax, fill=Condition), alpha=0.5) +
  geom_line(data=out_predictions[out_predictions$Experiment %in% c("Glucose", "TMAO"), ], aes(x=conc, y=p, col=Condition), lwd=0.8) +
  facet_wrap(~Experiment, scales="free_y")+
  #ylim(c(0, 75)) +
  xlab("Ribose concentration (mM)") + ylab("Tm of RbsK") +
  scale_fill_manual(values = my_palette[c("Glycerol", "Control", "TMAO", "Glucose")]) +
  scale_color_manual(values = my_palette[c("Glycerol", "Control", "TMAO", "Glucose")]) +
  geom_text(data =Labels_rbsk[c(1,3), ],  aes(label = Experiment), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size=5.5) +
  #theme_minimal() +
  theme(axis.title = element_text(size = 15, color="black"),
        text=element_text(color="black"),
        axis.text = element_text(size=15, color="black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.background = element_blank(), axis.line = element_line(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), legend.position = "none") +
  ylim(c(min(Output_Parameters$Tm_C), max(Output_Parameters$Tm_C+1))) +
  theme(panel.background = element_blank(), 
        axis.line = element_line())
ggsave("/Users/moni/Documents/Phd/Osmolyte_paper/Figure_4C_short.pdf", Figure_4C_short, height=3, widt=6)

Figure_4C_empty <- ggplot(out_data_raw[out_data_raw$Condition == "Control" & out_data_raw$Experiment %in% c("Glucose", "TMAO"),], aes(x = x, y = y)) +
  #geom_point(aes(col=Condition)) + 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", width=0.2, aes(color=Condition, group=Condition), col="black") +
  stat_summary(fun=mean, geom="point",aes(fill=Condition), size=2, col="black", pch=21)  +
  geom_ribbon(data=out_predictions[out_predictions$Condition == "Control" & out_predictions$Experiment %in% c("Glucose", "TMAO"), ], aes(x=conc, y=p, ymin=pmin, ymax=pmax, fill=Condition), alpha=0.5) +
  geom_line(data=out_predictions[out_predictions$Condition == "Control"& out_predictions$Experiment %in% c("Glucose", "TMAO"), ], aes(x=conc, y=p, col=Condition), lwd=0.8) +
  facet_wrap(~Experiment, scales="free_y")+
  #ylim(c(0, 75)) +
  xlab("Ribose concentration (mM)") + ylab("Tm of RbsK") +
  scale_fill_manual(values = my_palette[c("Glycerol", "Control", "TMAO", "Glucose")]) +
  scale_color_manual(values = my_palette[c("Glycerol", "Control", "TMAO", "Glucose")]) +
  #geom_text(data =Labels_rbsk,  aes(label = Experiment), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size=5.5) +
  #theme_minimal() +
  theme(axis.title = element_text(size = 15, color="black"),
        text=element_text(color="black"),
        axis.text = element_text(size=15, color="black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.background = element_blank(), axis.line = element_line(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), legend.position = "none") +
  ylim(c(min(Output_Parameters$Tm_C), max(Output_Parameters$Tm_C+1))) +
  theme(panel.background = element_blank(), 
        axis.line = element_line())
ggsave("/Users/moni/Documents/Phd/Osmolyte_paper/Figure_4C_empty.pdf", Figure_4C_empty, height=3, widt=6)


p <- ggplot(out_data_raw[out_data_raw$Experiment %in% c("Glucose", "TMAO"),], aes(x = x, y = y, fill=Condition)) +
  geom_point(pch=21, col="black", size=5) +
  #scale_fill_manual(values = my_palette[c( "Control", "Glucose", "Glycerol", "TMAO")]) +
  scale_fill_manual(values = my_palette[c( "Control", "Glucose", "TMAO")]) +
  theme(legend.text = element_text(size=10), 
        legend.position = "bottom", 
        legend.title=element_blank())

leg <- ggpubr::get_legend(p)

pdf("/Users/moni/Documents/Phd/Osmolyte_paper/Color_palette_RbsK_limited.pdf", width=4, height=0.25)
ggpubr::as_ggplot(leg)
dev.off()

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

tiff("/Users/moni/Documents/Phd/Experiments/Final_data/All_compounds_correlation_combined.tiff", units="in", width=4, height=4, res=300)
ggplot(Output_unique[Output_unique$Condition != "Control",], aes(x=Protein_stabilisation, y=dT)) +
  geom_smooth(method="lm", col="gray40", lty=20) +
  geom_point( aes(shape = Sample, col=Condition), size=3 ) +
  geom_errorbar(aes(ymin = (dT - error_propagated), ymax=(dT + error_propagated ), col = Condition)) +
  scale_color_manual(values=my_palette[1:6]) +
  geom_label(x=20, y=2.5,aes(label=paste("R2 = ", R2, sep=""))) +
  theme(panel.grid = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), panel.background = element_blank()) +
  ylab(paste("\U2206Tm",  "[", "\U00B0", "C]")) +
  xlab("Stability score")
dev.off()


tiff("/Users/moni/Documents/Phd/Experiments/Final_data/All_compounds_correlation.tiff", units="in", width=7, height=3, res=300)
ggplot(Output_unique[Output_unique$Condition != "Control",], aes(x=Protein_stabilisation, y=dT)) +
  geom_smooth(method="lm", col="gray40", lty=20) +
  geom_point( aes(shape = Sample, col=Condition), size=3 ) +
  geom_errorbar(aes(ymin = (dT - error_propagated), ymax=(dT + error_propagated ), col = Condition)) +
  scale_color_manual(values=my_palette[1:6]) +
  facet_wrap(~Sample) +
  geom_label(x=20, y=2.5, data = Output_unique_R, aes(label=paste("R2 = ", R2, sep=""))) +
  theme(panel.grid = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), panel.background = element_blank()) +
  ylab(paste("\U2206Tm",  "[", "\U00B0", "C]")) +
  xlab("Stability score")
dev.off()


ggplot(Output_unique[Output_unique$Condition != "Control",], aes(x=Protein_aggregation, y=dT)) +
  geom_smooth(method="lm", col="gray40", lty=20) +
  geom_point( aes(shape = Sample, col=Condition), size=3 ) +
  geom_errorbar(aes(ymin = (dT - error_propagated), ymax=(dT + error_propagated ), col = Condition)) +
  scale_color_manual(values=my_palette[1:6]) +
  geom_label(x=20, y=2.5,aes(label=paste("R2 = ", R2, sep=""))) +
  theme(panel.grid = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), panel.background = element_blank()) +
  ylab(paste("\U2206Tm",  "[", "\U00B0", "C]")) +
  xlab("Stability score")

ggplot(Output_unique[Output_unique$Condition != "Control",], aes(x=Coverage, y=dT)) +
  geom_smooth(method="lm", col="gray40", lty=20) +
  geom_point( aes(shape = Sample, col=Condition), size=3 ) +
  geom_errorbar(aes(ymin = (dT - error_propagated), ymax=(dT + error_propagated ), col = Condition)) +
  scale_color_manual(values=my_palette[1:6]) +
  geom_label(x=20, y=2.5,aes(label=paste("R2 = ", R2, sep=""))) +
  theme(panel.grid = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(), panel.background = element_blank()) +
  ylab(paste("\U2206Tm",  "[", "\U00B0", "C]")) +
  xlab("Stability score")
