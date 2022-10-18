# Script for fitting the DSF data, analysing the lysate melting temperatures and comparing it to LiP data
# The script was used to fit and plot the data for figures 2D, 2E and S2D of the osmolyte paper
# Author: Monika Pepelnjak

### Define the function for fitting the flourescence data  ####
trcFunc <- function(x,a,b,c,d,H,Tm){
  R = 8.31446261815324
  ((a+b*x) + ((c + d*x) * exp((log10(H)*(x-Tm))/x*R*Tm)))/(1 + exp((log10(H)*(x-Tm))/x*R*Tm)) 
  
  }

### load the libaries
library(minpack.lm)
library(tidyverse)
library(investr)
library(magrittr)
library(ggplot2)

# Define the general path
Path <- "/Users/moni/Documents/Phd/Collaborations/C6/Thermal_shift/"

# Define the experiment
experiment <- c("01", "02")

# Load in the data from RFU measurements
results_1 <- read.csv(paste(Path, experiment[1],"_RFU_data.csv"))[,-1]
results_2 <- read.csv(paste(Path, experiment[2],"_RFU_data.csv"))[,-1]

# Load in annotation data from RFU measurements
annotation_data_1 <- read.csv(paste(Path, experiment[1], "_Annotation.csv")) %>%
  mutate(Position = paste(Position, "1", sep="_"))
annotation_data_2 <- read.csv(paste(Path, experiment[2], "_Annotation.csv")) %>%
  mutate(Position = paste(Position, "2", sep="_"))
annotation_data <- rbind(annotation_data_1,annotation_data_2)

# Load in LiP-MS significance data
significant_TMAO <- readRDS("/Volumes/My Passport for Mac/data/Significant_both_quantiles.rds")

# Load in Savitski data
Savitski_data <- read.csv("/Users/moni/polybox/code_review/Monika/Savitski_melting_temp.csv")

# Palette for plotting
my_palette = c(Betaine = "#B5C3AD", 
               Trehalose = "#d4a993", 
               TMAO = "#87A9C6", 
               Glycerol = "#a997a2", 
               Proline = "#edd9bb", 
               Glucose = "#C78C8C")

# Define the osmolytes analysed
sm_list <- c("TMAO","Betaine", "Glycerol", "Proline", "Trehalose", "Glucose")

# Define the sample analysed - this should match the sample from the annotation file 
sample_analysed <- "lysate"

# Combine the two experiments - this part is not necessary if all the data is in the single file
y_all_1 <- results_1 %>%
  reshape2::melt(., id.var="Temperature") %>%
  `colnames<-`(c("Temperature", "Position", "RFU")) %>%
  mutate(Position = paste(Position, "1", sep="_")) %>%
  plyr::join(., annotation_data_1, by="Position") %>%
  filter(! Condition %in% c("TMAO", "Control", "Glucose"))

y_all_2 <- results_2 %>%
    reshape2::melt(., id.var="Temperature") %>%
    `colnames<-`(c("Temperature", "Position", "RFU")) %>%
    mutate(Position = paste(Position, "2", sep="_")) %>%
  plyr::join(., annotation_data_2, by="Position")   
  
y_all <- rbind(y_all_1, y_all_2)

Negative_control <- y_all[y_all$Condition == "Negative",]

# Transform the dataset:
y_all %<>%
group_by(Position) %>%
  mutate(Temperature = Temperature + 273.15) %>% # For the fitting equation to work, the temperature has to be in K
  filter(Temperature < (80 + 273.15) & Temperature > (37 + 273.15), # Remove the extreme temperatures - the flourescence gets a bit weird, so makes sense to look at it only at our temperatures of interest
         Sample == sample_analysed) %>% # only select the sample you want to analyse - this is needed in case you measure different samples in the same experiment
  filter(!Condition %in% "Negative") %>% # Remove the negative control
  mutate(y_scaled = ((RFU-min(RFU))/(max(RFU)-min(RFU)))) %>% # min-max scaling of the data 
  filter(Position != "E12_1") # remove a clear outlier - this is hardcoded for this data, should be taken out for others

all_positions <- y_all$Position %>% unique()

interval_out <- NULL
all_fits <- list()
for( i in all_positions){
  y <- y_all$y_scaled[y_all$Position %in% i] # get the x
  x <- y_all$Temperature[y_all$Position %in% i] # get the y for the model
  
  # fit the model
  # The starting parameters here were optimised for this specific case, so that the model converges
  # The starting parameters have to be changed for a different curve
  # This could be optimised more to find the appropriate starting parameters
  model <-nlsLM(y~trcFunc(x, a, b, c, d, H, Tm),start=list(a=max(y),b=-1,c=min(y), d=1, H =100, Tm = 330))
  
  all_fits[[i]] <- model
  
  # Design a new df
  new.data <- data.frame(x=seq((37 + 273.15), (80 + 273.15), by = 0.1))
  # predict the fit and the confidence intervals
  interval <- as_tibble(predFit(model, newdata = new.data, interval = "confidence", level= 0.99)) %>% 
    mutate(x = new.data$x, 
           Position = i)
  
  # Export the fits
  interval_out <- rbind(interval_out, interval)
  
  
}

# Join with the annotation data
interval_out %<>%
  plyr::join(., annotation_data, by="Position")
 
t <- lapply(all_fits, function(x) coef(x)) # get the coefficients from the model such that you can extract the Tm
Melting_points <- data.frame(Position = names(t),
                             Tm = unlist(lapply(t, function(x) x[6])),  # This is hard coded for this specific fitting model! Be careful if you change the model, this is also different!
                             H = unlist(lapply(t, function(x) x[5]))) %>% # This is hard coded for this specific fitting model! Be careful if you change the model, this is also different!
  mutate(Tm_C = Tm - 273.15) %>% # Get melting temperature in ??C
  plyr::join(., annotation_data, by="Position") %>%  # Combine with annotation
  mutate(Tm_control = mean(Tm_C[Condition == "Control"])) %>% # Get the mean temperature for control
  group_by(Condition) %>%
  mutate(dT = mean(Tm_C) -  Tm_control, 
         Tm_mean = mean(Tm_C))

# Calculate the mean stabilisation score from LiP-MS data
mean_all <- lapply(significant_TMAO, function(x) mean(x$Protein_level$Protein_stabilisation))

overview <- data.frame(Condition = names(mean_all), 
                       Mean_all = unlist(mean_all))

# Combine melting points with the LiP-MS data
Melting_points_sub <- Melting_points %>%
  dplyr::select(Condition, Tm_C, dT) %>%
  unique() %>%
  plyr::join(., overview, by="Condition") 

# Remove the TMAO at high concentration - we do not look at this here - again hard-coded for this experiment!
Melting_points_sub %<>%
  filter(! Condition %in% c("TMAO_2"))

# Assign 0 to Controls
Melting_points_sub[is.na(Melting_points_sub)] <- 0

fit <- lm(Melting_points_sub$Tm_C~Melting_points_sub$Mean_all) # fit a linear model between LiP-MS stabilisation score and DSF data
s <- round(summary(fit)$adj.r.squared, 2) # Calculate R2 from the model

# Plot the figure 2E of osmolyte paper
pdf("/Users/moni/Documents/Phd/Osmolyte_paper/Figure_2E.pdf", width=5, height=4)
ggplot(Melting_points_sub, aes(x=Mean_all, y=Tm_C)) +
  geom_smooth(method="lm", lty="dashed", col="gray65", fill="gray85") +
  stat_summary(col="black", fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", width=0.2) +
  stat_summary(aes(fill = Condition), col="black", pch=21, fun=mean, geom="point", size=5) +
  theme_classic() +
  theme(axis.title = element_text(size = 18, color="black"),
        text=element_text(color="black", size=20),
        axis.text = element_text(size=18, color="black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.background = element_blank(), axis.line = element_line(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  scale_fill_manual(values = c(Control = "azure3", my_palette[sm_list])) +
  ylab(paste("Lysate Tm ", "[", "\U00B0", "C]", sep="" )) + 
  xlab("Mean stabilisation score")
dev.off()

# Duplicate the point such that you can plot control and osmolyte in each plot
y_all_duplicated <- NULL
interval_out_2 <- NULL
for(sm in sm_list){
  y_combined <- y_all %>%
    filter(Condition %in% c("Control", sm)) %>%
    mutate(Con = sm)
  y_all_duplicated <- rbind(y_all_duplicated, y_combined)
  
  interval_out_test <- interval_out %>%
    filter(Condition %in% c("Control", sm)) %>%
    mutate(Con = sm)
  
  interval_out_2 <- rbind(interval_out_test, interval_out_2)
  
}

# Get a small df for Condition and dT to plot
Melting_points_2 <- Melting_points_sub %>%
  dplyr::select(Condition, dT) %>%
  filter(Condition %in% sm_list) %>%
  unique() %>%
  `colnames<-`(c("Con", "dT"))

# Order the osmolytes in the same order as always
y_all_duplicated$Con <- factor(y_all_duplicated$Con, levels=c(sm_list))
Melting_points_2$Con <- factor(Melting_points_2$Con, levels=c(sm_list))

# Plot figure 2D
p1 <- ggplot(data = y_all_duplicated) +  
  stat_summary(aes(group=Condition, x=Temperature-273.15, y=y_scaled, col=Condition), fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", width=0.5, alpha=0.5) +
  stat_summary(aes(col = Condition, x=Temperature-273.15, y=y_scaled), fun=mean, geom="point", size=0.5,  alpha=0.5) +
  
  xlab(paste("Temperature ", "[", "\U00B0", "C]", sep="" )) + 
  ylab("Normalised fluorescence") +
  facet_wrap(~Con, scales="free") +
  ylim(c(-0.1,1.1)) +
  geom_label(data = Melting_points_2, aes(label = paste("\U2206Tm = ", round(dT,2), "\U00B0", "C", sep=""), x=70, y= 0.15, group = Con,fill=Con))

# Add labels for plot names
Labels <- interval_out_2 %>%
  ungroup() %>%
  dplyr::select(Con) %>%
  unique()

# Calculate the ribbon values
interval_out_2 %<>%
  group_by(Condition, Con, x) %>% 
  mutate(Median_fit = median(fit)) %>%
  mutate(low_fit = min(fit), 
         high_fit = max(fit))

# Put all the osmolytes in the same order again
# If one of the data frames is not in the same order, the ordering does not work
interval_out_2$Con <- factor(interval_out_2$Con, levels=c(sm_list))
Labels$Con <- factor(Labels$Con, levels=c(sm_list))

# Plot and save figure 2D
cairo_pdf("/Users/moni/Documents/Phd/Osmolyte_paper/Figure_2D.pdf", width=8, height=3, family="sans")
Figure_2D <- p1+
  geom_line(data=interval_out_2, aes(x = (x-273.15), y = Median_fit,col=Condition, group = Position ), lwd=0.6)+
  geom_ribbon(data=interval_out_2, aes(x=(x-273.15), ymin=low_fit, ymax=high_fit, fill=Condition), alpha=0.5, inherit.aes=F)+
  scale_color_manual(values=c("Control" = "gray75", my_palette[sm_list])) +
  scale_fill_manual(values=c("Control" = "gray75", my_palette[sm_list])) +
  geom_text(data =Labels,  aes(label = Con), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size=4.5) +
  #theme_minimal() +
  theme(axis.title = element_text(size = 15, color="black"),
        text=element_text(color="black", family="sans"),
        axis.text = element_text(size=10, color="black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.background = element_blank(), axis.line = element_line(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), legend.position = "none")
print(Figure_2D)
dev.off()

### Plot an empty plot just with controls (For presentations) #### 
p1_empty <- ggplot(data = y_all_duplicated[y_all_duplicated$Condition == "Control",]) +  
  stat_summary(aes(group=Condition, x=Temperature-273.15, y=y_scaled, col=Condition), fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", width=0.5, alpha=0.5) +
  stat_summary(aes(col = Condition, x=Temperature-273.15, y=y_scaled), fun=mean, geom="point", size=0.5,  alpha=0.5) +
  
  xlab(paste("Temperature ", "[", "\U00B0", "C]", sep="" )) + 
  ylab("Normalised fluorescence") +
  facet_wrap(~Con, scales="free") +
  ylim(c(-0.1,1.1))

cairo_pdf("/Users/moni/Documents/Phd/Osmolyte_paper/Figure_2D_empty.pdf", width=8, height=3, family="sans")
Figure_2D_empty <- p1_empty+
  geom_line(data=interval_out_2[interval_out_2$Condition == "Control",], aes(x = (x-273.15), y = Median_fit,col=Condition, group = Position ), lwd=0.6)+
  geom_ribbon(data=interval_out_2[interval_out_2$Condition == "Control",], aes(x=(x-273.15), ymin=low_fit, ymax=high_fit, fill=Condition), alpha=0.5, inherit.aes=F)+
  scale_color_manual(values=c("Control" = "gray75", my_palette[sm_list])) +
  scale_fill_manual(values=c("Control" = "gray75", my_palette[sm_list])) +
  #geom_text(data =Labels,  aes(label = Con), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size=4.5) +
  #theme_minimal() +
  theme(axis.title = element_text(size = 15, color="black"),
        text=element_text(color="black", family="sans"),
        axis.text = element_text(size=10, color="black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.background = element_blank(), axis.line = element_line(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), legend.position = "none")
print(Figure_2D_empty)
dev.off()

## Analyse multiple control experiments for comparison 
sample_analysed <- "lysate"
# Define the general path
experiment_list <- c("02", "04", "03")

# Define the function to extract the data out of the three experiments
data_shape <- function(Path, Condition_wanted, experiment) {
  output <- list()
  results_1 <- read.csv(paste(Path,  experiment, "_RFU_data.csv", sep=""))[,-1]
  annotation_data_1 <- read.csv(paste(Path,  experiment, "_Annotation.csv", sep="")) %>%
    mutate(Position = paste(Position, experiment, sep="_"))
  
  y_all_1 <- results_1 %>%
    reshape2::melt(., id.var="Temperature") %>%
    `colnames<-`(c("Temperature", "Position", "RFU")) %>%
    mutate(Position = paste(Position, experiment, sep="_")) %>%
    plyr::join(., annotation_data_1, by="Position") %>%
    filter(Condition %in% Condition_wanted)
  
  output[["Data"]] <- y_all_1
  output[["Annotation"]] <- annotation_data_1
  
  return(output)
}

# Extract the data from the tree experiments
all_shift <- list()
for(experiment in experiment_list){
  
  all_shift[[experiment]] <- data_shape(Path, "Control", experiment)
  
}

# Combine all tree experiments in a single df
Data_all <- rbind(all_shift[[experiment_list[1]]]$Data, all_shift[[experiment_list[2]]]$Data) %>%
  rbind(.,all_shift[[experiment_list[3]]]$Data ) %>%
  filter(Condition == "Control") %>%
  group_by(Position) %>%
  mutate(Temperature = Temperature + 273.15) %>% # transform temperature into K
  filter(Temperature < (80 + 273.15) & Temperature > (37 + 273.15), # Filter for only the interesting range
         Sample == sample_analysed) %>% # only select the lysate
  mutate(y_scaled = ((RFU-min(RFU))/(max(RFU)-min(RFU)))) # scale the data

# Add the annotation
Annotation_all <- rbind(all_shift[[experiment_list[1]]]$Annotation, all_shift[[experiment_list[2]]]$Annotation) %>%
  rbind(.,all_shift[[experiment_list[3]]]$Annotation ) %>%
  filter(Condition == "Control")

# Get all positions for fitting
all_positions <- Data_all$Position %>% unique()

# Fit all the data as before
interval_out <- NULL
all_fits <- list()
for( i in all_positions){
  y <- Data_all$y_scaled[Data_all$Position %in% i]
  x <- Data_all$Temperature[Data_all$Position %in% i]
  
  model <-nlsLM(y~trcFunc(x, a, b, c, d, H, Tm),start=list(a=max(y),b=-1,c=min(y), d=1, H =100, Tm = 326))
  
  all_fits[[i]] <- model
  
new.data <- data.frame(x=seq((37 + 273.15), (80 + 273.15), by = 0.1))
interval <- as_tibble(predFit(model, newdata = new.data, interval = "confidence", level= 0.99)) %>% 
 mutate(x = new.data$x, 
           Position = i)
  
  interval_out <- rbind(interval_out, interval)
  
  
}

interval_out %<>%
  plyr::join(., Annotation_all, by="Position")
# Extract the melting temperatures
t <- lapply(all_fits, function(x) coef(x))

Melting_points_control <- data.frame(Position = names(t),
                             Tm = unlist(lapply(t, function(x) x[6])),
                             H = unlist(lapply(t, function(x) x[5]))) %>%
  mutate(Tm_C = Tm - 273.15)

model_fits <- lapply(all_fits, function(x) predict(x, newdata = new.data, interval = "confidence"))

# Extract the fitted data
repeated_new_data <- rep(new.data, length(model_fits)) %>% unlist()
repeated_names <- rep(names(model_fits), each = length(new.data$x))

# Compose a df
fitted_values <- data.frame(Position = repeated_names, 
                            x =repeated_new_data, 
                            fit = unlist(model_fits)) %>%
  plyr::join(., Annotation_all, by="Position")

# Group the positions into experiments
# This is again hard-coded for this experiment. This can be automated by designing another annotation file
conditions <- list("Experiment 1" = c("D3_03", "D4_03"), 
                   "Experiment 2" = c("D5_03", "D6_03"),
                   "Experiment 3" = c("D4_04", "D5_04", "D6_04"), 
                   "Experiment 4" = c("D7_04", "D8_04", "D9_04"),
                   "Experiment 5" = c("F1_02", "F2_02", "F3_02", "F4_02"))

# Awkward way of renaming the positions inot the conditions
for(i in  1:length(conditions)){
  Data_all$Condition[Data_all$Position %in% conditions[[i]]] <- names(conditions[i])
  fitted_values$Condition[fitted_values$Position %in% conditions[[i]]] <- names(conditions[i])
}
# Define colors for experiments
color_codes <- c("#7E9E64", "#649E91", "#64919E", "#647B9E","#64699E" )

# Plot the figure
Figure_S2D <- ggplot(data = Data_all) +  
  geom_point(aes(x=(Temperature -273.15) , y=y_scaled, col=Condition), size=0.5) +
  xlab("Temperature") + 
  ylab("Normalised fluorescence") +
  geom_line(data=fitted_values, aes(x = (x-273.15), y = fit,col=Condition, group = Position )) +
  scale_color_manual(values = color_codes) +
  theme_classic() +
  geom_vline(xintercept = mean(Melting_points_control$Tm_C), lty="dashed", col="gray65") +
  theme(axis.title = element_text(size = 15, color="black"),
        text=element_text(color="black", family="sans"),
        axis.text = element_text(size=15, color="black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.background = element_blank(), axis.line = element_line(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), legend.text = element_text(size=10), 
        legend.title = element_text(size=12)) +
  guides(colour = guide_legend(override.aes = list(size=2)))

# Calculate the mean and sd of the experiments  
mean(Melting_points_control$Tm_C)
sd(Melting_points_control$Tm_C)

# Calculate the mean and sd of the median Savitski temperature  
Melting_Savitski <- c(median(Savitski_data$meltPoint_Lysate_Mg_1[Savitski_data$location == "Cytosol"], na.rm=TRUE), median(Savitski_data$meltPoint_Lysate_Mg_2[Savitski_data$location == "Cytosol"], na.rm=TRUE),
                      median(Savitski_data$meltPoint_Lysate_Mg_3[Savitski_data$location == "Cytosol"], na.rm=TRUE))
mean(Melting_Savitski)
sd(Melting_Savitski)
