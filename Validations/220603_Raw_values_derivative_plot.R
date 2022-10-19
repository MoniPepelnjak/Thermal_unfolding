### Plot unfolding curves and derivatives ###
### Goal of the script: Plot the curve and derivative in instances where we have two melting points
### Author: Monika Pepelnjak

# !! If plyr is loaded, the code does not work well - so unload the plyr first!! #
if("package:plyr" %in% search()) detach("package:plyr", unload=TRUE, force = TRUE) 

# Load in the packages
library(ggplot2)
library(dplyr)
library(magrittr)

# Define the sample and experiment: 
sample <- c("DnaK")
experiment <- c("06")

# Define the color palette:
my_palette = c(Betaine = "#B5C3AD", 
               Trehalose = "#d4a993", 
               TMAO = "#87A9C6", 
               Glycerol = "#a997a2", 
               Proline = "#edd9bb", 
               Glucose = "#C78C8C",
               Betaine_2 = "#8E9C86",
               TMAO_2 = "#718CA2",
               Glycerol_2 = "#86697B",
               Proline_2 = "#D8BF9A")

#Load in the annotation 
i <- 1 # the script can be automated, but here we plot only one experiment and protein
annotation <- read.csv(paste("/Users/moni/Documents/Phd/Collaborations/C6/Thermal_shift/", experiment[i], "_Annotation.csv", sep="")) %>%
  filter(Sample %in% sample[i])

#Load in the raw values and select the relevant data
raw_values <- read.csv("/Users/moni/Documents/Phd/Collaborations/C6/Thermal_shift/06_RFU_data.csv") %>%
  .[,-1] %>%
  reshape2::melt(., id.var ="Temperature") %>%
  `colnames<-`(c("Temperature", "Position", "Y_Fitted")) %>%
  plyr::join(.,annotation, by="Position" ) %>% 
  filter(Sample %in% sample[i]) %>%
  filter(Condition != "Negative")

# Scale the data and calculate the derivative
raw_values %<>%
  group_by(Position) %>%
  mutate(scaled = (Y_Fitted-min(Y_Fitted))/(max(Y_Fitted)-min(Y_Fitted)) ) %>%
  mutate(diff_calc = c(diff(Y_Fitted)/diff(Temperature), NA)) %>%
  mutate(diff_calc_scaled = c(diff(scaled)/diff(Temperature), NA))

# Compose a df for plotting in a way that the control and osmolyte can be plotted in the single plot
final_plotting <- NULL
sm_selected <- c("TMAO", "Proline")
for(sm in sm_selected){
  
  raw_values_sub <- raw_values %>%
    filter(Condition %in% c(sm, "Control") ) %>% # always select the osmolyte and control condition
    mutate(Osmolyte = sm) %>% 
    filter(Temperature >= 37 & Temperature <= 76) # only select our relevant temperatures
  
  final_plotting <- rbind(final_plotting, raw_values_sub)
  
}

# Make a nicer df
final_plotting_long <- final_plotting %>%
  ungroup() %>%
  dplyr::select(Temperature, diff_calc_scaled, scaled, Osmolyte, Condition) %>%
  mutate(diff_calc_scaled = diff_calc_scaled) %>%
  reshape2::melt(id.var = c("Temperature", "Osmolyte", "Condition"))

# Rename the column values
final_plotting_long$Analysis <- ifelse(final_plotting_long$variable == "scaled", "Relative flourescene", "Derivative")
final_plotting_long$Analysis <- factor(final_plotting_long$Analysis, levels = c("Relative flourescene", "Derivative"))

# Compose a labels df - this is used for plot names
Labels_examples <- data.frame(Analysis = c("Relative flourescene", "Relative flourescene","Derivative", "Derivative" ), 
                     Osmolyte = c("TMAO", "Proline", "TMAO", "Proline" ), 
                     Label = c("TMAO", "Proline", "", ""), 
                     Condition = c("TMAO", "Proline", "TMAO", "Proline")) 

# Factorise it to have relative flourescence on top
Labels_examples$Analysis <- factor(Labels_examples$Analysis, levels = c("Relative flourescene", "Derivative"))

# Plot and save the plot
g <- ggplot(data=final_plotting_long, aes(x=Temperature, y=value, col=Condition)) +
  stat_summary(aes(col = Condition), fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", width=0.2) +
  stat_summary(aes(fill = Condition), fun=mean, geom="point", size=1, pch=21, col="black") +
  stat_summary(aes(col = Condition), fun=mean, geom="line", size=1) +
  geom_text(data =Labels_examples,  aes(label = Label), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size=5.5, col="black") +
  scale_color_manual(values=c("Control" = "gray60", my_palette[sm_selected])) +
  scale_fill_manual(values=c("Control" = "gray60", my_palette[sm_selected])) +
  facet_grid(Analysis~Osmolyte, scales = "free", space = "free", switch = "y") +
  theme(axis.title = element_text(size = 15, color="black"),
        text=element_text(color="black"),
        axis.text = element_text(size=15, color="black"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=15, color="black"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.y = element_text(angle = 270, size = 15),
        strip.text.x = element_blank(),
        panel.background = element_blank(), axis.line = element_line(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), legend.position = "none") +
  ylab("") +
  xlab("Temperature")

library(grid)
gt = ggplot_gtable(ggplot_build(g))
gt$heights[11] = 2.5*gt$heights[11] # This is hard-coded for this example !! if you have more plots, this might change
pdf("/Users/moni/Documents/Phd/Osmolyte_paper/Figure_6E.pdf", height = 4.4, width=5)
grid.draw(gt)
dev.off()

