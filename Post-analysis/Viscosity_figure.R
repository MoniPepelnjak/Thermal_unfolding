# Viscosity figure
# Author: Monika Pepelnjak

# The script loads in the data from the literature:
  # Osmolyte viscosity measured by us
  # fPSA and ???Gtr from the literature
# The script then fits the linear model between the mean LiP stabilisation score and the properties of osmolytes

path <- "/Volumes/My Passport for Mac/Final_data/"

# Load in the viscosity data
viscosity <- read.csv(paste(path, "/Validations/Viscosity/Viscosity.csv", sep = ""))

# Load in the LiP-MS significance data
significant_all <- readRDS(paste(path, "LiP/Significant_both_quantiles.rds", sep=""))

# Load in the libraries
library(dplyr)
library(magrittr)
library(ggplot2)
library("MASS")

#### Define the color palette
my_palette = c(Betaine = "#B5C3AD", 
               Trehalose = "#d4a993", 
               TMAO = "#87A9C6", 
               Glycerol = "#a997a2", 
               Proline = "#edd9bb", 
               Glucose = "#C78C8C")

# Calculate the mean values from LiP-MS data and make a df
mean_values <- lapply(significant_all, FUN = function(x) mean(x$Protein_level$Protein_stabilisation )) 
df <- data.frame(Condition = names(mean_values), 
                 Mean_stability = unlist(mean_values))

# Join the data with the literature data
df %<>% plyr::join(., viscosity, by="Condition") 

# Make one long format so it is easier to plot in ggplot
df_long <- df %>%
  reshape2::melt(., id.var = c("Condition", "Mean_stability"))

# Sugars have a very weird behaviour based on the fPSA, so we define which osmolytes
df_long$is_sugar <- ifelse(df$Condition%in% c("Trehalose"),  "Trehalose",
                           ifelse(df$Condition%in% c("Glucose"), "Glucose", "Not_Sugar"))
df$is_sugar <- ifelse(df$Condition%in% c("Trehalose"),  "Trehalose",
                      ifelse(df$Condition%in% c("Glucose"), "Glucose", "Not_Sugar"))

# fPSA
cairo_pdf(paste(path, "Figures/Figure2I.pdf", sep=""), width=3, height=3)
ggplot(df_long[df_long$variable == "fPSA",], aes(y=value, x=Mean_stability)) +
  geom_smooth(method="lm", lty="dashed", aes(group_by = is_sugar), col = "gray60", lty="dashed", alpha=0.2) + 
  geom_point(aes(fill=Condition), colour="black",pch=21, size=4) + 
  theme_classic() +
  scale_fill_manual(values = my_palette[sm_list]) +
  theme(axis.title = element_text(size = 15, color="black"),
        text=element_text(color="black"),
        axis.text = element_text(size=12, color="black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.background = element_blank(), axis.line = element_line(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), legend.position = "none") +
  theme(legend.position = "none") +
  #geom_label(aes( label=paste("R2 = ", df_R2_MS, sep=""), x=6, y=0.05 )) +
  ylab("fPSA") +
  xlab("Mean stabilisation score")
dev.off()

df_2 <- df[!df$Condition %in% c("Glucose", "Trehalose"), ]
summary(lm(Mean_stability~fPSA, data=df_2 ))

df_R2_Viscosity <- round(summary(lm(df_long$value[df_long$variable == "Viscosity"] ~ df_long$Mean_stability[df_long$variable == "Viscosity"]))$adj.r.squared,2)

# Viscosity plot
cairo_pdf(paste(path, "Figures/Figure2G.pdf", sep=""), width=3, height=3)
ggplot(df_long[df_long$variable == "Viscosity",], aes(y=value, x=Mean_stability)) +
  geom_smooth(method="rlm", col = "gray60", lty="dashed", alpha=0.2) + 
  geom_point(aes(fill=Condition), colour="black",pch=21, size=4) + 
  theme_classic() +
  scale_fill_manual(values = my_palette[sm_list]) +
  theme(axis.title = element_text(size = 15, color="black"),
        text=element_text(color="black"),
        axis.text = element_text(size=12, color="black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.background = element_blank(), axis.line = element_line(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), legend.position = "none") +
  theme(legend.position = "none") +
  #geom_label(aes( label=paste("R2 = ", df_R2_Viscosity, sep=""), x=6, y=0.8 )) +
  ylab("Viscosity [cP]") +
  theme(legend.position = "none") +
  xlab("Mean stabilisation score")
dev.off()

df_R2_MS <- round(summary(rlm(df_long$value[df_long$variable == "dG"] ~ df_long$Mean_stability[df_long$variable == "dG"]))$adj.r.squared,2)
summary(lm(df_long$value[df_long$variable == "TPSA" & df_long$Condition %in% c("TMAO", "Glycerol", "Proline", "Betaine")  ] ~ df_long$Mean_stability[df_long$variable == "TPSA" & df_long$Condition %in% c("TMAO", "Glycerol", "Proline", "Betaine")]))

# Mass concentration figure
cairo_pdf(paste(path, "Figures/Figure2F.pdf", sep=""), width=3, height=3)
ggplot(df_long[df_long$variable == "Mass_conc",], aes(y=value, x=Mean_stability)) +
  geom_smooth(method="rlm", col = "gray60", lty="dashed", alpha=0.2) + 
  geom_point(aes(fill=Condition), colour="black",pch=21, size=4) + 
  scale_fill_manual(values = my_palette[1:6]) +
  theme_classic() +
  scale_fill_manual(values = my_palette[sm_list]) +
  theme(axis.title = element_text(size = 15, color="black"),
        text=element_text(color="black"),
        axis.text = element_text(size=12, color="black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.background = element_blank(), axis.line = element_line(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), legend.position = "none") +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none") +
  #geom_label(aes( label=paste("R2 = ", df_R2_MS, sep=""), x=6, y=0.05 )) +
  ylab("Mass Concentration (w/v)") +
  xlab("Mean stabilisation score")
dev.off()

# dG figure
cairo_pdf(paste(path, "Figures/Figure2H.pdf", sep=""), width=3, height=3)
ggplot(df_long[df_long$variable == "dG",], aes(y=value, x=Mean_stability)) +
  geom_smooth(method="rlm", col = "gray60", lty="dashed", alpha=0.2) + 
  geom_point(aes(fill=Condition), colour="black",pch=21, size=4) + 
  theme_classic() +
  scale_fill_manual(values = my_palette[sm_list]) +
  theme(axis.title = element_text(size = 15, color="black"),
        text=element_text(color="black"),
        axis.text = element_text(size=12, color="black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.background = element_blank(), axis.line = element_line(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), legend.position = "none") +
  theme(legend.position = "none") +
  #geom_label(aes( label=paste("R2 = ", df_R2_MS, sep=""), x=6, y=0.05 )) +
  ylab(bquote( ''*Delta~g[tr]*' [cal/mol]'))+
  xlab("Mean stabilisation score")
  #xlab("Mean stabilisation score")
dev.off()

df_R2_MS <- round(summary(lm(df_long$value[df_long$variable == "dG"] ~ df_long$Mean_stability[df_long$variable == "dG"]))$adj.r.squared,2)

# Extract the legend. Returns a gtable
p <- ggplot(df, aes(y=Viscosity, x=Mean_stability)) +
  geom_smooth(method="lm", col = "gray60", lty="dashed", alpha=0.2) + 
  geom_point(aes(fill=Condition), colour="black",pch=21, size=4) + 
  theme_classic() +
  theme(legend.position = "bottom",
        text=element_text(size=10)) +
  scale_fill_manual(values = my_palette[sm_list]) +
  ylab("Viscosity") +
  xlab("Mean stabilisation score") +
  guides(fill = guide_legend(nrow = 1))

p <- p + 
  theme(text=element_text(size=15))
leg <- get_legend(p, position = "bottom") 

# Convert to a ggplot and print
cairo_pdf(paste(path, "Figures/Color_palette.pdf", sep=""), width=8, height=0.25)
as_ggplot(leg)
dev.off()



