## PCA analysis
## Quality check script ##
## Analyse the data for figure S1E - the effect of scaling on the PCA distribution
## Goal of the script: Show that with scaling we reduce the differences caused by changes in PK activity
## Author: Monika Pepelnjak

# Set working environment
setwd("/Users/moni/Documents/Phd/Experiments/015/")

#Import the files of DIA expeirment
DIA_report <- read.csv("DIA_report_HT_full.csv")

#Import the annotation table 
Annotation_table <- read.csv("Annotation_table.csv") %>% 
  filter(.$variable2!="trypsin") # remove the trypsin control samples

#Required packages
library(ggplot2)
library(magrittr)
library(plyr)
library(dplyr)
library(reshape2)
library(ggfortify)
library(UpSetR)
library(ComplexHeatmap)
library(gridExtra)

##Design the color scheme
temperature_colors <- c(rep(colorRampPalette(c("#FFF28B", "#D53616"))(10),6))
Annotation_table$colors <- temperature_colors

## Reshaping the df for PCA analysis
DIA_abundance <- DIA_report %>% merge(., Annotation_table, by="R.FileName") %>% 
  mutate(run = paste(.$variable2, .$variable1, .$replicate, sep="_")) %>% 
  filter(PEP.Quantity > 1000) %>% # remove the low peptide intensities - they are usually false positives
  dplyr::group_by(EG.StrippedSequence, variable2, replicate) %>% 
  dplyr::mutate(Scaled = (PEP.Quantity-min(PEP.Quantity))/(max(PEP.Quantity)- min(PEP.Quantity))) %>% # Scale the data
  ungroup() %>%
  reshape2::dcast(., EG.StrippedSequence ~ run, value.var="Scaled", fun.aggregate = function(x) mean(x, na.rm=T)) %>% #reshape the data
  `rownames<-`(.$EG.StrippedSequence) %>% .[,-1]

# Reformatting for PCA plotting
Annotation_table$run <- paste(Annotation_table$variable2, 
                              Annotation_table$variable1, 
                              Annotation_table$replicate, sep="_")

var2 <- levels(factor(Annotation_table$variable2) )
Annotation_table_var <-  Annotation_table 

DIA_for_pca <- DIA_abundance %>% 
  na.omit() %>% t() %>% as.matrix()

# Calculate and plot the scaled PCA
for_pca <- prcomp(DIA_for_pca, 
                  scale=TRUE)

Annotation_table_var$variable2 <- factor(Annotation_table_var$variable2)

DIA_for_pca_annotated <- DIA_for_pca %>% as.data.frame %>% mutate(run = rownames(.)) %>%
  merge(., Annotation_table_var, by= "run")

b <- autoplot(for_pca, data = DIA_for_pca_annotated, fill="variable1", shape = "variable2",
              size=5) 

b + scale_fill_gradient2(low="#6FD4D0", mid= "#FFF28B", high ="#D53616", midpoint = 54) + 
  theme_bw() +
  scale_shape_manual(values=c(21,22, 24)) +
  theme( 
    panel.grid = element_blank(), 
    text = element_text(size=20), 
    legend.position = "none") 

##Abundance correlation plot non-scaled
DIA_abundance<- DIA_report %>% merge(., Annotation_table, by="R.FileName") %>% 
  mutate(run = paste(.$variable2, .$variable1, .$replicate, sep="_")) %>% 
  reshape2::dcast(., EG.StrippedSequence ~ run, value.var="PEP.Quantity", fun.aggregate = function(x) mean(x, na.rm=T)) %>%
  `rownames<-`(.$EG.StrippedSequence) %>% .[,-1] %>% log2() 

#PCA
Annotation_table$run <- paste(Annotation_table$variable2, 
                              Annotation_table$variable1, 
                              Annotation_table$replicate, sep="_")

var2 <- levels(factor(Annotation_table$variable2) )
Annotation_table_var <-  Annotation_table #%>% filter(., .$variable2==var2[3])

DIA_for_pca <- DIA_abundance %>%  #.[,grepl(var2[3], colnames(.))] %>%
  na.omit() %>% t() %>% as.matrix()

for_pca <- prcomp(DIA_for_pca, 
                  scale=TRUE)

Annotation_table_var$variable2 <- factor(Annotation_table_var$variable2)

DIA_for_pca_annotated <- DIA_for_pca %>% as.data.frame %>% mutate(run = rownames(.)) %>%
  merge(., Annotation_table_var, by= "run")

a <- autoplot(for_pca, data = DIA_for_pca_annotated, fill="variable1", shape = "variable2",
              size=5) 

Non_Scaled_PCA <- a + scale_color_continuous(low="#FFF28B", high ="#D53616") + 
  scale_fill_gradient2(low="#6FD4D0", mid= "#FFF28B", high ="#D53616", midpoint = 54) + 
  theme_bw() +
  scale_shape_manual(values=c(21,22, 24)) +
  theme( 
    panel.grid = element_blank(), 
    text = element_text(size=20), 
    legend.position = "none") 

Scaled_PCA <- b + scale_color_continuous(low="#FFF28B", high ="#D53616") + 
  scale_fill_gradient2(low="#6FD4D0", mid= "#FFF28B", high ="#D53616", midpoint = 54) + 
  theme_bw() +
  scale_shape_manual(values=c(21,22, 24)) +
  theme( 
    panel.grid = element_blank(), 
    text = element_text(size=20), 
    legend.position = "right") +
  guides(fill=guide_legend(title="Temp"), 
         shape = guide_legend(title = "PK conc."))

pdf("/Users/moni/Documents/Phd/Thesis/Chapter1/PCA_before_after_scaling.pdf", height = 4, width=10.5)
grid.arrange(Non_Scaled_PCA, Scaled_PCA, 
              nrow = 1, 
             widths = c(1, 1.3))
dev.off()
