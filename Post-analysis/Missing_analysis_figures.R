### Osmolyte paper figure analysis ### 
## Author: Monika Pepelnjak
# The script performs the remaining analysis for the Osmolyte paper

# Required files:
  # Fitted data
  # Significantly changing protein list
  # TPP data (fitted)
  # TPP data significantly changing proteins

# Load in the files
all_results <- readRDS("/Volumes/My Passport for Mac/data/All_data_1M.rds")
significant_TMAO <- readRDS("/Volumes/My Passport for Mac/data/Significant_both_quantiles.rds")
TPP <- readRDS("/Users/moni/polybox/data/Python_TPP_fits_all.rds")
TPP_significant <-  readRDS("/Users/moni/polybox/data/Significant_TPP.rds")

## Libraries
library(ggplot2)
library(magrittr)
library(dplyr)
library(ggpubr)

## Calculate overlap between the different osmolytes ###
stabilised_list <- lapply(significant_TMAO, FUN = function(x) x$Protein_level$Protein[x$Protein_level$Protein_stabilisation > 0])
all_list <- lapply(significant_TMAO, FUN = function(x) x$Protein_level$Protein)
sm_list <- c("TMAO", "Betaine", "Glycerol", "Proline", "Trehalose", "Glucose")
values_overlap <- c()
list_overlap <- list()

# Calculate pairwise overlap as in figure 3A
for(sm1 in sm_list){
  
  for(i in 1:length(sm_list)){
    sm2 <- sm_list[i]
    
    all <- intersect(all_list[[sm1]], all_list[[sm2]])
    stabilised_A <- intersect(stabilised_list[[sm1]], all)
    stabilised_B <- intersect(stabilised_list[[sm2]], all)
    
    values_overlap[i] <- length(intersect(stabilised_A, stabilised_B))/length(union(stabilised_A, stabilised_B))

  }

  list_overlap[[sm1]] <- values_overlap
  
}

for(sm1 in sm_list){
  
  for(i in 1:length(sm_list)){
    sm2 <- sm_list[i]
    
    values_overlap[i] <- length(intersect(stabilised_list[[sm1]], stabilised_list[[sm2]]))/ length(intersect(stabilised_list[[sm1]], all_list[[sm2]]))
    
  }
  
  list_overlap[[sm1]] <- values_overlap
  
}


matrix_values <- list_overlap %>%
  unlist() %>%
  matrix(ncol=6) %>% # make the data into matrix
  `colnames<-`(names(list_overlap)) %>% 
  `rownames<-`(names(list_overlap)) %>%
  round(., 2) %>% # round the numbers
  t() 

#both <- intersect(c(stabilised_list[["TMAO"]], stabilised_list[["Glycerol"]]), c(all_list[["TMAO"]])) %>%
#  intersect(., all_list[["Glycerol"]])

# Plot figure 3A from the paper
cairo_pdf("/Users/moni/Documents/Phd/Osmolyte_paper/Figure_3A_red.pdf", height=4, width=5)

ComplexHeatmap::Heatmap(matrix_values, cluster_columns = FALSE, cluster_rows = FALSE, row_names_side = "left", 
                        rect_gp = grid::gpar(col = "white", lwd = 2), show_heatmap_legend = FALSE,
                        col = circlize::colorRamp2(c(0.24, 1), c("white", "#924A4A")),
                        column_title = "Compound B", 
                        column_title_side = "bottom",
                        row_names_gp = grid::gpar(fontsize=14),
                        column_names_gp = grid::gpar(fontsize=14),
                        row_title = "Compound A",
                        na_col = "white",
                        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                          color = ifelse(is.na(matrix_values[i, j]), "white", "black")
                          grid::grid.text(matrix_values[i, j], x, y, 
                                          gp = grid::gpar(fontsize=12, col=color))
            })
dev.off()

#### Get number of all stabilised conditions/protein
output <- NULL

# Extract protein-level information for all the osmolytes
# Combine the data into a single data frame
for(sm in sm_list){
  
  sm_df <- significant_TMAO[[sm]]$Protein_level %>%
    dplyr::select(Protein, Protein_stabilisation) %>%
    mutate(Condition = sm)
  
  output <- rbind(output, sm_df)
}

# Design a summary 
output_2 <- output %>%
  group_by(Protein) %>%
  mutate(n_condition = n()) %>% 
  filter(n_condition == 6) %>% # Filter only for proteins that are detected in all 6 conditions
  summarise(n_stabilised = sum(Protein_stabilisation > 0), # get the number of stabilising conditions
            # Get the most stabilising osmolyte per protein
            which_max = ifelse(any(Protein_stabilisation > 0), Condition[Protein_stabilisation == max(Protein_stabilisation)], "Not stabilised"))

# Plot figure 3B
Figure_3B <- ggplot(output_2, aes(x=n_stabilised)) +
  geom_bar(fill="#908B7F", col="black", lwd=0.8) +
  theme_classic() +
  scale_x_continuous(breaks=c(0:6)) +
  theme(axis.title = element_text(size = 15, color="black"),
        text=element_text(color="black"),
        axis.text = element_text(size=15, color="black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.background = element_blank(), axis.line = element_line(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), legend.position = "none") +
  xlab("Number of stabilising conditions") +
  ylab("Number of proteins")
ggsave("/Users/moni/Documents/Phd/Osmolyte_paper/Figure_3B.pdf", Figure_3B, height = 2.5, width = 5.2)

output_2$which_max <- factor(output_2$which_max, levels = c("TMAO", "Betaine", "Glycerol", "Proline", "Trehalose", "Glucose", "Not stabilised"))

## Identify which is the best stabiliser 
Figure_3C <- ggplot(output_2, aes(x=which_max, fill=which_max)) +
  geom_bar(aes(y = (..count..)/sum(..count..)), col="black") +
  scale_fill_manual(values = c(my_palette[sm_list], "Not stabilised" = "azure3")) +
  theme_classic() +
  theme(axis.title = element_text(size = 15, color="black"),
        text=element_text(color="black"),
        axis.text = element_text(size=15, color="black"),
        axis.text.x = element_text(angle=45, hjust=1),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.background = element_blank(), axis.line = element_line(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), legend.position = "none") +
  theme(legend.position = "none") +
  xlab("Best stabiliser") +
  ylab("Percentage of proteins")
ggsave("/Users/moni/Documents/Phd/Osmolyte_paper/Figure_3C.pdf", Figure_3C, height = 3.2, width = 5.2)

### Correlation analysis for figure 3D
# Select for proteins where at least 1 osmolytes is stabilising and it is detected at all conditions
output_3 <- output %>%
  group_by(Protein) %>%
  mutate(n_condition = n(), 
         n_stabilised = sum(Protein_stabilisation != 0)) %>%
  filter(n_condition == 6)

## Mean values for osmolytes accross all the measured proteins
mean_all <- lapply(significant_TMAO, function(x) mean(x$Protein_level$Protein_stabilisation))
df_mean <- data.frame(Condition = names(mean_all), 
                      Mean_LIP = unlist(mean_all))

# Combine the individual protein values with mean values df
output_3 %<>%
  plyr::join(., df_mean, by="Condition")

# Plot for examples shown in figure 4B
# Selection and renaming of the selected proteins
examples <- output_3 %>%
  filter(Protein %in% c("P0A6V8", "P23917", "P0A9J6")) %>%
  mutate(Protein = ifelse(Protein == "P0A6V8", "Glk", 
                          ifelse(Protein == "P0A9J6", "RbsK", "Mak")))

examples$Protein <- factor(examples$Protein, levels=c("Glk", "Mak", "RbsK"))
Labels_examples <- examples %>%
  dplyr::select(Protein) %>%
  unique()

Figure_4B <- ggplot(examples, aes(x=Mean_LIP, y=Protein_stabilisation) ) +
  geom_smooth(method="lm", col = "gray60", lty="dashed", alpha=0.2) + 
  geom_point(aes(fill=Condition), colour="black",pch=21, size=4) + 
  scale_fill_manual(values = my_palette[sm_list]) +
  facet_wrap(~Protein, scales="free") +
  theme_classic() +
  geom_text(data =Labels_examples,  aes(label = Protein), x = -Inf, y = Inf, hjust = -0.2, vjust = 1.5, size=5.5) +
  geom_text(data =Labels_examples,  aes(label = paste("Cor = ", cor, sep="")), x = Inf, y = -Inf, vjust = -0.5, hjust = 1, size=5.5) +
  theme(axis.title = element_text(size = 15, color="black"),
        text=element_text(color="black"),
        axis.text = element_text(size=15, color="black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.background = element_blank(), axis.line = element_line(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), legend.position = "none") +
  ylab("Protein stabilisation score") +
  xlab("Mean stabilisation score")
ggsave("/Users/moni/Documents/Phd/Osmolyte_paper/Figure_4B.pdf", height=3, widt=8)

# Set seed for randomised number selection
set.seed(126)
random_numbers <- sample(1:50000, 1000, replace = FALSE)
random_out <- NULL
# randomisation_loop
for(i in 1:length(random_numbers)){
  # set seed within the loop by choosing the numbers 
  set.seed(random_numbers[i])

  random_output <- output_3 %>%
    filter(n_stabilised >= 4) %>% # only focus on proteins with at least 4 stabilising conditions
    group_by(Protein,) %>%
    summarise(cor_random = cor.test(Mean_LIP[sample(6, replace=FALSE)], Protein_stabilisation, method = "spearman")$estimate )) %>% #calculate the correlation between Protein score and randomised mean value
    mutate(seed = random_numbers[i])
  random_out <- rbind(random_out, random_output)
  print(i)
}


output_3 %<>%
  group_by(Protein, n_stabilised) %>%
  summarise(pvalue = cor.test(Mean_LIP, Protein_stabilisation, method = "spearman")$p.value,
            cor = cor.test(Mean_LIP, Protein_stabilisation, method = "spearman")$estimate 
            )

#seeds <- unique(random_out$seed)
#ggplot(random_out[random_out$seed==seeds[7],], aes(x=cor_random)) +
#  geom_density()

# ggplot(output_3, aes(x=cor)) +
#   geom_histogram(bins=15, fill="gray70", col="gray20") +
#   theme_classic() +
#   #facet_wrap(~n_stabilised)+
#   scale_x_continuous(breaks=c(0:6)) +
#   xlab("Number of stabilising conditions") +
#   ylab("Number of proteins")

# ggplot(output_3, aes(x=cor_random)) +
#   geom_histogram(bins=15, fill="gray70", col="gray20") +
#   theme_classic() +
#   #scale_x_continuous(breaks=c(0:6)) +
#   xlab("Number of stabilising conditions") +
#   ylab("Number of proteins")
# 
# output_4 <- output_3 %>%
#   reshape2::melt()

#density_values <- ggplot(random_out, aes(x=cor_random,group=seed )) +
#  geom_density() +
#  theme_classic() +
#  #scale_x_continuous(breaks=c(0:6)) +
#  xlab("Number of stabilising conditions") +
#  ylab("Number of proteins") +
#  theme(legend.position = "none")

# Get the densities of the individual random correlations
frame_out <- NULL
for(i in 1:length(random_numbers)){
  test <- density( na.omit(random_out$cor_random[random_out$seed == random_numbers[i]]), from = -1, to=1, n=500)
  frame <- data.frame(seed = random_numbers[i], 
             x = test$x, 
             y= test$y)
  frame_out <- rbind(frame, frame_out)
  print(i)
}

# Calculate the density of the true correlations
density_true <-  density( na.omit(output_3$cor[output_3$n_stabilised >=4], from = -1, to=1, n=500))
frame_true <- data.frame(x = density_true$x, 
                    y= density_true$y)

ggplot(output_3, aes(x=cor)) +
  geom_density()
three_examples <- data.frame( x = output_3$cor[output_3$Protein %in% c("P00350", "P0A6P9", "P0A9J6")], 
                              y = 0)
Figure_3D <- ggplot(frame_out, aes(x=x, y=y)) +
  geom_line(aes(group = seed), col="#ADA293", alpha=0.01) +
  geom_smooth(col="#7E7260", lty="dashed") +
  geom_line(data = frame_true, aes(x=x, y=y), col="#852D2D", lwd=1.5)+
  geom_ribbon(data=frame_true ,aes(x=x,ymax=y),ymin=0,alpha=0.5, fill="#852D2D") +
  theme_classic() +
  theme(axis.title = element_text(size = 15, color="black"),
        text=element_text(color="black"),
        axis.text = element_text(size=15, color="black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.background = element_blank(), axis.line = element_line(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), legend.position = "none") +
  #geom_point(data = three_examples, aes(x=x, y=y)) +
  xlab("Correlation") +
  ylab("Density") +
  xlim(c(-1,1))

ggsave("/Users/moni/Documents/Phd/Osmolyte_paper/Figure_3D.pdf", Figure_3D, height=2.5, width=5.2)

ggplot(output_3, aes(x=cor)) +
  geom_density(lwd=1.5, fill="#A69797") +
  theme_classic() +
  xlab("Spearman correlation") +
  ylab("Density")


### Stabilised destabilised 
Stabilised <- lapply(significant_TMAO, function(x) sum(x$Protein_level$Protein_stabilisation > 0)/length(x$Protein_level$Protein_stabilisation) )
Destabilised <- lapply(significant_TMAO, function(x) sum(x$Protein_level$Protein_stabilisation < 0)/length(x$Protein_level$Protein_stabilisation) )
Stabilisation_df <- data.frame(Condition = names(Stabilised), 
                               Stabilised = unlist(Stabilised), 
                               Destabilised = unlist(Destabilised)) %>%
  reshape2::melt()

Stabilised_peptide <- lapply(significant_TMAO, function(x) sum(x$scores$Stabilisation > 0)/length(x$scores$Stabilisation) )
Destabilised_peptide <- lapply(significant_TMAO, function(x) sum(x$scores$Stabilisation < 0)/length(x$scores$Stabilisation) )
Stabilisation_df_peptide <- data.frame(Condition = names(Destabilised_peptide), 
                               Stabilised = unlist(Stabilised_peptide), 
                               Destabilised = unlist(Destabilised_peptide)) %>%
  reshape2::melt()

sm_list <- c("TMAO", "Betaine", "Glycerol", "Proline", "Trehalose", "Glucose")
Stabilisation_df$Condition <- factor(Stabilisation_df$Condition, levels=sm_list)
Stabilisation_df_peptide$Condition <- factor(Stabilisation_df_peptide$Condition, levels=sm_list)

Figure_2B <- ggplot(Stabilisation_df[Stabilisation_df$variable == "Stabilised",], aes(x=Condition, fill=Condition, y=value)) +
  geom_bar(  stat = "identity", position = "dodge", col="black") +
  scale_fill_manual(values = c(my_palette[sm_list])) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=15, angle=45, vjust=1, hjust=1, color="black"), 
        axis.text.y = element_text(size=15,  color="black"), 
        axis.title = element_text(size=15),
        axis.title.x = element_blank()) +  
  xlab("") +
  ylim(c(0,0.8)) +
  ylab("Fraction of stabilised proteins") +
  xlab("")

ggsave("/Users/moni/Documents/Phd/Osmolyte_paper/Figure2B.pdf", height = 3.5, width = 4)


all_stabilised <- NULL
for(sm in sm_list){
  
  s <- significant_TMAO[[sm]]$Protein_level %>%
    filter(Protein_stabilisation > 0) %>%
    mutate(Condition = sm) 
  
  if(any(colnames(s) %in% "mean_p_score")) {
    s %<>%
      dplyr::select(-mean_p_score)
  }

      all_stabilised <- rbind(all_stabilised, s)
}

all_stabilised$Condition <- factor(all_stabilised$Condition, levels=sm_list)
ggplot(all_stabilised, aes(x=Condition, fill=Condition, y=Protein_stabilisation)) +
  geom_boxplot(col = "gray20", outlier.colour = "gray20") +
  scale_fill_manual(values = c(my_palette[sm_list])) +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("Condition") +
  ylab("Percentage of stabilised proteins") +
  xlab("")

ggplot(Stabilisation_df[Stabilisation_df$variable == "Destabilised",], aes(x=Condition, fill=Condition, y=value)) +
  geom_bar(  stat = "identity", position = "dodge", col="gray20") +
  scale_fill_manual(values = c(my_palette[sm_list])) +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("Condition") +
  ylim(c(0,0.1)) +
  ylab("Percentage of destabilised proteins")

ggplot(Stabilisation_df_peptide[Stabilisation_df_peptide$variable == "Stabilised",], aes(x=Condition, fill=Condition, y=value)) +
  geom_bar(  stat = "identity", position = "dodge", col="gray20") +
  scale_fill_manual(values = c(my_palette[sm_list])) +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("Condition") +
  ylim(c(0,0.8)) +
  ylab("Percentage of stabilised proteins")

##### Binding story completed #### 
Binding_partners <- list()

Binding_partners[["Glucose"]] <- c(read.csv("/Users/moni/Documents/Phd/databases/Glucose_binding_uniprot.csv")$Entry, "P0AEE5")
Binding_partners[["Proline"]] <- c(read.csv("/Users/moni/Documents/Phd/databases/Proline_binding_uniprot.csv")$Entry)
Binding_partners[["Betaine"]] <- c("P17445", "P17444", "P33362")
Binding_partners[["Monosaccharides"]] <- read.csv("/Users/moni/Documents/Phd/databases/Monosaccharide_binding.csv")$Entry

bp_comparison <- c("Glucose", "Proline", "Betaine")

Combined_binding <- NULL
for(bp in bp_comparison){
  
  Comparison_df <- significant_TMAO[[bp]]$scores %>%
    mutate(Binds = Protein %in% Binding_partners[[bp]]) %>%
    filter(Stabilisation > 0) %>%
    mutate(Scaled_stabilisation = (Stabilisation - median(Stabilisation))/(quantile(Stabilisation, 0.75) - quantile(Stabilisation, 0.25)) ) %>%
    mutate(Median_scaled = (Stabilisation - median(Stabilisation))) %>%
    dplyr::select(Protein, Scaled_stabilisation, Binds, Stabilisation, Median_scaled) %>%
    mutate(Condition = bp)
    
  Combined_binding <- rbind(Comparison_df, Combined_binding)
  
}
#max(Combined_binding[Combined_binding$Stabilisation >0,]$Scaled_stabilisation)
Combined_binding$Binds <- ifelse(Combined_binding$Binds  == TRUE, "Binding", "Not binding")
Figure4A_left <- ggplot(Combined_binding[Combined_binding$Stabilisation >0,], aes(x=Binds, y=Scaled_stabilisation, fill=Binds)) +
  geom_violin(alpha=0.5) +
  geom_boxplot(width=0.2, outlier.alpha = 0) +
  #facet_wrap(~Condition) +
  #geom_jitter(width=0.2) +
  scale_fill_manual(values=c("Not binding" = "#CAC6BA", "Binding" = "#494439", "Monosaccharide"="#7A7362"))+
  theme_classic() +
  ylim(c(-1.21, 6)) +
  theme(legend.position = "none",
        axis.text.x = element_text(size=15, angle=45, vjust=1, hjust=1, color="black"), 
        axis.text.y = element_text(size=15,  color="black"), 
        axis.title = element_text(size=15),
        axis.title.x = element_blank()) +
  stat_compare_means(comparisons = list(c("Binding", "Not binding")), label = "p.signif", size=8) +
  ylab("Scaled stabilisation")
ggsave("/Users/moni/Documents/Phd/Osmolyte_paper/Figure_4A_left.pdf", height=4, width = 2.5)


MS_df <- significant_TMAO[["Glucose"]]$scores %>%
  mutate(Binds = ifelse(Protein %in% Binding_partners[["Glucose"]], "Glucose", 
                        ifelse(Protein %in% Binding_partners[["Monosaccharides"]], "Not binding", "Not binding"))) %>%
  mutate(Scaled_stabilisation = (Stabilisation - median(Stabilisation))/(quantile(Stabilisation, 0.75) - quantile(Stabilisation, 0.25)) ) %>%
  dplyr::select(Protein, Scaled_stabilisation, Binds, Stabilisation) %>%
  mutate(Condition = "Glucose")
MS_df$Outlier <- ifelse(MS_df$Protein == "P0A9J6", "TRUE", "FALSE")
point <- MS_df %>%
  filter(Outlier == TRUE) %>%
  filter(Stabilisation == max(Stabilisation))
MS_df$Binds[MS_df$Binds == "Monosaccharide"] <- "Monosacch."

Figure_4A_right <- ggplot(MS_df[MS_df$Stabilisation > 0,], aes(x=Binds, y=Stabilisation, fill=Binds)) +
  geom_violin(alpha=0.5, col=NA) +
  geom_boxplot(width=0.2, outlier.alpha = 0.8) +
  #geom_point(data = point, aes(x=Binds, y=Stabilisation), col="red", size = 3, pch=21, fill=NA, stroke = 2) +
  scale_fill_manual(values=c("Not binding" = "#CAC6BA", "Glucose" = "#7A7362", "Monosacch."="#7A7362"))+
  #geom_jitter() +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=15, angle=45, vjust=1, hjust=1, color="black"), 
        axis.text.y = element_text(size=15,  color="black"), 
        axis.title = element_text(size=15),
        axis.title.x = element_blank()) +
  stat_compare_means(comparisons = list(c("Glucose", "Not binding")), , label = "p.signif", size=8) +
  ylim(c(0,48)) +
  ylab("Stabilisation score")
ggsave("/Users/moni/Documents/Phd/Osmolyte_paper/Figure_4A_right.pdf", Figure_4A_right, height=4, width = 3)

output_binding <- list()
sm_list <- c("TMAO", "Betaine", "Glycerol", "Proline", "Trehalose",  "Glucose")
for(sm in sm_list){
  Paths <- read.csv("/Users/moni/Documents/Phd/Experiments/TD_analysis/Path_directions.csv")
  
  DIA_import <- read.delim(Paths$path[Paths$sm == sm & Paths$file == "Spectronaut_export"], header=T)
  
  Annotation <- read.csv(Paths$path[Paths$sm == sm & Paths$file == "annotation_file"]) %>%
    filter(variable2 %in% c(sm, "Control"), replicate %in% c(1,2)) %>%
    mutate(run = paste(.$variable2, .$variable1, .$replicate, sep="_"))
  
  #numb_rep <- max(Annotation$replicate)
  
  Binding_annotated <- DIA_import %>% 
    plyr::join(., Annotation, by="R.FileName") %>%
    filter(.$variable2 %in% c("Control", sm) &
             .$replicate %in% c(1,2))  %>%
    filter(.$variable1 %in% c(37, 40.5)) %>%
    filter(grepl("ecoli", PG.FastaFiles), 
           PEP.IsProteotypic == "True")
  
  Binding <- Binding_annotated %>%
    mutate(Abundance = as.numeric(as.character(.$PEP.Quantity))) %>%
    reshape2::dcast(., PEP.StrippedSequence ~ run, value.var = "Abundance",
                    fun.aggregate = function(x) sum(x, na.rm = T))
  
  Binding[Binding < 100] <- NA
  
  #Remove values that have more than 1 missing value
  Binding_check <- Binding %>%
    mutate(control_missing = apply(.[,grepl("Control", colnames(.))], MARGIN=1, FUN=function(x) sum(is.na(x))),
           glucose_missing = apply(.[,grepl(sm, colnames(.))], MARGIN=1, FUN=function(x) sum(is.na(x)))) %>%
    filter(.$control_missing < 1 & .$glucose_missing < 1) %>%
    .[,!(names(.) %in% c("control_missing", "glucose_missing"))]
  
  pvalues_binding <- Binding_annotated[Binding_annotated$variable1 %in% c(37, 40.5),] %>%
    filter(PEP.StrippedSequence %in% Binding_check$PEP.StrippedSequence) %>%
    group_by(PEP.StrippedSequence, PG.ProteinAccessions, PEP.DigestType....Trypsin.P., PG.ProteinDescriptions) %>%
    #  mutate(n = n()) %>%
    #  filter(n >= 8 ) %>%
    summarise(pvalue = t.test(PEP.Quantity[variable2 == sm], PEP.Quantity[variable2 == "Control"])$p.value,
              FC = log2(mean(PEP.Quantity[variable2 == sm])/ mean(PEP.Quantity[variable2 == "Control"])))
  
  Binding_annotated$variable1 %>% unique()
  
  pvalues_binding %<>%
    ungroup() %>%
    mutate(qvalue = p.adjust(pvalue, method="BH"))
  
  pvalues_binding %<>%
    mutate(is_significant = ifelse(abs(FC) > 1.5 & qvalue < 0.01, "Significant", "Not_significant"))

  significant <- pvalues_binding %>%
    filter(abs(FC) > 1.5 & qvalue < 0.01 )
  
  output_binding[[sm]][["Significant"]] <- significant
  
  pvalues_binding_count <- pvalues_binding %>%
    group_by(PG.ProteinAccessions, PG.ProteinDescriptions) %>%
    summarise(n=n())
  
  significant_count <- significant %>%
    group_by(PG.ProteinAccessions, PG.ProteinDescriptions) %>%
    summarise(n=n())
  
  output_binding[[sm]][["significant_count"]] <- significant_count
  
  output_binding[[sm]][["Binding_analysis"]] <- pvalues_binding
  
  
}


summarised_value_all <- NULL
median_value <- list()
for(sm in sm_list){
  
  summarised_value <- output_binding[[sm]]$Binding_analysis %>%
    group_by(PG.ProteinAccessions) %>%
    summarise( n_all = n() , 
               n_changing = sum(qvalue < 0.01 & abs(FC) > 1.5)) %>%
    mutate(percentage_changing = n_changing/n_all) %>%
    mutate(Condition = sm)
  
  summarised_value_all <- rbind(summarised_value, summarised_value_all)
  
  median_value[[sm]] <- median(summarised_value$percentage_changing[summarised_value$percentage_changing > 0 ])

  
}

summarised_value_all$Condition <- factor(summarised_value_all$Condition, levels=sm_list)
ggplot(summarised_value_all[summarised_value_all$percentage_changing > 0 ,], aes(x= Condition, y=percentage_changing, fill=Condition)) +
  geom_boxplot() +
  #scale_fill_manual(values=my_palette) +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("") +
  ylab("Percentage of changing peptides/protein")

summarised_value_all %<>%
  group_by(Condition) %>%
  mutate(median_changing = mean(percentage_changing[percentage_changing>0]))

summarised_value_all_changing <- summarised_value_all %>%
  group_by(PG.ProteinAccessions) %>%
  mutate(n_appears = n()) %>% 
  filter( n_appears == 6) %>%
  summarise(n_os_changing= sum( percentage_changing > 0))

ggplot(summarised_value_all_changing, aes(x=n_os_changing)) +
  geom_bar(fill="gray70", col="gray20") +
  theme_classic() +
  scale_x_continuous(breaks=c(0:6)) +
  xlab("Number of changing conditions") +
  ylab("Number of proteins")


#### Binding proteins

Changing_proteins <- lapply(output_binding, function(x) length(unique(x$Significant$PG.ProteinAccessions)))
All_proteins <- lapply(output_binding, function(x) length(unique(x$Binding_analysis$PG.ProteinAccessions)))
#length(unique(output_binding$Betaine$Binding_analysis$PG.ProteinAccessions))
df_protein <- data.frame(Condition = names(Changing_proteins),
                            Changing_percentage = unlist(Changing_proteins), 
                         All_proteins = unlist(All_proteins)) %>%
  mutate(Percentage_changing =Changing_percentage/All_proteins )
df_protein$Condition <- factor(df_protein$Condition, levels = sm_list)

ggplot(df_protein, aes(x=Condition, y=Percentage_changing, fill=Condition)) +
  geom_bar(stat="identity", col="gray20") +
  scale_fill_manual(values = c(my_palette[sm_list])) +
  ylim(c(0,0.5)) +
  theme_classic() +
  theme(legend.position = "none") +
  ylab("Percentage of proteins") +
  xlab("")

cc <- significant_TMAO$Betaine$AA_level %>%
  group_by(protein_position_id) %>%
  summarise(n = n()) %>%
  group_by(n) %>%
  summarise(n_combined=n()) %>%
  ungroup() %>%
  mutate(percentage = n_combined/sum(n_combined))


### Compare significant binding
sm_list <- c( "TMAO", "Betaine", "Glycerol", "Proline", "Trehalose", "Glucose")
percentage_stabilised <- c()
output_test_binding <- list()
all_pvalues <- NULL
ggplot_frame_combined <- NULL
for(sm in sm_list){
group1 <- output_binding[[sm]]$Binding_analysis$PG.ProteinAccessions[output_binding[[sm]]$Binding_analysis$qvalue < 0.01 & abs(output_binding[[sm]]$Binding_analysis$FC) > 2 ]
#group2 <- (output_3$Protein[output_3$cor < 0.8 & output_3$cor > 0.5]) 
group2 <- output_binding[[sm]]$Binding_analysis$PG.ProteinAccessions %>% unique() 
group2 <- group2[!group2 %in% group1 ]
group_stabilised <- significant_TMAO[[sm]]$Protein_level$Protein[significant_TMAO[[sm]]$Protein_level$Protein_stabilisation>0]
group_not_stabilised <- significant_TMAO[[sm]]$Protein_level$Protein[significant_TMAO[[sm]]$Protein_level$Protein_stabilisation<=0]

#group1 <- group1[!group1 %in% ribosomal_proteins$X...Entry]
#group2 <- group2[!group2 %in% ribosomal_proteins$X...Entry]

#output_test <- Test_two_groups(All_combined_predictions, group1, group2, "Binding", "Not_binding")

#output_test_binding[[sm]] <- output_test
percentage_stabilised[sm] <- length(intersect(group_stabilised, group1))/length(intersect(group_stabilised, c(group1, group2)))

#ggplot_frame <- output_test$All_values %>%
#  reshape2::melt(., id.vars = c("Protein", "group")) %>%
  #filter(variable == "") %>%
#  mutate(Condition = sm)

#ggplot_frame_combined <- rbind(ggplot_frame_combined, ggplot_frame)

#pvalues <- output_test$Significance_test %>%
#  mutate(Condition = sm)
#all_pvalues <- rbind(all_pvalues, pvalues)
}

df_binding_stabilised <- data.frame(Condition = names(percentage_stabilised), 
                                    Percentage = unlist( percentage_stabilised) )
df_binding_stabilised$Condition <- factor(df_binding_stabilised$Condition, levels = df_binding_stabilised$Condition)
all_figures[["Figure_5H"]] <- Figure_5H
saveRDS(all_figures, "/Volumes/My Passport for Mac/data/All_figures.RDS")
Figure_4D <- ggplot(df_binding_stabilised, aes(x=Condition, y=Percentage, fill=Condition)) +
  geom_bar(stat="identity", col="black") +
  scale_fill_manual(values = c(my_palette[sm_list])) +
  ylim(c(0,0.5)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=15, angle=45, vjust=1, hjust=1, color="black"), 
        axis.text.y = element_text(size=15,  color="black"), 
        axis.title = element_text(size=15),
        axis.title.x = element_blank()) +
  theme(legend.position = "none") +
  ylab("Fraction of stabilised Proteins") +
  xlab("")
ggsave("/Users/moni/Documents/Phd/Osmolyte_paper/Figure_4D.pdf", Figure_4D, height=3.6, width=6)
#names(significant_TMAO)
sm_list <- c("TMAO", "Trehalose", "Glycerol", "Proline", "Betaine", "Glucose")
MS_df_all <- NULL
for(sm in sm_list){
MS_df <- significant_TMAO[[sm]]$Protein_level %>%
  mutate(LIP_binding = ifelse(Protein %in% output_binding[[sm]]$Significant$PG.ProteinAccessions, "Changing", 
                              ifelse(Protein %in% output_binding[[sm]]$Binding_analysis$PG.ProteinAccessions, "Not Changing", NA))) %>%
  filter(Protein_stabilisation > 0) %>%
  #mutate(coverage = n_peptide >= 5) %>%
  mutate(Scaled_stabilisation = (Protein_stabilisation - median(Protein_stabilisation))) %>%
  dplyr::select(Protein, Scaled_stabilisation, LIP_binding, Protein_stabilisation) %>%
  mutate(Condition = sm)

MS_df_all <- rbind(MS_df_all, MS_df)
}

MS_df_plot <- MS_df_all %>%
  filter(LIP_binding %in% c("Changing", "Not Changing"), 
         Protein_stabilisation > 0) %>%
  filter(Condition %in% c("Betaine", "Trehalose", "TMAO", "Proline"))

ggplot(MS_df_plot, aes(x=LIP_binding, y=Protein_stabilisation)) +
  geom_violin(aes(fill=Condition), alpha=0.5) +
  geom_boxplot(width=0.5, outlier.alpha = 0, aes(fill=Condition)) +
  facet_wrap(~Condition, scales="free", nrow=1) +
  stat_compare_means(label = "p.signif", comparisons =  list(c("Changing", "Not Changing")), vjust = 0.5) +
  scale_fill_manual(values=my_palette[sm_list]) +
  #ylim(c(0,10))+
  #geom_jitter() +
  theme_classic() +
  ylab("Stabilisation Score") +
  xlab("Native structure") +
  theme(legend.position = "none")


#### Heatmap binding
volcano <- NULL

for(sm in sm_list){
  sub <- output_binding[[sm]]$Binding_analysis %>%
    mutate(Condition = sm)
  
  volcano <- rbind(volcano, sub)
}
binding_partners <- unlist(Binding_partners[c("Proline", "Glucose", "Betaine")])

second_option <- volcano %>%
  filter(PG.ProteinAccessions %in%  binding_partners) %>%
  group_by(Condition, PG.ProteinAccessions) %>%
  summarise(n_significant = sum(qvalue < 0.05 & abs(FC) > 1.5)) %>%
  group_by(PG.ProteinAccessions) %>%
  filter(any(n_significant > 0)) %>%
  `colnames<-`(c("Condition", "Protein", "n_significant"))

gene_annotation <- read.csv("/Users/moni/Documents/Phd/databases/Gene_name_ids.csv")

second_option %<>%
  plyr::join(., gene_annotation, by="Protein")

Binding_Compounds <- Binding_partners[c("Proline", "Glucose", "Betaine")]
n_int <- lapply(Binding_Compounds, function(x) length(x)) %>% unlist()
rep_n <- rep(names(n_int), n_int)
binding_partners <- unlist(Binding_Compounds)
binding_df <- data.frame(binding_partners, rep_n) %>%
  `colnames<-`(c("Protein", "Condition"))

binding_df %<>%
  plyr::join(., gene_annotation, by="Protein")%>%
  dplyr::select(Condition, Gene_names)

wide_format_2 <- second_option %>%
  mutate(Condition = factor(Condition, levels = sm_list)) %>%
  reshape2::dcast(Gene_names ~ Condition, vaue.var = "n_significant") %>%
  na.omit() %>%
  mutate(Gene_names = factor(Gene_names, levels = binding_df$Gene_names))

wide_format_2 <- wide_format_2[order(wide_format_2$Gene_names),] %>%
  `rownames<-`(.$Gene_names) %>%
  .[,-1] %>%
  as.matrix() 

binding_df %<>%
  filter(Gene_names %in% rownames(wide_format_2)) %>%
  `rownames<-`(.$Gene_names) %>%
  .[,-2]

row_ha <- ComplexHeatmap::rowAnnotation(foo = binding_df, col=list(foo = c("Betaine" = "#B5C3AD", "Glucose" = "#C78C8C", "Proline" = "#edd9bb")),
                                        annotation_label = "Binding to")
f2 = circlize::colorRamp2(seq(0, max(wide_format_2), length = 3), c("grey98",  "#8E897E", "#635D52"))

library(ComplexHeatmap)
hm <- ComplexHeatmap::Heatmap(wide_format_2, cluster_rows = FALSE, cluster_columns = FALSE, col = f2,
                              row_names_gp = grid::gpar(fontsize = 10), show_row_names = TRUE, raster_device = "CairoPNG", 
                              left_annotation = row_ha, 
                              rect_gp = grid::gpar(col = "gray20", lwd = 1),
                              heatmap_legend_param = list(
                                title = "#Peptides", direction = "horizontal"))

ComplexHeatmap::draw(hm, heatmap_legend_side = "bottom", 
                     annotation_legend_side = "bottom")

############ 
plot_peptide_grid(all_results, names(all_results), "Control", "ASPSLLDGIVVEYYGTPTPLR") +
  theme(panel.grid = element_blank())

plot_peptide_grid(all_results, names(all_results), "Control", "IAAVAEDGEPCVTYIGADGAGHYVK") +
  theme(panel.grid = element_blank())

####
aggregation_percentage <- lapply(significant_TMAO, function(x) sum(x$scores$Aggregation > 0)/length(x$scores$Aggregation))
aggregation_df <- data.frame(Condition = names(aggregation_percentage),
                             Aggregation = unlist(aggregation_percentage))
aggregation_df$Condition <- factor(aggregation_df$Condition, levels = sm_list)
ggplot(aggregation_df, aes(x=Condition, fill=Condition, y=Aggregation)) +
  geom_bar(  stat = "identity", position = "dodge", col="gray20") +
  scale_fill_manual(values = c(my_palette[sm_list])) +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("Condition") +
  ylim(c(0,0.2)) +
  ylab("Fraction of Aggregation profiles") +
  xlab("")

#### 
selection <- TPP_significant[["proline"]]$Combined %>%
  na.omit() %>%
  filter(define_peak == "aggregation")


peptides <- selection$peptide[selection$peak_sum > 0]
plot_peptide_grid(TPP, names(TPP), peptide=peptides[6])   

plot_peptide_grid(all_results, names(all_results), control="Control", peptide="FAAQAVMGSAK")                                             


###

Comparison_validations <- data.frame(Condition = c("TMAO", "Glucose", "Trehalose", "Proline", "Betaine", "Glycerol"), 
                                     DSF = c(5.05, 4.7, 3.29, 2.23, 3.32, 1.63), 
                                     TPP = c(5.11, 4.15, 3.8, 0.52, 2.34, 0.65))

ggplot(Comparison_validations, aes(x=DSF, y=TPP) ) +
  geom_smooth(method="rlm", col = "gray60", lty="dashed", alpha=0.2) + 
  geom_point(aes(fill=Condition), colour="black",pch=21, size=4) + 
  scale_fill_manual(values = my_palette[1:6]) +
  #facet_wrap(~Protein, scales="free") +
  theme(axis.line = element_line(), panel.background = element_blank(),
        panel.grid = element_blank(), legend.position = "none")+
  ylab("Melting temperature TPP") +
  xlab("dT for DSF measurements")

summary(lm(Comparison_validations$DSF~Comparison_validations$TPP))

### Boxplots 
sm_list <- c("TMAO", "Betaine", "Glycerol", "Proline", "Trehalose", "Glucose")
all_boxplot <- NULL
for(sm in sm_list) {
  
  df <- significant_TMAO[[sm]]$Protein_level[significant_TMAO[[sm]]$Protein_level$Protein_stabilisation > 0,] %>%
    mutate(Condition = sm)
  
  if(any(colnames(df) %in% "mean_p_score")) {
    df %<>%
      dplyr::select(-mean_p_score)
  }
  
  all_boxplot <- rbind(all_boxplot, df)
    
}
all_boxplot$Condition <- factor(all_boxplot$Condition, levels = sm_list)
# all_boxplot$Condition  <- 
#   recode_factor(all_boxplot$Condition, "" = "Beta", "Glycerol" = "Glyc", "Proline" = "Pro", "Trehalose" = "Treh", "Glucose" = "Gluco")
#all_boxplot$Condition <- factor(all_boxplot$Condition, levels=c("TMAO", "Beta", "Glyc", "Pro", "Treh", "Gluco"))
Figure_2C <- ggplot(all_boxplot, aes(x=Condition, y=Protein_stabilisation, fill=Condition)) +
  geom_violin() +
  geom_boxplot(width=0.2, outlier.alpha = 0) +
  scale_fill_manual(values=my_palette) +
  #stat_compare_means(label = "p.signif", method = "wilcox.test",
  #                   ref.group = ".all.", size=8, vjust=1) +
  theme_classic() +
  #ylim(c(0,45))+
  #geom_hline(yintercept = mean(all_boxplot$Protein_stabilisation), linetype = 2, alpha=0.5, lwd=1)+ # Add horizontal line at base mean
  theme(legend.position = "none",
        axis.text.x = element_text(size=15, angle=45, vjust=1, hjust=1, color="black"), 
        axis.text.y = element_text(size=15,  color="black"), 
        axis.title = element_text(size=15),
        axis.title.x = element_blank()) +  
  xlab("") +
  ylab("Protein stabilisation score")

ggsave("/Users/moni/Documents/Phd/Osmolyte_paper/Figure2C.pdf",Figure_2C, height = 3.5, width = 4)

ggplot(all_boxplot, aes(x=Condition, y=Protein_stabilisation, fill=Condition)) +
  #geom_violin() +
  geom_boxplot() +
  scale_fill_manual(values=my_palette) +
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     
                     ref.group = ".all.") +
  theme_classic() +
  geom_hline(yintercept = mean(all_boxplot$Protein_stabilisation), linetype = 2, alpha=0.5, lwd=1)+ # Add horizontal line at base mean
  theme(legend.position = "none",
        axis.text = element_text(size=11)) +
  xlab("") +
  ylab("Protein stabilisation score")
ggsave("/Users/moni/Documents/Phd/Experiments/Final_data/Boxplot_all_stabilisation_pvalue_ind.pdf")

all_barplot <- lapply(significant_all, function(x) sum(x$Protein_level$Protein_stabilisation > 0)/ length(x$Protein_level$Protein_stabilisation))
df_barplot <- data.frame(Condition = names(all_barplot), 
                         Percentage = unlist(all_barplot))

df_barplot$Condition <-factor(df_barplot$Condition, levels=sm_list)
ggplot(df_barplot, aes(x=Condition, y=Percentage, fill=Condition)) +
  #geom_violin() +
  geom_bar(stat="identity", col="black") +
  scale_fill_manual(values=my_palette) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size=15), 
        axis.title = element_text(size=15)) +  
  xlab("") +
  ylab("Fraction stabilised proteins")
ggsave("/Users/moni/Documents/Phd/Experiments/Final_data/Barplot_all_stabilisation.pdf")

### Correlations for Glucokinase, Ribokinase, Fructokinase
selection_df <- data.frame(Protein = c("P0A6V8", "P23917", "P0A9J6"), 
                           names_proteins = c("Glucokinase", "Fructokinase", "Ribokinase"))

mean_values <- lapply(significant_all, FUN = function(x) mean(x$Protein_level$Protein_stabilisation )) 
df_mean <- data.frame(Condition = names(mean_values), 
                 Mean_stability = unlist(mean_values))

selection_all <- NULL
for(sm in sm_list){
  
  selected_df <- significant_all[[sm]]$Protein_level[significant_all[[sm]]$Protein_level$Protein %in% selection_df$Protein,] %>%
    dplyr::select(Protein, Protein_stabilisation) %>%
    mutate(Condition = sm)
  
  selection_all <- rbind(selection_all, selected_df)
}

selection_all %<>%
  plyr::join(., selection_df, by="Protein") %>%
  plyr::join(., df_mean, by="Condition")
selection_all$names_proteins <- factor(selection_all$names_proteins, levels=selection_df$names_proteins)


ggplot(selection_all, aes(y=Protein_stabilisation, x=Mean_stability)) +
  geom_smooth(method="lm", col = "gray60", lty="dashed", alpha=0.2) + 
  geom_point(aes(fill=Condition), colour="black",pch=21, size=5) + 
  scale_fill_manual(values = my_palette[1:6]) +
  ylab("Protein stabilisation score") +
  facet_wrap(~names_proteins, scales = "free_y") +
  xlab("Mean stabilisation score") +
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(), 
        axis.text = element_text(size=11),
        strip.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.position = "none")
ggsave("/Users/moni/Documents/Phd/Experiments/Final_data/Gluco_Fructo_Ribo_correlation.pdf")

