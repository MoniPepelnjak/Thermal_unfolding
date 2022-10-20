## Author: Monika Pepelnjak
## Date: 14.08.2020
## Post-python analysis
# python post analysis functions

library(seqinr)
library(magrittr)
library(dplyr)
library(Biostrings)
library(seqinr)
library(ggplot2)
library(ggpubr)
#Function that loads the two python files and makes them usable
load_data <- function(path_dir, small_molecule) {
  
  for(i in 1:nrow(path_dir)) {
    if (file.exists(path_dir$path[i]) == FALSE) {
      print(paste("File ", path_dir$file[i], " does not exist!", sep="" ))
    }
    }
  python_fit_HT <- read.csv(path_dir$path[path_dir$small.molecule == small_molecule & path_dir$file == "python_fit_HT"])
  python_fit_FT <- read.csv(path_dir$path[path_dir$small.molecule == small_molecule & path_dir$file == "python_fit_FT"])
  
  python_fit <- rbind(python_fit_HT, python_fit_FT)
  
  python_score_FT <- read.csv(path_dir$path[path_dir$small.molecule == small_molecule & path_dir$file == "python_score_FT"])
  python_score_FT$trypticity <- "FT"
  python_score_HT <- read.csv(path_dir$path[path_dir$small.molecule == small_molecule & path_dir$file == "python_score_HT"])
  python_score_HT$trypticity <- "HT"
  
  output <- list()
  output[["python_fit"]] <- rbind(python_fit_HT, python_fit_FT)
  output[["python_score"]] <- rbind(python_score_FT, python_score_HT)
  
  library <- read.csv(path_dir$path[path_dir$small.molecule == small_molecule & path_dir$file == "library"]) %>% 
    dplyr::select(Sequence, Master.Protein.Accessions) %>%
    `colnames<-`(c("peptide", "protein"))
  
  output[["python_score"]] <- output[["python_score"]] %>% 
    merge(., library, by="peptide")
  
  output[["python_fit"]] <- output[["python_fit"]] %>% 
    merge(., library, by="peptide")
  
  output
}

load_data_FT <- function(path_dir, small_molecule) {
  
  for(i in 1:nrow(path_dir)) {
    if (file.exists(path_dir$path[i]) == FALSE) {
      print(paste("File ", path_dir$file[i], " does not exist!", sep="" ))
    }
  }
  python_fit <- read.csv(path_dir$path[path_dir$small.molecule == small_molecule & path_dir$file == "python_fit_FT"])
  python_score <- read.csv(path_dir$path[path_dir$small.molecule == small_molecule & path_dir$file == "python_score_FT"])
  output <- list()
  
  output[["python_fit"]] <- python_fit
  output[["python_score"]] <- python_score
  
  library <- read.csv(path_dir$path[path_dir$small.molecule == small_molecule & path_dir$file == "library"]) %>% 
    dplyr::select(Sequence, Master.Protein.Accessions) %>%
    `colnames<-`(c("peptide", "protein"))
  
  output[["python_score"]] <- output[["python_score"]] %>% 
    merge(., library, by="peptide")
  
  output[["python_fit"]] <- output[["python_fit"]] %>% 
    merge(., library, by="peptide")
  
  output
}


#Eucl. calculation function
eucl_calculation <- function(files, small_molecule, control = "control") {
  
  Eucl_calc <- files[["python_fit"]] %>% 
    na.omit()
  
  seq_freq <- seq(1, 101, by = 5)
  reduced_temp <- Eucl_calc$t %>% unique() %>% .[seq_freq]
  
  Eucl_calc %<>% 
    filter(.$condition %in% c(small_molecule, control),
           .$t %in% c(reduced_temp)) %>%
    ungroup() %>%
    unique()
  
  Eucl_calc%<>%
    group_by(peptide, t) %>%
    mutate(sq_diff = (y[condition == small_molecule] - y[condition == control])^2) %>% 
    ungroup()
  
  Eucl_calc %<>% 
    dplyr::select(peptide, sq_diff) %>% 
    unique() %>% 
    group_by(peptide) %>% 
    mutate(eucl = sqrt(sum(sq_diff))) %>% 
    ungroup()
  
  Eucl_calc %<>% 
    dplyr::select(peptide, eucl) %>% 
    unique()
  
  files[["python_score"]] %<>% 
    plyr::join(., Eucl_calc, by="peptide")
  
  files[["python_score"]]$score <- files[["python_score"]]$BF * files[["python_score"]]$eucl
  files
}

# calculate p-values by using the control document 
pvalue_calculations <- function(files) {
  control_dataset <- read.csv("/Users/moni/Documents/Phd/Experiments/020/scores_control.csv")
  control_mean <- mean(log2(control_dataset$score[control_dataset$c1>0.1 & control_dataset$c2>0.1]))
  control_sd <- sd(log2(control_dataset$score[control_dataset$c1>0.1 & control_dataset$c2>0.1]))
  
  pvalue <- rep(NA, length(files[["python_score"]]$score))
  
  for(i in 1:length(files[["python_score"]]$score)){
    pvalue[i] <- 1-pnorm(log2(files[["python_score"]]$score[i]), control_mean, control_sd)}
  
  files[["python_score"]]$pvalue <- pvalue
  files[["python_score"]]$padj <- p.adjust(files[["python_score"]]$pvalue, "BH")
  files
}

## Calculate peptide positions

ranges_calculations <- function(files, fasta_file) {
  Fasta_ecoli <- fasta_file
  proteins_ecoli <- names(Fasta_ecoli) 
  proteins_fasta <- lapply(strsplit(proteins_ecoli, split="\\|"), FUN = function(x) "[[" (x,2)) %>% unlist() %>% #from the FASTA get just the UniprotIDs
    intersect(., files[["python_score"]]$protein) 
  
  #comment: Like this we make sure that we map the proteins that are also in the FASTA and not for example iRT peptides we use for DIA analysis
  
  
  #Prepare empty data frame 
  position_export <- data.frame(matrix(ncol = 7, nrow=0))
  colnames(position_export) <- c("protein","start_match", "end_match", "peptide", "protein_length")
  
  
  for(i in 1:length(proteins_fasta)) {
    #extract protein sequence from the FASTA file
    protein_sequence <- Fasta_ecoli[[grep(proteins_fasta[i], names(Fasta_ecoli))[1]]] %>%
      paste(., sep="", collapse="") %>% AAString()
    
    #get the peptides from that protein
    peptides <- files[["python_score"]] %>% 
      filter(.$protein == proteins_fasta[i])
    
    peptide_sequences <- peptides$peptide %>% AAStringSet()
    
    #match the peptides with the sequence
    matches <- Biostrings::matchPDict(peptide_sequences, protein_sequence) %>% unlist() 
    
    if(length(ranges(matches)) == length(peptides$peptide)){ #some proteins have peptides that map to more than one sequence, those are excluded in the first round
      protein_length <- length(protein_sequence)
      start_match <- Biostrings::start(matches) #Give out the vector of start positions of a match
      end_match <- Biostrings::end(matches) #Give out the vector of end positions of a match
      
      #comment: if there is a region with both significant peptides, missing peptides... that region will count for both, so the sum of all the region will not be 100% 
      #This can be adjusted if needed
      
      intervals <- data.frame(protein = proteins_fasta[i],
                              start_match, 
                              end_match,
                              peptide = as.character(peptides$peptide),
                              protein_length = protein_length)
      
      intervals <- intervals[order(intervals$start_match),]
      
      #Order the levels of the peptides based on the position
      intervals$peptide <- factor(intervals$peptide, levels=unique(intervals$peptide))
      position_export <- rbind(position_export, intervals)
      
      
    }
    
  else {
    print(proteins_fasta[i])
  }
    
    
  }
  position_export
  
}

# Calculate ranges

determine_range <- function(files, small_molecule, control = "control") {
  
  small_molecule <- small_molecule
  control <- control
  
  output_ranges <- data.frame(matrix(nrow = 0, ncol=(ncol(files[["peptide_matching"]]) + 1)))
  colnames(output_ranges) <- c(colnames(files[["peptide_matching"]]), "ranges")
  
  significant_proteins <- files[["peptide_matching"]]$protein %>% 
    plyr::count(.) %>%
    filter(.$freq > 1) %>%
    .$x %>%
    as.character()
  
  for(j in 1:length(significant_proteins)) {
    one_protein <- files[["peptide_matching"]] %>% 
      filter(.$protein == significant_proteins[j])
    
    range <- 1
    output <- c()
    for(k in 1:(length(one_protein$peptide)-1)) {
      intersect <- IntervalSurgeon::intersects(x=as.matrix(one_protein[k, c(2,3)]), y=as.matrix(one_protein[k+1, c(2,3)]))
      #uni <- unions(x=as.matrix(one_protein[k, c(2,3)]), y=as.matrix(one_protein[k+1, c(2,3)]))
      intersect_length <-  ifelse(length(intersect) > 0,  intersect[2] - intersect[1] + 1, 0)
      shorter_peptide <-  min(nchar(as.character(one_protein$peptide[k])), 
                              nchar(as.character(one_protein$peptide[k+1])))
      percentage_overlap <- intersect_length/shorter_peptide
      is_overlap <- (percentage_overlap < 0.5)
      range <- range + is_overlap
      output[k] <- range
    }
    output <- c(1, output)
    one_protein$ranges <- output 
    
    output_ranges <- rbind(output_ranges, one_protein)}

  output_ranges %<>%
    dplyr::select(-c(protein)) %>%
    plyr::join(., files[["python_score"]], by="peptide")
  
  output_ranges %<>%
    group_by(protein, ranges) %>%
    mutate(start_range = min(start_match), 
           end_range = max(end_match)) %>% 
    ungroup() %>%
    mutate(gof = abs(!!ensym(small_molecule) + !!ensym(control))) %>% 
    data.frame() 
  
  output_ranges %<>%
    group_by(protein, ranges) %>%
    mutate(score_mean = mean(score),
           #  score_best_fit = score[gof == max(gof)],
           score_weighted = stats::weighted.mean(score, gof)) %>%
    mutate(overlaps = end_match[start_match == min(start_match)] > start_match ) %>%
    ungroup() %>% 
    group_by(protein) %>% 
    mutate(eucl_weighted = stats::weighted.mean(eucl, gof)) %>%
    ungroup()
  
  only_ranges <- output_ranges %>%
    dplyr::select(protein, protein_length, ranges, score_weighted, start_range, end_range, eucl_weighted ) %>% 
    unique() %>% 
    ungroup() %>% 
    mutate(range_len = end_range - start_range ) %>% 
    mutate(is_significant = score_weighted > 4) %>% 
    group_by(protein) %>% 
    mutate(covered = sum(range_len)/protein_length, 
           significant = sum(range_len[is_significant == TRUE])/protein_length) %>% 
    ungroup() %>%
    mutate(signif_covered = significant/covered)
  
  
  }

protein_level_calculation <- function(all_results, small_molecule, control = "control"){
  proteins <- all_results[[small_molecule]][["ranges"]] %>% 
    group_by(protein) %>% 
    mutate(mean_score = mean(score_weighted),
           max_score = max(score_weighted)) %>% 
    dplyr::select(protein, protein_length, covered, significant, signif_covered, mean_score, max_score) %>% 
    unique()
  
  p_value <- rep(NA, length(proteins$protein))
  for(i in 1:length(proteins$protein)) {
    one <- all_results[[small_molecule]][["ranges"]] %>%
      filter(.$protein == proteins$protein[i]) 
    
    k_test <- ks.test( log2(all_results[["control"]][["ranges"]]$score_weighted) , log2(one$score_weighted))
    p_value[i] <- k_test$p.value
  }

  proteins$prot_pvalue <- p_value
  proteins$prot_pvalue_adj <- p.adjust(proteins$prot_pvalue, "BH")
  proteins
}


protein_from_peptides <- function(all_results, small_molecule, control = "control"){
  proteins <- all_results[[small_molecule]][["python_score"]] %>% 
    group_by(protein) %>% 
    mutate(mean_score = mean(score_weighted),
           max_score = max(score_weighted)) %>% 
    dplyr::select(protein, protein_length, covered, significant, signif_covered, mean_score, max_score) %>% 
    unique()
  
  p_value <- rep(NA, length(proteins$protein))
  for(i in 1:length(proteins$protein)) {
    one <- all_results[[small_molecule]][["ranges"]] %>%
      filter(.$protein == proteins$protein[i]) 
    
    k_test <- ks.test( log2(all_results[["control"]][["ranges"]]$score_weighted) , log2(one$score_weighted))
    p_value[i] <- k_test$p.value
  }
  
  proteins$prot_pvalue <- p_value
  proteins$prot_pvalue_adj <- p.adjust(proteins$prot_pvalue, "BH")
  proteins
}


plot_individual_peptides <- function(all_results, small_molecule, control = "control", peptide_set) {
  
  peptide_set <- all_results[[small_molecule]][["python_score"]]$peptide[all_results[[small_molecule]][["python_score"]]$peptide %in% peptide_set] %>% unique()
  
  split_peptides <- ceiling(length(peptide_set)/12)

    for(i in 1:split_peptides) {
    test <-  all_results[[small_molecule]][["python_fit"]][all_results[[small_molecule]][["python_fit"]]$peptide %in% peptide_set[((12*(i-1))+1) : ((12*(i-1))+12)], ] %>%
      filter(condition != "full model")
    
    test$condition <- factor(test$condition, levels = c(small_molecule, control))
    levels(test$peptide) <- peptide_set
    
    sm_color <- c(Control = "#8D8D92",
                  Proline = "#3C4285", 
                  Betaine = "#2C938B",
                  Trehalose = "#A84343",
                  Glucose = "#6B1818",
                  TMAO_2 = "#6B1818",
                  PK_0.01 = "#A84343",
                  TMAO = "#2C8AA4", 
                  TMAO1 = "#2C8AA4", 
                  ATP = "#3C4285", 
                  Low_PK = "#4AB698", 
                  Glycerol_2 = "#87667F")

    plot3 <- ggplot(test) +  
    geom_ribbon(data = test[test$type == "fitted",], aes(ymin=conflik_lower, ymax=conflik_upper, x= t, fill=condition), alpha=0.2) +
    geom_point(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), size=2) + 
    geom_line(data = test[test$type == "fitted",], aes(y=y, x=t, col=condition), lwd=1.2) +
    scale_color_manual(values = c(sm_color[names(sm_color) == small_molecule], "#8F8F8F")) +
    scale_fill_manual(values = c(sm_color[names(sm_color) == small_molecule], "#8F8F8F"))+
    #  geom_smooth(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), method="loess") +
    facet_wrap(~peptide) +
      theme_minimal() +
      theme(axis.line = element_blank(), panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(), legend.position = "none") +
      annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
      annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) + 
      xlab("Temperature") +
      ylab("Scaled intensity")
  return(plot3)}

}

plot_peptides_one_condition <- function(all_results, small_molecule, condition_plot = "Control", group1, group2, color_selected = c("gray20", "gray50")) {
  
  peptide_set <- c(group1, group2)
  peptide_set <- all_results[[small_molecule]][["python_score"]]$peptide[all_results[[small_molecule]][["python_score"]]$peptide %in% peptide_set] %>% unique()
  
  split_peptides <- ceiling(length(peptide_set)/12)
  
  for(i in 1:split_peptides) {
    test <-  all_results[[small_molecule]][["python_fit"]][all_results[[small_molecule]][["python_fit"]]$peptide %in% peptide_set[((12*(i-1))+1) : ((12*(i-1))+12)], ] %>%
      filter(condition == condition_plot)
    
    #test$condition <- factor(test$condition, levels = c(small_molecule, control))
    levels(test$peptide) <- peptide_set
    test$color_peptide <- ifelse(test$peptide %in% group1, color_selected[1], color_selected[2])
    plot3 <- ggplot(test) +  
      geom_ribbon(data = test[test$type == "fitted",], aes(ymin=conflik_lower, ymax=conflik_upper, x= t, group=peptide, fill=color_peptide), alpha=0.2) +
      geom_point(data = test[test$type == "measured",], aes(y=y, x=t, group=peptide, col=color_peptide), size=2) + 
      geom_line(data = test[test$type == "fitted",], aes(y=y, x=t, group=peptide, col=color_peptide), lwd=1.2) +
      scale_color_manual(values = color_selected) +
      scale_fill_manual(values = color_selected)+
      #  geom_smooth(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), method="loess") +
      #facet_wrap(~peptide) +
      theme_minimal() +
      theme(axis.line = element_blank(), panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(), legend.position = "none") +
      annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
      annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) + 
      xlab("Temperature") +
      ylab("Scaled intensity")
    return(plot3)}
  
}


plot_individual_peptides_conf <- function(all_results, small_molecule, control = "control", peptide_set) {
  
  peptide_set <- all_results[[small_molecule]][["python_score"]]$peptide[all_results[[small_molecule]][["python_score"]]$peptide %in% peptide_set] %>% unique()
  
  split_peptides <- ceiling(length(peptide_set)/12)
  
  for(i in 1:split_peptides) {
    test <-  all_results[[small_molecule]][["python_fit"]][all_results[[small_molecule]][["python_fit"]]$peptide %in% peptide_set[((12*(i-1))+1) : ((12*(i-1))+12)], ] %>%
      filter(condition != "full model")
    
    test$condition <- factor(test$condition, levels = c(small_molecule, control))
    levels(test$peptide) <- peptide_set
    
    sm_color <- c(control = "#8D8D92",
                  proline = "#3C4285", 
                  betaine = "#2C938B",
                  trehalose = "#A84343",
                  glucose = "#6B1818",
                  TMAO_2 = "#6B1818",
                  PK_0.01 = "#A84343",
                  TMAO = "#2C8AA4", 
                  TMAO1 = "#2C8AA4", 
                  ATP = "#3C4285", 
                  low_PK = "#A84343")
    
    plot3 <- ggplot(test) +  
      geom_ribbon(data = test[test$type == "fitted",], aes(ymin=conf_lower, ymax=conf_upper, x= t, fill=condition), alpha=0.2) +
      geom_point(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), size=2) + 
      geom_line(data = test[test$type == "fitted",], aes(y=y, x=t, col=condition), lwd=1.2) +
      scale_color_manual(values = c(sm_color[names(sm_color) == small_molecule], "#8F8F8F")) +
      scale_fill_manual(values = c(sm_color[names(sm_color) == small_molecule], "#8F8F8F"))+
      #  geom_smooth(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), method="loess") +
      facet_wrap(~peptide) +
      ylab("Scaled abundance") +
      xlab("Temperature [C]")+
      theme_bw() +
      theme(legend.position="none")
    return(plot3)}
  
}

plot_individual_peptides2 <- function(all_results, small_molecule1, small_molecule2, control = "control", peptide_set) {
  
  peptide_set <- all_results[[small_molecule1]][["python_score"]]$peptide[all_results[[small_molecule1]][["python_score"]]$peptide %in% peptide_set] %>% unique()

  split_peptides <- ceiling(length(peptide_set)/12)
  
  for(i in 1:split_peptides) {
    test1 <-  all_results[[small_molecule1]][["python_fit"]][all_results[[small_molecule1]][["python_fit"]]$peptide %in% peptide_set[((12*(i-1))+1) : ((12*(i-1))+12)], ] %>%
      filter(condition != "full model") %>%
      mutate(exp="a") %>% 
      unique()
    
    test2 <-  all_results[[small_molecule2]][["python_fit"]][all_results[[small_molecule2]][["python_fit"]]$peptide %in% peptide_set[((12*(i-1))+1) : ((12*(i-1))+12)], ] %>%
      filter(condition != "full model") %>%
      mutate(exp="b") %>% 
      unique()
    
    test <- rbind(test1, test2)
    
    test$condition <- factor(test$condition, levels = c(small_molecule1, small_molecule2, control))
    levels(test$peptide) <- peptide_set
    
    sm_color <- c(control = "#8D8D92",
                  proline = "#3C4285", 
                  betaine = "#2C938B",
                  trehalose = "#A84343",
                  glucose = "#6B1818",
                  TMAO_2 = "#3C4285",
                  TMAO2 = "#A84343",
                  TMAO1 = "#3C4285",
                  glycerol_2 ="#A84343",
                  TMAO = "#2C8AA4", 
                  proline0 = "#2C8AA4", 
                  proline_2 = "#0093BC")
    
    plot3 <- ggplot(test) +  
      geom_ribbon(data = test[test$type == "fitted",], aes(ymin=conflik_lower, ymax=conflik_upper, x= t, fill=condition, group=interaction(exp, condition)), alpha=0.2) +
      geom_point(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), size=2) + 
      geom_line(data = test[test$type == "fitted",], aes(y=y, x=t, col=condition, group=interaction(exp,condition)), lwd=1.2) +
      scale_color_manual(values = c(sm_color[names(sm_color) == small_molecule1],sm_color[names(sm_color) == small_molecule2], "#8F8F8F")) +
      scale_fill_manual(values = c(sm_color[names(sm_color) == small_molecule1], sm_color[names(sm_color) == small_molecule2], "#8F8F8F"))+
      #  geom_smooth(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), method="loess") +
      facet_wrap(~peptide) +
      ylab("Scaled abundance") +
      xlab("Temperature [C]")+
      theme_bw() +
      theme(legend.position="none")
    print(plot3)}
  
}


plot_peptide_grid <- function(all_results, sm_all, control = "control", peptide){
 # sm_all <- c("betaine")
  test <- data.frame(matrix(nrow=0, ncol = (ncol(all_results[[sm_all[1]]]$python_fit)+1)))
  colnames(test) <- c(colnames(all_results[[sm_all[1]]]$python_fit), "exp")
  
  for(sm in sm_all){
    
    test1 <-  all_results[[sm]][["python_fit"]][all_results[[sm]][["python_fit"]]$peptide %in% peptide, ] %>%
      filter(condition != "full model") %>%
      mutate(exp=sm) %>% 
      unique()
    
    test <- rbind(test, test1)
  }
  
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
  
  
  plot3 <- ggplot(test) +  
    geom_ribbon(data = test[test$type == "fitted",], aes(ymin=conflik_lower, ymax=conflik_upper, x= t, fill=condition, group=interaction(exp, condition)), alpha=0.2) +
    geom_point(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), size=2) + 
    geom_line(data = test[test$type == "fitted",], aes(y=y, x=t, col=condition, group=interaction(exp,condition)), lwd=1.2) +
    scale_color_manual(values = my_palette) +
    scale_fill_manual(values = my_palette )+
    #  geom_smooth(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), method="loess") +
    facet_wrap(~exp) +
    ylab("Scaled abundance") +
    xlab("Temperature [C]")+
    theme_bw() +
    theme(legend.position="none") +
    labs(title = peptide)
  return(plot3)
  
}


plot_individual_peptides3 <- function(all_results, small_molecule1, small_molecule2, small_molecule3, control = "control", peptide_set) {
  
  peptide_set <- all_results[[small_molecule1]][["python_score"]]$peptide[all_results[[small_molecule1]][["python_score"]]$peptide %in% peptide_set] %>% unique()
  
  split_peptides <- ceiling(length(peptide_set)/12)
  
  for(i in 1:split_peptides) {
    test1 <-  all_results[[small_molecule1]][["python_fit"]][all_results[[small_molecule1]][["python_fit"]]$peptide %in% peptide_set[((12*(i-1))+1) : ((12*(i-1))+12)], ] %>%
      filter(condition != "full model") %>%
      mutate(exp="a") %>% 
      unique()
    
    test2 <-  all_results[[small_molecule2]][["python_fit"]][all_results[[small_molecule2]][["python_fit"]]$peptide %in% peptide_set[((12*(i-1))+1) : ((12*(i-1))+12)], ] %>%
      filter(condition != "full model") %>%
      mutate(exp="b") %>% 
      unique()
    
    test3 <-  all_results[[small_molecule3]][["python_fit"]][all_results[[small_molecule3]][["python_fit"]]$peptide %in% peptide_set[((12*(i-1))+1) : ((12*(i-1))+12)], ] %>%
      filter(condition != "full model") %>%
      mutate(exp="c") %>% 
      unique()
    
    test <- rbind(test1, test2) %>% rbind(., test3)
    
    test$condition <- factor(test$condition, levels = c(small_molecule1, small_molecule2, small_molecule3, control))
    levels(test$peptide) <- peptide_set
    
    sm_color <- c(control = "#8D8D92",
                  proline = "#3C4285", 
                  betaine = "#2C938B",
                  trehalose = "#A84343",
                  glucose = "#6B1818",
                  TMAO_2 = "#3C4285",
                  TMAO2 = "#A84343",
                  TMAO1 = "#3C4285",
                  glycerol_2 ="#A84343",
                  TMAO = "#2C8AA4", 
                  proline0 = "#2C8AA4", 
                  proline_2 = "#0093BC")
    
    plot3 <- ggplot(test) +  
      geom_ribbon(data = test[test$type == "fitted",], aes(ymin=conflik_lower, ymax=conflik_upper, x= t, fill=condition, group=interaction(exp, condition)), alpha=0.2) +
      geom_point(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), size=2) + 
      geom_line(data = test[test$type == "fitted",], aes(y=y, x=t, col=condition, group=interaction(exp,condition)), lwd=1.2) +
      scale_color_manual(values = c(sm_color[names(sm_color) == small_molecule1],sm_color[names(sm_color) == small_molecule2], sm_color[names(sm_color) == small_molecule3], "#8F8F8F")) +
      scale_fill_manual(values = c(sm_color[names(sm_color) == small_molecule1], sm_color[names(sm_color) == small_molecule2], sm_color[names(sm_color) == small_molecule3], "#8F8F8F"))+
      #  geom_smooth(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), method="loess") +
      facet_wrap(~peptide, nrow = 1) +
      ylab("Scaled abundance") +
      xlab("Temperature [C]")+
      theme_bw() +
      theme(legend.position="none")
    return(plot3)}
}

plot_peptides_from_protein <- function(all_results, small_molecule, uniprot_code, n_split = 12) {
#  uniprot_code <- "P0A6Y8"
#  small_molecule <- "ATP"
  peptide_set <- all_results[[small_molecule]][["peptide_matching"]][all_results[[small_molecule]][["peptide_matching"]]$protein == uniprot_code,] %>%
    .[order(.$start_match),] %>%
    mutate(peptide = factor(peptide, levels=unique(.$peptide)))
  
  split_peptides <- ceiling(length(peptide_set$peptide)/n_split)

    for(i in 1:split_peptides) {
    test <-  all_results[[small_molecule]][["python_fit"]][all_results[[small_molecule]][["python_fit"]]$peptide %in% peptide_set$peptide[((n_split*(i-1))+1) : ((n_split*(i-1))+n_split)], ]
    
    plot3 <- ggplot(test) +  
      geom_ribbon(data = test[test$type == "fitted",], aes(ymin=conflik_lower, ymax=conflik_upper, x= t, fill=condition), alpha=0.2) +
      geom_point(data = test[test$type == "measured",], aes(y=y, x=t, col=condition)) + 
      geom_line(data = test[test$type == "fitted",], aes(y=y, x=t, col=condition)) +
      scale_color_manual(values = c("#00059C", "#8F8F8F", "#767791")) +
      scale_fill_manual(values = c("#00059C", "#8F8F8F", "#767791"))+
      #  geom_smooth(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), method="loess") +
      facet_wrap(~peptide)
    print(plot3)}
  
}

correlation_analysis <- function(all_results, small_molecule, control = "control", protein, domain_borders = c()) {
  peptides <- all_results[[small_molecule]][["peptide_matching"]]$peptide[all_results[[small_molecule]][["peptide_matching"]]$protein == protein 
                                                               #           & all_results[[small_molecule]][["peptide_matching"]]$end_match > 179 
                                                                          ]
   
  proteins_sm <- all_results[[small_molecule]][["python_fit"]] %>% 
    filter(condition == small_molecule, 
           peptide %in% peptides,
           type == "fitted" 
           #        peptide %in% second_domain_proteins$peptide
    ) %>% 
    reshape2::dcast(peptide~t, value.var = "y", fun.aggregate = mean)
  
  proteins_control <- all_results[[small_molecule]][["python_fit"]] %>% 
    filter(condition == "control", 
           peptide %in% peptides,
           type == "fitted" 
           #        peptide %in% second_domain_proteins$peptide
    ) %>% 
    reshape2::dcast(peptide~t, value.var = "y", fun.aggregate = mean)
  
  
  corrplot::corrplot(cor((proteins_sm[,-1])))
  corrplot::corrplot(cor((proteins_control[,-1])))
  
  cor_values_control <- cor((proteins_control[,-1]))
  cor_values_sm <- cor((proteins_sm[,-1])) 
  
  plot(as.numeric(colnames(proteins_control[-1])), abs(cor_values_control[1,]))
  plot(as.numeric(colnames(proteins_control[-1])), abs(cor_values_sm[1,]))
  
  plot(as.numeric(colnames(proteins_control[-1])), abs(cor_values[min(abs(cor_values[1,])) == abs(cor_values[1,]),]))
  plot(as.numeric(colnames(proteins_control[-1])), abs(cor_values[50,]))
  
  min(abs(cor_values[1,])) == abs(cor_values[1,])

  tree <- hclust(dist(cor_values_control), method="complete")
  plot(tree)
  cut_tree <- cutree(tree, k = NULL, h = 1.5)
  plot(cut_tree/15 ~ names(cut_tree), cex=0.5)
}

library(drc)
summarise_melting_temperature <- function(all_results, small_molecule, sm, control, protein, k, 
                                          domain_start = 1, domain_end = 10000) {

peptides <- all_results[[small_molecule]][["peptide_matching"]]$peptide[all_results[[small_molecule]][["peptide_matching"]]$protein == protein 
                                                                       & all_results[[small_molecule]][["peptide_matching"]]$end_match <= domain_end 
                                                                       & all_results[[small_molecule]][["peptide_matching"]]$start_match >= domain_start 
                                                                        ]
#peptides <- intersect(peptides, all_results[[small_molecule]][["python_score"]]$peptide[all_results[[small_molecule]][["python_score"]]$control>0.3])

proteins_sm <- all_results[[small_molecule]][["python_fit"]] %>% 
  filter(condition == sm, 
         peptide %in% peptides,
         type == "fitted" 
         #        peptide %in% second_domain_proteins$peptide
  ) %>% 
  reshape2::dcast(peptide~t, value.var = "y", fun.aggregate = mean)

proteins_control <- all_results[[small_molecule]][["python_fit"]] %>% 
  filter(condition == control, 
         peptide %in% peptides,
         type == "fitted" 
         #        peptide %in% second_domain_proteins$peptide
  ) %>% 
  reshape2::dcast(peptide~t, value.var = "y", fun.aggregate = mean)

cor_values_control <- cor((proteins_control[,-1]))
cor_values_sm <- cor((proteins_sm[,-1])) 

tree_control <- hclust(dist(cor_values_control), method="complete")
cut_tree_control <- cutree(tree_control, k = k, h = NULL)

tree_sm <- hclust(dist(cor_values_sm), method="complete")
cut_tree_sm <- cutree(tree_sm, k = k, h = NULL)

frame_control <- data.frame(values = (cut_tree_control-1)/(k-1), temp = as.numeric(names(cut_tree_control)), condition = control)
frame_sm <- data.frame(values = (cut_tree_sm-1)/(k-1), temp = as.numeric(names(cut_tree_sm)), condition = sm)
frame1 <- rbind(frame_control, frame_sm)

newdata_sm <- expand.grid(conc=seq(37, 76, length=100))
newdata_control <- expand.grid(conc=seq(37, 76, length=100))

m_control <- drm(frame_control$values~frame_control$temp, fct=LL.3())
m_sm <- drm(frame_sm$values~frame_sm$temp, fct=LL.3())

proteins_control %<>%
  `colnames<-`(c(round(as.numeric(colnames(.[,-1])), 2)))

proteins_sm %<>%
  `colnames<-`(c(round(as.numeric(colnames(.[,-1])), 2)))

(corrplot::corrplot(cor((proteins_sm[,-1])), tl.cex=0.5, tl.col="gray40"))
(corrplot::corrplot(cor((proteins_control[,-1])),tl.cex=0.5, tl.col="gray40"))

s_control <- summary(m_control)
s_control$coefficients[3]

s_sm <- summary(m_sm)
s_sm$coefficients[3]

pm_sm <- predict(m_sm, newdata=newdata_sm, interval="confidence")
newdata_sm$p <- pm_sm[,1]
newdata_sm$pmin <- pm_sm[,2]
newdata_sm$pmax <- pm_sm[,3]
newdata_sm$condition <- sm

pm_control <- predict(m_control, newdata=newdata_control, interval="confidence")
newdata_control$p <- pm_control[,1]
newdata_control$pmin <- pm_control[,2]
newdata_control$pmax <- pm_control[,3]
newdata_control$condition <- control

newdata <- rbind(newdata_control, newdata_sm)

plot1 <- ggplot(frame1, aes(x = temp, y = values)) +
  geom_point( aes(col=condition)) +
  geom_ribbon(data=newdata, aes(x=conc, y=p, ymin=pmin, ymax=pmax, col=condition, fill=condition), alpha=0.2) +
  geom_line(data=newdata, aes(x=conc, y=p, col=condition)) + 
  coord_trans(x="log") +
  xlab("Temperature [C]") + ylab("cluster split") + 
  geom_label(x=40, y=0.85, 
            label=paste("Tm control:",  round(s_control$coefficients[3], 2), "\n", 
                                      "Tm sm:", round(s_sm$coefficients[3], 2), sep=""))+
  labs(title = protein)
print(plot1)
}

melting_T_calculation <- function(all_results, experiment, condition1, condition2){
  all_proteins <- plyr::count(all_results[[experiment]][["peptide_matching"]]$protein) %>% 
    filter(.$freq > 4) %>% 
    .$x 
  
  Tm <- rep(NA, length(all_proteins))
  Tm_error <- rep(NA, length(all_proteins))
    for(i in 1:length(all_proteins)) {
    peptides <- all_results[[experiment]][["peptide_matching"]]$peptide[all_results[[experiment]][["peptide_matching"]]$protein == all_proteins[i]]
    
    proteins_control <- all_results[[experiment]][["python_fit"]] %>% 
      filter(condition == condition1,
             peptide %in% peptides,
             type == "fitted"
             #        peptide %in% second_domain_proteins$peptide
      ) %>% 
      reshape2::dcast(peptide~t, value.var = "y", fun.aggregate = mean)
    
    cor_values_control <- cor((proteins_control[,-1]))
    
    tree_control <- hclust(dist(cor_values_control), method="complete")
    cut_tree_control <- cutree(tree_control, k = k, h = NULL)
    
    frame_control <- data.frame(values = (cut_tree_control-1)/(k-1), temp = as.numeric(names(cut_tree_control)), condition = control)
    frame_sm <- data.frame(values = (cut_tree_sm-1)/(k-1), temp = as.numeric(names(cut_tree_sm)), condition = sm)
    frame1 <- rbind(frame_control, frame_sm)
    
    newdata_sm <- expand.grid(conc=seq(37, 76, length=100))
    newdata_control <- expand.grid(conc=seq(37, 76, length=100))
    
    m_control <- drm(frame_control$values~frame_control$temp, fct=LL.3())
    
   sumar <-  summary(m_control)
   Tm[i] <- sumar$coefficients[3]
   Tm_error[i] <-   sumar$coefficients[6]
  
  }
  
  export_table <- data.frame(
    protein = all_proteins, 
    Tm = Tm, 
    Tm_error = Tm_error
  )
  
  export_table
  
}


protein_correlations <- function(all_results, small_molecule, sm, control, protein,
                                                                  domain_start = 1, domain_end = 10000, title) {
  
  peptides <- all_results[[small_molecule]][["peptide_matching"]]$peptide[all_results[[small_molecule]][["peptide_matching"]]$protein == protein 
                                                                          & all_results[[small_molecule]][["peptide_matching"]]$end_match <= domain_end 
                                                                          & all_results[[small_molecule]][["peptide_matching"]]$start_match >= domain_start 
                                                                          ]
  #peptides <- intersect(peptides, all_results[[small_molecule]][["python_score"]]$peptide[all_results[[small_molecule]][["python_score"]]$control>0.3])
  
  temp_reduced <- all_results[[small_molecule]][["python_fit"]]$t %>% 
    unique() %>%  
    .[order(.)] %>% .[seq(1, 100, by=2)]
  

  
  proteins_sm <- all_results[[small_molecule]][["python_fit"]] %>% 
    filter(condition == sm, 
           peptide %in% peptides,
           type == "fitted",
           t %in% temp_reduced
           #        peptide %in% second_domain_proteins$peptide
    ) 
  
  for_ordering <- proteins_sm %>%
    filter(t %in% 37) %>%
    .[order(.$y),] %>%
    unique()
  
  proteins_sm %<>%
    mutate(peptide = factor(peptide, levels = for_ordering$peptide)) %>%
    mutate(t = round(t, 2)) %>%
    mutate(peptide_examples = ifelse( peptide %in% for_ordering$peptide[1:6], "beginning", 
                                      ifelse(peptide %in% for_ordering$peptide[(round(length(for_ordering$peptide)/2,0)-5):(round(length(for_ordering$peptide)/2,0))], "middle",
                                             ifelse(peptide %in% for_ordering$peptide[(length(for_ordering$peptide)-5):length(for_ordering$peptide)], "end", "not shown"))) )
  proteins_sm$peptide_examples <- factor(proteins_sm$peptide_examples, levels = c("beginning", "middle", "end", "not shown"))
  col_palette_moni <- c("#fb9092", "#c03d56", "#782d4d", "#8D8D92")  
  ggplot(proteins_sm, aes(x=peptide, y=y, col=peptide_examples)) +
    geom_point() +
    facet_wrap(~t) +
    scale_color_manual(values = col_palette_moni) +
    labs(title = title)
}

plot_just_corplot <- function(all_results, small_molecule, sm, control, protein,
                              domain_start = 1, domain_end = 10000) {

  peptides <- all_results[[small_molecule]][["peptide_matching"]]$peptide[all_results[[small_molecule]][["peptide_matching"]]$protein == protein 
                                                                          & all_results[[small_molecule]][["peptide_matching"]]$end_match <= domain_end 
                                                                          & all_results[[small_molecule]][["peptide_matching"]]$start_match >= domain_start 
                                                                          ]
  
  proteins_sm <- all_results[[small_molecule]][["python_fit"]] %>% 
    filter(condition == sm, 
           peptide %in% peptides,
           type == "fitted" 
           #        peptide %in% second_domain_proteins$peptide
    ) %>% 
    reshape2::dcast(peptide~t, value.var = "y", fun.aggregate = mean)
  
  proteins_control <- all_results[[small_molecule]][["python_fit"]] %>% 
    filter(condition == control, 
           peptide %in% peptides,
           type == "fitted" 
           #        peptide %in% second_domain_proteins$peptide
    ) %>% 
    reshape2::dcast(peptide~t, value.var = "y", fun.aggregate = mean)
  
  proteins_control %<>%
    `colnames<-`(c(round(as.numeric(colnames(.[,-1])), 2)))
  
  proteins_sm %<>%
    `colnames<-`(c(round(as.numeric(colnames(.[,-1])), 2)))
  
  corrplot::corrplot(cor((proteins_control[,-1])),tl.cex=0.5, tl.col="gray40", title=protein)
  corrplot::corrplot(cor((proteins_sm[,-1])), tl.cex=0.5, tl.col="gray40", title=protein)
  
}

plot_minvalues <- function(all_results, small_molecule, sm, control, protein,
                           domain_start = 1, domain_end = 10000, title) {
  
  peptides <- all_results[[small_molecule]][["peptide_matching"]]$peptide[all_results[[small_molecule]][["peptide_matching"]]$protein == protein 
                                                                          & all_results[[small_molecule]][["peptide_matching"]]$end_match <= domain_end 
                                                                          & all_results[[small_molecule]][["peptide_matching"]]$start_match >= domain_start 
                                                                          ]

  proteins_sm <- all_results[[small_molecule]][["python_fit"]] %>% 
    filter(condition == sm, 
           peptide %in% peptides,
           type == "fitted" 
           #        peptide %in% second_domain_proteins$peptide
    ) %>% 
    reshape2::dcast(peptide~t, value.var = "y", fun.aggregate = mean)
  
  proteins_control <- all_results[[small_molecule]][["python_fit"]] %>% 
    filter(condition == control, 
           peptide %in% peptides,
           type == "fitted" 
           #        peptide %in% second_domain_proteins$peptide
    ) %>% 
    reshape2::dcast(peptide~t, value.var = "y", fun.aggregate = mean)
  
  proteins_control %<>%
    `colnames<-`(c("peptide", round(as.numeric(colnames(.[,-1])), 2)))
  
  proteins_sm %<>%
    `colnames<-`(c("peptide", round(as.numeric(colnames(.[,-1])), 2)))
  
  values_control <- data_frame(condition = "control", absolute_correlation = abs(cor((proteins_control[,-1]))[1,]), temperature = as.numeric(names(proteins_control[,-1]))) 
  values_sm <- data_frame(condition = small_molecule, absolute_correlation = abs(cor((proteins_sm[,-1]))[1,]), temperature = as.numeric(names(proteins_sm[,-1]))) 
  
  values <- rbind(values_control, values_sm)
  
  min_value_control <- mean(values_control$temperature[values_control$absolute_correlation<0.2])
  min_value_sm <- mean(values_sm$temperature[values_sm$absolute_correlation<0.2])
  values$condition <- factor(values$condition, levels=c(sm, control))
  
  ggplot(values, aes(x=temperature, y=absolute_correlation)) +
  geom_point(aes(col=condition)) +
  labs(title = title) +
  geom_label(x=40, y=0.85, 
             label=paste("Tm control:",  round(min_value_control, 2), "\n", 
                         "Tm sm:", round(min_value_sm, 2), sep=""))
}


plot_example_peptides <- function(all_results, small_molecule, sm, control, protein,
                                 domain_start = 1, domain_end = 10000) {

  
  peptides <- all_results[[small_molecule]][["peptide_matching"]]$peptide[all_results[[small_molecule]][["peptide_matching"]]$protein == protein 
                                                                          & all_results[[small_molecule]][["peptide_matching"]]$end_match <= domain_end 
                                                                          & all_results[[small_molecule]][["peptide_matching"]]$start_match >= domain_start 
                                                                          ]
  #peptides <- intersect(peptides, all_results[[small_molecule]][["python_score"]]$peptide[all_results[[small_molecule]][["python_score"]]$control>0.3])
  
  
  proteins_sm <- all_results[[small_molecule]][["python_fit"]] %>% 
    filter(condition == sm, 
           peptide %in% peptides,
           type == "fitted",
           t == 37
           #        peptide %in% second_domain_proteins$peptide
    ) 
  
  for_ordering <- proteins_sm %>%
    filter(t %in% 37) %>%
    .[order(.$y),] %>%
    unique() %>% 
    .$peptide %>%
    unique()
  
  test_start <-  all_results[[small_molecule]][["python_fit"]][all_results[[small_molecule]][["python_fit"]]$peptide %in% for_ordering[1:6], ]
  test_middle <-  all_results[[small_molecule]][["python_fit"]][all_results[[small_molecule]][["python_fit"]]$peptide %in% for_ordering[round((length(for_ordering)/2)-5, 0):round(length(for_ordering)/2, 0)], ]
  test_end <-  all_results[[small_molecule]][["python_fit"]][all_results[[small_molecule]][["python_fit"]]$peptide %in% for_ordering[(length(for_ordering)-5):length(for_ordering)], ]
  
  plot3 <- ggplot(test_start) +  
    geom_ribbon(data = test_start[test_start$type == "fitted",], aes(ymin=conflik_lower, ymax=conflik_upper, x= t, fill=condition), alpha=0.2) +
    geom_point(data = test_start[test_start$type == "measured",], aes(y=y, x=t, col=condition)) + 
    geom_line(data = test_start[test_start$type == "fitted",], aes(y=y, x=t, col=condition)) +
    scale_color_manual(values = c("#fb9092", "#8F8F8F", "#767791")) +
    scale_fill_manual(values = c("#fb9092", "#8F8F8F", "#767791"))+
    #  geom_smooth(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), method="loess") +
    facet_wrap(~peptide) + 
    labs(title = "Peptides - beginning")
  print(plot3)
  
  plot1 <- ggplot(test_middle) +  
    geom_ribbon(data = test_middle[test_middle$type == "fitted",], aes(ymin=conflik_lower, ymax=conflik_upper, x= t, fill=condition), alpha=0.2) +
    geom_point(data = test_middle[test_middle$type == "measured",], aes(y=y, x=t, col=condition)) + 
    geom_line(data = test_middle[test_middle$type == "fitted",], aes(y=y, x=t, col=condition)) +
    scale_color_manual(values = c("#c03d56", "#8F8F8F", "#767791")) +
    scale_fill_manual(values = c("#c03d56", "#8F8F8F", "#767791"))+
    #  geom_smooth(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), method="loess") +
    facet_wrap(~peptide)+ 
    labs(title = "Peptides - middle")
  print(plot1)

  plot2 <- ggplot(test_end) +  
    geom_ribbon(data = test_end[test_end$type == "fitted",], aes(ymin=conflik_lower, ymax=conflik_upper, x= t, fill=condition), alpha=0.2) +
    geom_point(data = test_end[test_end$type == "measured",], aes(y=y, x=t, col=condition)) + 
    geom_line(data = test_end[test_end$type == "fitted",], aes(y=y, x=t, col=condition)) +
    scale_color_manual(values = c("#782d4d", "#8F8F8F", "#767791")) +
    scale_fill_manual(values = c("#782d4d", "#8F8F8F", "#767791"))+
    #  geom_smooth(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), method="loess") +
    facet_wrap(~peptide)+ 
    labs(title = "Peptides - end")
  print(plot2) 
}

minvalues_calc <- function(all_results, small_molecule, sm, control, protein,
                           domain_start = 1, domain_end = 10000) {
  
  peptides <- all_results[[small_molecule]][["peptide_matching"]]$peptide[all_results[[small_molecule]][["peptide_matching"]]$protein == protein 
                                                                          & all_results[[small_molecule]][["peptide_matching"]]$end_match <= domain_end 
                                                                          & all_results[[small_molecule]][["peptide_matching"]]$start_match >= domain_start 
                                                                          ]
  
  proteins_sm <- all_results[[small_molecule]][["python_fit"]] %>% 
    filter(condition == sm, 
           peptide %in% peptides,
           type == "fitted" 
           #        peptide %in% second_domain_proteins$peptide
    ) %>% 
    reshape2::dcast(peptide~t, value.var = "y", fun.aggregate = mean)
  
  proteins_control <- all_results[[small_molecule]][["python_fit"]] %>% 
    filter(condition == control, 
           peptide %in% peptides,
           type == "fitted" 
           #        peptide %in% second_domain_proteins$peptide
    ) %>% 
    reshape2::dcast(peptide~t, value.var = "y", fun.aggregate = mean)
  
  proteins_control %<>%
    `colnames<-`(c("peptide", round(as.numeric(colnames(.[,-1])), 2)))
  
  proteins_sm %<>%
    `colnames<-`(c("peptide", round(as.numeric(colnames(.[,-1])), 2)))
  
  values_control <- data_frame(condition = "control", absolute_correlation = abs(cor((proteins_control[,-1]))[1,]), temperature = as.numeric(names(proteins_control[,-1]))) 
  values_sm <- data_frame(condition = small_molecule, absolute_correlation = abs(cor((proteins_sm[,-1]))[1,]), temperature = as.numeric(names(proteins_sm[,-1]))) 
  
  values <- rbind(values_control, values_sm)
  
  min_value_control <- mean(values_control$temperature[values_control$absolute_correlation<0.2])
  min_value_sm <- mean(values_sm$temperature[values_sm$absolute_correlation<0.2])
  values$condition <- factor(values$condition, levels=c(sm, control))

  print(min_value_control, min_value_sm)
}

extract_melt_temp <- function(all_results, small_molecule, sm, control) {
  proteins_to_investigate <- all_results[[small_molecule]]$python_score %>% 
    filter(!!ensym(control) > 0.1 &
             !!ensym(sm) > 0.1) %>% 
    dplyr::count(protein) %>%
    filter(.$n > 4) %>% 
    .$protein 
  
  melt_temp_sm <- rep(NA, length(proteins_to_investigate))
  sd_temp_sm <- rep(NA, length(proteins_to_investigate))
  
  melt_temp_control <- rep(NA, length(proteins_to_investigate))
  sd_temp_control <- rep(NA, length(proteins_to_investigate))
  
  #i <- 1
  proteins_to_investigate <- intersect( ribosomal_proteins$X...Entry, ribos$glucose$python_score$protein) %>% unique()
  #small_molecule <- "glucose"
  #control <- "control"
  #sm <- "glucose"
  for(i in 1:length(proteins_to_investigate)){
    
    peptides <- all_results[[small_molecule]][["python_score"]] %>%
      filter(protein == proteins_to_investigate[i], 
             !!ensym(control) > 0.1,
             !!ensym(sm) > 0.1) %>%
      .$peptide
    
    proteins_sm <- all_results[[small_molecule]][["python_fit"]] %>% 
      filter(condition == small_molecule, 
             peptide %in% peptides,
             type == "fitted" 
             #        peptide %in% second_domain_proteins$peptide
      ) %>% 
      reshape2::dcast(peptide~t, value.var = "y", fun.aggregate = mean)
    
    proteins_control <- all_results[[small_molecule]][["python_fit"]] %>% 
      filter(condition == "control", 
             peptide %in% peptides,
             type == "fitted" 
             #        peptide %in% second_domain_proteins$peptide
      ) %>% 
      reshape2::dcast(peptide~t, value.var = "y", fun.aggregate = mean)
    
    cor_matrix <- (cor((proteins_sm[,-1])))
    cor_matrix2 <- (cor((proteins_control[,-1])))
    
    tem <- rownames(cor_matrix) %>% as.numeric()
    
    f <- rowMeans(abs(cor_matrix[,1:5]))
    temp <- tem[f<0.2]
    
    f_control <- rowMeans(abs(cor_matrix2[,1:5]))
    temp_control <- tem[f_control<0.2]
    
    if(length(temp)>1) {
      tree <- hclust(dist(temp), method="complete")
      cut_tree <- cutree(tree, h = 3)
      
      #temp[cut_tree == 1]
      melt_temp_sm[i] <- mean(temp[cut_tree == 1])
      sd_temp_sm[i] <- sd(temp[cut_tree == 1])}
    
    else {
      melt_temp_sm[i] <- ifelse(length(temp)>0, temp, NA)
      sd_temp_sm[i] <- 0
    }
    
    if(length(temp_control)>1) {
      tree_control <- hclust(dist(temp_control), method="complete")
      cut_tree_control <- cutree(tree_control, h = 3)
      
      #melt_temp_control[i] <- temp_control[cut_tree_control == 1]
      melt_temp_control[i] <- mean(temp_control[cut_tree_control == 1])
      sd_temp_control[i] <- sd(temp_control[cut_tree_control == 1])}
    
    else {
      melt_temp_sm[i] <- ifelse(length(temp_control)>0, temp_control, NA)
      sd_temp_sm[i] <- 0
    }
    
  }
  
  output_proteins <- data.frame(protein = proteins_to_investigate,
                                melt_control = melt_temp_control, 
                                sd_control = sd_temp_control,
                                melt_sm = melt_temp_sm,
                                sd_sm = sd_temp_sm) %>%
    mutate(dT = melt_sm - melt_control) }



sigmoid_curve_plot <- function(all_results, small_molecule, sm, control, protein,
                           domain_start = 1, domain_end = 10000) {

peptides <- all_results[[small_molecule]][["peptide_matching"]]$peptide[all_results[[small_molecule]][["peptide_matching"]]$protein == protein 
                                                                          & all_results[[small_molecule]][["peptide_matching"]]$end_match <= domain_end 
                                                                          & all_results[[small_molecule]][["peptide_matching"]]$start_match >= domain_start 
                                                                          ]
  
  proteins_sm <- all_results[[small_molecule]][["python_fit"]] %>% 
    filter(condition == sm, 
           peptide %in% peptides,
           type == "fitted" 
           #        peptide %in% second_domain_proteins$peptide
    ) %>% 
    reshape2::dcast(peptide~t, value.var = "y", fun.aggregate = mean)
  
  proteins_control <- all_results[[small_molecule]][["python_fit"]] %>% 
    filter(condition == control, 
           peptide %in% peptides,
           type == "fitted" 
           #        peptide %in% second_domain_proteins$peptide
    ) %>% 
    reshape2::dcast(peptide~t, value.var = "y", fun.aggregate = mean)
  
  proteins_control %<>%
    `colnames<-`(c("peptide", round(as.numeric(colnames(.[,-1])), 2)))
  
  proteins_sm %<>%
    `colnames<-`(c("peptide", round(as.numeric(colnames(.[,-1])), 2)))
  
  values_control <- data_frame(condition = "control", absolute_correlation = (cor((proteins_control[,-1]))[1,]), temperature = as.numeric(names(proteins_control[,-1]))) 
  values_sm <- data_frame(condition = small_molecule, absolute_correlation = (cor((proteins_sm[,-1]))[1,]), temperature = as.numeric(names(proteins_sm[,-1]))) 
  
  values <- rbind(values_control, values_sm)
  
  ggplot(values, aes(x=temperature, y=absolute_correlation, col=condition)) +
    geom_point()
}

plot_ranges <- function(all_results, small_molecule, protein) {
  protein <- "P0A6Y8"
  small_molecule <- "trehalose"
  one_protein_range <- all_results[[small_molecule]]$ranges[all_results[[small_molecule]]$ranges$protein %in% protein,]
  
  sm_color <- c(control = "#8D8D92",
                ATP = "#3C4285", 
                betaine = "#2C938B",
                trehalose = "#A84343",
                glucose = "#6B1818", 
                TMAO_2 = "#6B1818", 
                PK_0.01 = "#6B1818")
  
  p <- ggplot(one_protein_range, aes(y=protein, yend=protein, col = log2(score_weighted))) +
    geom_segment(aes(x=start_range, xend=end_range), size=10) +
    guides(guide = guide_legend(reverse=TRUE)) +
    scale_x_continuous(breaks=seq(0,1000,100)) +
    # scale_color_manual(values=setNames(c("darkblue", "indianred3"), c(FALSE, TRUE)), na.value="gray80", drop = FALSE)+
    theme_classic(base_size=5) +
    scale_color_gradient2(low="white", mid = "gray80", high=sm_color[names(sm_color) == small_molecule], midpoint = log2(1))+
    # xlim(c(0,protein_length)) +
    theme(axis.line.y=element_blank(),
          axis.ticks.y=element_blank()) + 
   geom_vline(xintercept = 383, lty="dashed", col="gray20", lwd=1)
  
  p+theme(axis.text.x = element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          legend.position = "none") +
    ylab("")+
    xlab("")
  
  print(p)

}


protein_domain <- data.frame(protein = c("P0A6Y8", "P0A6Y8"),
                              domains = c("ATP binding", "Substrate binding"), 
                             start_domain = c(1, 384),
                             end_domain = c(383, 603))

p <- ggplot(protein_domain, aes(y=protein, yend=protein,col = domains)) +
  geom_segment(aes(x=start_domain, xend=end_domain), size=10) +
  guides(guide = guide_legend(reverse=TRUE)) +
  scale_x_continuous(breaks=seq(0,1000,100)) +
  theme_classic(base_size=5) +
  geom_label(aes(label=domains, x=(start_domain + 100))) +
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank())# + 
# geom_vline(xintercept = c(one_protein_range$start_range,one_protein_range$end_range ), lty="dashed", col="gray90")

p+theme(axis.text.x = element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none") +
  ylab("")+
  xlab("")

#print(p)



map_cleavage_site <- function(files, fasta_file) {
  Fasta_ecoli <- fasta_file
  proteins_ecoli <- names(Fasta_ecoli) 
  proteins_fasta <- lapply(strsplit(proteins_ecoli, split="\\|"), FUN = function(x) "[[" (x,2)) %>% unlist() %>% #from the FASTA get just the UniprotIDs
    intersect(., files[["python_score"]]$protein) 
  
  #comment: Like this we make sure that we map the proteins that are also in the FASTA and not for example iRT peptides we use for DIA analysis
  
  
  #Prepare empty data frame 
  position_export <- data.frame(matrix(ncol = 11, nrow=0))
  colnames(position_export) <- c("protein","start_match", "end_match", "peptide", "protein_length", "start1", "start2", "end1", "end2")
  
  
  for(i in 1:length(proteins_fasta)) {
    #extract protein sequence from the FASTA file
    protein_sequence <- Fasta_ecoli[[grep(proteins_fasta[i], names(Fasta_ecoli))[1]]] %>%
      paste(., sep="", collapse="") %>% AAString()
    
    prot_string <- Fasta_ecoli[[grep(proteins_fasta[i], names(Fasta_ecoli))[1]]] %>%
      paste(., sep="", collapse="") %>% unlist() %>% as.character()

    #get the peptides from that protein
    peptides <- files[["python_score"]] %>% 
      filter(.$protein == proteins_fasta[i])
    
    peptide_sequences <- peptides$peptide %>% AAStringSet()
    
    #match the peptides with the sequence
    matches <- Biostrings::matchPDict(peptide_sequences, protein_sequence) %>% unlist() 
    
    if(length(ranges(matches)) == length(peptides$peptide)){ #some proteins have peptides that map to more than one sequence, those are excluded in the first round
      protein_length <- length(protein_sequence)
      start_match <- Biostrings::start(matches) #Give out the vector of start positions of a match
      end_match <- Biostrings::end(matches) #Give out the vector of end positions of a match
      
      #comment: if there is a region with both significant peptides, missing peptides... that region will count for both, so the sum of all the region will not be 100% 
      #This can be adjusted if needed
      
      intervals <- data.frame(protein = proteins_fasta[i],
                              start_match, 
                              end_match,
                              peptide = as.character(peptides$peptide),
                              protein_length = protein_length)
      
      intervals %<>%
        group_by(peptide) %>%
        mutate(start1 = substr(prot_string, (start_match-1), (start_match-1)),
               start2 = substr(prot_string, start_match, start_match),
               end1 = substr(prot_string, end_match, end_match),
               end2 = substr(prot_string, (end_match+1), (end_match+1))) %>%
        ungroup()
      
      
      intervals <- intervals[order(intervals$start_match),]
      
      #Order the levels of the peptides based on the position
      intervals$peptide <- factor(intervals$peptide, levels=unique(intervals$peptide))
      position_export <- rbind(position_export, intervals)
      
      
    }
    
    else {
      print(proteins_fasta[i])
    }
    
    
  }
  position_export
  
}


minvalues_easy <- function(all_results, small_molecule, sm, control, proteins_input,
                           domain_start = 1, domain_end = 10000, title) {
  
  proteins <- all_results[[small_molecule]]$python_score %>%
    filter(protein %in% proteins_input) %>%
    filter(control > 0.1) %>%
    dplyr::count(protein) %>%
    filter(n > 5) %>%
    .$protein %>% as.character()
  
  min_value_control <- rep(NA, length(proteins))
  min_value_sm <- rep(NA, length(proteins))
  for(i in 1:length(proteins)) {
    
    peptides <- all_results[[small_molecule]][["peptide_matching"]]$peptide[all_results[[small_molecule]][["peptide_matching"]]$protein == proteins[i]
                                                                            & all_results[[small_molecule]][["peptide_matching"]]$end_match <= domain_end 
                                                                            & all_results[[small_molecule]][["peptide_matching"]]$start_match >= domain_start 
                                                                            ]
    
    proteins_sm <- all_results[[small_molecule]][["python_fit"]] %>% 
      filter(condition == sm, 
             peptide %in% peptides,
             type == "fitted" 
             #        peptide %in% second_domain_proteins$peptide
      ) %>% 
      reshape2::dcast(peptide~t, value.var = "y", fun.aggregate = mean)
    
    proteins_control <- all_results[[small_molecule]][["python_fit"]] %>% 
      filter(condition == control, 
             peptide %in% peptides,
             type == "fitted" 
             #        peptide %in% second_domain_proteins$peptide
      ) %>% 
      reshape2::dcast(peptide~t, value.var = "y", fun.aggregate = mean)
    
    proteins_control %<>%
      `colnames<-`(c("peptide", round(as.numeric(colnames(.[,-1])), 2)))
    
    proteins_sm %<>%
      `colnames<-`(c("peptide", round(as.numeric(colnames(.[,-1])), 2)))
    
    values_control <- data_frame(condition = "control", absolute_correlation = abs(cor((proteins_control[,-1]))[1,]), temperature = as.numeric(names(proteins_control[,-1]))) 
    values_sm <- data_frame(condition = small_molecule, absolute_correlation = abs(cor((proteins_sm[,-1]))[1,]), temperature = as.numeric(names(proteins_sm[,-1]))) 
    
    values <- rbind(values_control, values_sm)
    
    min_value_control[i] <- mean(values_control$temperature[values_control$absolute_correlation<0.2])
    min_value_sm[i] <- mean(values_sm$temperature[values_sm$absolute_correlation<0.2])
    
    
  }
  
  output <- data.frame(proteins, min_value_control, min_value_sm)

}

extract_correlation <- function(all_results, small_molecule, control = "control") {
#  proteins <- "P0A6Y8"
#  small_molecule <- "glucose"
  
  all_prot_peptides <- all_results[[small_molecule]][["python_score"]] %>% 
    filter(control > 0.5)
  
  
  proteins <- plyr::count(all_prot_peptides$protein) %>% 
    filter(freq > 4) %>%
    .$x
  
  
  list_proteins <- vector("list", length(proteins))
  names(list_proteins) <- proteins 
  
  for(i in 1:length(proteins)){
    
    peptides <- all_prot_peptides$peptide[all_prot_peptides$protein == proteins[i]]
    
    proteins_control <- all_results[[small_molecule]][["python_fit"]] %>% 
      filter(condition == "control", 
             peptide %in% peptides,
             type == "fitted" 
             #        peptide %in% second_domain_proteins$peptide
      ) %>% 
      reshape2::dcast(peptide~t, value.var = "y", fun.aggregate = mean)
    
    list_proteins[[proteins[i]]] <- cor((proteins_control[,-1]))[,1]
  }
  list_proteins
}
