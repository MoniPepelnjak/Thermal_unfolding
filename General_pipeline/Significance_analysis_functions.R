### Function for significance analysis, combining the clusters

find_significant_peptides <- function(all_results, clusters, sm, control = "Control") {
  
  # Calculate the distance between the curves when the confidence intervals do not overlap
  area_calculated <- area_calculation_2(all_results,small_molecule = sm, cutoff = 0.5, control = control)
  area_calculated <- asign_peak_apex(area_calculated)
  area_calculated$define_peak <- ifelse(area_calculated$define_new_peak == "not_aggregation", "stabilisation", area_calculated$define_peak)
  
  #aggregated_area <- aggregate_effects(area_calculated)
  ## Create a wide format just for the peaks
  area_wide <- area_calculated %>%
    ungroup() %>%
    mutate(t_apex = ifelse(apex, t, NA) ) %>%
    dplyr::select(peptide, peak_sum, define_peak, goal, t_apex) %>%
    reshape2::dcast(peptide + goal + t_apex ~define_peak, value.var = "peak_sum", fun.aggregate = mean) 

  colnames(area_wide)[1] <- "Peptide"
  
  # Join the data with clusters
  area_wide %<>% plyr::join(., clusters[[sm]]$data[clusters[[sm]]$data$Condition == sm & clusters[[sm]]$data$is_flip == FALSE,], by="Peptide" )

  # Based on cluster shape, determine whether it is stabilisation or destabilisation
  area_wide %<>% 
    mutate(effect = ifelse(Cluster %in% c(1, 2, 3) & stabilisation > 0, "stabilisation", 
                           ifelse ( Cluster %in% c(1, 2, 3) & stabilisation < 0, "destabilisation", 
                                    ifelse ( Cluster %in% c(4, 5, 6) & stabilisation > 0, "destabilisation", 
                                             ifelse ( Cluster %in% c(4, 5, 6) & stabilisation < 0, "stabilisation",
                                                      ifelse ( Cluster %in% c(7) & stabilisation > 0, "stabilisation",
                                                               ifelse ( Cluster %in% c(7) & stabilisation < 0, "destabilisation",
                                                                        ifelse ( Cluster %in% c(8) & stabilisation < 0, "stabilisation",
                                                                                 ifelse ( Cluster %in% c(8) & stabilisation > 0, "destabilisation", NA))))))))) %>%
    mutate(effect_agg = ifelse(Cluster %in% c(1, 2, 3) & aggregation > 0, "aggregation_start", 
                           ifelse ( Cluster %in% c(1, 2, 3) & aggregation < 0, "aggregation_stop", 
                                    ifelse ( Cluster %in% c(4, 5, 6) & aggregation > 0, "aggregation_stop", 
                                             ifelse ( Cluster %in% c(4, 5, 6) & aggregation < 0, "aggregation_start",
                                                      ifelse ( Cluster %in% c(7) & aggregation > 0, "aggregation_start",
                                                               ifelse ( Cluster %in% c(7) & aggregation < 0, "aggregation_stop",
                                                                        ifelse ( Cluster %in% c(8) & aggregation < 0, "aggregation_start",
                                                                                 ifelse ( Cluster %in% c(8) & aggregation > 0, "aggregation_stop", NA)))))))))
  
  ### Different assessment for aggregation / binding / stabilisation
  ### Make a data frame of all peptides that are included in clustering and calculations
  all_peptides <- area_wide %>%
    filter(Cluster %in% c(1:8)) %>%
    dplyr::select(Peptide) %>%
    unique()
  
  ### Aggregation: 
  aagregation_df <- area_wide %>%
    filter(!is.na(effect_agg)) %>%
    dplyr::select(-t_apex) %>%
    unique() %>%
    group_by(Peptide) %>%
    mutate(aggregation = ifelse(effect_agg == "aggregation_start", abs(aggregation), -abs(aggregation))) %>%
    summarise(Aggregation = aggregation[abs(aggregation) == max(abs(aggregation))]) %>%
    dplyr::select(Peptide, Aggregation)

  all_peptides %<>%
    plyr::join(., aagregation_df, by="Peptide")
  
  area_wide_sigmoid <- area_wide %>%
    filter(Cluster %in% c(1:6)) %>%
    filter(!is.na(stabilisation),
           !is.na(t_apex)) %>% 
    mutate(effect_value = ifelse(effect == "stabilisation", abs(stabilisation), -abs(stabilisation)) ) %>%
    unique() %>%
    group_by(Peptide) %>%
    summarise(max_value = effect_value[abs(effect_value) == max(abs(effect_value))]) %>%
    dplyr::select(Peptide, max_value) %>%
    `colnames<-`(c("Peptide", "Stabilisation"))
  
  test_assumptions <- all_results[[sm]]$python_fit[all_results[[sm]]$python_fit$peptide %in% area_wide$Peptide[area_wide$Cluster %in% c(7,8)],] %>%
    .[.$type == "fitted" & .$condition == sm,] %>%
    group_by(peptide ) %>% 
    summarise( t_extreme = ifelse( peptide %in% area_wide$Peptide[area_wide$Cluster %in% c(7)], t[y == min(y)], t[y == max(y)])) %>%
    `colnames<-`(c("Peptide", "t_extreme")) %>%
    unique()
  
  area_wide_check <-area_wide  %>%
    filter(Cluster %in% c(7,8)) %>%
    plyr::join(., test_assumptions, by="Peptide") %>%
    mutate(is_lower = t_apex < t_extreme) %>%
    group_by(Peptide) %>%
    filter(!is.na(stabilisation),
           !is.na(t_apex)) %>%
    unique() %>%
    mutate(real_effect = ifelse(is_lower, effect, 
                                ifelse(!is_lower & effect == "destabilisation", "stabilisation", "destabilisation"))  ) %>%
    mutate(effect_value = ifelse(real_effect == "stabilisation", abs(stabilisation), -abs(stabilisation)) )
  
  area_wide_check %<>%
    group_by(Peptide) %>%
    summarise(aggregated_effect = sum(effect_value, na.rm=TRUE), 
              max_value = effect_value[abs(effect_value) == max(abs(effect_value))]) %>%
    dplyr::select(Peptide, max_value) %>%
    `colnames<-`(c("Peptide", "Stabilisation"))
    
  area_stabilisation <- rbind(area_wide_check, area_wide_sigmoid)
  
  all_peptides %<>% plyr::join(., area_stabilisation, by="Peptide")
  
  area_binding <- area_wide %>%
    filter(!is.na(binding)) %>%
    dplyr::select(Peptide, binding) %>%
    unique() %>%
    mutate(binding = abs(binding))
  
  all_peptides %<>% plyr::join(., area_binding, by="Peptide")
  all_peptides[is.na(all_peptides)] <- 0
  
  gof <-  all_results[[sm]]$python_score %>%
          unique() %>%
    group_by(peptide) %>%
          mutate(gof = sum(!!ensym(control), !!ensym(sm))) %>%
          dplyr::select(peptide, gof) %>%
    filter(peptide %in% all_peptides$Peptide) %>%
    `colnames<-`(c("Peptide", "gof"))
  
  all_peptides %<>%
    plyr::join(.,gof, by="Peptide" ) 
  
  output <- list()
  output$scores <- all_peptides
  output$long <- area_calculated
  return(output) }

################# Peptide matching #############
library(seqinr)

fasta_ecoli <- read.fasta("/Users/moni/Documents/Phd/databases/official_161018_uniprot_ecoli_K12.fasta", "AA")

peptide_fasta_matching <- function(peptide_list, fasta_ecoli){
  
  protein_sequence <- fasta_ecoli %>%
    lapply(., FUN = function(x) paste(x, sep="", collapse=""))
  proteins_ecoli <- names(protein_sequence)
  names(protein_sequence) <-  lapply(strsplit(proteins_ecoli, split="\\|"), FUN = function(x) "[[" (x,2))
  
  df <- data.frame(matrix(unlist(protein_sequence), nrow=length(protein_sequence), byrow=TRUE))
  df$proten <- names(protein_sequence)
  colnames(df) <- c("Sequence", "Protein")
  
  matches <- sapply(peptide_list, FUN = function(x) names(protein_sequence)[grep(x, protein_sequence)])

  is_unique <- sapply(matches, length)
  matches_2 <- matches[is_unique == 1]
  df_matched <- data.frame(matrix(unlist(matches_2), nrow=length(matches_2), byrow=TRUE))
  df_matched$Peptide <- names(matches_2)
  colnames(df_matched) <- c("Protein", "Peptide")
  
  joined_df <- plyr::join(df_matched, df, by="Protein")
  
  joined_df %<>%
    mutate(start = stringr::str_locate(Sequence, Peptide)) %>%
    mutate(end = start + nchar(Peptide) - 1) %>%
    mutate(prot_len = nchar(Sequence)) %>%
    dplyr::select(-Sequence) %>%
    group_by(Peptide) %>%
    mutate(start = start[1], 
           end = end[1])
  
  
  return(joined_df)
  
}

##################### AA level score calculation ########################

score_aggregation <- function(significant_all, sm, matched_df, pep_level_quantile = 0.75, prot_level_quantile = 0.75) {
  all_peptides <- significant_all[[sm]]$scores %>%
    unique() %>%
    plyr::join(., matched_df, by="Peptide") %>%
    na.omit() %>%
    mutate(gof = ifelse(gof < 0, 0.0001, gof))  %>%
    data.frame()
  
  scoring_full_proteome <- all_peptides %>%
    group_by(Peptide) %>%
    mutate(start = unlist(start)[1], 
           end = unlist(end)[1]) %>%
    mutate(positions = list(start[1]:end[1])) %>%
    tidyr::unnest(positions) %>%
    #now mutate a column which is unique for each protein
    mutate(protein_position_id = paste(Protein, positions, sep = "_")) %>%
    #groupby the protein position id to extract the mean value per position
    group_by(protein_position_id) %>%
    mutate(AA_score = as.numeric(wtd.quantile(x = Stabilisation, q=pep_level_quantile, weight = gof,  na.rm=TRUE)), 
           median_gof = median(gof),
           AA_aggregation = wtd.quantile(x=Aggregation, weight=gof, q=pep_level_quantile),
           AA_binding = wtd.quantile(x=binding, weight=gof, q=pep_level_quantile)
           )
  
  #scoring_full_proteome$AA_score <- as.numeric(scoring_full_proteome$AA_score)
  scoring_full_proteome$median_gof <- as.integer(scoring_full_proteome$median_gof)
  scoring_full_proteome %<>% unique() 
  
  #prot_level_quantile <- 0.75
  scoring_protein_level <- scoring_full_proteome %>%
    group_by(Protein) %>%
    mutate(n_peptide = length(unique(Peptide))) %>%
    dplyr::select(Protein, positions, AA_score,median_gof, AA_aggregation, AA_binding, n_peptide, prot_len ) %>%
    unique() %>%
    group_by(Protein, n_peptide, prot_len) %>%
    summarise(Protein_stabilisation = reldist::wtd.quantile(x=as.numeric(AA_score), weight=median_gof, q=prot_level_quantile),
              Percentage_stabilised = sum(AA_score>0)/n(),
              Protein_aggregation = reldist::wtd.quantile(x=AA_aggregation, weight=median_gof, q=prot_level_quantile),
              Protein_binding = reldist::wtd.quantile(x=AA_binding, weight=median_gof, q=prot_level_quantile),
              Percentage_aggregated = sum(abs(AA_aggregation)>0)/n(),
              Percentage_binding = sum(AA_binding>0)/n(),
              Coverage = n()/prot_len
    )

  
  scoring_protein_level %<>% unique() %>%
    ungroup() %>%
    mutate(prot_rank = rank(-Protein_stabilisation))
  
  significant_all[[sm]]$AA_level <- scoring_full_proteome
  significant_all[[sm]]$Protein_level <- scoring_protein_level
  significant_all[[sm]]$scores <- all_peptides
  return(significant_all)
  
}

score_aggregation_both_quantiles <- function(significant_all, sm, pep_level_quantile = 0.75, prot_level_quantile = 0.75) {
  all_peptides <- significant_all[[sm]]$scores %>%
    unique() %>%
    na.omit() %>%
    mutate(gof = ifelse(gof < 0, 0.0001, gof))
  
  scoring_full_proteome <- all_peptides %>%
    group_by(Peptide) %>%
    mutate(start = start[1], 
           end = end[1]) %>%
    mutate(positions = list(start:end)) %>%
    tidyr::unnest(positions) %>%
    #now mutate a column which is unique for each protein
    mutate(protein_position_id = paste(Protein, positions, sep = "_")) %>%
    #groupby the protein position id to extract the mean value per position
    group_by(protein_position_id) %>%
    mutate(mean_score = weighted.mean(Stabilisation, gof)) %>%
    mutate(AA_score = ifelse(mean_score >=0, wtd.quantile(x=Stabilisation, weight=gof, q=pep_level_quantile), wtd.quantile(x=Stabilisation, weight=gof, q=(1-pep_level_quantile))) , 
           median_gof = median(gof),
           AA_aggregation = wtd.quantile(x=Aggregation, weight=gof, q=pep_level_quantile),
           AA_binding = wtd.quantile(x=binding, weight=gof, q=pep_level_quantile)) %>%
    dplyr::select(-mean_score)
  
  scoring_full_proteome %<>% unique()
  
  scoring_protein_level <- scoring_full_proteome %>%
    group_by(Protein) %>%
    mutate(n_peptide = length(unique(Peptide))) %>%
    dplyr::select(Protein, positions, AA_score,median_gof, AA_aggregation, AA_binding, n_peptide, prot_len ) %>%
    unique() %>%
    group_by(Protein) %>%
    mutate(mean_p_score = weighted.mean(AA_score, median_gof)) %>%
    group_by(Protein, n_peptide, prot_len,mean_p_score ) %>%
    summarise(Protein_stabilisation = ifelse(mean_p_score >= 0, wtd.quantile(x=AA_score, weight=median_gof,q= prot_level_quantile), wtd.quantile(x=AA_score, weight=median_gof, q=(1-prot_level_quantile))),
              Percentage_stabilised = sum(AA_score>0)/n(),
              Protein_aggregation = wtd.quantile(x=AA_aggregation, weight=median_gof, q=prot_level_quantile),
              Protein_binding = wtd.quantile(x=AA_binding, weight=median_gof, q=prot_level_quantile),
              Percentage_aggregated = sum(abs(AA_aggregation)>0)/n(),
              Percentage_binding = sum(AA_binding>0)/n(),
              Coverage = n()/prot_len
    ) %>%
    dplyr::select(-mean_p_score)
  
  scoring_protein_level %<>% unique() %>%
    ungroup() %>%
    mutate(prot_rank = rank(-Protein_stabilisation))
  
  significant_all[[sm]]$AA_level <- scoring_full_proteome
  significant_all[[sm]]$Protein_level <- scoring_protein_level
  significant_all[[sm]]$scores <- all_peptides
  return(significant_all)
  
}


