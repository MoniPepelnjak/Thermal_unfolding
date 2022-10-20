#### Function ####
score_aggregation <- function(significant_all, sm, matched_df, pep_level_quantile = 0.75, prot_level_quantile = 0.75) {
  
  all_peptides <- significant_all[[sm]]$scores %>%
    unique() %>%
    plyr::join(., matched_df, by="Peptide") %>%
    na.omit() %>%
    mutate(gof = ifelse(gof < 0, 0.0001, gof))
  pep_level_quantile <- 0.75
  
  scoring_full_proteome <- all_peptides %>%
    group_by(Peptide) %>%
    mutate(start = start[1], 
           end = end[1]) %>%
    mutate(positions = list(start[1]:end[1])) %>%
    tidyr::unnest(positions) %>%
    #now mutate a column which is unique for each protein
    mutate(protein_position_id = paste(Protein, positions, sep = "_")) %>%
    #groupby the protein position id to extract the mean value per position 
    group_by(protein_position_id) %>%
    mutate(mean_score = wtd.mean(x = Stabilisation,weight =  gof)) %>%
    mutate(AA_score = ifelse(mean_score >=0, wtd.quantile(x = Stabilisation, weight = gof, q =  pep_level_quantile), wtd.quantile(x = Stabilisation, weight = gof, q = (1-pep_level_quantile))) , 
           median_gof = median(gof),
           AA_aggregation = wtd.quantile(x = Aggregation, weight = gof, q = pep_level_quantile),
           AA_binding = wtd.quantile(x = binding, weight = gof,q =  pep_level_quantile))
  #head(scoring_full_proteome)  
  scoring_full_proteome %<>% unique()
  
  scoring_protein_level <- scoring_full_proteome %>%
    group_by(Protein) %>%
    mutate(n_peptide = length(unique(Peptide))) %>%
    dplyr::select(Protein, positions, AA_score,median_gof, AA_aggregation, AA_binding, n_peptide, prot_len ) %>%
    unique() %>%
    group_by(Protein) %>%
    mutate(mean_p_score = wtd.mean(x = AA_score, weight = median_gof)) %>%
    group_by(Protein, n_peptide, prot_len) %>%
    summarise(Protein_stabilisation = ifelse(mean_p_score >= 0, wtd.quantile(x = AA_score, weight = median_gof, q = prot_level_quantile), wtd.quantile(x = AA_score,weight =  median_gof, q = (1-prot_level_quantile))),
              Percentage_stabilised = sum(AA_score>0)/n(),
              Protein_aggregation = wtd.quantile(x = AA_aggregation, weight = median_gof, q = prot_level_quantile),
              Protein_binding = wtd.quantile(x = AA_binding, weight = median_gof, q = prot_level_quantile),
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
