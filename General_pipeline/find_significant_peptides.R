#### Functions ####
find_significant_peptides <- function(area_calculated, clusters, sm, control = "Control") {
    #aggregated_area <- aggregate_effects(area_calculated)
  ## Create a wide format just for the peaks

  # Calculate the distance between the curves when the confidence intervals do not overlap
  area_calculated <- area_calculation_2(all_results,small_molecule = sm, cutoff = 0.5, control = control)
  area_calculated <- asign_peak_apex(area_calculated)
  area_calculated$define_peak <- ifelse(area_calculated$define_new_peak == "not_aggregation", "stabilisation", area_calculated$define_peak)
  
  area_wide <- area_calculated %>%
    ungroup() %>%
    dplyr::mutate(t_apex = ifelse(apex, t, NA) ) %>%
    dplyr::select(peptide, peak_sum, define_peak, goal, t_apex) %>%
    reshape2::dcast(peptide + goal + t_apex ~define_peak, value.var = "peak_sum", fun.aggregate = mean) 
  
  colnames(area_wide)[1] <- "Peptide"
  
  # Join the data with clusters
  area_wide %<>% plyr::join(., clusters[[sm]]$data[clusters[[sm]]$data$Condition == sm & clusters[[sm]]$data$is_flip == FALSE,], by="Peptide" )
  
  # Based on cluster shape, determine whether it is stabilisation or destabilisation
  area_wide %<>% 
    mutate(effect = ifelse( Cluster %in% c(1,3) & stabilisation > 0,"stabilisation",
                           ifelse( Cluster %in% c(1,3) & stabilisation < 0,"destabilisation",
                                  ifelse( Cluster %in% c(2,4) & stabilisation > 0,"destabilisation",
                                         ifelse( Cluster %in% c(2,4) & stabilisation < 0,"stabilisation",NA))))) %>%
    
    mutate(effect_agg = ifelse(Cluster %in% c(1,3) & aggregation > 0, "aggregation_start",
                               ifelse ( Cluster %in% c(1,3) & aggregation < 0, "aggregation_stop",
                                        ifelse ( Cluster %in% c(2,4) & aggregation > 0, "aggregation_stop",
                                                 ifelse ( Cluster %in% c(2,4) & aggregation < 0, "aggregation_start",NA)))))

  
  
  ### Different assessment for aggregation / binding / stabilisation
  ### Make a data frame of all peptides that are included in clustering and calculations
  all_peptides <- area_wide %>%
    filter(Cluster %in% c(1:4)) %>%
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
    filter(Cluster %in% c(1:2)) %>%
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
    filter(Cluster %in% c(3,4)) %>%
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
  
  gof_t <-  all_results[[sm]]$python_score %>%
    unique()

    #group_by(peptide) %>%
    #mutate(gof = sum(!!ensym(control), !!ensym(sm))) %>%
    #dplyr::select(peptide, gof) %>%
    #filter(peptide %in% all_peptides$Peptide) %>%
    #`colnames<-`(c("Peptide", "gof"))
  gof_t$gof <- rowSums(gof_t[,2:3])
  gof <- gof_t %>% dplyr::select(peptide,gof) %>%
    filter(peptide %in% all_peptides$Peptide) %>%
    `colnames<-`(c("Peptide","gof"))

  all_peptides %<>%
    plyr::join(.,gof, by="Peptide" ) 
  
  output <- list()
  output$scores <- all_peptides
  output$long <- area_calculated
  return(output) }
