## Clustering functions
## Author: Monika Pepelnjak
## Date 22.01.2021

## Fuzzy k-means clustering

# Peptide level analysis of aggregation
### Clustering script for one condition:

# Goal of the script:
  # Take one of the fits for each peptide (control/sm)
  # Sample points (not to take all 100)
  # Perform fuzzy clustering on peptides
  # Define clusters based on their characteristics
  # Return df with fixed cluster numbers 

# Parameters specified and their meaning
  # all results - List of result files, fits after pyhthon analysis of scaled values
  # sm - name of the experiment (for example small molecule used)
  # condition analysed - whether control or small molecule condition within the sm is analysed -> it is always either the same as sm or "control" (also default)
  # k = number of clusters. Minimal number is 8


fuzzy_clusters_two_conditions <- function(all_results, sm, control="control", k = 8, show_clusters = TRUE, cutoff= 1.65) {
  set.seed(123) 

  high_quality <- all_results[[sm]]$python_score %>%   # extract python "score" data
    filter(!!ensym(sm) > 0.5 &
             !!ensym(control) > 0.5) %>%      # only take the ones that have a good enough fit
    dplyr::select(peptide, trypticity) %>%
    `colnames<-`(c("Peptide", "trypticity"))

  # Select the fitted peptides for control
  res <- all_results[[sm]]$python_fit[all_results[[sm]]$python_fit$condition %in% c(sm, control) & all_results[[sm]]$python_fit$type == "fitted" & all_results[[sm]]$python_fit$peptide %in% high_quality$Peptide,] %>%
     mutate(peptide = paste(peptide, condition, sep="%%%"))
  
  # only select few points from the fit
  temp <- res$t %>% sort() %>% unique() %>% .[seq(from = 1, to = length(.), by=5 )]
  
  # make a matrix for clustering
  res %<>%
    filter( t %in% temp) %>%
    reshape2::dcast(peptide ~ t, value.var = "y", fun.aggregate = mean)

  rownames(res) <- res$peptide
  res <- res[,-1]
  
  # Perform fuzzy k-means clustering

  fuzzy <- fclust::FKM.ent(res, k=k, stand = 1, ent=2)
  
  fuzzy_plot <- fuzzy$H %>%
      reshape2::melt() %>%
      as.data.frame()
  
  # Take the median values of clusters and describe their characteristics:
    # increasing trend (up_down), whether it has two step profile ("is_agg")
  median_clusters <- fuzzy_plot %>%
    `colnames<-`(c("cluster", "temperature", "median_value")) %>%
    group_by(cluster) %>%
    mutate(diff = median_value[temperature == 37] - median_value[temperature == 76]) %>% 
    ungroup() %>%
    mutate(is_agg = ifelse(abs(diff) > max(abs(diff))/cutoff , FALSE, TRUE)) %>%
   # mutate(is_agg = ifelse(abs(diff) > 0.75 , FALSE, TRUE)) %>%
    group_by(cluster) %>%
    mutate(up_down = ifelse(diff>0, "down", "up") ) %>% 
    mutate(interesting_temp = ifelse(temperature > 45 & temperature < 65, TRUE, FALSE)) %>%
    group_by(cluster, interesting_temp) %>%
    mutate(agg_parameter = ifelse(sign(median_value[abs(median_value) == max(abs(median_value))]) > 0, TRUE, FALSE)) %>%
    ungroup() %>%
    group_by(cluster) %>%
    mutate(max_diff = (max(median_value)+2) - (min(median_value)+2)) %>%
    mutate(is_seven = ifelse(any(interesting_temp == TRUE & agg_parameter == FALSE ), TRUE, FALSE) )

  median_clusters %<>%
    group_by(cluster) %>%
    mutate(sum_median = sum(median_value))
  
  cluster_order <- median_clusters %>%
    ungroup() %>%
    # sum median - easy way to get the area under the curve --> differentiate between thermo sensitive/stable proteins
    dplyr::select(cluster, is_agg, up_down, sum_median, is_seven, max_diff,diff) %>%
    unique()
  
  cluster_order %<>%
    # Based on sum-median --> group clusters in low-stability, max-stability or medium stability.
    mutate(up_down = ifelse(is_agg == FALSE, up_down, "up")) %>%
    group_by(is_agg, up_down) %>%
    mutate(value_order = ifelse(sum_median == min(sum_median), "min_value",
                                ifelse(sum_median == max(sum_median), "max_value", "medium"))) 
  
  # Define cluster numbers based on previously defined characteristics: up_down, is_agg, value_order
  # If Number of clusters is 8, for most datasets we will get one cluster for each number. If number is higher, 
  # we get more clusters in the mean range
  
  # Whole goal of this: Always get the same number for clusters with similar profile
  cluster_order %<>%
    ungroup() %>%
    mutate(manual_order = ifelse(max_diff < 0.8, "9",
                                 ifelse(
                                    is_agg == FALSE & 
                                   up_down == "down" &
                                   value_order == "min_value", "1", 
                                 
                                 ifelse(is_agg == FALSE & 
                                          up_down == "down" &
                                          value_order == "medium", "2", 
                                        
                                        ifelse(is_agg == FALSE & 
                                                 up_down == "down" &
                                                 value_order == "max_value", "3",
                                               
                                               ifelse(is_agg == FALSE & 
                                                        up_down == "up" &
                                                        value_order == "max_value", "4",
                                                      
                                                      ifelse(is_agg == FALSE & 
                                                               up_down == "up" &
                                                               value_order == "medium", "5",
                                                             
                                                             ifelse(is_agg == FALSE & 
                                                                      up_down == "up" &
                                                                      value_order == "min_value", "6",
                                                                    
                                                                    ifelse(is_agg == TRUE & 
                                                                             value_order == "max_value" |
                                                                             is_seven, "7", "8") ) ) ) ) ) ) ) )
  
  
  cluster_order %<>% 
    dplyr::select(cluster, manual_order)
  
  output <- list()
  output$full_set <- fuzzy$U  %>%
      data.frame() %>%
      mutate(Peptide = rownames(.)) %>%
#      mutate(cluster = paste("Clus", as.character(.$Cluster), sep=" ")) %>%
      mutate(Condition = strsplit(Peptide, "%%%") %>% lapply(., function(x) '[['(x,2)) %>% unlist() ) %>%
      mutate(Peptide = strsplit(Peptide, "%%%") %>% lapply(., function(x) '[['(x,1)) %>% unlist() ) %>%
      reshape2::melt(., id.var=c("Peptide", "Condition")) %>%
      mutate(Peptide = strsplit(Peptide, "%%%") %>% lapply(., function(x) '[['(x,1)) %>% unlist() ) %>%
      mutate(variable = gsub(pattern = "\\.", replacement = " ",  x = variable)) %>%
      `colnames<-`(c("Peptide", "Condition", "cluster", "Membership")) %>%
      plyr::join(., cluster_order, by="cluster") %>%
      plyr::join(., high_quality, by="Peptide")
  
    
  output$manual_order <- cluster_order
  
  rm(fuzzy)
  rm(cluster_order)
  rm(fuzzy_plot)
  rm(median_clusters)
  rm(fuzzy_clusters)

  return(output)

}

cluster_analysis_all <- function(LIP, sm_all, k=15, diff_cutoff=0.85, control="Control") {
  
  Cluster_output <- data.frame(matrix(ncol = 9, nrow=0))
  colnames(Cluster_output) <- c("Peptide","trypticity","Control_shape","Sm_shape","Is_changing","Extra_filter","Change", "Experiment")
  
  k_num <- k
  for(sm in sm_all){
  cluster_LIP_all <- fuzzy_clusters_two_conditions(LIP, sm, control, k=k_num, show_clusters = FALSE)
  Cluster_LIP_full <- cluster_LIP_all$full_set
  rm(cluster_LIP_all)
  Cluster_LIP_full$description <- ifelse(Cluster_LIP_full$manual_order %in% c(1, 2, 3), "down", 
                                         ifelse(Cluster_LIP_full$manual_order %in% c(4, 5, 6), "up", "two-step"))
  
  Cluster_LIP_full %<>%
    group_by(Peptide, Condition, description, trypticity) %>%
    dplyr::summarise(Sum_membership = sum(Membership))
  
  Cluster_LIP_full %<>%
    ungroup( ) %>%
    group_by(Peptide, description) %>%
    mutate(diff = Sum_membership[Condition == sm] - Sum_membership[Condition == control])
  
  Cluster_diff <- Cluster_LIP_full %>%
    group_by(Peptide) %>%
    mutate(Is_changing = ifelse(max(abs(diff)) > diff_cutoff, TRUE, FALSE )) %>%
    mutate(Control_shape = description[ Sum_membership == max(Sum_membership[Condition == control]) & Condition == control ]) %>%
    mutate(Sm_shape = description[ Sum_membership == max(Sum_membership[Condition == sm]) & Condition == sm ]) %>%
    mutate(Extra_filter = ifelse(Sm_shape != Control_shape & Is_changing == FALSE, TRUE, FALSE))
  
  Cluster_diff %<>%
    dplyr::select(Peptide, trypticity, Control_shape, Sm_shape, Is_changing, Extra_filter) %>%
    unique()
  
  Cluster_diff %<>%
    mutate(Change = ifelse(Is_changing == FALSE, "No change", 
                           ifelse((Control_shape == "up" & Sm_shape == "down") |  (Control_shape == "down" & Sm_shape == "up"), "Flip", 
                                  ifelse(Control_shape %in% c("up", "down") & Sm_shape == "two-step", "Two step start", 
                                         ifelse(Sm_shape %in% c("up", "down") & Control_shape == "two-step", "Two step stop", NA)))))
  
  
  Cluster_diff$Experiment <- sm

  Cluster_output <- rbind(Cluster_output, Cluster_diff)
  
  }
  

  return(Cluster_output)
}



 