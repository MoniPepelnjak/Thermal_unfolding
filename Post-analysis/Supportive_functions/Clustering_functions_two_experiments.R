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


fuzzy_clusters_two_experiments <- function(LIP, sm_list, control="control", k = 8, show_clusters = TRUE) {

  set.seed(123) 
#  all_results <- list()
#  all_results$betaine <- LIP$betaine
#  sm <- "betaine"
#  control <- "control"
#  k <- 10
  
  # only select few points from the fit
  temp <- intersect(LIP[[sm_list[1]]]$python_fit$t[LIP[[sm_list[1]]]$python_fit$type == "fitted"], LIP[[sm_list[2]]]$python_fit$t[LIP[[sm_list[2]]]$python_fit$type == "fitted"])   %>% sort() %>% unique() %>% .[seq(from = 1, to = length(.), by=5 )]
  
  #res_output <- data.frame(matrix(ncol=10, nrow=0))
  #colnames(res_output) <- colnames(LIP[[sm_list[1]]]$python_fit)
#  sm_list <- c("trehalose", "betaine")
  res_output <- NULL
  for(sm in sm_list){
    
    high_quality <- LIP[[sm]]$python_score %>%   # extract python "score" data
      filter(!!ensym(control) > 0.5) %>%      # only take the ones that have a good enough fit
      dplyr::select(peptide, trypticity, protein) %>%
      `colnames<-`(c("Peptide", "trypticity", "Protein"))
    
    # Select the fitted peptides for control
    res <- LIP[[sm]]$python_fit %>% 
      filter(condition %in% c(control), 
             type == "fitted", 
             peptide %in% high_quality$Peptide) %>%
      mutate(peptide = paste(peptide, sm, sep="%%%"))

        res_output <- rbind(res_output, res)
        
        print(sm)
  
  }

  # make a matrix for clustering
  res_output %<>%
    filter( t %in% temp) %>%
    reshape2::dcast(peptide ~ t, value.var = "y", fun.aggregate = function(x) mean(x))

  rownames(res_output) <- res_output$peptide
  res_output <- res_output[,-1]
  
  # Perform fuzzy k-means clustering
  
#  k <- 15

  fuzzy <- fclust::FKM.ent(res_output, k=k, stand=1, ent=1)
  
  fuzzy_plot <- fuzzy$H %>%
      reshape2::melt() %>%
      as.data.frame()
  
  # all_plot <- fuzzy$X %>%
  #   reshape2::melt() %>%
  #   `colnames<-`(c("Peptide_con", "temp", "y"))
  # 
  # clus <- fuzzy$clus %>%
  #   as.data.frame() %>%
  #   mutate(Peptide_con = rownames(.))
  # 
  # plot_all <- plyr::join(all_plot,clus, by="Peptide_con")
  # 
  # ggplot(plot_all, aes(x=temp, y=y, col=Cluster, alpha=`Membership degree`, group = Peptide_con)) +
  #   geom_line() +
  #   facet_wrap(~Cluster)
  
  # Take the median values of clusters and describe their characteristics:
    # increasing trend (up_down), whether it has two step profile ("is_agg")
  median_clusters <- fuzzy_plot %>%
    `colnames<-`(c("cluster", "temperature", "median_value")) %>%
    group_by(cluster) %>%
    mutate(diff = median_value[temperature == 37] - median_value[temperature == 76]) %>% 
    mutate(is_agg = ifelse(abs(diff)> 0.75 , FALSE, TRUE)) %>%
    mutate(up_down = ifelse(diff>0, "down", "up") ) 
  
  median_clusters %<>%
    group_by(cluster) %>%
    mutate(sum_median = sum(median_value))
  
  cluster_order <- median_clusters %>%
    ungroup() %>%
    # sum median - easy way to get the area under the curve --> differentiate between thermo sensitive/stable proteins
    dplyr::select(cluster, is_agg, up_down, sum_median) %>%
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
    mutate(manual_order = ifelse(is_agg == FALSE & 
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
                                                                             value_order == "max_value", "7", 8) ) ) ) ) ) ) )
  
  
  cluster_order %<>% 
    dplyr::select(cluster, manual_order)
  
  median_clusters <- fuzzy_plot %>%
    `colnames<-`(c("cluster", "temperature", "median_value")) %>%
    group_by(cluster) %>%
    mutate(diff = median_value[temperature == 37] - median_value[temperature == 76]) %>%
    mutate(is_agg = ifelse(abs(diff)>0.75, FALSE, TRUE)) %>%
    mutate(up_down = ifelse(diff>0, "down", "up") )  %>%
    plyr::join(., cluster_order, by="cluster")
  
  
  fuzzy_clusters <- fuzzy$clus %>% 
    as.data.frame() %>%
    mutate(Peptide = rownames(.)) %>%
    mutate(cluster = paste("Clus", as.character(.$Cluster), sep=" ")) %>%
    mutate(Condition = strsplit(Peptide, "%%%") %>% lapply(., function(x) '[['(x,2)) %>% unlist() ) %>%
    mutate(Peptide = strsplit(Peptide, "%%%") %>% lapply(., function(x) '[['(x,1)) %>% unlist() ) 
  
  fuzzy_clusters %<>%
    plyr::join(., cluster_order, by="cluster")
  
  fuzzy_clusters %<>%
    dplyr::select(-Cluster, -cluster) %>%
    `colnames<-`(c("Membership_degree", "Peptide", "Condition", "Cluster")) %>%
    plyr::join(., high_quality, by="Peptide")
  
  output <- list()
#  show_clusters=TRUE
  if(show_clusters){
    cluster_plot <- ggplot(median_clusters, aes(x=temperature, y=median_value, col=cluster)) +
      geom_line() +
      facet_wrap(~manual_order) +
      theme_minimal() +
      theme(legend.position = "none") 
    
    output$plot <- cluster_plot
    
  }
  

  output$data <- fuzzy_clusters
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
  
  return(output)

}

cluster_analysis_experiments <- function(LIP, sm_all, k=15) {
  #sm_all <- c("glucose", "betaine")
  #k=15
  
  Cluster_output <- data.frame(matrix(ncol = 9, nrow=0))
  colnames(Cluster_output) <- c("Protein", "Peptide","trypticity","Control_shape","Sm_shape","Is_changing","Extra_filter","Change", "Experiment")
  
  k_num <- k
  sm_all <- sm_all
  cluster_LIP_all <- fuzzy_clusters_two_experiments(LIP, sm_list=sm_all, control = "Control", k=k_num, show_clusters = FALSE)
  
  Cluster_LIP_full <- cluster_LIP_all$full_set
  Cluster_LIP_full$description <- ifelse(Cluster_LIP_full$manual_order %in% c(1, 2, 3), "down", 
                                         ifelse(Cluster_LIP_full$manual_order %in% c(4, 5, 6), "up", "two-step"))
  
  count_peptide <- Cluster_LIP_full %>%
    group_by(Peptide) %>%
    unique() %>%
    summarise(n = n()) %>%
    filter(n == 2*k)
  
  Cluster_LIP_full %<>%
    filter(Peptide %in% count_peptide$Peptide) %>%
    group_by(Peptide, Condition, description, Protein, trypticity) %>%
    dplyr::summarise(Sum_membership = sum(Membership))
  
  Cluster_LIP_full %<>%
    ungroup( ) %>%
    group_by(Peptide, description) %>%
    mutate(diff = Sum_membership[Condition == sm_all[1]] - Sum_membership[Condition == sm_all[2]])
  
  Cluster_diff <- Cluster_LIP_full %>%
    group_by(Peptide) %>%
    mutate(Is_changing = ifelse(max(abs(diff)) > 0.5, TRUE, FALSE )) %>%
    mutate(Control_shape = description[ Sum_membership == max(Sum_membership[Condition == sm_all[1]]) & Condition == sm_all[1] ]) %>%
    mutate(Sm_shape = description[ Sum_membership == max(Sum_membership[Condition == sm_all[2]]) & Condition == sm_all[2] ]) %>%
    mutate(Extra_filter = ifelse(Sm_shape != Control_shape & Is_changing == FALSE, TRUE, FALSE))
  
  Cluster_diff %<>%
    dplyr::select(Protein, Peptide, trypticity, Control_shape, Sm_shape, Is_changing, Extra_filter) %>%
    unique()
  
  Cluster_diff %<>%
    mutate(Change = ifelse(Is_changing == FALSE, "No change", 
                           ifelse((Control_shape == "up" & Sm_shape == "down") |  (Control_shape == "down" & Sm_shape == "up"), "Flip", 
                                  ifelse(Control_shape %in% c("up", "down") & Sm_shape == "two-step", "Two step start", 
                                         ifelse(Sm_shape %in% c("up", "down") & Control_shape == "two-step", "Two step stop", NA)))))
  
  
  Cluster_diff$Experiment <- paste(sm_all[1], sm_all[2], sep="_")

  Cluster_output <- rbind(Cluster_output, Cluster_diff)
  
  
  

  return(Cluster_output)
}


fuzzy_clusters_one <- function(LIP, sm, plot_which = "control", control="control", k = 8, show_clusters = TRUE) {
  #sm <- "Trehalose"
  #plot_which <- "Control"
  #control <- "Control"
  
  set.seed(123) 
  #  all_results <- list()
  #  all_results$betaine <- LIP$betaine
  #  sm <- "betaine"
  #  control <- "control"
  #  k <- 10
  
  # only select few points from the fit
  temp <- LIP[[sm]]$python_fit$t[LIP[[sm]]$python_fit$type == "fitted"]  %>% sort() %>% unique() %>% .[seq(from = 1, to = length(.), by=5 )]
  
  #  sm_list <- c("trehalose", "betaine")
  high_quality <- LIP[[sm]]$python_score %>%   # extract python "score" data
      filter(!!ensym(plot_which) > 0.75) %>%      # only take the ones that have a good enough fit
      dplyr::select(peptide, trypticity, protein) %>%
      `colnames<-`(c("Peptide", "trypticity", "Protein"))
    
    # Select the fitted peptides for control
  res <- LIP[[sm]]$python_fit %>% 
      filter(condition %in% c(plot_which), 
             type == "fitted", 
             peptide %in% high_quality$Peptide) %>%
      mutate(peptide = paste(peptide, sm, sep="%%%"))
  
  # make a matrix for clustering
  res %<>%
    filter( t %in% temp) %>%
    reshape2::dcast(peptide ~ t, value.var = "y", fun.aggregate = function(x) mean(x))
  
  rownames(res) <- res$peptide
  res <- res[,-1]
  
  # Perform fuzzy k-means clustering
  
  #  k <- 15
  fuzzy <- fclust::FKM.ent(res, k=k, stand=1, ent=2)
  
  fuzzy_plot <- fuzzy$H %>%
    reshape2::melt() %>%
    as.data.frame()
  
  all_plot <- fuzzy$X %>%
     reshape2::melt() %>%
    `colnames<-`(c("Peptide_con", "temp", "y"))

  clus <- fuzzy$clus %>%
    as.data.frame() %>%
    mutate(Peptide_con = rownames(.))

  
  # Take the median values of clusters and describe their characteristics:
  # increasing trend (up_down), whether it has two step profile ("is_agg")
  median_clusters <- fuzzy_plot %>%
    `colnames<-`(c("cluster", "temperature", "median_value")) %>%
    group_by(cluster) %>%
    mutate(diff = median_value[temperature == 37] - median_value[temperature == 76]) %>% 
    mutate(is_agg = ifelse(abs(diff)> 1 , FALSE, TRUE)) %>%
    mutate(up_down = ifelse(diff>0, "down", "up") ) 
  
  median_clusters %<>%
    group_by(cluster) %>%
    mutate(sum_median = sum(median_value))
  
  cluster_order <- median_clusters %>%
    ungroup() %>%
    # sum median - easy way to get the area under the curve --> differentiate between thermo sensitive/stable proteins
    dplyr::select(cluster, is_agg, up_down, sum_median) %>%
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
    mutate(manual_order = ifelse(is_agg == FALSE & 
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
                                                                             value_order == "max_value", "7", 8) ) ) ) ) ) ) )
  
  
  cluster_order %<>% 
    dplyr::select(cluster, manual_order)
  
  median_clusters <- fuzzy_plot %>%
    `colnames<-`(c("cluster", "temperature", "median_value")) %>%
    group_by(cluster) %>%
    mutate(diff = median_value[temperature == 37] - median_value[temperature == 76]) %>%
    mutate(is_agg = ifelse(abs(diff)> 1, FALSE, TRUE)) %>%
    mutate(up_down = ifelse(diff>0, "down", "up") )  %>%
    plyr::join(., cluster_order, by="cluster")
  
  
  fuzzy_clusters <- fuzzy$clus %>% 
    as.data.frame() %>%
    mutate(Peptide = rownames(.)) %>%
    mutate(cluster = paste("Clus", as.character(.$Cluster), sep=" ")) %>%
    mutate(Condition = strsplit(Peptide, "%%%") %>% lapply(., function(x) '[['(x,2)) %>% unlist() ) %>%
    mutate(Peptide = strsplit(Peptide, "%%%") %>% lapply(., function(x) '[['(x,1)) %>% unlist() ) 
  
  fuzzy_clusters %<>%
    plyr::join(., cluster_order, by="cluster")
  
  fuzzy_clusters %<>%
    dplyr::select(-Cluster, -cluster) %>%
    `colnames<-`(c("Membership_degree", "Peptide", "Condition", "Cluster")) %>%
    plyr::join(., high_quality, by="Peptide")
  
  output <- list()
  #  show_clusters=TRUE
  if(show_clusters){
    cluster_plot <- ggplot(median_clusters, aes(x=temperature, y=median_value, col=cluster)) +
      geom_line() +
      facet_wrap(~manual_order) +
      theme_minimal() +
      theme(legend.position = "none") 
    
    output$plot <- cluster_plot
    
  }
  
  
  output$data <- fuzzy_clusters
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
  output$whole_plot <-   plyr::join(all_plot,clus, by="Peptide_con") %>%
    mutate(cluster = paste("Clus", Cluster, sep=" ")) %>%
    plyr::join(., cluster_order, by="cluster")
  
  return(output)
  
}


