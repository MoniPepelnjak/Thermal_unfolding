#### Functions ####
fuzzy_clusters_two_conditions <- function(all_results, sm, control="control", k = 8, show_clusters = FALSE, cutoff= 1.65) {
  set.seed(123) 
  
  high_quality <- all_results[[sm]]$python_score %>%   # extract python "score" data
    filter(!!ensym(sm) > 0.5 &
             !!ensym(control) > 0.5) %>%      # only take the ones that have a good enough fit
    dplyr::select(peptide, trypticity) %>%
    `colnames<-`(c("Peptide", "trypticity"))
  
  # Select the fitted peptides for control
  res <- all_results[[sm]]$python_fit %>% 
    filter(condition %in% c(sm, control), 
           type == "fitted", 
           peptide %in% high_quality$Peptide) %>%
    mutate(peptide = paste(peptide, condition, sep="%%%"))
  
  # only select every 5th point from the fit
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
  
  all_plot <- fuzzy$X %>%
    reshape2::melt() %>%
    `colnames<-`(c("Peptide_con", "temp", "y"))
  
  clus <- fuzzy$clus %>%
    as.data.frame() %>%
    mutate(Peptide_con = rownames(.))
  
  plot_all <- plyr::join(all_plot,clus, by="Peptide_con")
  
  # 
  # Take the median values of clusters and describe their characteristics:
  # increasing trend (up_down), whether it has two step profile ("is_agg")
  median_clusters <- fuzzy_plot %>%
    `colnames<-`(c("cluster", "temperature", "median_value")) %>%
    group_by(cluster) %>%
    mutate(diff = median_value[temperature == 37] - median_value[temperature == 76]) %>% # difference at 1st and last Temperature
    ungroup() %>%
    mutate(is_agg = ifelse(abs(diff) > max(abs(diff))/cutoff , FALSE, TRUE)) %>%
    group_by(cluster) %>%
    mutate(up_down = ifelse(diff>0, "down", "up"),
           max_diff = (max(median_value)+2) - (min(median_value+2)),
           inverted_bellcurve = median_value[temperature == 37] + median_value[temperature == 76])

  cluster_order <- median_clusters %>%
    ungroup() %>%
    dplyr::select(cluster, is_agg, up_down, max_diff,diff,inverted_bellcurve) %>%
    unique() 

  
  # Define cluster numbers based on previously defined characteristics: up_down, is_agg, value_order
  # If Number of clusters is 8, for most datasets we will get one cluster for each number. If number is higher, 
  # we get more clusters in the mean range
  
  
  # Whole goal of this: Always get the same number for clusters with similar profile
  # cluster into 4 groups: up, down and 2 bimodal profile. Unmatched profile go into a fifth group
  cluster_order %<>%
    ungroup() %>%
    mutate(manual_order = ifelse(max_diff < 0.8, 5, # noise
                                 ifelse(
                                    is_agg == FALSE & 
                                   up_down == "down", 1, # sigmoidal curve going down
                                               ifelse(is_agg == FALSE & 
                                                        up_down == "up", 2, # sigmoidal curve going up
                                                                    ifelse(inverted_bellcurve > 0, 3, 4))))) # bimodal 3:^, 4:u

  
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
  show_clusters=TRUE
  if(show_clusters){
    cluster_plot <- ggplot(median_clusters, aes(x=temperature, y=median_value, col=cluster)) +
      geom_line() +
      facet_wrap(~manual_order) +
      theme_minimal() #+
    theme(legend.position = "none") 
    
    output$plot <- cluster_plot
    
  }
  
  output$plot
  output$data <- fuzzy_clusters %>%
    group_by(Peptide) %>%
    mutate(is_flip = ifelse( any(Cluster %in% c(1)) & any(Cluster %in% c(2)) | any(Cluster %in% c(3)) & any(Cluster %in% c(4)), TRUE, FALSE )  )
  

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
  rm(fuzzy_clusters)
  rm(median_clusters)
  rm(res)
  rm(temp)
  rm(cluster_order)
  rm(clus)
  rm(high_quality)
  rm(fuzzy_plot)
  rm(all_plot)
  rm(plot_all)

  
  return(output)
}
