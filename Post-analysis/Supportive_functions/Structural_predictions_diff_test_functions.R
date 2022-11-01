## Test significant prediction features functions

#All_combined_predictions <- read.csv("...")

Test_two_groups <-  function(All_combined_predictions, group1, group2, name_group_1,name_group_2 ){
  
  All_combined_predictions$group <- ifelse(All_combined_predictions$Protein %in% group1, 
                                           name_group_1, ifelse(All_combined_predictions$Protein %in% group2, name_group_2, "not_interesting") )
  
  
  All_combined_predictions_sub <- All_combined_predictions[All_combined_predictions$Protein %in% c(group1, group2),]
  
  
  All_combined_predictions_sub %<>%
    reshape2::melt(., id.vars=c("Protein", "group")) %>%
    group_by(variable) %>%
    summarise(p_value = wilcox.test(value~group)$p.value)
  
  output <- list()
  output[["Significance_test"]] <- All_combined_predictions_sub
  output[["All_values"]] <- All_combined_predictions[All_combined_predictions$Protein %in% c(group1, group2),]
  return(output)
  
}

Test_three_groups_pairwise <-  function(All_combined_predictions, groupA, groupB, groupC, name_group_A,name_group_B, name_group_C ){
  
  All_combined_predictions$group <- ifelse(All_combined_predictions$Protein %in% groupA, 
                                           name_group_A, ifelse(All_combined_predictions$Protein %in% groupB, name_group_B,
                                                                ifelse(All_combined_predictions$Protein %in% groupC, name_group_C, "not_interesting") ))
  
  
  All_combined_predictions_sub <- All_combined_predictions[All_combined_predictions$Protein %in% c(groupA, groupB, groupC),]
  
  
  All_combined_predictions_sub %<>%
    reshape2::melt(., id.vars=c("Protein", "group")) %>%
    group_by(variable) %>%
    summarise(p_value_AB = wilcox.test(value[group == name_group_A], value[group == name_group_B])$p.value, 
              p_value_BC = wilcox.test(value[group == name_group_C], value[group == name_group_B])$p.value,
              p_value_AC = wilcox.test(value[group == name_group_A], value[group == name_group_C])$p.value) 
  
  output <- list()
  output[["Significance_test"]] <- All_combined_predictions_sub
  output[["All_values"]] <- All_combined_predictions[All_combined_predictions$Protein %in% c(groupA, groupB, groupC),]
  return(output)
  
}
