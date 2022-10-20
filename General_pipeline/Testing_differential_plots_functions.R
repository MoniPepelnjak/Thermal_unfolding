# Functions for area and score calculations 
# Goal: Calculate the area between the curves for the control/small molecule and determine which one is aggregation/stabilisation
# Next step: Combine this with the cluster profile, to determine whether protein is stabilised, destabilised
# Author: Monika Pepelnjak
# Date: 25.02.2021

### Function for calculating the area differences
my_palette = c(Betaine = "#B5C3AD", 
               Trehalose = "#d4a993", 
               TMAO = "#87A9C6", 
               Glycerol = "#a997a2", 
               Proline = "#edd9bb", 
               Glucose = "#C78C8C",
               Betaine_2 = "#8E9C86",
               TMAO_2 = "#718CA2",
               Glycerol_2 = "#86697B",
               Proline_2 = "#D8BF9A")

find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

#Area calculation: Takes the python fits and calculates the area between the two curves when confidence intervals do not overlap
# Otherwise the area is set to 0
area_calculation <- function(all_results, small_molecule, cutoff = 0.5, control="Control") {
  
  selected_peptides <- all_results[[small_molecule]][["python_fit"]] %>%
    .[.$condition %in% c(control, small_molecule) & .$type == "fitted",] %>%
     .[!kit::fduplicated(.), ] %>%
    .[.$peptide %in% all_results[[small_molecule]][["python_score"]]$peptide[all_results[[small_molecule]][["python_score"]]["control"] > cutoff &all_results[[small_molecule]][["python_score"]][small_molecule] > cutoff ],]

  selected_peptides  %<>%
    group_by(peptide, t) %>%
    mutate(differences = ifelse((conflik_lower[condition == small_molecule] > conflik_lower[condition == control] &
                                     conflik_lower[condition == small_molecule] < conflik_upper[condition == control]) | 
                                    (conflik_upper[condition == small_molecule] > conflik_lower[condition == control] &
                                       conflik_upper[condition == small_molecule] < conflik_upper[condition == control] | 
                                       
                                       conflik_lower[condition == control] > conflik_lower[condition == small_molecule] &
                                       conflik_lower[condition == control] < conflik_upper[condition == small_molecule]) | 
                                    (conflik_upper[condition == control] > conflik_lower[condition == small_molecule] &
                                       conflik_upper[condition == control] < conflik_upper[condition == small_molecule]    
                                     
                                     
                                    ), 
                                  0, (y[condition == small_molecule] - y[condition == control]) )) %>%
    
    
    .[.$condition %in% c(small_molecule),]
}


area_calculation_2 <- function(all_results, small_molecule, cutoff=0.5, control="Control") {
  
  selected_peptides <- all_results[[small_molecule]][["python_fit"]] %>%
    .[.$condition %in% c(control, small_molecule) & .$type == "fitted",] %>%
    .[!kit::fduplicated(.), ]%>%
    .[.$peptide %in% all_results[[small_molecule]][["python_fit"]]$peptide[all_results[[small_molecule]][["python_score"]][control] > cutoff &all_results[[small_molecule]][["python_score"]][small_molecule] > cutoff ],]
  
  
  selected_peptides  %<>%
    group_by(peptide, t) %>%
    mutate(differences = ifelse((conflik_lower[condition == small_molecule] > conflik_lower[condition == control] &
                                   conflik_lower[condition == small_molecule] < conflik_upper[condition == control]) | 
                                  (conflik_upper[condition == small_molecule] > conflik_lower[condition == control] &
                                     conflik_upper[condition == small_molecule] < conflik_upper[condition == control] | 
                                     
                                     conflik_lower[condition == control] > conflik_lower[condition == small_molecule] &
                                     conflik_lower[condition == control] < conflik_upper[condition == small_molecule]) | 
                                  (conflik_upper[condition == control] > conflik_lower[condition == small_molecule] &
                                     conflik_upper[condition == control] < conflik_upper[condition == small_molecule]    
                                   
                                   
                                  ), 
                                0, (y[condition == small_molecule] - y[condition == control]) )) %>%
    
    mutate(strict_condition = ifelse((y[condition == small_molecule] > conf_lower[condition == control] &
                                   y[condition == small_molecule] < conf_upper[condition == control]) | 
                                     y[condition == control] > conf_lower[condition == small_molecule] &
                                     y[condition == control] < conf_upper[condition == small_molecule], 
                                0, (y[condition == small_molecule] - y[condition == control]) )) %>%
    
    .[.$condition %in% c(small_molecule),]
  
  return(selected_peptides)
}

### Find where is the apex of each peak - for some weird curves you will get more than one peak and assign aggregation/stabilisation

asign_peak_apex <-   function( selected_peptides ) {
    
    selected_peptides %<>% 
      group_by(peptide) %>%
      arrange(peptide, t) %>%
      mutate(apex = ifelse(t %in% t[find_peaks(abs(differences))], TRUE, FALSE))
    
    selected_peptides %<>%
      ungroup() %>%
      group_by(peptide) %>%
      arrange(t) %>%
      mutate(peak = ifelse(abs(differences) > 0, TRUE, FALSE), 
             temp1 = cumsum(!peak)) %>%
      group_by(peptide, temp1) %>%
      mutate(goal =  +(row_number() == which.max(peak) & any(peak))) %>%
      group_by(peptide) %>%
      mutate(goal = ifelse(peak, cumsum(goal), NA)) %>%
      dplyr::select(-peak, -temp1) %>%
      group_by(goal, peptide) %>%
      mutate(peak_sum = sum(differences)) %>%
      mutate(define_peak = ifelse( abs(peak_sum) == 0, NA, 
                                   ifelse( any(t < 39), "binding", 
                                           ifelse(any(t > 74), "aggregation", "stabilisation")))) %>%
      group_by(peptide) %>%
      mutate(define_new_peak = ifelse(define_peak == "aggregation" & any(strict_condition[t>74] == 0), "not_aggregation",  define_peak))
    
    return(selected_peptides)
}

##### Option with only area between the confidence intervals
area_between_CI <- function(all_results, small_molecule, cutoff=0.5, control="Control") {
  
  selected_peptides <- all_results[[small_molecule]][["python_fit"]] %>%
    .[.$condition %in% c(control, small_molecule) & .$type == "fitted",] %>%
    .[!kit::fduplicated(.), ]%>%
    .[.$peptide %in% all_results[[small_molecule]][["python_fit"]]$peptide[all_results[[small_molecule]][["python_score"]][control] > cutoff &all_results[[small_molecule]][["python_score"]][small_molecule] > cutoff ],]
  
  
  selected_peptides  %<>%
    group_by(peptide, t) %>%
    mutate(differences = ifelse((conflik_lower[condition == small_molecule] > conflik_lower[condition == control] &
                                   conflik_lower[condition == small_molecule] < conflik_upper[condition == control]) | 
                                  (conflik_upper[condition == small_molecule] > conflik_lower[condition == control] &
                                     conflik_upper[condition == small_molecule] < conflik_upper[condition == control] | 
                                     
                                     conflik_lower[condition == control] > conflik_lower[condition == small_molecule] &
                                     conflik_lower[condition == control] < conflik_upper[condition == small_molecule]) | 
                                  (conflik_upper[condition == control] > conflik_lower[condition == small_molecule] &
                                     conflik_upper[condition == control] < conflik_upper[condition == small_molecule]    
                                   
                                   
                                  ), 
                                0, 1)) %>%
    
    mutate(test_1 = ifelse(differences == 0, 0, 
                           (conflik_lower[condition == small_molecule] - conflik_upper[condition == control]))) %>%
    mutate(test_2 = ifelse(differences == 0, 0, 
                           (conflik_upper[condition == small_molecule] - conflik_lower[condition == control]))) %>%
    
    mutate(strict_condition = ifelse((y[condition == small_molecule] > conf_lower[condition == control] &
                                        y[condition == small_molecule] < conf_upper[condition == control]) | 
                                       y[condition == control] > conf_lower[condition == small_molecule] &
                                       y[condition == control] < conf_upper[condition == small_molecule], 
                                     0, (y[condition == small_molecule] - y[condition == control]) )) %>%
    
    .[.$condition %in% c(small_molecule),]
  
  selected_peptides %<>%
  mutate(differences = ifelse(abs(test_1) <= abs(test_2), test_1, test_2)) %>%
    dplyr::select(-test_1, -test_2)
  
  return(selected_peptides)
}

### Find where is the apex of each peak - for some weird curves you will get more than one peak and assign aggregation/stabilisation

asign_peak_apex <-   function( selected_peptides ) {
  
  selected_peptides %<>% 
    group_by(peptide) %>%
    arrange(peptide, t) %>%
    mutate(apex = ifelse(t %in% t[find_peaks(abs(differences))], TRUE, FALSE))
  
  selected_peptides %<>%
    ungroup() %>%
    group_by(peptide) %>%
    arrange(t) %>%
    mutate(peak = ifelse(abs(differences) > 0, TRUE, FALSE), 
           temp1 = cumsum(!peak)) %>%
    group_by(peptide, temp1) %>%
    mutate(goal =  +(row_number() == which.max(peak) & any(peak))) %>%
    group_by(peptide) %>%
    mutate(goal = ifelse(peak, cumsum(goal), NA)) %>%
    dplyr::select(-peak, -temp1) %>%
    group_by(goal, peptide) %>%
    mutate(peak_sum = sum(differences)) %>%
    mutate(define_peak = ifelse( abs(peak_sum) == 0, NA, 
                                 ifelse( any(t < 39), "binding", 
                                         ifelse(any(t > 74), "aggregation", "stabilisation")))) %>%
    group_by(peptide) %>%
    mutate(define_new_peak = ifelse(define_peak == "aggregation" & any(strict_condition[t>74] == 0), "not_aggregation",  define_peak))
  
  return(selected_peptides)
}



## Aggregate the information into peptide-level info 
# Option one: take the sum of all effects

aggregate_effects <- function(selected_peptides){
  #selected_peptides <- area_calculated
  
  aggregated_data <- selected_peptides %>%
    dplyr::select(peptide, protein, peak_sum, define_peak, goal) %>%
    .[!kit::fduplicated(.), ] #%>%
    #na.omit()
  
  wide_aggregation <- aggregated_data %>%
    reshape2::dcast(peptide + protein ~ define_peak, value.var="peak_sum", fun.aggregate = function(x) max(abs(x)))
  
  return(wide_aggregation)
}

plot_regions_colored <- function(all_results, selected_peptides, peptide_set, small_molecule, control="Control") {

  #selected_peptides <- area_calculated_2
  #peptide_set <- "NQIADLVGADPR"
  #small_molecule="Glucose"
  #control="Control"
  subset_peptides <-selected_peptides  %>%
    dplyr::select(peptide, define_peak, peak_sum, t, differences, goal) %>%
    filter(peptide %in% peptide_set)
  
  all_results_peptide <- all_results[[small_molecule]]$python_fit[all_results[[small_molecule]]$python_fit$peptide %in% peptide_set,]
  
  all_results_peptide %<>% 
    plyr::join(., subset_peptides, by=c("peptide", "t")) %>%
    filter(condition %in% c(control, small_molecule))
  
  only_peaks <- all_results_peptide %>%
    filter(type == "fitted") %>%
    filter(type == "fitted") %>%
    
    group_by(peptide,goal, t) %>%
    mutate( lower_point = ifelse(differences == 0, NA, 
                                 ifelse(differences < 0,  y[condition == small_molecule],  y[condition == control])),
            upper_point = ifelse(differences == 0, NA, 
                                 ifelse(differences > 0,  y[condition == small_molecule],  y[condition == control]))) %>%
    filter(!is.na(goal), 
           condition == small_molecule)
  
  test <- all_results_peptide
  small_molecule <- small_molecule
  test$condition <- factor(test$condition, levels=c(small_molecule, control))
  plot3 <- ggplot() +  
    geom_ribbon(data = only_peaks, aes(ymin=lower_point, ymax=upper_point, x=t, alpha=factor(goal), fill=define_peak)) +
    geom_point(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), size=2) + 
    geom_line(data = test[test$type == "fitted",], aes(y=y, x=t, col=condition), lwd=1.2) +
    scale_color_manual(values = c(my_palette[small_molecule], Control = "azure3")) +
    scale_fill_manual(values = c(binding = "#5D72A9", stabilisation = "#5DA993", aggregation="azure3"))+
    #  geom_smooth(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), method="loess") +
    facet_wrap(~peptide) +
    ylab("Scaled abundance") +
    xlab("Temperature [C]")+
    theme_bw() +
    scale_alpha_discrete(range = c(0.3, 0.300001)) +
    guides(alpha = FALSE)
  
  return(plot3)
  
}
  
  

plot_regions_colored_2 <- function(all_results, selected_peptides, peptide_set, small_molecule) {
  
  subset_peptides <-selected_peptides  %>%
    dplyr::select(peptide, define_peak, peak_sum, t, differences, st_dest) %>%
    filter(peptide %in% peptide_set)
  
  all_results_peptide <- all_results[[small_molecule]]$python_fit[all_results[[small_molecule]]$python_fit$peptide %in% peptide_set,]
  
  all_results_peptide %<>% 
    plyr::join(., subset_peptides, by=c("peptide", "t")) %>%
    filter(condition %in% c(control, small_molecule))
  
  only_peaks <- all_results_peptide %>%
    filter(type == "fitted") %>%
    filter(type == "fitted") %>%
    
    group_by(peptide,goal, t) %>%
    mutate( lower_point = ifelse(differences == 0, NA, 
                                 ifelse(differences < 0,  y[condition == small_molecule],  y[condition == control])),
            upper_point = ifelse(differences == 0, NA, 
                                 ifelse(differences > 0,  y[condition == small_molecule],  y[condition == control]))) %>%
    filter(!is.na(goal), 
           condition == small_molecule)
  
  only_peaks$st_dest <- factor(only_peaks$st_dest, levels = c("stabilisation", "destabilisation", "needs_revision"))
  
  test <- all_results_peptide
  test$condition <- factor(test$condition, levels=c(small_molecule, control))
  plot3 <- ggplot() +  
    geom_ribbon(data = only_peaks, aes(ymin=lower_point, ymax=upper_point, x=t, alpha=factor(goal), fill=st_dest)) +
    geom_point(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), size=2) + 
    geom_line(data = test[test$type == "fitted",], aes(y=y, x=t, col=condition), lwd=1.2) +
    scale_color_manual(values = c("indianred3", "#8F8F8F")) +
    scale_fill_manual(values = c("#00A071", "#A20000", "#FFF559"), drop=FALSE)+
    #  geom_smooth(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), method="loess") +
    facet_wrap(~peptide) +
    ylab("Scaled abundance") +
    xlab("Temperature [C]")+
    theme_bw() +
    scale_alpha_discrete(range = c(0.5, 0.500001)) +
    guides(alpha = FALSE)
  
  return(plot3)
  
}

plot_regions_colored_CI <- function(all_results, selected_peptides, peptide_set, small_molecule, control="Control") {
  
  #selected_peptides <- area_calculated_2
  #peptide_set <- "NQIADLVGADPR"
  #small_molecule="Glucose"
  #control="Control"
  subset_peptides <-selected_peptides  %>%
    dplyr::select(peptide, define_peak, peak_sum, t, differences, goal) %>%
    filter(peptide %in% peptide_set)
  
  all_results_peptide <- all_results[[small_molecule]]$python_fit[all_results[[small_molecule]]$python_fit$peptide %in% peptide_set,]
  
  all_results_peptide %<>% 
    plyr::join(., subset_peptides, by=c("peptide", "t")) %>%
    filter(condition %in% c(control, small_molecule))
  
  only_peaks <- all_results_peptide %>%
    filter(type == "fitted") %>%
    filter(type == "fitted") %>%
    
    group_by(peptide,goal, t) %>%
    mutate( lower_point = ifelse(differences == 0, NA, 
                                 ifelse(differences < 0,  conflik_upper[condition == small_molecule],  conflik_upper[condition == control])),
            upper_point = ifelse(differences == 0, NA, 
                                 ifelse(differences > 0,  conflik_lower[condition == small_molecule],  conflik_lower[condition == control]))) %>%
    filter(!is.na(goal), 
           condition == small_molecule)
  
  test <- all_results_peptide
  small_molecule <- small_molecule
  test$condition <- factor(test$condition, levels=c(small_molecule, control))
  plot3 <- ggplot() +  
    geom_ribbon(data = only_peaks, aes(ymin=lower_point, ymax=upper_point, x=t, alpha=factor(goal), fill=define_peak)) +
    geom_point(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), size=2) + 
    geom_line(data = test[test$type == "fitted",], aes(y=y, x=t, col=condition), lwd=1.2) +
    scale_color_manual(values = c(my_palette[small_molecule], Control = "azure3")) +
    scale_fill_manual(values = c(binding = "#5D72A9", stabilisation = "#5DA993", aggregation="azure3"))+
    #  geom_smooth(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), method="loess") +
    facet_wrap(~peptide) +
    ylab("Scaled abundance") +
    xlab("Temperature [C]")+
    theme_bw() +
    scale_alpha_discrete(range = c(0.3, 0.300001)) +
    guides(alpha = FALSE)
  
  return(plot3)
  
}
