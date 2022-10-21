#### Functions ####
### Area calculation: Takes the python fits and calculates the area between the two curves when confidence intervals do not overlap
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
                                     conflik_upper[condition == control] < conflik_upper[condition == small_molecule]), 
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

assign_peak_apex <-   function( selected_peptides ) {
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
    mutate(goal =  +(row_number() == which.max(peak) & any(peak))) %>% # identify peak regions to sum up area
    group_by(peptide) %>%
    mutate(goal = ifelse(peak, cumsum(goal), NA)) %>%
    dplyr::select(-peak, -temp1) %>%
    group_by(goal, peptide) %>%
    mutate(peak_sum = sum(differences)) %>%
    mutate(define_peak = ifelse( abs(peak_sum) == 0, NA, # define peak type by region,eg at beginnig binding
                                 ifelse( any(t < 39), "binding", 
                                         ifelse(any(t > 74), "aggregation", "stabilisation")))) %>%
    group_by(peptide) %>% # no aggregation if strict condition towards the end is 0?
    mutate(define_new_peak = ifelse(define_peak == "aggregation" & any(strict_condition[t>74] == 0), "not_aggregation",  define_peak))
  
  return(selected_peptides)
}
