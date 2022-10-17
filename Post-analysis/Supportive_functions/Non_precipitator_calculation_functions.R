## Define non-precipitators ##
## Author: Monika Pepelnjak
## Date: 16.04.2021

## Goal: define whether the protein is a "non-precipitator" 

define_precipitators_control <- function(TPP, Savitski = c(), both = FALSE, cutoff = 0.5){
  
  both=both
  # here the data is hard-coded to glucose because only one TPP experiment was performed, therefore all control samples are the same
  Non_precipitators <- TPP$glucose$python_fit %>% 
    filter(condition == "control") %>% # select control condition
    filter(peptide %in% TPP$glucose$python_score$peptide[TPP$glucose$python_score$control>0.5]) %>% 
    filter(t %in% c(37, 76)) %>% # Only select the first and last temperature
    group_by(peptide, t, type) %>%
    summarise(y = mean(y), 
              measured_all = n(), 
              t=t) %>%
    ungroup() %>%
    group_by(peptide) %>%
    summarise(FC = ifelse( 37 %in% t[type == "measured"] & 76 %in%  t[type == "measured"], 
                           y[type == "fitted" & t == 76]/y[type == "fitted" & t == 37], NA ) )
  
  
  Non_precipitators %<>%
    `colnames<-`(c("Protein", "FC_TPP")) 
  Non_precipitators$FC_TPP[is.na(Non_precipitators$FC_TPP)] <- 0
  
  Precipitator_output <- Non_precipitators
  
  # Based on the designed cutoff, determine whether the protein is a precipitator or not
  Precipitator_output %<>%
    mutate(Precipitator_TPP = ifelse(FC_TPP > cutoff, "Non-Precipitator", "Precipitator"))
  
  # Important: for this part the Savitski dataset is required
  # Non precipitator: there is more protein left in the supernatant than the designed cutoff
  if(both == TRUE){
    
    Savitski %<>%
      dplyr::select(Uniprot_id, contains("rel_fc_131L_Lysate_Mg")) %>%
      reshape2::melt(., id.var="Uniprot_id") %>%
      group_by(Uniprot_id) %>%
      summarise(FC = mean(value)) %>%
      `colnames<-`(c("Protein", "FC_Savitski")) %>%
      mutate(Precipitator_Savitski = ifelse(FC_Savitski > cutoff, "Non-Precipitator", "Precipitator" ))
    
    Precipitator_output <- plyr::join(Precipitator_output, Savitski, by="Protein", type="full")
    
  }
  
  return(Precipitator_output)

}

# The same function, only small molecules are analysed instead of the control
define_precipitators_sm <- function(TPP, small_molecule, cutoff = 0.5){
  
  Non_precipitators <- TPP[[small_molecule]]$python_fit %>%
    filter(condition == small_molecule) %>% 
    filter(peptide %in% TPP[[small_molecule]]$python_score$peptide[TPP[[small_molecule]]$python_score$control>0.5]) %>% 
    filter(t %in% c(37, 76)) %>% 
    group_by(peptide, t, type) %>%
    summarise(y = mean(y), 
              measured_all = n(), 
              t=t) %>%
    ungroup() %>%
    group_by(peptide) %>%
    summarise(FC = ifelse( 37 %in% t[type == "measured"] & 76 %in%  t[type == "measured"], 
                           y[type == "fitted" & t == 76]/y[type == "fitted" & t == 37], NA ) )
  
  
  Non_precipitators %<>%
    `colnames<-`(c("Protein", "FC_TPP")) 
  Non_precipitators$FC_TPP[is.na(Non_precipitators$FC_TPP)] <- 0
  
  Precipitator_output <- Non_precipitators
  
  Precipitator_output %<>%
    mutate(Precipitator_TPP = ifelse(FC_TPP > cutoff, "Non-Precipitator", "Precipitator")) %>%
    `colnames<-`(c("Protein", paste("FC_TPP_", small_molecule, sep=""),  paste("Precipitator_", small_molecule, sep=""))) 
  
  
  
  return(Precipitator_output)
  
}
