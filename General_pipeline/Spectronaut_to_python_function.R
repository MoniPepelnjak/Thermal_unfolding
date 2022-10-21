# Python import function
# Author: Monika Pepelnjak
# Description: The function takes the import from Spectronout and transforms the data for python fitting approach
# The function also saves the files in the described path

#Python analysis of thermal stabilization
# The script firstly filters the data - for now only complete cases, proteotypic peptides, removes contaminants and iRT peptides
# In the next step peptides are being scaled.

# sm = small molecule you want to analyse
# common_path = "common path where you would want to save your data"
# include replicates = which replicates would you want to keep in your data
# fasta = partial name of your fasta file of your specific organism, to remove the iRT peptides and contaminants

Spectronaut_to_python <- function(sm, path_file, common_path =  "/Users/moni/Documents/Phd/Experiments", include_replicates = c(1,2,3), fasta = "ecoli") {
  Import <- read.delim(path_file$path[path_file$file == "Spectronaut_export" & path_file$sm == sm], header=T) 

  Annotation <- read.csv(path_file$path[path_file$file == "annotation_file" & path_file$sm == sm])  %>% 
    mutate(run = paste(.$condition, .$temperature, .$replicate, sep="_"))
  
  colnames(Import)[grepl("StrippedSequence", colnames(Import))] <- "EG.StrippedSequence" # renaming such that column names are equalised from different possible spectronaut imports
  colnames(Import)[grepl("IsProteotypic", colnames(Import))] <- "EG.IsProteotypic"
  
  FT_peptides <- Import$EG.StrippedSequence[Import$PEP.DigestType....Trypsin.P. == "Specific"] %>% unique()
  HT_peptides <- Import$EG.StrippedSequence[Import$PEP.DigestType....Trypsin.P. %in% c( "Specific-N", "Specific-C")] %>% unique()
  
  Annotated <- Import %>% 
    #join with the annotation file
    plyr::join(., Annotation, by="R.FileName") %>%
    filter(PEP.Quantity > 100) %>%
    filter(EG.IsProteotypic %in% "True", # Only select the proteotypic peptides
           replicate %in% include_replicates) %>%  # Define which replicates to take - this is included so that you always select the same number of replicates for all conditions
    .[grepl(pattern = fasta, x = .$PG.FastaFiles, ignore.case = TRUE),] # only keep the proteins from your organism 
  
  # Allow only for up to 5 missing values
  Annotated_full <- Annotated %>%
    reshape2::dcast(EG.StrippedSequence ~ run, value.var = "PEP.Quantity", fun.aggregate = function(x) mean(x, na.rm=T)) %>% 
    reshape2::melt() %>% 
    `colnames<-`(c("EG.StrippedSequence", "run", "PEP.Quantity")) %>%
    plyr::join(., Annotation, by="run")
  
  Annotated_full %<>%
    group_by(EG.StrippedSequence) %>%
    mutate(number_missing = sum(is.na(PEP.Quantity)))
  
  Annotated_full %<>%
    mutate(which_missing = ifelse(is.na(PEP.Quantity), R.FileName, NA))
  
  python_export <- Annotated_full %>%
    ungroup() %>%
    filter(number_missing < 5) 

  python_export <- python_export %>% 
    ungroup() %>% 
    dplyr::select(-which_missing) %>%
    group_by(EG.StrippedSequence, condition, replicate) %>%
    na.omit() %>%
    # scale the data between 0 and 1
    mutate( pcr_minmax = (PEP.Quantity - min(PEP.Quantity))/(max(PEP.Quantity)-min(PEP.Quantity))) %>%
    ungroup() 
  
  # python export file
  export_table <- python_export %>%
    dplyr::select(EG.StrippedSequence, condition, temperature, pcr_minmax) %>% 
    `colnames<-`(c("uniqueID", "condition", "x", "y"))
  

  # create output directory
  path_output <- paste(common_path, paste("0", unique(path_file$Experiment[path_file$sm == sm]), sep=""), sm, "python_import", sep="/" )
  
  if (!dir.exists(path_output)) {
    dir.create(path_output,showWarnings = FALSE)
  }
  
  # export HT and FT as .csv
  write.csv(export_table[export_table$uniqueID %in% FT_peptides, ], paste(path_output,  "/", sm, "_FT", ".csv", sep=""))
  write.csv(export_table[export_table$uniqueID %in% HT_peptides, ], paste(path_output,  "/", sm, "_HT", ".csv", sep=""))

  print(paste("Succesfully exported", sm, "files. The files can be found in", path_output))
  
}








