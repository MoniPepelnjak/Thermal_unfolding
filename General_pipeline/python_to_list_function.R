### Function to load the python data and organise it with a list
### Peptide matching functions (by position)

python_to_list <- function(defined_list, sm, path_file, common_path, fasta = "ecoli") {
  #defined_list <- list()
  sm_list <- list()
  Import <- read.delim(path_file$path[path_file$file == "Spectronaut_export" & path_file$sm == sm], header=T) %>% 
    dplyr::select(PG.FastaFiles, PG.ProteinAccessions, PEP.StrippedSequence) %>%
    unique() %>%
    .[grepl(pattern = fasta, .$PG.FastaFiles),]
  
  Import %<>%
    dplyr::select(PG.ProteinAccessions, PEP.StrippedSequence) %>%
    `colnames<-`(c("protein", "peptide"))
  path <- paste(common_path, "0" , unique(path_file$Experiment[path_file$sm == sm]), "/", sm, "/", "python_output/", sep="" )
  
  # import HT and FT fit and combine into list
  HT_import_fit <- read.csv(paste(path, paste("solution_", sm, "_HT.csv", sep=""), sep=""))
  HT_import_fit %<>% plyr::join(., Import, by="peptide") %>%
    mutate(trypticity = "HT")
  
  FT_import_fit <- read.csv(paste(path, paste("solution_", sm, "_FT.csv", sep=""), sep=""))
  FT_import_fit %<>% plyr::join(., Import, by="peptide") %>%
    mutate(trypticity = "FT")
  
  sm_list[["python_fit"]] <- rbind(HT_import_fit, FT_import_fit)
  
  # import HT and FT MLL and combine into list
  HT_import_score <- read.csv(paste(path, paste("MLL_", sm, "_HT.csv", sep=""), sep=""))
  HT_import_score %<>% plyr::join(., Import, by="peptide") %>%
    mutate(trypticity = "HT")
  
  FT_import_score <- read.csv(paste(path, paste("MLL_", sm, "_FT.csv", sep=""), sep=""))
  FT_import_score %<>% plyr::join(., Import, by="peptide") %>%
    mutate(trypticity = "FT")
  
  sm_list[["python_score"]] <- rbind(HT_import_score, FT_import_score)
  
  defined_list[[sm]] <- sm_list
  
  return(defined_list)
}
