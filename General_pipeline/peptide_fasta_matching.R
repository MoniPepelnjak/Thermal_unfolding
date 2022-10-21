#### Function ####
peptide_fasta_matching <- function(peptide_list, fasta_ecoli){
  
  protein_sequence <- fasta_ecoli %>%
    lapply(., FUN = function(x) paste(x, sep="", collapse=""))
  proteins_ecoli <- names(protein_sequence)
  names(protein_sequence) <-  lapply(strsplit(proteins_ecoli, split="\\|"), FUN = function(x) "[[" (x,2))
  
  df <- data.frame(matrix(unlist(protein_sequence), nrow=length(protein_sequence), byrow=TRUE))
  df$proten <- names(protein_sequence)
  colnames(df) <- c("Sequence", "Protein")
  
  # match peptides to proteins
  matches <- sapply(peptide_list, FUN = function(x) names(protein_sequence)[grep(x, protein_sequence)])
  
  is_unique <- sapply(matches, length)
  matches_2 <- matches[is_unique == 1]
  df_matched <- data.frame(matrix(unlist(matches_2), nrow=length(matches_2), byrow=TRUE))
  df_matched$Peptide <- names(matches_2)
  colnames(df_matched) <- c("Protein", "Peptide")
  
  joined_df <- plyr::join(df_matched, df, by="Protein")
  
  joined_df %<>%
    mutate(start = stringr::str_locate(Sequence, Peptide)) %>%
    mutate(end = start + nchar(Peptide) - 1) %>%
    mutate(prot_len = nchar(Sequence)) %>%
    dplyr::select(-Sequence) %>%
    group_by(Peptide) %>%
    mutate(start = start[1], 
           end = end[1])
  
  
  return(joined_df)
  
}
