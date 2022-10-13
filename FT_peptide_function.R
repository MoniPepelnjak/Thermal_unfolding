## FT_HT_matching function
## The script defines the functions used to produce the supplementary figures
  # Matching HT and FT pairs
  # Plotting the matched sequences
  # Plotting the matched pairs

# Author: Monika Pepelnjak

# Cutoff = cutoff for goodness of the fit. Rationale = the noisy peptides contain less information and are more difficult to compare
matching_FT_HT <- function(LIP, sm, control = "Control", cutoff = 0.8) {
  
  # Select all peptides with good "goodness of fit"
  All_peptides <- LIP[[sm]]$python_score$peptide[LIP[[sm]]$python_score[control] > cutoff ] %>% as.character()
  FT_peptides <- LIP[[sm]]$python_score$peptide[LIP[[sm]]$python_score[control] > cutoff & LIP[[sm]]$python_score$trypticity == "FT" ] %>% as.character()
  
  # Define an empty vector of appropriate length
  FT_peptide_matched <- as.character(c(1:length(All_peptides)))
  
  # Match the FT peptide with HT peptide
  for( i in 1:length(All_peptides) ) {
    FT_peptide_matched[i] <- FT_peptides[grepl(All_peptides[i], FT_peptides)][1]
  }
  
  # Compose a df of the matches
  matched_peptides <- data.frame(FT_peptide = as.character(FT_peptide_matched), 
                                 HT_peptide = as.character(All_peptides)) %>% 
    na.omit() %>%  # remove rows when HT peptide was not matched
    mutate(matched = ifelse(FT_peptide == HT_peptide, "FT_peptide", "HT_match")) %>% # Figure out whether the peptide is matched or just the same FT peptide
    group_by(FT_peptide) %>% 
    mutate(n_matched = n()) # Calculate the number of HT peptides for each FT peptide
}

# Function to plot the matches - produces the figure S1A of the Osmolyte paper
plot_matches <- function(matched_peptides, FT_peptide_match, size1=12, size2=5){
  
  # Select the subset of peptides to plot by selecting the "parent" FT peptide 
  example <- matched_peptides %>%
    filter(FT_peptide == FT_peptide_match) %>%
    unique() %>%
    group_by(HT_peptide) %>%
    # Get the location of the start/end of the peptide
    mutate(match_start = gregexpr(HT_peptide, FT_peptide)[[1]]) %>%
    mutate(match_end = match_start + nchar(HT_peptide)) 
  # Order peptides by their start
  example <- example[order(-(example$match_start), (example$match_end)),]
  example$HT_peptide <- factor(example$HT_peptide, levels=c(example$HT_peptide))

  # Plot the matches
  plot4 <- ggplot(example, aes(x = match_start, xend=match_end, col=matched,y=HT_peptide, yend=HT_peptide)) +
    geom_text(aes(y=HT_peptide,x = 1, hjust=0, label=FT_peptide), col="gray80",size=size2, family="Courier") +
    geom_segment(size=size1) +
    geom_text(aes(y=HT_peptide,x = (example$match_start), hjust=0, label=HT_peptide), col="gray25",size=size2, family="Courier") +
    #scale_color_manual(values=c(FT_peptide = "#8BB9B4", HT_match = "#C79C8B")) +
    scale_color_manual(values=c(FT_peptide = "#A1948D", HT_match = "#DDD8D5")) +
    theme_bw() +
    #scale_y_discrete(name ="", 
    #                 limits=c("FT","HT","HT", "HT","HT")) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_text(size=15),
          axis.ticks.y=element_blank(),
          panel.grid = element_blank(), 
          panel.border = element_blank())+
    theme(legend.position = "")
  
  return(plot4)
    
  # gregexpr(example$HT_peptide[1], example$FT_peptide[1])
  
}


# Plot the peptide melting profiles for the matched peptides
# Select which condition to plot - in this case Control is plotted

plot_HT_pairs <- function(LIP, sm, plot_which = "Control", matched_peptides, FT_peptide_sel) {

  example <- matched_peptides %>%
    unique() %>%
    filter(FT_peptide == FT_peptide_sel) %>%
    group_by(HT_peptide) %>%
    mutate(match_start = gregexpr(HT_peptide, FT_peptide)[[1]]) %>%
    mutate(match_end = match_start + nchar(HT_peptide)) %>%
    filter(FT_peptide != HT_peptide)
  
  example <- example[order((example$match_start), -(example$match_end)),]
  example$HT_peptide <- factor(example$HT_peptide, levels=c(example$HT_peptide))
  
 
  test <- NULL
  
  for( i in 1:nrow(example) ){
    
    test_out <- LIP[[sm]]$python_fit[LIP[[sm]]$python_fit$peptide %in% c(as.character(example$HT_peptide[i]), example$FT_peptide[i]),] %>%
      filter(condition == "Control") %>% 
      mutate(is_FT_peptide = ifelse(peptide %in% example$FT_peptide[i], "FT_peptide", "HT_match")) %>%
      mutate(FT_match = example$HT_peptide[i])
    
    test <- rbind(test, test_out)
    
  }
  
  test$FT_match <- factor(test$FT_match, levels=example$HT_peptide)

    plot3 <- ggplot(test) +  
      geom_ribbon(data = test[test$type == "fitted",], aes(ymin=conflik_lower, ymax=conflik_upper, x= t, fill=is_FT_peptide, group=peptide), alpha=0.8) +
      geom_line(data = test[test$type == "fitted",], aes(y=y, x=t, col=is_FT_peptide, group=peptide), lwd=1.2) +
      geom_line(data = test[test$type == "fitted",], aes(y=y, x=t,group=peptide), lwd=0.5, col="black") +
      geom_point(data = test[test$type == "measured",], aes(y=y, x=t, fill=is_FT_peptide, group=peptide), size=2, pch=21, col="black") + 
      
      #scale_color_manual(values=c(FT_peptide = "#8BB9B4", HT_match = "#C79C8B")) +
      #scale_fill_manual(values=c(FT_peptide = "#8BB9B4", HT_match = "#C79C8B"))+
      scale_color_manual(values=c(FT_peptide = "#A1948D", HT_match = "#DDD8D5")) +
      scale_fill_manual(values=c(FT_peptide = "#A1948D", HT_match = "#DDD8D5")) +
      #  geom_smooth(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), method="loess") +
      facet_wrap(~FT_match, scales = "free") +
      theme_minimal() +
      theme(axis.line = element_blank(), panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(), legend.position = "none", 
            text = element_text(size=15), 
            panel.border = element_rect(fill=NA), 
            strip.text = element_text(size=12)) +
      annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
      annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) + 
      xlab("Temperature") +
      ylab("Scaled intensity")
    
    plot(plot3)
    return(plot3)}
  


# Plot the peptide melting profiles for the matched peptides
# Select which condition to plot - in this case Control is plotted
plot_HT_pairs_individual <- function(LIP, sm, plot_which = "Control", matched_peptides, FT_peptide_sel) {
  
  example <- matched_peptides %>%
    unique() %>%
    filter(FT_peptide == FT_peptide_sel) %>% # Selected FT peptides
    group_by(HT_peptide) %>%
    mutate(match_start = gregexpr(HT_peptide, FT_peptide)[[1]]) %>%
    mutate(match_end = match_start + nchar(HT_peptide)) %>%
    filter(FT_peptide != HT_peptide)
  
  example <- example[order((example$match_start), -(example$match_end)),]
  example$HT_peptide <- factor(example$HT_peptide, levels=c(example$HT_peptide))

  test <- LIP[[sm]]$python_fit[LIP[[sm]]$python_fit$peptide %in% c(as.character(example$HT_peptide), example$FT_peptide),] %>%
      filter(condition == "Control") %>% 
      mutate(is_FT_peptide = ifelse(peptide %in% example$FT_peptide, "FT_peptide", "HT_match"))
  plot_all <- list()
  all_peptides <- test$peptide %>% unique()
  
  for( i in 1:length(all_peptides)) {
  
    test_plot <- test[test$peptide %in% all_peptides[i],]
    
  plot3 <- ggplot(test_plot) +  
    geom_ribbon(data = test_plot[test_plot$type == "fitted",], aes(ymin=conflik_lower, ymax=conflik_upper, x= t, fill=is_FT_peptide, group=peptide), alpha=0.2) +
    geom_point(data = test_plot[test_plot$type == "measured",], aes(y=y, x=t, col=is_FT_peptide, group=peptide), size=2) + 
    geom_line(data = test_plot[test_plot$type == "fitted",], aes(y=y, x=t, col=is_FT_peptide, group=peptide), lwd=1.2) +
    scale_color_manual(values=c(FT_peptide = "#8BB9B4", HT_match = "#C79C8B")) +
    scale_fill_manual(values=c(FT_peptide = "#8BB9B4", HT_match = "#C79C8B"))+
    #  geom_smooth(data = test[test$type == "measured",], aes(y=y, x=t, col=condition), method="loess") +
    facet_wrap(~peptide, scales = "free") +
    theme_minimal() +
    theme(axis.line = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(), legend.position = "none") +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) + 
    xlab("Temperature") +
    ylab("Scaled intensity")
  
  plot_all[[all_peptides[i]]] <- plot3
  
  
  }
  
  
  return(plot_all)}


