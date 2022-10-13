# Script to generate the supplementary data for the Osmolyte paper
# Author: Monika Pepelnjak

# Exporting the important list
# Load in the important libraries
library(xlsx)
library(magrittr)
library(dplyr)
library(protti)

# Import all the necessary datasets
output_binding <- readRDS("/Volumes/My Passport for Mac/data/Binding_all_osmoytes.rds")
significant_all <- readRDS("/Volumes/My Passport for Mac/data/Significant_both_quantiles.rds")
TPP_significant <- readRDS("/Volumes/My Passport for Mac/data/Significant_TPP.rds")
comparisons_LIP_all <- readRDS("/Volumes/My Passport for Mac/data/Susceptibility_agg.rds")
comparisons <- readRDS("/Volumes/My Passport for Mac/data/Promoted_reduced_agg.rds")
output_3_export <- readRDS("/Volumes/My Passport for Mac/data/Promoted_reduced_agg.rds")

# List the osmolytes you want to include
sm_list <- c("TMAO", "Betaine", "Glycerol", "Proline", "Trehalose", "Glucose")

# Fetch uniprot annotation
Uni_anno <- protti::fetch_uniprot_proteome(organism_id = "83333", columns=c("id", "protein names"))
colnames(Uni_anno) <- c("Protein", "Protein_name")

# Create workbook
wb<-createWorkbook(type="xlsx")

TABLE_ROWNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE)
TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE) +
  Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") +
  Border(color="black", position=c("TOP", "BOTTOM"), 
         pen=c("BORDER_THIN", "BORDER_THICK")) 

TITLE_STYLE <- CellStyle(wb)+ Font(wb,  heightInPoints=16, 
                                  isBold=TRUE, underline=1)
SUB_TITLE_STYLE <- CellStyle(wb) + 
  Font(wb,  heightInPoints=14, 
       isItalic=TRUE, isBold=FALSE)

# Styles for the data table row/column names
TABLE_ROWNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE)
TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE) +
  Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") +
  Border(color="black", position=c("TOP", "BOTTOM"), 
         pen=c("BORDER_THIN", "BORDER_THICK")) 

xlsx.addTitle<-function(sheet, rowIndex, title, titleStyle){
  rows <-createRow(sheet,rowIndex=rowIndex)
  sheetTitle <-createCell(rows, colIndex=1)
  setCellValue(sheetTitle[[1,1]], title)
  setCellStyle(sheetTitle[[1,1]], titleStyle)
}

# Create stabilisation workbook for all osmolytes
for(sm in sm_list[1:6]){
  
  sheet <- createSheet(wb, sheetName = sm)
  xlsx.addTitle(sheet, rowIndex=1, title=paste(sm, " stabilisation", sep=""),
                titleStyle = TITLE_STYLE)
  
  Protein_df <- significant_all[[sm]]$Protein_level %>%
    plyr::join(., Uni_anno, by="Protein") %>%
    dplyr::select(Protein, Protein_name, Protein_stabilisation, n_peptide, Coverage) %>%
    filter(Protein_stabilisation != 0)
  
  colnames(Protein_df) <- c("UniprotID", "Protein_name", "Stabilisation_score", "N_peptides", "Coverage")
  
  addDataFrame(Protein_df, sheet, startRow=3, startColumn=2, row.names = FALSE, 
               colnamesStyle = TABLE_COLNAMES_STYLE)
  setColumnWidth(sheet, colIndex=c(1:(ncol(Protein_df)+1)), colWidth=16)
  
}

# Add correlation data
output_3_export <- output_3 %>% 
  filter(n_stabilised >= 4) %>%
  plyr::join(., Uni_anno, by="Protein") %>%
  dplyr::select(Protein, Protein_name, cor, n_stabilised) %>%
  `colnames<-`(c("Protein", "Protein_name", "Spearman_correlation", "N_stabilised")) %>%
  data.frame()

sheet <- createSheet(wb, sheetName = "Correlation")
xlsx.addTitle(sheet, rowIndex=1, title="Spearman correlation to proteome",
              titleStyle = TITLE_STYLE)

addDataFrame(output_3_export, sheet, startRow=3, startColumn=2, row.names = FALSE, 
             colnamesStyle = TABLE_COLNAMES_STYLE)
setColumnWidth(sheet, colIndex=c(1:(ncol(output_3_export)+1)), colWidth=16)

# Supplementary table 1
saveWorkbook(wb, "/Users/moni/Documents/Phd/Osmolyte_paper/Sup/Stabilisation_scores.xlsx")

# Create workbook 2
wb2 <- createWorkbook(type="xlsx")

TABLE_ROWNAMES_STYLE <- CellStyle(wb2) + Font(wb2, isBold=TRUE)
TABLE_COLNAMES_STYLE <- CellStyle(wb2) + Font(wb2, isBold=TRUE) +
  Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") +
  Border(color="black", position=c("TOP", "BOTTOM"), 
         pen=c("BORDER_THIN", "BORDER_THICK")) 
TITLE_STYLE <- CellStyle(wb2)+ Font(wb2,  heightInPoints=16, 
                                   isBold=TRUE, underline=1)
SUB_TITLE_STYLE <- CellStyle(wb2) + 
  Font(wb2,  heightInPoints=14, 
       isItalic=TRUE, isBold=FALSE)
# Styles for the data table row/column names
TABLE_ROWNAMES_STYLE <- CellStyle(wb2) + Font(wb2, isBold=TRUE)
TABLE_COLNAMES_STYLE <- CellStyle(wb2) + Font(wb2, isBold=TRUE) +
  Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") +
  Border(color="black", position=c("TOP", "BOTTOM"), 
         pen=c("BORDER_THIN", "BORDER_THICK")) 

for(sm in sm_list[1:6]){
  
  sheet <- createSheet(wb2, sheetName = paste(sm, "_binding", sep=""))
  xlsx.addTitle(sheet, rowIndex=1, title=paste(sm, " binding", sep=""),
                titleStyle = TITLE_STYLE)
  
  Protein_df <- output_binding[[sm]]$Significant %>%
    dplyr::select(PG.ProteinAccessions, PG.ProteinDescriptions, PEP.StrippedSequence, qvalue, FC) %>%
    data.frame()

  colnames(Protein_df) <- c("UniprotID", "Protein_name", "Peptide", "adj. p-value", "log2(FC)")
  
  addDataFrame(Protein_df, sheet, startRow=3, startColumn=2, row.names = FALSE, 
               colnamesStyle = TABLE_COLNAMES_STYLE)
  setColumnWidth(sheet, colIndex=c(1:(ncol(Protein_df)+1)), colWidth=16)
  
}
saveWorkbook(wb2, "/Users/moni/Documents/Phd/Osmolyte_paper/Sup/Binding_analysis.xlsx")

# Create workbook 3
wb3 <- createWorkbook(type="xlsx")

TABLE_ROWNAMES_STYLE <- CellStyle(wb3) + Font(wb3, isBold=TRUE)
TABLE_COLNAMES_STYLE <- CellStyle(wb3) + Font(wb3, isBold=TRUE) +
  Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") +
  Border(color="black", position=c("TOP", "BOTTOM"), 
         pen=c("BORDER_THIN", "BORDER_THICK")) 
TITLE_STYLE <- CellStyle(wb3)+ Font(wb3,  heightInPoints=16, 
                                    isBold=TRUE, underline=1)
SUB_TITLE_STYLE <- CellStyle(wb3) + 
  Font(wb3,  heightInPoints=14, 
       isItalic=TRUE, isBold=FALSE)
# Styles for the data table row/column names
TABLE_ROWNAMES_STYLE <- CellStyle(wb3) + Font(wb3, isBold=TRUE)
TABLE_COLNAMES_STYLE <- CellStyle(wb3) + Font(wb3, isBold=TRUE) +
  Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") +
  Border(color="black", position=c("TOP", "BOTTOM"), 
         pen=c("BORDER_THIN", "BORDER_THICK")) 
sm_TPP <- c("TMAO", "glucose", "proline")
for(sm in sm_TPP){
  name <- ifelse(sm == "glucose", "Glucose",
                 ifelse(sm == "proline", "Proline", 
                        ifelse(sm == "TMAO", "TMAO", NA)))
  
  sheet <- createSheet(wb3, sheetName = paste(name, "_TPP_stabilisation", sep=""))
  xlsx.addTitle(sheet, rowIndex=1, title=paste(name, " TPP stabilisation", sep=""),
                titleStyle = TITLE_STYLE)
  
  Protein_df <- TPP_significant[[sm]]$Stabilisation %>%
    `colnames<-`(c("Protein", "Stabilisation_score")) %>%
    plyr::join(., Uni_anno, by="Protein") %>%
    dplyr::select(Protein, Protein_name, Stabilisation_score) %>%
    filter(Stabilisation_score != 0) %>%
    na.omit() %>%
    data.frame()
  
  addDataFrame(Protein_df, sheet, startRow=3, startColumn=2, row.names = FALSE, 
               colnamesStyle = TABLE_COLNAMES_STYLE)
  setColumnWidth(sheet, colIndex=c(2:(ncol(Protein_df)+1)), colWidth=16)

  #write.xlsx(Protein_df, file="/Users/moni/Documents/Phd/Osmolyte_paper/Sup/Stabilisation_file.xlsx", sheetName=paste(sm, "_stabilisaiton", sep=""),append=TRUE, row.names=FALSE)
  
}

sheet <- createSheet(wb3, sheetName = "Aggregation_TPP")
xlsx.addTitle(sheet, rowIndex=1, title="Aggregation from TPP",
              titleStyle = TITLE_STYLE)
addDataFrame(comparisons, sheet, startRow=3, startColumn=2, row.names = FALSE, 
             colnamesStyle = TABLE_COLNAMES_STYLE)
setColumnWidth(sheet, colIndex=c(2:(ncol(Protein_df)+1)), colWidth=16)

for(sm in sm_list[1:6]){
  
  sheet <- createSheet(wb3, sheetName = paste(sm, "_susceptibility", sep=""))
  xlsx.addTitle(sheet, rowIndex=1, title=paste(sm, " susceptibility", sep=""),
                titleStyle = TITLE_STYLE)
  
  Protein_df <- comparisons_LIP_all[[sm]] %>%
    plyr::join(., Uni_anno, by="Protein") %>%
    dplyr::select(Protein, Protein_name, Peptide, qvalue, log_FC) %>%
    na.omit() %>%
    data.frame()
  
  colnames(Protein_df) <- c("UniprotID", "Protein_name", "Peptide", "adj. p-value", "log2(FC)")
  
  addDataFrame(Protein_df, sheet, startRow=3, startColumn=2, row.names = FALSE, 
               colnamesStyle = TABLE_COLNAMES_STYLE)
  setColumnWidth(sheet, colIndex=c(2:(ncol(Protein_df)+1)), colWidth=16)
    #write.xlsx(Protein_df, file="/Users/moni/Documents/Phd/Osmolyte_paper/Sup/Stabilisation_file.xlsx", sheetName=paste(sm, "_stabilisaiton", sep=""),append=TRUE, row.names=FALSE)
  
}
saveWorkbook(wb3, "/Users/moni/Documents/Phd/Osmolyte_paper/Sup/Aggregation_analysis.xlsx")

