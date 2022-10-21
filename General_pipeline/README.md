# General pipeline

This folder contains the scripts that are required to analyse the LiP-MS thermal unfolding data from Spectronaut file to calculating the stabilisation scores on peptide and protein level.

[Main_script_loading.R](Main_script_loading.R) is the main script for the analysis of the thermal unfolding data, combining the different steps of the analysis.

## Required files for the analysis:
  - Path file: file that defines the experiments number, small molecule analysed and the paths to most important files
    - The important files to include in the path file:
      - Spectronaut export file - exported in long format, as .tsv
      - Annotation file - in .csv format
      
  - Spectronaut export file with the following columns:
    - R.FileName
    - PG.FastaFile
    - PG.ProteinAccessions
    - PEP.IsProteotypic
    - PEP.StrippedSequence
    - EG.ModifiedSequence
    - PEP.DigestType....TrypsinP.
    - PEP.Quantity
    
  - Annotation file with the following information:
    - R.FileName: Exactly the same file name as the raw file imported in Spectronaut
    - variable1: The temperature of experiment
    - variable2: The condition analysed --> the variable2 name should match the small molecule and sm name in the path file!!
    - replicate: The replicate number of the sample
  
  - Potential formatting issues:
    - Currently only analysis with spectronaut is supported, as some column names are hard-coded in the functions. This could be further automated if needed.
    - All the csv files were created with the european standards (values separated by , and not by ;).
    - The code was writen with macOS, therefore some incompatibility with Windows might be an issue.
    
    
  
  
    

  -[]
