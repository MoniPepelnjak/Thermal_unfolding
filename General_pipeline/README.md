# General pipeline

This folder contains the scripts that are required to analyse the LiP-MS thermal unfolding data from Spectronaut file to calculating the stabilisation scores on peptide and protein level.

[Main_script_loading.R](Main_script_loading.R) is the main script for the analysis of the thermal unfolding data, combining the different steps of the analysis.

## Required files for the analysis:
  - **Path file:** file that defines the experiments number, small molecule analysed and the paths to most important files
    - The important files to include in the path file:
      - Spectronaut export file - exported in long format, as .tsv
      - Annotation file - in .csv format
      
  - **Spectronaut export file** with the following columns:
    - R.FileName
    - PG.FastaFile
    - PG.ProteinAccessions
    - PEP.IsProteotypic
    - PEP.StrippedSequence
    - EG.ModifiedSequence
    - PEP.DigestType....TrypsinP.
    - PEP.Quantity
    
  - **Annotation file** with the following information:
    - R.FileName: Exactly the same file name as the raw file imported in Spectronaut
    - variable1: The temperature of experiment
    - variable2: The condition analysed --> the variable2 name should match the small molecule and sm name in the path file!!
    - replicate: The replicate number of the sample
  - **Fasta file of your organism**
  
  - Potential formatting issues:
    - Currently only analysis with spectronaut is supported, as some column names are hard-coded in the functions. This could be further automated if needed.
    - All the csv files were created with the european standards (values separated by , and not by ;).
    - The code was writen with macOS, therefore some incompatibility with Windows might be an issue.
    
  ## Required scripts and analysis steps
  
The analysis can be performed by running the script [Main_script_loading.R](Main_script_loading.R) and uploading the path file, defining the small molecules that you wish to analyze and a common path where you want to save your data.

The script includes the following steps and supporting functions that are required:
  - [Spectronaut_to_python_function.R](Spectronaut_to_python_function.R) Loads in the data exported from Spectronaut, based on the defined path in the Path_file, filters and scales the data to prepare the data for GP fitting.
  - After the data has been successfully exported, the data has to be fit with python, as decribed in [python_fitting](python_fitting) folder. The [Spectronaut_to_python_function.R](Spectronaut_to_python_function.R) already creates the right folder structure such that running the scripts on Euler requires little optimisation.
  - After the data has been fit, the [python_to_list_function.R](python_to_list_function.R) imports the data from the python fitting and creates a single list from all the conditions analysed.
  - Data is then clustered using fuzzy clustering (function: fuzzy_clusters_two_conditions)that is defined in [fuzzy_clustering.R](fuzzy_clustering.R) script. The fuzzy clustering is performed in order to capture a general shape of the curve. This information is important to determine whether the observed change is stabilising or destabilising. 
  - In the next step, we determine whether the peptides are significantly stabilised or destabilised. This is achieved by (find_significant_peptides.R)[find_significant_peptides.R]. Within the script, we firstly calculate the area between two curves using (area_calculation.R)[area_calculation.R] script and in the second step determine whether the change is associated with binding, aggregation or stabilisation/destabilisation. 
  - Peptides are then mapped to proteins, to find the exact position in the protein using the [peptide_fasta_matching.R](peptide_fasta_matching.R). 
  - Lastly, the AA-level stabilisation scores and Protein-level stabilisation scores are calculated using score_aggregation function from [aggregate_scores.R](aggregate_scores.R). 

## Files exported from the analysis

The following two lists are exported in a format .RDS format.

- **all_results**: 
  - Combined list of all python fitting outputs, containing the fits, the goodness of fit and confidence interval information
- **significant_all**
  - Combined list of significant scores on peptide level, AA-level and protein level for all analysed small molecules
