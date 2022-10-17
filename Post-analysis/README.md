Post-analysis:
  -

Post-analysis contains scripts that were used to analyse the data after the significantly changing proteins and peptides were already identified.

- [Clustering_analysis_figures.R](https://github.com/MoniPepelnjak/Thermal_unfolding/blob/master/Post-analysis/Clustering_analysis_figures.R) is the script used in the paper to analyse the shape profiles of LiP-MS thermal unfolding data, analyse the distribution of profiles for HT and FT peptides and for precipitators and non-precipitators and to check the changes of shape profiles upon addition of osmolytes. 
  - The following figures are produced with this script:
    - 1C, 1D, 1F
    - S1G
    - S2A
  - The clustering for all osmolytes takes a long time, therefore the script is adapted in a way that either the clustering is performed or the files previously produced can be uploaded.
  - Several parts of the script are hard coded for the specific dataset used: Trehalose. This is the example used in the paper, however the script could be adapted in a way such that it is more automated.
  - The following supporting scripts are needed in the file:


- [FT_HT_correlations.R](https://github.com/MoniPepelnjak/Thermal_unfolding/blob/master/Post-analysis/FT_HT_correlations.R) script matches the HT peptides with its FT "parent peptide". The script requires the supportive script [FT_peptide_function.R](https://github.com/MoniPepelnjak/Thermal_unfolding/blob/master/Post-analysis/Supportive_functions/FT_peptide_function.R) with defined functions. 
  - The following figures are produced with this script:
    - S1A, S1B
  - The script applies filtering criteria for analysis, only the peptides with good enough fit are considered. This was applied after seeing that otherwise the data is too noisy. 
  - Script can further be automated to include several experiments. Currently only one osmolyte experiment at the time is considered.
- [Missing_analysis_figures.R](https://github.com/MoniPepelnjak/Thermal_unfolding/blob/master/Post-analysis/Missing_analysis_figures.R) script contains the remaining analysis for the figures 2 and 3 of the paper. Besides the plotting functions, the following analysis is found:
  - The following figures are produced by this script:
    - 2B, 2C
    - 3A, 3B, 3C, 3D, 3G
    - 4A, 4B, 4C, 4E
    - S4A, S4B, S4C
  - Randomised correlation analysis between individual proteins and the mean LiP stabilisation score
  - Calculation of the "best stabilisers" and number of stabilising conditions
  - Overlap calculation for figure 3A
- [221007_Aggregation_figure_analysis.R](https://github.com/MoniPepelnjak/Thermal_unfolding/blob/master/Post-analysis/221007_Aggregation_figure_analysis.R) performs the analysis of protein aggregation data. Significant changes for LiP and for TPP data should be uploaded for analysis.
  - The following figures are produced by this script:
    - 5B, 5C, 5D, 5E, 5F, 5H
    - S5A, S5B, S5C, S5D, S5E
  - The remaining figures are connected to validations and the scripts for those figures can be found in [Validations](https://github.com/MoniPepelnjak/Thermal_unfolding/tree/master/Validations)
