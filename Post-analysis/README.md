Post-analysis:

Post-analysis contains scripts that were used to analyse the data after the significantly changing proteins and peptides were already identified.

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
