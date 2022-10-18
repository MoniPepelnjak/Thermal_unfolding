Validations:
  -

This repository contains the collection of scripts used to analyse and plot the validation experiments of the Osmolyte paper.

- [CD_melting_transition_fit.py](https://github.com/MoniPepelnjak/Thermal_unfolding/blob/master/Validations/CD_melting_transition_fit.py) uses a thermodynamic model to fit the unfolding curve to the measurements obtained by CD. The function defined for fitting can be used to fit any sigmoidal transition with flexible ends (linear model at the start and end). For different proteins different starting parameters might have to be chosen, however the script at the moment is not yet fully automated.

- [CD_Frr_plotting.R](https://github.com/MoniPepelnjak/Thermal_unfolding/blob/master/Validations/CD_Frr_plotting.R) plots the CD data. It requires the raw values and fit data, created with [CD_melting_transition_fit.py](https://github.com/MoniPepelnjak/Thermal_unfolding/blob/master/Validations/CD_melting_transition_fit.py)
  - The following figure was produced with this script:
    - S5F
    
- [Frr_Bradford_quantification.R](https://github.com/MoniPepelnjak/Thermal_unfolding/blob/master/Validations/Frr_Bradford_quantification.R) plots and displays the Frr and Gnd aggregation validation data used in the Osmolyte paper in Figure 5 and Figure S5
  - The following figures were produced with this script:
    - 5G, 5I
    - S5G
    
- [Single_melting_transition_fit.py](https://github.com/MoniPepelnjak/Thermal_unfolding/blob/master/Validations/Single_melting_transition_fit.py) is essentially the same script as [CD_melting_transition_fit.py](https://github.com/MoniPepelnjak/Thermal_unfolding/blob/master/Validations/CD_melting_transition_fit.py), slightly more automated to directly fit the data obtained from the Differential scanning flourimetry.

- [RbsK_single_melting_point_analysis.R](https://github.com/MoniPepelnjak/Thermal_unfolding/blob/master/Validations/RbsK_single_melting_point_analysis.R) requires the data from  [Single_melting_transition_fit.py](https://github.com/MoniPepelnjak/Thermal_unfolding/blob/master/Validations/Single_melting_transition_fit.py) and plots the data for RbsK (figure 4C of the Osmolyte paper). The script also includes the function for fitting a Michaelis-Menten model curve.

- [Lysate_validations.R](https://github.com/MoniPepelnjak/Thermal_unfolding/blob/master/Validations/Lysate_validations.R) fits the sigmoidal model with flexible ends in R and plots the data for DSF data for lysates. The script was used only for the lysates, as model converged for lysate fits, but not for purified proteins. For fitting the model for different curves, the starting parameters need to be adjusted. Alternatively, I would suggest to rather use the [Single_melting_transition_fit.py](https://github.com/MoniPepelnjak/Thermal_unfolding/blob/master/Validations/Single_melting_transition_fit.py).
  - The following figures were produced with this script:
    - 2D, 2E
    - S2D
