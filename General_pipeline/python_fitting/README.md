Python fitting
  -
This folder contains scripts and files required to fit the thermal unfolding data in python. The script requires the scaled data that are produced by the [Main_script_loading.R] script.

In order to learn the temperature profiles for each peptide in different conditions, we used Gaussian processes (GP), which provide a probabilistic framework to learn a non-parametric relationship of peptide abundance to temperature. In detail, we used gpytorch version 1.4.2 with an ExactGP model choosing a constant mean function, a squared exponential kernel and a Gaussian likelihood. For each peptide, separate GP models for the peptide intensities in absence (control condition) and presence of an osmolyte (osmolyte condition)  as well as a joint model were defined and model hyperparameters were found by maximizing the sum marginal log-likelihood across all models using Adam optimizer with a learning rate of 0.1 and 1000 iterations. Based on the resulting posterior of the fit, predicted mean abundance profiles and confidence intervals based on 2 standard deviations around the mean were found for each peptide and condition. The residual sum of squares between the observed peptide intensities and the predicted intensities are calculated for each peptide and condition to assess the goodness of the fit.

The output of the script are the fitted curves, the confidence intervals and the information about the goodness of the fit. 

The script can be run in three different ways:
- Directly in R using a [Main_script_loading.R] script. 
  - Running this in R directly can be very inefficient. Therefore I would strongly recommend exporting the files from R and importing and running the script in one of the other ways.
- Running python script [] locally
  - This can be done when running a small number of conditions.
  - Running the script is easy, by defining the files and paths directly in the script and requires little adjustments.
  - Running one LiP experiment can take up to one day of running, so it is not ideal to run it on your local computer.
- Running the scripts on Euler:
  - This is by far the fastest approach when fitting multiple conditions, as the approach can be parallelised.
  - For this approach, you will need to have a single folder with:
    - The bash script [Running_all_compound.sh](https://github.com/MoniPepelnjak/Thermal_unfolding/blob/master/General_pipeline/python_fitting/Running_all_compound.sh) to run the files on euler
    - [Path_directions.csv](https://github.com/MoniPepelnjak/Thermal_unfolding/blob/master/General_pipeline/python_fitting/Path_directions.csv) file with defined experiment number and small molecule
    - [Main_script_fitGPs_2.py](https://github.com/MoniPepelnjak/Thermal_unfolding/blob/master/General_pipeline/python_fitting/Main_script_fitGPs_2.py) - the main script running script
    - [fitGPs2_function.py](https://github.com/MoniPepelnjak/Thermal_unfolding/blob/master/General_pipeline/python_fitting/fitGPs2_function.py) the script with defined fitting function
    - Location of the input files in a structured folder:
      - General_path/Experiment_number/Small_molecule/python_import/Small_molecule_trypticity.csv
      - An example for TMAO: General_path/023/TMAO/python_import/TMAO_FT.csv
      - Empty General_path/Experiment_number/Small_molecule/python_output/ folder has to be provided
