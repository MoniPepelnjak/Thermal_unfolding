# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 18:17:54 2021

@author: pmonika
"""

import torch
import gpytorch
import pandas as pd
import numpy as np
from gpytorch.mlls import SumMarginalLogLikelihood
import matplotlib.pyplot as plt
# from matplotlib import cm
# import matplotlib.colors as colors
import scipy.stats as ss
# import seaborn as sns
from fitGPs2_function import FitGPs_function
import sys
# TODO make an interface to R package here

###################
### Data import ###
###################

# replace with the path to your data
# the input needs to be a table with columns:
#   - uniqueID : peptide-ID
#   - x : temperature
#   - y : fraction non-denaturated
#   - condition : treatment/control
# data = pd.read_csv("Example_peptides_for_fit.csv")
#data = pd.read_csv("~/Documents/Phd/Experiments/018/python/HT_ATP_all.csv")    # example data from the TPP NPARC package


path_list = pd.read_csv(sys.argv[1] + "Path_directions.csv")
sm_list = [sys.argv[3]]
for sm in sm_list:
    FitGPs_function(sm, path_list, sys.argv[4])
