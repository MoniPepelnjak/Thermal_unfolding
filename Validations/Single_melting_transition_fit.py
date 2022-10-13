#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 10:19:13 2021

@author: moni
"""

# Script to fit the data of thermal unfolding from DSF
from scipy.optimize import curve_fit
import pandas as pd
import math
import numpy as np
from matplotlib import pyplot
from mpmath import *

# Define the experiment number and sample name
experiment = "11_04"
path = "/Users/moni/Documents/Phd/Collaborations/C6/Thermal_shift/"
sample = "Eno"

# The datasets need to contain the experiment name 
dataset = pd.read_csv(path + experiment + "_RFU_data.csv")
annotation = pd.read_csv(path + experiment + "_annotation.csv")
keep = np.logical_and((annotation["Condition"] != "Negative"), (annotation["Sample"] ==  sample))

interesting_positions = annotation["Position"][keep]

# Transform to kelvins
temperature = dataset["Temperature"]
keep_temperature = np.logical_and(temperature >= 37, temperature <= 76)
x = temperature[keep_temperature] + 273.15
R = 8.31446261815324

# Define the function
def objective(x, a, b, c, d, H, Tm):
	return  ((a+b*x) + ((c + d*x) * np.exp(((H)*(x-Tm))/x*R*Tm)))/(1 + np.exp(((H)*(x-Tm))/x*R*Tm))

def objective_2(x, a, b, c, d, e, g, h, i, k):
   return (a+b*x+np.exp((((i)/(8.3145*h))*((x- h)/h)))*(k+(c+d*x)* np.exp(((g)/(8.3145*e))*((x-e)/e))))/(1+np.exp(((i)/(8.3145*h))*((x-h)/h))* (1+np.exp(((g)/( 8.3145*e))*((x-e)/e))))

myfunction2= np.vectorize(objective)
myfunction3= np.vectorize(objective_2)

Tm_out=[]
output_fit_frame = []
output_y_frame = []
H_out = []
#interesting_positions = interesting_positions[:2]
for i in interesting_positions:
    y = dataset[i][keep_temperature]
    y = (y-min(y))/(max(y)-min(y))
    pyplot.scatter(x, y)
    
    y_condition = np.logical_and(y<0.55,y>0.45)
    border_x = np.mean(x[y == max(y)]).item()
    x_condition = x < border_x
    mid_temp = np.mean(x[np.logical_and(y_condition, x_condition)])

    # fit curve
    popt, _ = curve_fit(myfunction2, x, y, p0 = [0,-1,1,1,1,mid_temp], maxfev = 30000)
#    popt_2, _ = curve_fit(myfunction3, x, y, p0 = [max(y),-1,min(y),0,335,100000,320,10000,40], maxfev=500000)

    a, b, c, d, H, Tm = popt

#    a, b, c, d, e, g, h, i, k = popt_2
    second_try, _ = curve_fit(myfunction2, x, y, p0 = [a, b, c, d, H, Tm])
    
    a_2, b_2, c_2, d_2, H_2, Tm_2 = second_try


    x_line = np.arange(min(x), max(x), 0.1)
    y_line = myfunction2(x_line,a_2, b_2, c_2, d_2, H_2, Tm_2)
#    y_line_2 = myfunction3(x_line, a, b, c, d, e, g, h, i, k)

    #pyplot.scatter(x, y)
    #pyplot.plot(x_line, y_line_2, '--', color='red')
    pyplot.plot(x_line, y_line, '--', color='orange')
    #pyplot.show()
    
    y_line_df = pd.DataFrame(y_line)
    y_line_df = y_line_df.rename(columns = {0: i})
    y_line_df = y_line_df.set_index(x_line)
    
    scaled_y = pd.DataFrame(y)
    scaled_y = scaled_y.rename(columns = {0: i})
    scaled_y = scaled_y.set_index(x)

    
    Tm_out.append(Tm_2)
    H_out.append(H)
    output_fit_frame.append(y_line_df)
    output_y_frame.append(scaled_y)

    
output_fit_frame_conc = pd.concat(output_fit_frame, axis=1)
output_y_conc = pd.concat(output_y_frame, axis=1)
output_parameters = {"Position":interesting_positions, "Tm":Tm_out, "H":H_out}
output_parameters = pd.DataFrame(data =output_parameters )
output_fit_frame_conc.to_csv(path + "Fitted_data/" + experiment + "_" + sample + "_Fitted_values.csv" )
output_y_conc.to_csv(path + "Fitted_data/" + experiment + "_" + sample + "_scaled_y_values.csv" )
output_parameters.to_csv(path + "Fitted_data/" + experiment + "_" + sample + "_output_parameters.csv" )

sub_int = interesting_positions[0:9]
for i in sub_int:
    y = dataset[i][keep_temperature]
    #y = (y-min(y))/(max(y)-min(y))
    pyplot.scatter(x, y)
