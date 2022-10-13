#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 10:19:13 2021

@author: Monika Pepelnjak
"""
from scipy.optimize import curve_fit
import pandas as pd
import math
import numpy as np
from matplotlib import pyplot
from mpmath import *

# Path with the raw data from CD experiment
path = "/Users/moni/Documents/Phd/Collaborations/C6/CD/20210903_Frr_thermaltrans.csv"

dataset = pd.read_csv(path)

# The two conditions analysed
interesting_positions = ["Control", "TMAO"]

temperature = dataset["Temperature"]
# Only select specific temperature range
keep_temperature = np.logical_and(temperature >= 37, temperature <= 76)
# Transform the temperature into K
x = temperature[keep_temperature] + 273.15
R = 8.31446261815324

# Define the fitting function
def objective(x, a, b, c, d, H, Tm):
#	return  ((a+b*x) + ((c + d*x) * np.exp((math.log10(H)*(x-Tm))/x*R*Tm)))/(1 + np.exp((math.log10(H)*(x-Tm))/x*R*Tm))
	return  ((a+b*x) + ((c + d*x) * np.exp(((H)*(x-Tm))/x*R*Tm)))/(1 + np.exp(((H)*(x-Tm))/x*R*Tm))

def objective_2(x, a, b, c, d, e, g, h, i, k):
 #  return (a+b*x+math.exp(((math.log10(i)/(8.3145*h))*((x- h)/h)))*(k+(c+d*x)* math.exp((math.log10(g)/(8.3145*e))*((x-e)/e))))/(1+math.exp((math.log10(i)/(8.3145*h))*((x-h)/h))* (1+math.exp((math.log10(g)/( 8.3145*e))*((x-e)/e))))
   return (a+b*x+np.exp((((i)/(8.3145*h))*((x- h)/h)))*(k+(c+d*x)* np.exp(((g)/(8.3145*e))*((x-e)/e))))/(1+np.exp(((i)/(8.3145*h))*((x-h)/h))* (1+np.exp(((g)/( 8.3145*e))*((x-e)/e))))

myfunction2= np.vectorize(objective)
myfunction3= np.vectorize(objective_2)

Tm_out=[]
output_fit_frame = []
output_y_frame = []
H_out = []

for i in interesting_positions:
    y = dataset[i][keep_temperature]
    y = (y-min(y))/(max(y)-min(y)) # scale the data between 0 and 1
    pyplot.scatter(x, y)
    
    y_condition = np.logical_and(y<0.55,y>0.45)
    border_x = np.mean(x[y == max(y)]).item()
    x_condition = x < border_x
    mid_temp = np.mean(x[np.logical_and(y_condition, x_condition)])

    # fit curve
    popt, _ = curve_fit(myfunction2, x, y, p0 = [0,-1,1,1,1,mid_temp], maxfev = 30000)

    a, b, c, d, H, Tm = popt

    second_try, _ = curve_fit(myfunction2, x, y, p0 = [a, b, c, d, H, Tm])
    
    a_2, b_2, c_2, d_2, H_2, Tm_2 = second_try


    x_line = np.arange(min(x), max(x), 0.1)
    y_line = myfunction2(x_line,a_2, b_2, c_2, d_2, H_2, Tm_2)
    pyplot.plot(x_line, y_line, '--', color='orange')    #pyplot.show()
    
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

      
    #y-predictions
output_fit_frame_conc = pd.concat(output_fit_frame, axis=1)
output_y_conc = pd.concat(output_y_frame, axis=1)
output_parameters = {"Position":interesting_positions, "Tm":Tm_out, "H":H_out}
output_parameters = pd.DataFrame(data =output_parameters )
output_fit_frame_conc.to_csv("/Users/moni/Documents/Phd/Collaborations/C6/CD/" +"Frr_Fitted_values.csv" )
output_y_conc.to_csv("/Users/moni/Documents/Phd/Collaborations/C6/CD/" + "Frr_scaled_y_values.csv" )
output_parameters.to_csv("/Users/moni/Documents/Phd/Collaborations/C6/CD/" + "Frr_output_parameters.csv" )

