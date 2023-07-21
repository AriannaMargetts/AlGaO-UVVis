# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 15:51:36 2023

@author: Ari
"""

#imports
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import stats


# retrieving data from csv file
# INPUT DATA FILE NAME IN LINE BELOW
df = pd.read_csv("GSW8_1UVVis.csv").values


#variables
absorbance = np.zeros(len(df))
hv = np.zeros(len(df))
absorption_coefficient = np.zeros(len(df))
ahv2 = np.zeros(len(df))
scaled_ahv2 = np.zeros(len(df))
film_thickness_input = float(input("What is the film thickness"
                                   " in nanometres? "))
film_thickness = film_thickness_input / 10000000
sample_name = input("What is the sample name? ")

### DATA COLLECTING


# calculating values needed for Tauc plot
for i in range(len(df)):
    if df[i,1] <= 0 :
        absorbance[i] == 0
    else:
        absorbance[i] = -math.log10(df[i, 1]/100)
    hv[i] = 1240 / df[i, 0]     
    absorption_coefficient[i] = 2.303 * (absorbance[i] / film_thickness)
    ahv2[i] = (absorption_coefficient[i] * hv[i]) **2
    scaled_ahv2[i] = ahv2[i] / (10**9)

# filtering out zero values for Tauc plot
index = 0 
idx = []
for scaled_ahv2[i] in scaled_ahv2:
    if scaled_ahv2[i] == 0:
        idx.append(index)
    index += 1

# deleting 0 values    
scaled_ahv2 = np.delete(scaled_ahv2, idx)
hv = np.delete(hv, idx)

# deleting NaN values
scaled_ahv2 = scaled_ahv2[~np.isnan(scaled_ahv2)]
hv = hv[~np.isnan(hv)]

# plotting Tauc plot
fig, ax = plt.subplots(figsize = (10,8)) 
plt.scatter(hv, scaled_ahv2,color='#950757', marker = '.')
plt.xlim(left=4.5, right=hv[len(hv) - 1] + 0.25)
plt.ylim(bottom=-0.05*np.max(scaled_ahv2), top=np.max(scaled_ahv2) + 0.05*np.max(scaled_ahv2))
plt.title("Tauc plot of {}" .format(sample_name))
plt.xlabel("hv, eV")
plt.ylabel(r"$\alpha hv^2 /10^9$")
plt.show() 


### DATA ANALYSIS


# creating baseline measurements
baseline_array = np.zeros(1)

for i in range(len(hv)): # baseline defined as anything below 4.7 eV
    if hv[i] <= 4.7:
        baseline_array = np.append(baseline_array, scaled_ahv2[i])
        
baseline_array = np.delete(baseline_array, 0)
twenty_standard_deviations = 20 * np.std(baseline_array)
baseline_average = np.average(baseline_array)

working_scaled_ahv2_data = np.zeros(len(scaled_ahv2))
working_hv_data = np.zeros(len(scaled_ahv2))

# finding data above baseline + 20 sigma, and between 4.7 to 5.5 eV
for i in range(len(working_scaled_ahv2_data)): 
    if 4.7 <= hv[i] <= 5.5 and scaled_ahv2[i] >= baseline_average + twenty_standard_deviations:
        working_scaled_ahv2_data[i] = scaled_ahv2[i]
        working_hv_data[i] = hv[i]

# filtering out zero values from each side of working data arrays
index = 0 
idx = []
for working_scaled_ahv2_data[i] in working_scaled_ahv2_data:
    if working_scaled_ahv2_data[i] == 0:
        idx.append(index)
    index += 1
    
working_scaled_ahv2_data = np.delete(working_scaled_ahv2_data, idx)
working_hv_data = np.delete(working_hv_data, idx)

graph_ahv_data = working_scaled_ahv2_data
graph_hv_data = working_hv_data

#functions

def take_from_right(working_scaled_ahv2_data, working_hv_data, k):
    '''This function takes the data currently being analysed in both hv and 
    ahv2 arrays, and a step k to reduce the number of elements by. The last k
    points of data are removed, and the new data set replaces the current data 
    set if Pearson's r-value is larger than that of the previous data set.'''
    for i in range(len(hv)):
        if len(working_hv_data) > 20:
            working_r_value = stats.pearsonr(working_hv_data, working_scaled_ahv2_data)[0]
            k = k
            test_ahv2_data = working_scaled_ahv2_data[: len(working_scaled_ahv2_data) - k]
            test_hv_data = working_hv_data[: len(working_hv_data) - k]
            if len(test_hv_data) >= 2:
                test_r_value = stats.pearsonr(test_hv_data, test_ahv2_data)[0]
                if test_r_value > working_r_value:
                    working_hv_data = test_hv_data
                    working_scaled_ahv2_data = test_ahv2_data
                else:
                    break
                break
            break

    return working_hv_data, working_scaled_ahv2_data

def take_from_left(working_scaled_ahv2_data, working_hv_data, k):
    '''This function takes the data currently being analysed in both hv and 
    ahv2 arrays, and a step k to reduce the number of elements by. The first k
    points of data are removed, and the new data set replaces the current data 
    set if Pearson's r-value is larger than that of the previous data set.'''
    for i in range(len(hv)):
        if len(working_hv_data) > 20:
            working_r_value = stats.pearsonr(working_hv_data, working_scaled_ahv2_data)[0]
            k = k
            test_ahv2_data = working_scaled_ahv2_data[k:]
            test_hv_data = working_hv_data[k:]
            if len(test_hv_data) >= 2:
                test_r_value = stats.pearsonr(test_hv_data, test_ahv2_data)[0]
                if test_r_value > working_r_value and len(test_hv_data) >= k+1:
                    working_hv_data = test_hv_data
                    working_scaled_ahv2_data = test_ahv2_data
                else:
                    break
                break
            break
    
    return working_hv_data, working_scaled_ahv2_data

def take_from_both(working_scaled_ahv2_data, working_hv_data, k):
    '''This function takes the data currently being analysed in both hv and 
    ahv2 arrays, and a step k to reduce the number of elements by. The first
    and last k points of data are removed, and the new data set replaces the 
    current data set if Pearson's r-value is larger than that of the previous
    data set.'''
    for i in range(len(hv)):
        if len(working_hv_data) > 20:
            working_r_value = stats.pearsonr(working_hv_data, working_scaled_ahv2_data)[0]
            k = k
            test_ahv2_data = working_scaled_ahv2_data[k:len(working_scaled_ahv2_data) - k]
            test_hv_data = working_hv_data[k:len(working_scaled_ahv2_data) - k]
            if len(test_hv_data) >= 2:
                test_r_value = stats.pearsonr(test_hv_data, test_ahv2_data)[0]
                if test_r_value > working_r_value and len(test_hv_data) >= (2*k)+2:
                    working_hv_data = test_hv_data
                    working_scaled_ahv2_data = test_ahv2_data
                else:
                    break
                break
            break
    
    return working_hv_data, working_scaled_ahv2_data
    

# reducing data points used for linearisation using functions take_from_right() and take_from_left()
# data removed in chunks of 20, then 15, 10, then 5
# if data isn't reducing how you desire, try uncommenting the 2 lines below

#working_hv_data, working_scaled_ahv2_data = take_from_right(working_scaled_ahv2_data, working_hv_data, 40)
#working_hv_data, working_scaled_ahv2_data = take_from_left(working_scaled_ahv2_data, working_hv_data, 40)
working_hv_data, working_scaled_ahv2_data = take_from_right(working_scaled_ahv2_data, working_hv_data, 20)
working_hv_data, working_scaled_ahv2_data = take_from_left(working_scaled_ahv2_data, working_hv_data, 20)
working_hv_data, working_scaled_ahv2_data = take_from_right(working_scaled_ahv2_data, working_hv_data, 15)
working_hv_data, working_scaled_ahv2_data = take_from_left(working_scaled_ahv2_data, working_hv_data, 15)
working_hv_data, working_scaled_ahv2_data = take_from_right(working_scaled_ahv2_data, working_hv_data, 10)
working_hv_data, working_scaled_ahv2_data = take_from_left(working_scaled_ahv2_data, working_hv_data, 10)
working_hv_data, working_scaled_ahv2_data = take_from_right(working_scaled_ahv2_data, working_hv_data, 5)
working_hv_data, working_scaled_ahv2_data = take_from_left(working_scaled_ahv2_data, working_hv_data, 5)

# reducing data points used for linearisation using function take_from_both() 
working_hv_data, working_scaled_ahv2_data = take_from_both(working_scaled_ahv2_data, working_hv_data, 15)
working_hv_data, working_scaled_ahv2_data = take_from_both(working_scaled_ahv2_data, working_hv_data, 5)

# calculating bandgap
gradient, intercept = np.polyfit(working_hv_data, working_scaled_ahv2_data, 1)
bandgap = round(-intercept / gradient, 3)
print("The bandgap of {0} is {1} eV. The linearised data used to obtain"
      " this value is shown in dark purple in the second graph." .format(sample_name, bandgap))

# calculating aluminium concentration
al_conc_miller = round((bandgap - 4.857) / 1.54, 4)*100
print("The aluminium concentration of {0} is {1}% using Miller et al.'s formula." .format(sample_name, al_conc_miller))
al_conc_agnitron = round((bandgap - 4.857) / 1.7, 4)*100
print("The aluminium concentration of {0} is {1}% using Agnitron's formula." .format(sample_name, al_conc_agnitron))

# plotting Tauc plot showing linear data
fig, ax = plt.subplots(figsize = (10,8)) 
plt.scatter(hv, scaled_ahv2, color='hotpink', marker = '.')
plt.scatter(working_hv_data, working_scaled_ahv2_data, color='#950757', marker = '.')
plt.axline((0,intercept), slope=gradient, color='orangered')
plt.xlim(left=(((-baseline_average*twenty_standard_deviations)-intercept)/gradient)-0.5, right=5.5)
plt.ylim(bottom=-5, top=np.max(working_scaled_ahv2_data)+10)
# -baseline_average*twenty_standard_deviations/4
plt.title("Tauc plot of {}" .format(sample_name))
plt.xlabel("hv, eV")
plt.ylabel(r"$\alpha hv^2 /10^9$")
plt.show() 

