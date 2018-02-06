# -*- coding: utf-8 -*-
"""
Created on Sat Dec 16 17:12:43 2017

@author: Lana
"""

import glob as g
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import lmfit
import random
import scipy.stats

data = countData[countData['Dilution'] == 1]
dataA = data[data['Name']== '0mM']


x = list(dataA['Day'])
y = list(dataA['Count'])

x = x[0:2]
y = y[0:2]

m, b, r, p, stdderr = scipy.stats.linregress(x, y=y)

# Plot the data and fitted model
plt.figure(0)
plt.plot(x, y, 'k.')
plt.plot(x, m*np.array(x)+b, 'r-')
plt.xlabel(r'Day')
plt.ylabel(r'Cell Count/ml')
plt.title(r'Linear regression of Cell count Timeseries')
plt.show()

print 'equation is', m,'*x+',b
print 'R is', r
print 'p-value is', p

print 'doubling time is ', b/m
