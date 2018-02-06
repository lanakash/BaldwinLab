# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 16:19:07 2017

@author: Lana
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 16:15:52 2017

@author: Lana
"""

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import lmfit

def toFloat(pd):
    pd[' Ordinate (A)'] = np.array(pd[' Ordinate (A)'], dtype = float)



dataZ4 = pd.read_csv('C:\Users\Lana\Desktop\zinc_aug_25th.csv')
dataZ = pd.read_csv('C:\Users\Lana\Desktop\Lana Nitrate Zinc August-25-17 4_56 PM Pacific Daylight Time\Results Table.csv')
dataC = pd.read_csv('C:\Users\Lana\Desktop\Lana Nitrate Cad August-25-17 5_06 PM Pacific Daylight Time\Results Table.csv')


StdZ4 = dataZ4.loc[[0,1,2,3],'Abs']
StdZ = dataZ.loc[[0,1,2,3],' Ordinate (A)']
StdC = dataC.loc[[0,1,2,3,4,5],' Ordinate (A)']

Std = [StdZ, StdC]
StdName = ['C rep1', 'C rep2', 'C rep3', 'C sel', 'Z']

    
plt.figure(0)
plt.suptitle('Standards Aug 25th', size = 15)
plt.subplot(2,1,1)
x = [0,0.05,0.1,0.2]
y1 = np.array(StdZ4)
y1 = y1 - y1[0]
fit1 = np.polyfit(x, y1, 1)
fit1_fn = np.poly1d(fit1) 
core = np.corrcoef(x,y1)**2
plt.plot(x,y1, 'ro', x, fit1_fn(x), '-r', label = 'ZStd4, p = %.3f' % core[0,1])


x = [0,0.05,0.1,0.2]
y1 = np.array(StdZ)
y1 = y1 - y1[0]
fit1 = np.polyfit(x, y1, 1)
fit1_fn = np.poly1d(fit1) 
core = np.corrcoef(x,y1)**2
plt.plot(x,y1, 'bo', x, fit1_fn(x), '-b', label = 'ZStd, p = %.3f' % core[0,1])
plt.legend(loc = 'lower right')


plt.subplot(2,1,2)
x = [0,0.05,0.1,0.2,0.5,1.0]
y1 = np.array(StdC)
y1 = y1 - y1[0]
fit1 = np.polyfit(x, y1, 1)
fit1_fn = np.poly1d(fit1)
core = np.corrcoef(x,y1)**2 
plt.plot(x,y1, 'bo', x, fit1_fn(x), '-b', label = 'CStd, p = %.3f' % core[0,1])



plt.legend(loc = 'lower right')

#plt.figure(1)
#plt.suptitle('Nitrate Standards Plots')
#i = 1
#for stdType in Std:
#    plt.subplot(3,2,i)
#    x = stdType[' Concentration']
#    y = stdType[' Ordinate (A)']
#    fit = np.polyfit(x, y, 1)
#    fit_fn = np.poly1d(fit) 
#    plt.plot(x,y, 'ro', x, fit_fn(x), '-r', label = StdName[i-1])
#    plt.xlabel(r'Concentration Nitrate (mg/L)')
#    plt.ylabel(r'Absorbance (543nm)')
#    plt.legend(loc ='lower right')
#    i += 1

# plot concentration vs std deviation (make a residual plot and also a avg+std deviation)

#    for a in stdlist:
#        #plt.subplot(2,len(StdR),i)
#        fitavg = np.polyfit(conc, avg, 1)
#        fitavg_fn = np.poly1d(fitavg) 
#        plt.plot(conc, np.array(a[' Ordinate (A)']), rbg[l], conc, fitavg_fn(conc), rbg[l+1], label = StdRName[i-1])
#        #plt.plot(conc, np.array(a[' Ordinate (A)']), '.b')
#        plt.errorbar(conc, avg, yerr=stddev, fmt='.', color='k')    
#        plt.ylabel(r'Absorbance (543nm)')
#        plt.xlabel(r'Concentration Nitrate (mg/L)')
#        plt.title(StdRName[i-1])
#        plt.legend(loc ='lower right')
#        plt.axis([-0.1, 1.1, 0, 3])
#        l = l+2
#        i += 1
#    #make residual plot or std dev plot
#    plt.subplot(2,len(StdR),i+len(StdR))
#    plt.plot(conc, stddev, 'bo',conc, stddev, 'b-')
#    plt.ylabel(r'std deviation')
#    plt.xlabel(r'Concentration Nitrate (mg/L)')
#    plt.legend(loc ='upper right')
#    plt.axis([-0.1, 1.1, 0, 0.5])
   

