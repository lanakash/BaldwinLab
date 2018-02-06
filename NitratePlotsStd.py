# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 11:07:29 2017

@author: Lana
"""

import glob as g
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import lmfit

#This code will analyse Data of one sample, select enrichment in the following:
EnrichDate = 'July11'    #input('Enter string of enrichment batch desired with no space ei. "May16":')

SampName = ['Neg0', 'Neg8', '0mM', '6mM', '8mM'] #sample names
SampNameR = ['eg0', 'eg8', '0', '6', '8'] #sample names

dataplot = pd.DataFrame(columns = ['Name','Nitrate', 'Nitrate_Error', 'Nitrite', 'Nitrite_Error', 'Nit+Nit','Nit+Nit_Error','Day'])
    
#####MAKE FILE LOCATION MORE GENERAL?############
Files = g.glob('C:\Users\Lana\Desktop\Nitrate Analysis\*\*'+ EnrichDate + '*.csv')


j = 0
i = 0
days = 0

for fileName in Files:
    
    # Fitted Model
    def model(parameters, xi):
    	a = parameters['a'].value
    	b = parameters['b'].value     
    	return a*xi+b
    # Residuals to be Minimized
    def residuals(fit_parameters, x_data, y_data, y_errors):
    	return (y_data - model(fit_parameters, x_data)) / y_errors
    
    
    # Load Data
    data = pd.read_csv(fileName)
    
    Std = data[data['Name'].str.contains('Standard')]
    Samp = data[~data['Name'].str.contains('Standard')]
    x = Std['Concentration (mg/L)']
    y = Std['Ordinate(A)']
    y_err = Std['Error']
    
    
    # Initialize fit parameters object
    p = lmfit.Parameters()
    #           (Name,   Value,   Vary,   Min,   Max,   Expr)
    p.add_many(('a', 1.0,     True,   None,  None,  None),
               ('b', 0.0,     True,   None,  None,  None))
    
    # Perform least-squares optimization
    fit_results = lmfit.minimize(residuals, p, args=(x, y, y_err))
    
    # Generate points on the fitted line for plotting
    num_points = 1000
    x_fit = np.linspace(min(x), max(x), num=num_points)
    y_fit = model(p, x_fit)
    
    # Output Results
    print '\nFit Parameters:'
    print 'a = %.3f, Standard Error in a = %.3f' % (p['a'].value, p['a'].stderr)
    print 'b = %.3f, Standard Error in b = %.3f' % (p['b'].value, p['b'].stderr)
    
    print '\n\'Goodness\' of Fit:'
    print 'Chi^2 = ', fit_results.chisqr
    print 'Reduced Chi^2 = ', fit_results.redchi
    
    
    # Plot the data and fitted model
    plt.figure(0)
    plt.errorbar(x, y, yerr=y_err, fmt='.', color='k')
    plt.plot(x_fit, y_fit, 'r-')
    plt.xlabel(r'Concentration Nitrate (mg/L)')
    plt.ylabel(r'Absorbance (543nm)')
    plt.title(r'Nitrate Analysis')
    plt.axis([-0.1, 1.1, -0.1, 0.5])
    plt.show()
    
    date= fileName.lower()
    date = date.split(".")[0]
    date = int(date.split('day')[1])
    
    Nit = Samp[Samp['Name'].str.contains('ni')]
    Samp = Samp[~Samp['Name'].str.contains('ni')]
    
    
    for n in SampNameR:
        for a in Samp.index:
            if n.lower() in Samp.loc[a,'Name']:
                dataplot.loc[i,'Nit+Nit'] = (Samp.loc[a,'Ordinate(A)']-p['b'].value)/p['a'].value
                dataplot.loc[i,'Nit+Nit_Error'] = dataplot.loc[i,'Nit+Nit']*np.sqrt((np.sqrt(Samp.loc[a,'Error']**2+p['b'].stderr**2)/(Samp.loc[a,'Ordinate(A)']-p['b'].value))**2+(p['a'].stderr/p['a'].value)**2)
                Samp = Samp.drop(a)     
        for a in Nit.index:
            if n.lower() in Nit.loc[a,'Name']:
                dataplot.loc[i,'Nitrite'] = (Nit.loc[a,'Ordinate(A)']-p['b'].value)/p['a'].value
                dataplot.loc[i,'Nitrite_Error'] = dataplot.loc[i,'Nitrite']*np.sqrt((np.sqrt(Nit.loc[a,'Error']**2+p['b'].stderr**2)/(Nit.loc[a,'Ordinate(A)']-p['b'].value))**2+(p['a'].stderr/p['a'].value)**2)
                Nit = Nit.drop(a)
        i += 1
    
    
    count = 0          
    while count < 5:
        dataplot.loc[j,'Nitrate'] = dataplot.loc[j,'Nit+Nit'] - dataplot.loc[j,'Nitrite']
        dataplot.loc[j,'Nitrate_Error'] = np.sqrt(dataplot.loc[j,'Nitrite_Error']**2 + dataplot.loc[j,'Nit+Nit_Error']**2)
        dataplot.loc[j,'Day'] =  date       
        count += 1    
        j += 1
    
    days += 1
    
i=0
for k in range(0,days):
    for n in SampName:
        dataplot.loc[i,'Name'] = n
        i += 1
########################################################################################
#make graph (with many subplots) of the standard of each day.
for samp in SampName:
    
    xy = dataplot[dataplot['Name'] == samp]
    plt.subplot(3,len(SampName), i)
    plt.errorbar(xy['Day'], xy['Nitrate'], yerr=xy['Nitrate_Error'], fmt='.', color='k')
    plt.axis([-0.1, days-0.9, 0, 1])
    plt.title(r'%s'% samp)
    i += 1
#calculate uncertainty ect, ect, (error bars)

#fit to curve


#insert equation for absorbance -> conc

# make a table of sample, day, nitrate conc

#plot time series