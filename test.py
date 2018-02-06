# -*- coding: utf-8 -*-
"""
Created on Wed May 24 21:34:31 2017

@author: Lana
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import random
import lmfit

def toFloat(pd):
    pd[' Concentration'] = np.array(pd[' Concentration'], dtype = float)
    pd[' Ordinate (A)'] = np.array(pd[' Ordinate (A)'], dtype = float)

dataZ = pd.read_csv('C:\Users\Lana\Desktop\zinc_aug_23ed.csv')
dataC = pd.read_csv('C:\Users\Lana\Desktop\cadmium_aug_23ed.csv')

stdZ = dataZ[dataZ['Type'] == 'Std']
stdC = dataC[dataC['Type'] == 'Std']

Sample = pd.DataFrame(columns = ['Name','Type',' Ordinate (A)', ' Concentration', 'Error'])
i = 0
k = 0

for data in [dataZ[dataZ['Type'] != 'Std'], dataC[dataC['Type'] != 'Std']]:

    for group in np.unique(data['Type']):
        glist = data[data['Type']== group]
        avg = np.average(glist['Abs'])
        stddev = np.average((glist['Abs'] - avg)**2)
        Sample.loc[i, ' Ordinate (A)'] = avg
        Sample.loc[i, 'Error'] = stddev
        Sample.loc[i, 'Name'] = group
        if k == 0:
            Sample.loc[i, 'Type'] = 'Nitrite'
        else:
            Sample.loc[i, 'Type'] = 'Nit+Nit'
        i += 1
    k += 1
    


# Fitted Model
def model(parameters, xi):
	a = parameters['a'].value
	b = parameters['b'].value     
	return a*xi+b
# Residuals to be Minimized
def residuals(fit_parameters, x_data, y_data, y_errors):
	return (y_data - model(fit_parameters, x_data)) / y_errors

x = stdZ['Sample']
x = np.array(map(float,x))
y = np.array(stdZ['Abs'])
y_err =  np.array([0.1]*6)

# Initialize fit parameters object
p = lmfit.Parameters()
#           (Name,   Value,   Vary,   Min,   Max,   Expr)
p.add_many(('a', 1.0,     True,   None,  None,  None),
           ('b', 5.0,     True,   None,  None,  None))

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

print '\n\'Goodness\' of Fit Zinc:'
print 'Chi^2 = ', fit_results.chisqr
print 'Reduced Chi^2 = ', fit_results.redchi
core = np.corrcoef(x,y)
print 'R^2 = ', core[0,1]**2

# Plot the data and fitted model
plt.figure(0)
plt.subplot(2,1,1)
plt.errorbar(x, y, yerr=y_err, fmt='o', color='r', label = 'Zinc std 4thfloor')
plt.plot(x_fit, y_fit, 'r-')
plt.xlabel(r'Nitrate Concentration (mg/L)')
plt.ylabel(r'Absorbance (543nm)')
plt.title(r'Zinc Standard')
plt.axis([-0.1,1.1,0,2])
plt.legend()
plt.show()

sampZ = Sample[Sample['Type'] == 'Nitrite']
for a in sampZ.index:
    Sample.loc[a, ' Concentration'] = 1000*(sampZ.loc[a,' Ordinate (A)']-p['b'].value)/p['a'].value

# Fitted Model
def model(parameters, xi):
	a = parameters['a'].value
	b = parameters['b'].value     
	return a*xi+b
# Residuals to be Minimized
def residuals(fit_parameters, x_data, y_data, y_errors):
	return (y_data - model(fit_parameters, x_data)) / y_errors

x = stdC['Sample']
x = np.array(map(float,x))
y = np.array(stdC['Abs'])
y_err =  np.array([0.1]*6)

# Initialize fit parameters object
p = lmfit.Parameters()
#           (Name,   Value,   Vary,   Min,   Max,   Expr)
p.add_many(('a', 1.0,     True,   None,  None,  None),
           ('b', 5.0,     True,   None,  None,  None))

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

print '\n\'Goodness\' of Fit Cad:'
print 'Chi^2 = ', fit_results.chisqr
print 'Reduced Chi^2 = ', fit_results.redchi
core = np.corrcoef(x,y)
print 'R^2 = ', core[0,1]**2

# Plot the data and fitted model
plt.subplot(2,1,2)
plt.errorbar(x, y, yerr=y_err, fmt='o', color='r', label = 'Cad Std 4thfloor')
plt.plot(x_fit, y_fit, 'r-')
plt.xlabel(r'Nitrate Concentration (mg/L)')
plt.ylabel(r'Absorbance (543nm)')
plt.title(r'Cadmium standard')
plt.axis([-0.5,10.5,0,1.2])
plt.legend()
plt.show()

sampC = Sample[Sample['Type'] == 'Nit+Nit']
for a in sampC.index:
    Sample.loc[a, ' Concentration'] = 100*(sampC.loc[a,' Ordinate (A)']-p['b'].value)/p['a'].value


PlotData = pd.DataFrame(columns = ['Name','Nitri', 'Nit+Nit', 'Nitra'])

i = 0

PlotData.loc[i,'Name'] = 'Neg8 Zinc'
PlotData.loc[i,'Nit+Nit'] = Sample.loc[5,' Concentration']
PlotData.loc[i,'Nitri'] = Sample.loc[4,' Concentration']
if PlotData.loc[i,'Nit+Nit'] < 0:
    PlotData.loc[i,'Nit+Nit'] = 0
if PlotData.loc[i,'Nitri'] < 0:
    PlotData.loc[i,'Nitri'] = 0
PlotData.loc[i,'Nitra'] = PlotData.loc[i,'Nit+Nit'] - PlotData.loc[i,'Nitri']
i += 1
PlotData.loc[i,'Name'] = '8mm Zinc'
PlotData.loc[i,'Nit+Nit'] = Sample.loc[2,' Concentration']
PlotData.loc[i,'Nitri'] = Sample.loc[1,' Concentration']
if PlotData.loc[i,'Nit+Nit'] < 0:
    PlotData.loc[i,'Nit+Nit'] = 0
if PlotData.loc[i,'Nitri'] < 0:
    PlotData.loc[i,'Nitri'] = 0
PlotData.loc[i,'Nitra'] = PlotData.loc[i,'Nit+Nit'] - PlotData.loc[i,'Nitri']
if PlotData.loc[i,'Nitra'] < 0:
    PlotData.loc[i,'Nitra'] = 0
i += 1
Sample = Sample.drop(2)
Sample = Sample.drop(5)


for name in np.unique(Sample['Name']):
    PlotData.loc[i,'Name'] = name
    dup = Sample[Sample['Name'] == name]
    PlotData.loc[i,'Nit+Nit'] = np.array(dup[dup['Type'] == 'Nit+Nit'][' Concentration'])[0]
    PlotData.loc[i,'Nitri'] = np.array(dup[dup['Type'] == 'Nitrite'][' Concentration'])[0]
    if PlotData.loc[i,'Nit+Nit'] < 0:
        PlotData.loc[i,'Nit+Nit'] = 0
    if PlotData.loc[i,'Nitri'] < 0:
        PlotData.loc[i,'Nitri'] = 0
    PlotData.loc[i,'Nitra'] = PlotData.loc[i,'Nit+Nit'] - PlotData.loc[i,'Nitri']
    i += 1


#######################################################################


dataZ = pd.read_csv('C:\Users\Lana\Desktop\Lana Nitrate Zinc August-23-17 1_21 PM Pacific Daylight Time\Results Table.csv')
dataC = pd.read_csv('C:\Users\Lana\Desktop\Lana Nitrate Cad August-23-17 1_31 PM Pacific Daylight Time\Results Table.csv')

stdZ = dataZ.drop(range(6,23))
stdC = dataC.drop(range(6,17))

Sample = pd.DataFrame(columns = ['Name','Type',' Ordinate (A)', ' Concentration', 'Error'])
i = 0
k = 0

sampZ = dataZ.drop(range(0,6))
sampC = dataC.drop(range(0,6))


    
for a in [sampZ, sampC]:
    for i in a.index:
        name = a.loc[i,'Sample ID']
        if ' ' in name:
            a.loc[i,'Sample ID'] = name.split(' ')[0]
        else:
            a.loc[i,'Sample ID'] = name + ' Zinc'

for data in [sampZ, sampC]:
    for group in np.unique(data['Sample ID']):
        glist = data[data['Sample ID']== group]
        avg = np.average(glist[' Ordinate (A)'])
        stddev = np.average((glist[' Ordinate (A)'] - avg)**2)
        Sample.loc[i, ' Ordinate (A)'] = avg
        Sample.loc[i, 'Error'] = stddev
        Sample.loc[i, 'Name'] = group
        if k == 0:
            Sample.loc[i, 'Type'] = 'Nitrite'
        else:
            Sample.loc[i, 'Type'] = 'Nit+Nit'
        i += 1
    k += 1
   
Sample[' Ordinate (A)'] = np.array(Sample[' Ordinate (A)']) - Sample.loc[19, ' Ordinate (A)']
Sample = Sample.drop(19)


# Fitted Model
def model(parameters, xi):
	a = parameters['a'].value
	b = parameters['b'].value     
	return a*xi+b
# Residuals to be Minimized
def residuals(fit_parameters, x_data, y_data, y_errors):
	return (y_data - model(fit_parameters, x_data)) / y_errors

toFloat(stdZ)

x = stdZ[' Concentration']
y = np.array(stdZ[' Ordinate (A)']) - stdZ.loc[0, ' Ordinate (A)']
y_err =  np.array([0.1]*6)

# Initialize fit parameters object
p = lmfit.Parameters()
#           (Name,   Value,   Vary,   Min,   Max,   Expr)
p.add_many(('a', 1.0,     True,   None,  None,  None),
           ('b', 5.0,     True,   None,  None,  None))

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

print '\n\'Goodness\' of Fit Zinc5:'
print 'Chi^2 = ', fit_results.chisqr
print 'Reduced Chi^2 = ', fit_results.redchi
core = np.corrcoef(x,y)
print 'R^2 = ', core[0,1]**2

# Plot the data and fitted model
plt.figure(0)
plt.subplot(2,1,1)
plt.errorbar(x, y, yerr=y_err, fmt='o', color='b', label = 'Zinc std 5thfloor')
plt.plot(x_fit, y_fit, 'b-')
plt.xlabel(r'Nitrate Concentration (mg/L)')
plt.ylabel(r'Absorbance (543nm)')
plt.title(r'Zinc Standard')
plt.legend(loc = 'lower right')
plt.axis([-0.1,1.1,0,2])
plt.show()

sampZ = Sample[Sample['Type'] == 'Nitrite']
for a in sampZ.index:
    Sample.loc[a, ' Concentration'] = 1000*(sampZ.loc[a,' Ordinate (A)']-p['b'].value)/p['a'].value

# Fitted Model
def model(parameters, xi):
	a = parameters['a'].value
	b = parameters['b'].value     
	return a*xi+b
# Residuals to be Minimized
def residuals(fit_parameters, x_data, y_data, y_errors):
	return (y_data - model(fit_parameters, x_data)) / y_errors

toFloat(stdC)
x = stdC[' Concentration']
y =  np.array(stdC[' Ordinate (A)']) - stdC.loc[0, ' Ordinate (A)']
y_err =  np.array([0.1]*6)

# Initialize fit parameters object
p = lmfit.Parameters()
#           (Name,   Value,   Vary,   Min,   Max,   Expr)
p.add_many(('a', 1.0,     True,   None,  None,  None),
           ('b', 5.0,     True,   None,  None,  None))

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

print '\n\'Goodness\' of Fit Cad5:'
print 'Chi^2 = ', fit_results.chisqr
print 'Reduced Chi^2 = ', fit_results.redchi
core = np.corrcoef(x,y)
print 'R^2 = ', core[0,1]**2

# Plot the data and fitted model
plt.subplot(2,1,2)
plt.errorbar(x, y, yerr=y_err, fmt='o', color='b', label = 'Cad Std 5thfloor')
plt.plot(x_fit, y_fit, 'b-')
plt.axis([-0.5,10.5,0,1.2])
plt.xlabel(r'Nitrate Concentration (mg/L)')
plt.ylabel(r'Absorbance (543nm)')
plt.title(r'Cadmium standard')
plt.legend(loc = 'lower right')

plt.show()

sampC = Sample[Sample['Type'] == 'Nit+Nit']
for a in sampC.index:
    Sample.loc[a, ' Concentration'] = 100*(sampC.loc[a,' Ordinate (A)']-p['b'].value)/p['a'].value


PlotData5 = pd.DataFrame(columns = ['Name','Nitri', 'Nit+Nit', 'Nitra'])

i = 0

PlotData5.loc[i,'Name'] = 'Neg8 Zinc'
PlotData5.loc[i,'Nit+Nit'] = Sample.loc[22,' Concentration']
PlotData5.loc[i,'Nitri'] = Sample.loc[21,' Concentration']
if PlotData5.loc[i,'Nit+Nit'] < 0:
    PlotData5.loc[i,'Nit+Nit'] = 0
if PlotData5.loc[i,'Nitri'] < 0:
    PlotData5.loc[i,'Nitri'] = 0
PlotData5.loc[i,'Nitra'] = PlotData5.loc[i,'Nit+Nit'] - PlotData5.loc[i,'Nitri']
i += 1
PlotData5.loc[i,'Name'] = '8mm Zinc'
PlotData5.loc[i,'Nit+Nit'] = Sample.loc[18,' Concentration']
PlotData5.loc[i,'Nitri'] = Sample.loc[17,' Concentration']
if PlotData5.loc[i,'Nit+Nit'] < 0:
    PlotData5.loc[i,'Nit+Nit'] = 0
if PlotData5.loc[i,'Nitri'] < 0:
    PlotData5.loc[i,'Nitri'] = 0
PlotData5.loc[i,'Nitra'] = PlotData5.loc[i,'Nit+Nit'] - PlotData5.loc[i,'Nitri']
if PlotData5.loc[i,'Nitra'] < 0:
    PlotData5.loc[i,'Nitra'] = 0
i += 1
Sample = Sample.drop(18)
Sample = Sample.drop(22)


for name in np.unique(Sample['Name']):
    PlotData5.loc[i,'Name'] = name
    dup = Sample[Sample['Name'] == name]
    PlotData5.loc[i,'Nit+Nit'] = np.array(dup[dup['Type'] == 'Nit+Nit'][' Concentration'])[0]
    PlotData5.loc[i,'Nitri'] = np.array(dup[dup['Type'] == 'Nitrite'][' Concentration'])[0]
    if PlotData5.loc[i,'Nit+Nit'] < 0:
        PlotData5.loc[i,'Nit+Nit'] = 0
    if PlotData5.loc[i,'Nitri'] < 0:
        PlotData5.loc[i,'Nitri'] = 0
    PlotData5.loc[i,'Nitra'] = PlotData5.loc[i,'Nit+Nit'] - PlotData5.loc[i,'Nitri']
    i += 1