# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 17:07:15 2017

@author: Lana
"""


import glob as g
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import lmfit
import random

def capitalize(line):
    return ' '.join(s[0].upper() + s[1:] for s in line.split(' '))

def toFloat(pd):
    pd[' Concentration'] = np.array(pd[' Concentration'], dtype = float)
    pd[' Ordinate (A)'] = np.array(pd[' Ordinate (A)'], dtype = float)

# For this code ensure the following:
#- the sample name is identical accross the same sample over days and dilutions 
#  (including Caps and all details), make sure there is no underscore in the 
#  sample name (ie. 0mM20b_ not 0mM_20_b_)
#- seperate information with an underscore ie. OmM_1-100_PI not 0mM1-100PI
#- do not seperate the same information ie. do not do _Day_9_ do _Day9_

#This code will analyse Data of one sample, select enrichment in the following:
EnrichDate = 'October 16'    #input('Enter string of enrichment batch desired with no space ei. "May16":')

#####MAKE FILE LOCATION MORE GENERAL?############
#Files = g.glob('C:\Users\Lana\Desktop\Nitrate Analysis\*'+ EnrichDate)

Folders = g.glob('C:\Users\Lana\Desktop\Nitrate AutoSpec\*'+ EnrichDate +'\*')

files = pd.DataFrame(columns = ['File','Type','Day'])


replicate = 3

i = 0
day0 = int(EnrichDate.split(' ')[1])
#fix for 31 -> 1 ect ect

for folder in Folders:
    files.loc[i,'File'] = folder

    folder = folder.lower()
    
    month = EnrichDate.split(' ')[0]
    date = folder.split(month.lower()+'-',1)[1]
    date = date.split("-")[0]
    files.loc[i,'Day'] = int(date) - day0
    
    if 'cad' in folder:
        files.loc[i,'Type'] = 'Nitrate+Nitrite'
    else:
        files.loc[i,'Type'] = 'Nitrite'
        
    i += 1


SampName = ['Neg0', 'Neg8', '0mM', '6mM', '8mM'] #sample names
SampNameR = ['eg0', 'eg8', '0', '6', '8'] #sample names

dataplot = pd.DataFrame(columns = ['Name','Nitrate', 'Nitrate_Error', 'Nitrite', 'Nitrite_Error', 'Nit+Nit','Nit+Nit_Error','Day'])


n = 1

for k in files.index:
    
    # Fitted Model
    def model(parameters, xi):
    	a = parameters['a'].value
    	#b = parameters['b'].value     
    	return a*xi
    # Residuals to be Minimized
    def residuals(fit_parameters, x_data, y_data, y_errors):
    	return (y_data - model(fit_parameters, x_data)) / y_errors
    
    
    # Load Data
    Std = pd.read_csv(files.loc[k,'File'] + '\Standards Table.csv' )
    toFloat(Std)
    Samp = pd.read_csv(files.loc[k,'File'] + '\Sample Table.csv' )
    Samp[' Ordinate (A)'] = np.array(Samp[' Ordinate (A)'], dtype = float)
      
    
    if files.loc[k,'Type'] == 'Nitrite':
        blank = Samp.loc[0, ' Ordinate (A)']
        S1 = Samp.drop(0)
        S1 = S1.drop(range(7,21))
        S2 = Samp.drop(range(0,7))
        S2 = S2.drop(range(13,21))
        S3 = Samp.drop(range(0,13))
        S3 = S3.drop(range(19,21))
        
        Extra = Samp.drop(range(0,19))
        Extra['Sample ID'] = ['neg8 Zinc','8mm Zinc']
        Extra['Error'] = np.array([0.1]*len(Extra.index))
    else:
        blank = Samp.loc[0, ' Ordinate (A)']
        S1 = Samp.drop(0)
        S1 = S1.drop(range(6,16))
        S2 = Samp.drop(range(0,6))
        S2 = S2.drop(range(11,16))
        S3 = Samp.drop(range(0,11))
          
    samplist = [S1, S2, S3]

    avg = [0] * len(samplist[0])
    stddev = [0] * len(samplist[0])
        
    for a in samplist:
        avg += np.array(a[' Ordinate (A)'])
    avg = avg/len(samplist)
    
    for a in samplist:
        stddev += (np.array(a[' Ordinate (A)']) - avg)**2
    stddev = np.sqrt(stddev/(len(samplist)-1))
    
    stddev = np.append(blank, stddev)
    avg = np.append(blank, avg)
    
    if files.loc[k,'Type'] == 'Nitrite':
        Samp = Samp.drop(range(7, 21))
        Samp[' Ordinate (A)'] = avg
        Samp['Error'] = stddev
    else:
        Samp = Samp.drop(range(6,16))
        Samp[' Ordinate (A)'] = avg
        Samp['Error'] = stddev
    

    x = [0]
    x.extend(Std[' Concentration'])
    x = np.array(x)    
    x = x[1:6]    
    
    y = [Samp.loc[0, ' Ordinate (A)']]
    y.extend(Std[' Ordinate (A)'])
    y = y - Samp.loc[0, ' Ordinate (A)']
    y = np.array(y)
    y = y[1:6]
    y_err = np.array([0.1]*len(y))
    
    if files.loc[k,'Type'] == 'Nitrite':
        Samp[' Ordinate (A)'] = Samp[' Ordinate (A)'] - Samp.loc[6, ' Ordinate (A)'] 
        Samp = Samp.drop(6)
    else:
        Samp[' Ordinate (A)'] = Samp[' Ordinate (A)'] - Samp.loc[0, ' Ordinate (A)'] 
    
    Samp = Samp.drop(0)
    
    # Initialize fit parameters object
    p = lmfit.Parameters()
    #           (Name,   Value,   Vary,   Min,   Max,   Expr)
    p.add_many(('a', 1.0,     True,   None,  None,  None))
    
    # Perform least-squares optimization
    fit_results = lmfit.minimize(residuals, p, args=(x, y, y_err))
    
    # Generate points on the fitted line for plotting
    num_points = 1000
    x_fit = np.linspace(0-max(x)/10, max(x) + max(x)/10, num=num_points)
    y_fit = model(p, x_fit)
    
    # Output Results
    print '\nFit Parameters:'
    print 'a = %.3f, Standard Error in a = %.3f' % (p['a'].value, p['a'].stderr)
    #print 'b = %.3f, Standard Error in b = %.3f' % (p['b'].value, p['b'].stderr)
    
    
    # R VALUE ###################
    print '\n\'Goodness\' of Fit for day',files.loc[k,'Day'], files.loc[k,'Type'],':'
    print 'Chi^2 = ', fit_results.chisqr
    print 'Reduced Chi^2 = ', fit_results.redchi
    core = np.corrcoef(x,y)
    print 'R^2 = ', core[0,1]**2
    # pearson correlations between x and y    
    
    plt.figure(0)
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized()    
    
    # Plot the data and fitted model
    fig = plt.figure(0)
    plt.suptitle(r'Standard Curves for Nitrate Quantification', fontsize = 20)

    #files.loc[k,'Type']
    
    ax = plt.subplot(2,np.ceil(len(Folders)/2) ,n)
    plt.errorbar(x, y, yerr=y_err, fmt='.', color='k')
    plt.plot(x_fit, y_fit, 'r-')
    
    if files.loc[k,'Type'] == 'Nitrite':
        plt.xlabel(r'Nitrite-N (mg/L)')
        plt.ylabel(r'Absorbance (543nm)')
        maxy = 1.5
    else:
        plt.xlabel(r'Nitrate-N + Nitrite-N (mg/L)')
        plt.ylabel(r'Absorbance (400nm)')
        maxy = 1.65
        
    ax.text(0, maxy-0.25, r'R^2 = %s'%core[0,1]**2, fontsize=11)
    plt.axis([0-max(x)/10, max(x)+ max(x)/10, -0.15, maxy])
    plt.title(r'%s Std Day %s'% (files.loc[k,'Type'], files.loc[k,'Day']))    # Change to date
    plt.show()
    
    n += 1

    plt.tight_layout(rect = [0.02,0,0.97,0.95])   
        
    #fig.text(0.02, 0.80, 'Absorbance (400nm)', ha='right', va='top', rotation=90, size = 15)
    #fig.text(0.02, 0.20, 'Absorbance (543nm)', ha='right', va='bottom', rotation=90, size = 15)


    
    i = 0
    
    for a in Samp.index:
        dataplot.loc[5*files.loc[k,'Day']  + i,'Name'] = Samp.loc[a,'Sample ID'] 
        dataplot.loc[5*files.loc[k,'Day']  + i,'Day'] =  files.loc[k,'Day']      

        if files.loc[k,'Type'] == 'Nitrate+Nitrite':
            dataplot.loc[5*files.loc[k,'Day']  + i,'Nit+Nit'] = 100*Samp.loc[a,' Ordinate (A)']/p['a'].value
            dataplot.loc[5*files.loc[k,'Day']  + i,'Nit+Nit_Error'] = dataplot.loc[5*files.loc[k,'Day']  + i,'Nit+Nit']*np.sqrt((Samp.loc[a,'Error']/Samp.loc[a,' Ordinate (A)'])**2+(np.sqrt(len(y)-1)*p['a'].stderr/p['a'].value)**2)/np.sqrt(replicate-1)     
        else:
            dataplot.loc[5*files.loc[k,'Day']  + i,'Nitrite'] = 1000*Samp.loc[a,' Ordinate (A)']/p['a'].value
            dataplot.loc[5*files.loc[k,'Day']  + i,'Nitrite_Error'] = dataplot.loc[5*files.loc[k,'Day']  + i,'Nitrite']*np.sqrt((Samp.loc[a,'Error']/Samp.loc[a,' Ordinate (A)'])**2+(np.sqrt(len(y)-1)*p['a'].stderr/p['a'].value)**2)/np.sqrt(replicate-1)
        i += 1
        
###############################################
#    num = 0
#    i = 0
#    m = 25
#    if files.loc[k,'Type'] == 'Nitrite':
#        for a in Extra.index:
#            dataplot.loc[5*files.loc[k,'Day']  + i + m,'Name'] = Extra.loc[a,'Sample ID'] 
#            dataplot.loc[5*files.loc[k,'Day']  + i + m,'Day'] =  files.loc[k,'Day']      
#    
#            dataplot.loc[5*files.loc[k,'Day']  + i + m,'Nit+Nit'] = 1000*(Extra.loc[a,' Ordinate (A)']-p['b'].value)/p['a'].value
#            dataplot.loc[5*files.loc[k,'Day']  + i + m,'Nit+Nit_Error'] = 100*dataplot.loc[i,'Nit+Nit']*np.sqrt((np.sqrt(Extra.loc[a,'Error']**2+p['b'].stderr**2)/(Extra.loc[a,' Ordinate (A)']-p['b'].value))**2+(p['a'].stderr/p['a'].value)**2)     
#            if (num % 2 == 0): #even 
#                dataplot.loc[5*files.loc[k,'Day']  + i + m,'Nitrite'] = 1000*(Samp.loc[2,' Ordinate (A)']-p['b'].value)/p['a'].value
#                dataplot.loc[5*files.loc[k,'Day']  + i + m,'Nitrite_Error'] = 1000*dataplot.loc[i,'Nitrite']*np.sqrt((np.sqrt(Samp.loc[2,'Error']**2+p['b'].stderr**2)/(Samp.loc[2,' Ordinate (A)']-p['b'].value))**2+(p['a'].stderr/p['a'].value)**2)
#            else: #odd        
#                dataplot.loc[5*files.loc[k,'Day']  + i + m,'Nitrite'] = 1000*(Samp.loc[5,' Ordinate (A)']-p['b'].value)/p['a'].value
#                dataplot.loc[5*files.loc[k,'Day']  + i + m,'Nitrite_Error'] = 1000*dataplot.loc[i,'Nitrite']*np.sqrt((np.sqrt(Samp.loc[5,'Error']**2+p['b'].stderr**2)/(Samp.loc[2,' Ordinate (A)']-p['b'].value))**2+(p['a'].stderr/p['a'].value)**2)
#            num += 1
#            i += 1
#        
 #####################################################   

         
for j in dataplot.index:
    if dataplot.loc[j,'Nit+Nit'] < 0:
        dataplot.loc[j,'Nit+Nit'] = 0
    if dataplot.loc[j,'Nitrite'] < 0:
        dataplot.loc[j,'Nitrite'] = 0
        
    dataplot.loc[j,'Nitrate'] = (dataplot.loc[j,'Nit+Nit'] - dataplot.loc[j,'Nitrite'])
    dataplot.loc[j,'Nitrate_Error'] = np.sqrt(dataplot.loc[j,'Nitrite_Error']**2 + dataplot.loc[j,'Nit+Nit_Error']**2)
    
 

 ####################################


for k in dataplot.index:
    if 'a' in dataplot.loc[k,'Name']:
        name = dataplot.loc[k,'Name']
        dataplot.loc[k,'Name'] = name.split(' ')[0]
########################################
#########################################################################################
#make graph (with many subplots) of the standard of each day.
nameList = np.unique(dataplot['Name'])
nameList.sort()        

#plt.figure(1)
#mng = plt.get_current_fig_manager()
#mng.window.showMaximized() 
#
#i=1
#plt.figure(1)
#plt.clf()
#plt.suptitle('Time Series of Nitrate Quantification', fontsize = 20)
#for samp in nameList:
#    xy = dataplot[dataplot['Name'] == samp]
#    plt.subplot(1,len(nameList), i)
#    plt.plot(xy['Day'], xy['Nitrate'], 'ok')
#    plt.plot(xy['Day'], xy['Nitrate'], '-k')
#    plt.errorbar(xy['Day'], xy['Nitrate'], yerr=xy['Nitrate_Error'], fmt='.', color='k')
#    plt.axis([-0.1, max(files['Day'])+0.1, -max(dataplot['Nitrate'])/20, max(dataplot['Nitrate'])+max(dataplot['Nitrate'])/20]) 
#    plt.xlabel(r'Day', fontsize = 16)
#    plt.title(r'%s'% capitalize(samp), fontsize = 20)
#    i += 1
#plt.subplot(1,len(nameList), 1)
#plt.ylabel(r'Nitrate-N (mg/L)', fontsize = 16)
#plt.tight_layout(rect = [0.0,0,0.99,0.95])    
#   
#
#plt.figure(2)
#mng = plt.get_current_fig_manager()
#mng.window.showMaximized() 
#
#i=1
#plt.figure(2)
#plt.clf()
#plt.suptitle('Time Series of Nitrite Quantification', fontsize = 20)
#for samp in nameList:
#    xy = dataplot[dataplot['Name'] == samp]
#    plt.subplot(1,len(nameList), i)
#    plt.plot(xy['Day'], xy['Nitrite'], 'ok')
#    plt.plot(xy['Day'], xy['Nitrite'], '-k')
#    plt.errorbar(xy['Day'], xy['Nitrite'], yerr=xy['Nitrite_Error'], fmt='.', color='k')
#    plt.axis([-0.1, max(files['Day'])+0.1, -max(dataplot['Nitrite'])/20, max(dataplot['Nitrite'])+max(dataplot['Nitrite'])/20])
#    plt.xlabel(r'Day', fontsize = 16)
#    plt.title(r'%s'% capitalize(samp), fontsize = 20)
#    i += 1
#plt.subplot(1,len(nameList), 1)
#plt.ylabel(r'Nitrite-N (mg/L)', fontsize = 16)
#plt.tight_layout(rect = [0.0,0,0.99,0.95])     
#    
#
#plt.figure(3)
#mng = plt.get_current_fig_manager()
#mng.window.showMaximized() 
#
#i=1    
#plt.figure(3)
#plt.clf()
#plt.suptitle('Time Series of Nitrate and Nitrite Quantification', fontsize = 20)
#for samp in nameList:
#    xy = dataplot[dataplot['Name'] == samp]
#    plt.subplot(1,len(nameList), i)
#    plt.plot(xy['Day'], xy['Nit+Nit'], 'ok')
#    plt.plot(xy['Day'], xy['Nit+Nit'], '-k')
#    plt.errorbar(xy['Day'], xy['Nit+Nit'], yerr=xy['Nit+Nit_Error'], fmt='.', color='k')
#    plt.axis([-0.1, max(files['Day'])+0.1, -max(dataplot['Nit+Nit'])/20, max(dataplot['Nit+Nit'])+max(dataplot['Nit+Nit'])/20])
#    #plt.axis([-0.1, max(files['Day'])+0.1 + 1, 0, 500])
#    plt.xlabel(r'Day', fontsize = 16)
#    plt.title(r'%s'% capitalize(samp), fontsize = 20)
#    i += 1
#plt.subplot(1,len(nameList), 1)
#plt.ylabel(r'Nitrate-N + Nitrite-N (mg/L)', fontsize = 16)
#plt.tight_layout(rect = [0.0,0,0.99,0.95]) 
#

plt.figure(4)
mng = plt.get_current_fig_manager()
mng.window.showMaximized() 

fig = plt.figure(4)
plt.clf()
i=1
plt.suptitle('Time Series of Nitrate Quantification Overview Plot', fontsize = 20)
for samp in nameList:
    xy = dataplot[dataplot['Name'] == samp]
    xy = xy.sort_values('Day')
    plt.subplot(3,len(nameList), i)
    plt.plot(xy['Day'], xy['Nitrate'], 'ok')
    plt.plot(xy['Day'], xy['Nitrate'], '-k')
    plt.errorbar(xy['Day'], xy['Nitrate'], yerr=xy['Nitrate_Error'], fmt='.', color='k')
    ind = dataplot['Nitrate'].idxmax()
    indm = dataplot['Nitrate'].idxmin()
    plt.axis([-0.1, max(files['Day'])+0.1, dataplot.loc[indm, 'Nitrate']-1.5*abs(dataplot.loc[indm, 'Nitrate_Error']), dataplot.loc[ind, 'Nitrate']+1.5*abs(dataplot.loc[ind, 'Nitrate_Error'])])    
    plt.title(r'%s'% capitalize(samp))
    i += 1
for samp in nameList:
    xy = dataplot[dataplot['Name'] == samp]
    xy = xy.sort_values('Day')
    plt.subplot(3,len(nameList), i)
    plt.plot(xy['Day'], xy['Nitrite'], 'ok')
    plt.plot(xy['Day'], xy['Nitrite'], '-k')
    plt.errorbar(xy['Day'], xy['Nitrite'], yerr=xy['Nitrite_Error'], fmt='.', color='k')
    ind = dataplot['Nitrite'].idxmax()
    indm = dataplot['Nitrite'].idxmin()
    plt.axis([-0.1, max(files['Day'])+0.1, dataplot.loc[indm, 'Nitrite']-1.5*abs(dataplot.loc[indm, 'Nitrite_Error']), dataplot.loc[ind, 'Nitrite']+1.5*abs(dataplot.loc[ind, 'Nitrite_Error'])])    
    i += 1        
for samp in nameList:
    xy = dataplot[dataplot['Name'] == samp]
    xy = xy.sort_values('Day')
    plt.subplot(3,len(nameList), i)
    plt.plot(xy['Day'], xy['Nit+Nit'], 'ok')
    plt.plot(xy['Day'], xy['Nit+Nit'], '-k')
    plt.errorbar(xy['Day'], xy['Nit+Nit'], yerr=xy['Nit+Nit_Error'], fmt='.', color='k')
    ind = dataplot['Nit+Nit'].idxmax()
    indm = dataplot['Nit+Nit'].idxmin()
    plt.axis([-0.1, max(files['Day'])+0.1, dataplot.loc[indm, 'Nit+Nit']-1.5*abs(dataplot.loc[indm, 'Nit+Nit_Error']), dataplot.loc[ind, 'Nit+Nit']+1.5*abs(dataplot.loc[ind, 'Nit+Nit_Error'])])
    plt.xlabel(r'Day')
    i += 1




plt.tight_layout(rect = [0.02,0,0.99,0.95])   

fig.text(0.02, 0.85, 'Nitrate-N (mg/L)', ha='right', va='top', rotation=90, size = 11)
fig.text(0.02, 0.4875, 'Nitrite-N (mg/L)', ha='right', va='center', rotation=90, size = 11)
fig.text(0.02, 0.1, 'Nitrate-N + Nitrite-N (mg/L)', ha='right', va='bottom', rotation=90, size = 11)


plt.show()


#plt.subplot(1,len(SampName)+2, 2)
#plt.plot(2,PlotData.loc[1,'Nit+Nit'],'ok') 
#plt.subplot(1,len(SampName)+2, 3)
#plt.plot(2,PlotData.loc[3,'Nit+Nit'],'ok') 
#plt.subplot(1,len(SampName)+2, 4)
#plt.plot(2,PlotData.loc[4,'Nit+Nit'],'ok') 
#plt.subplot(1,len(SampName)+2, 5)
#plt.plot(2,PlotData.loc[0,'Nit+Nit'],'ok') 
#plt.subplot(1,len(SampName)+2, 6)
#plt.plot(2,PlotData.loc[5,'Nit+Nit'],'ok') 
#plt.subplot(1,len(SampName)+2, 0)
#plt.plot(2,PlotData.loc[2,'Nit+Nit'],'ok') 



#calculate uncertainty ect, ect, (error bars)

#fit to curve


#insert equation for absorbance -> conc

# make a table of sample, day, nitrate conc

#plot time series
