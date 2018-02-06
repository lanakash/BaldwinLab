# -*- coding: utf-8 -*-
"""
Created on Tue Jun 06 13:58:16 2017

@author: Lana
"""
import glob as g
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def capitalize(line):
    return ' '.join(s[0].upper() + s[1:] for s in line.split(' '))
    
# For this code ensure the following:
#- the sample name is identical accross the same sample over days and dilutions 
#  (including Caps and all details), make sure there is no underscore in the 
#  sample name (ie. 0mM20b_ not 0mM_20_b_)
#- seperate information with an underscore ie. OmM_1-100_PI not 0mM1-100PI
#- do not seperate the same information ie. do not do _Day_9_ do _Day9_

#This code will analyse Data of one sample, select enrichment in the following:
EnrichDate = 'July11'    #input('Enter string of enrichment batch desired with no space ei. "May16":')

#####MAKE FILE LOCATION MORE GENERAL?############
Files = g.glob('C:\Users\Lana\Desktop\FACS\*\*'+ EnrichDate + '*.csv')

data = pd.DataFrame(columns = ['File','Name','Dilution','Day', 'Stain', 'Filtered'])

# Collect file info into data frame
i = 0 
for csv in Files:
    data.loc[i,'File'] = csv
    name = csv.split("\\")[-1]    
    data.loc[i,'Name'] =  capitalize(name.split("_")[0])
    
    
    
    ##################### There is ones that say Day0_DAy1_
    csv = csv.lower()
    date= csv.split("\\")[-1]
    date= date.split('day',1)[1]
    date = date.split('_',1)[0]
    data.loc[i,'Day'] = int(date.split('.',1)[0])
    
    if 'fil' in csv:
        data.loc[i,'Filtered'] = True
    else:
        data.loc[i,'Filtered'] = False

#####ALLOW CHOOSING FOR DYE?#####
    if 'pi' in csv:
        data.loc[i,'Stain'] = True
    else:
        data.loc[i,'Stain'] = False
    
    if '100' in csv:
        data.loc[i,'Dilution'] = 100
    elif '10' in csv:
        
        #################this includes those random numbers at the end, cut it off
        data.loc[i,'Dilution'] = 10
    else:
        data.loc[i,'Dilution'] = 1
        
    i += 1
        
        
# Dataframe to plot time series
countData = pd.DataFrame(columns = ['Name','Dilution','Day', 'Count'])
j = 0   

for fil in [False,True]:
    #Filter for dataframe for this filter status
    dataFil = data[data['Filtered'] == fil]
    
    for day in list(set(dataFil['Day'])):   
        #Filter for dataframe for this day 
        dataDay = dataFil[dataFil['Day'] == day]
        
        #Select and record medium controls then remove from data frame
        for i in dataDay.index:
            name = dataDay.loc[i,'Name']
            if 'med' in name.lower():
                if dataDay.loc[i,'Stain']:
                    medium_stain_file = dataDay.loc[i,'File']
                else:
                    medium_nostain_file = dataDay.loc[i,'File']
                dataDay = dataDay.drop(i)
       
        for dilu in [1, 10, 100]:
            #filter dataframe for dilutions 
            dataDayDilu = dataDay[dataDay['Dilution'] == dilu]
            
            #Select and record exp sample control then remove from data frame
            for i in dataDayDilu.index:
                name = dataDayDilu.loc[i,'Name']
                if 'samp' in name.lower() or 'exp' in name.lower():
                    if dataDayDilu.loc[i,'Stain'] == False:
                        samp_nostain_file = dataDayDilu.loc[i,'File']
                    dataDayDilu = dataDayDilu.drop(i)        
            
            #For each experimental sample with stain calculate cell count and record in list 
            for i in dataDayDilu.index:  
                samp_stain_file = dataDayDilu.loc[i,'File']
                
                #FILTER###############################################
                #logging makes this code super slow
                s_dye = pd.read_csv(samp_stain_file)
                #log_s_dye = np.log(s_dye+500)
                s = s_dye[s_dye['FITC-A']>np.exp(6.2)]
                s = s[s['PI/PE-Cy5.5-A']>np.exp(5.8)]
                #######################################################
                
                
                #store values for time series
                if dataDayDilu.loc[i,'Filtered']:
                    countData.loc[j,'Name'] = dataDayDilu.loc[i,'Name'] + ' Filtered'
                else:
                    countData.loc[j,'Name'] = dataDayDilu.loc[i,'Name']
                countData.loc[j,'Dilution'] = dilu
                countData.loc[j,'Day'] = day
                countData.loc[j,'Count'] = len(s.index)
                
                j += 1

           

#plot time series  
matplotlib.style.use('ggplot')                 
fig = plt.figure(0)
plt.suptitle(r'Time series of Flow cytometry Data filtered by PI and sybr', fontsize = 18)

nameList = list(set(countData['Name']))
i = 1
MaxData1 = countData[countData['Dilution'] == 1]
Max1 = max(MaxData1['Count'])/82*1000
MaxDay1 = max(MaxData1['Day'])
MaxData10 = countData[countData['Dilution'] == 10]
Max10 = max(MaxData10['Count'])/82*1000
MaxDay10 = max(MaxData10['Day'])
MaxData100 = countData[countData['Dilution'] == 100]
Max100 = max(MaxData100['Count'])/82*1000
MaxDay100 = max(MaxData100['Day'])

#######get rid of numbers in middle plots (axis), change backgroud to white
#######make it pretty

for samp in nameList:
    
    plotData = countData[countData['Name'] == samp]
    plotData1 = plotData[plotData['Dilution'] == 1]
    plt.subplot(3,len(nameList), i)
    plt.plot(plotData1['Day'], plotData1['Count']/82*1000, color = 'k')
    plt.axis([0, MaxDay1, 0, Max1+(Max1/10)])
    plt.title(r'%s'% samp)

    plotData10 = plotData[plotData['Dilution'] == 10]
    plt.subplot(3,len(nameList), i+len(nameList))
    plt.plot(plotData10['Day'], plotData10['Count']/82*1000, color = 'k')
    plt.axis([0, MaxDay10, 0, Max10+(Max10/10)])
    
    plotData100 = plotData[plotData['Dilution'] == 100]
    plt.subplot(3,len(nameList), i+len(nameList)*2)
    plt.plot(plotData100['Day'], plotData100['Count']/82*1000, color = 'k')
    plt.axis([0, MaxDay100, 0, Max100+(Max100/10)])
    plt.xlabel('Day')
    
    i += 1
   
    
plt.subplot(3,len(nameList), 1) 
plt.ylabel('Cell count') 
plt.subplot(3,len(nameList), 1+len(nameList))
plt.ylabel('Cell count') 
plt.subplot(3,len(nameList), 1+len(nameList)*2)
plt.ylabel('Cell count')    
    
fig.text(0.93, 0.80, '1:1', ha='right', va='top', rotation=-90)
fig.text(0.93, 0.5, '1:10', ha='right', va='center', rotation=-90)
fig.text(0.93, 0.20, '1:100', ha='right', va='bottom', rotation=-90)
  
    
    