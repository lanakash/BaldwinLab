# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 18:12:54 2017

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



#This code will analyese Data of one sample, select enrichment in the following:
EnrichDates = ['July11', 'July17', 'July24', 'Aug8', 'Aug14', 'Aug28']    #input('Enter string of enrichment batch desired with no space ei. "May16":')

countData = pd.DataFrame(columns = ['Name','Dilution','Day', 'Count'])     
j = 0 

for EnrichDate in EnrichDates:
    #####MAKE FILE LOCATION MORE GENERAL?############
    Files = g.glob('C:\Users\Lana\Desktop\FACS\*\*'+ EnrichDate + '*.csv')
    
    data = pd.DataFrame(columns = ['File','Name','Dilution','Day', 'Stain', 'Filtered'])
    
    # Collect file info into data frame
    i = 0 
    for csv in Files:
        data.loc[i,'File'] = csv
        name = csv.split("\\")[-1]    
        name =  name.split("_")[0]
        data.loc[i,'Name'] =  capitalize(name.split(" ")[0])
        
        
        
        ##################### There is ones that say Day0_DAy1_
        csv = csv.lower()
        date= csv.split("\\")[-1]
        date= date.split('day',1)[1]
        date = date.split('_',1)[0]
        data.loc[i,'Day'] = int(date.split('.',1)[0])
    
        #####ALLOW CHOOSING FOR DYE?#####
        if 'sybr' in csv:
            data.loc[i,'Stain'] = True
        else:
            data.loc[i,'Stain'] = False    
        
        csv = csv.split("\\")[-1]
        csv = csv.split("sybr")[0]
        
        
        if 'fil' in csv:
            data.loc[i,'Filtered'] = True
        else:
            data.loc[i,'Filtered'] = False
    
    
        
        if '100' in csv:
            data.loc[i,'Dilution'] = 100
        elif '10' in csv:
            
            #################this includes those random numbers at the end, cut it off
            data.loc[i,'Dilution'] = 10
        else:
            data.loc[i,'Dilution'] = 1
            
        i += 1
    
    #data we dont want in the plot        
    data = data[~data['Name'].str.contains('sip')]
    data = data[~data['Name'].str.contains('live')]
    data = data[~data['Name'].str.contains('test')]
    data = data[~data['Name'].str.contains('Sip')]
    data = data[~data['Name'].str.contains('Live')]
    data = data[~data['Name'].str.contains('Test')]
    #data = data[~data['Name'].str.contains('-')]
        
    # Dataframe to plot time series
        
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
                    s = s_dye[s_dye['SYBR Green-A']>np.exp(8)]
                    #s = s[s['PI_2-A']>np.exp(7)]
                    #######################################################
    #                m_dye = pd.read_csv(medium_stain_file)
    #                m = m_dye[m_dye['SYBR Green-A']>np.exp(6.2)]
                    #m = m[m['PI-A']>np.exp(5.8)]
    
    #                s = s.drop(['BL2-A', 'BL2-H','BL2-W','Time', 'PI-H', 'PI-W', 'SYBR Green-H', 'SYBR Green-W', 'FSC-H', 'FSC-W', 'SSC-H', 'SSC-W'], axis=1)
    #                m = m.drop(['BL2-A', 'BL2-H','BL2-W','Time', 'PI-H', 'PI-W', 'SYBR Green-H', 'SYBR Green-W', 'FSC-H', 'FSC-W', 'SSC-H', 'SSC-W'], axis=1)
    #                
    #                log_s_dye = np.log(s)
    #                log_m_dye = np.log(m)
    #                #go through each point left in med, look to see if there is anything within a 
    #                #small range from it, take down index # in list. 
    #                # delete collected index numbers, 
    #                overlapArea = 0.1
    #                
    #                for particle in m.index:
    #                    done = False
    #                    for k in s.index:
    #                        if done==False and (m.loc[particle,'FSC-A']-s.loc[k,'FSC-A']) < overlapArea:
    #                            if (m.loc[particle,'SSC-A']-s.loc[k,'SSC-A']) < overlapArea:
    #                                s = s.drop(k)
    #                                m = m.drop(particle)   #append to m_bg instead?
    #                                done = True
                            
                    
                    #store values for time series
                    if dataDayDilu.loc[i,'Filtered']:
                        countData.loc[j,'Name'] = dataDayDilu.loc[i,'Name'] + ' Filtered'
                    else:
                        countData.loc[j,'Name'] = dataDayDilu.loc[i,'Name']
                    countData.loc[j,'Dilution'] = dilu
                    countData.loc[j,'Day'] = day
                    countData.loc[j,'Count'] = len(s.index)
                    
                    j += 1



tripCountData = pd.DataFrame(columns = ['Name','Dilution','Day', 'Count', 'Count Error'])
j = 0 

for name in list(set(countData['Name'])):
    nameData = countData[countData['Name'] == name]
    
    for day in list(set(nameData['Day'])):
        dayData = nameData[nameData['Day'] == day]
        
        for dilu in list(set(dayData['Dilution'])):
            diluData = dayData[dayData['Dilution'] == dilu]
            
            tripCountData.loc[j,'Name'] = name  
            tripCountData.loc[j,'Day'] = day
            tripCountData.loc[j,'Dilution'] = dilu
            tripCountData.loc[j,'Count'] = sum(diluData['Count'])/len(diluData['Count'])
            
            if len(diluData['Count']) > 1:
                stddev = (diluData['Count'] - tripCountData.loc[j,'Count'])**2
                tripCountData.loc[j,'Count Error'] = np.sqrt(sum(stddev)/(len(diluData['Count'])-1))           
            else:
                tripCountData.loc[j,'Count Error'] = 0
            j += 1
            
            
        
countData = tripCountData



plt.figure(0)
mng = plt.get_current_fig_manager()
mng.window.showMaximized()
        

#plot time series  
matplotlib.style.use('ggplot')                 
fig = plt.figure(0)
plt.clf()
plt.suptitle(r'Time Series of Flow Cytometry Data Filtered by Sybrsafe', fontsize = 20)

nameList = list(set(countData['Name']))

countData = countData.drop(116)

nameList.sort()
i = 1
MaxData1 = countData[countData['Dilution'] == 1]
Max1 = max(MaxData1['Count'])/50*1000
MaxDay1 = max(MaxData1['Day'])
MaxData10 = countData[countData['Dilution'] == 10]
Max10 = max(MaxData10['Count'])/50*1000
MaxDay10 = max(MaxData10['Day'])
MaxData100 = countData[countData['Dilution'] == 100]
Max100 = max(MaxData100['Count'])/50*1000
MaxDay100 = max(MaxData100['Day'])

#######get rid of numbers in middle plots (axis), change backgroud to white
#######make it pretty

for samp in nameList:
    
    plotData = countData[countData['Name'] == samp]
    plotData1 = plotData[plotData['Dilution'] == 1]
    ax = plt.subplot(3,len(nameList), i)
    plt.plot(plotData1['Day'], plotData1['Count']/50*1000, color = 'k')
    plt.errorbar(plotData1['Day'], plotData1['Count']/50*1000, yerr=plotData1['Count Error'], fmt='.', color='k')
    plt.axis([0, MaxDay1, 0, Max1+(Max1/10)])
    plt.title(r'%s'% samp)
    plt.rc('font', size = 12)
    plt.xticks(fontsize = 10)
    plt.yticks(fontsize = 10)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useMathText=True)
    ax.yaxis.major.formatter._useMathText = True    
    #t = ax.yaxis.get_offset_text()
    #t.set_x(-0.15)    

    plotData10 = plotData[plotData['Dilution'] == 10]
    ax = plt.subplot(3,len(nameList), i+len(nameList))
    plt.plot(plotData10['Day'], plotData10['Count']/50*1000, color = 'k')
    plt.errorbar(plotData10['Day'], plotData10['Count']/50*1000, yerr=plotData10['Count Error'], fmt='.', color='k')
    plt.axis([0, MaxDay10, 0, Max10+(Max10/10)])
    plt.xticks(fontsize = 10)
    plt.yticks(fontsize = 10)    
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax.yaxis.major.formatter._useMathText = True    
    
    plotData100 = plotData[plotData['Dilution'] == 100]
    ax = plt.subplot(3,len(nameList), i+len(nameList)*2)
    plt.plot(plotData100['Day'], plotData100['Count']/50*1000, color = 'k')
    plt.errorbar(plotData100['Day'], plotData100['Count']/50*1000, yerr=plotData100['Count Error'], fmt='.', color='k')
    plt.axis([0, MaxDay100, 0, Max100+(Max100/10)])
    plt.xlabel('Day', fontsize = 11)
    plt.xticks(fontsize = 10)
    plt.yticks(fontsize = 10)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax.yaxis.major.formatter._useMathText = True    
    
    i += 1
   
    
#plt.subplot(3,len(nameList), 1) 
#plt.ylabel('Cell count/ml') 
#plt.subplot(3,len(nameList), 1+len(nameList))
#plt.ylabel('Cell count/ml') 
#plt.subplot(3,len(nameList), 1+len(nameList)*2)
#plt.ylabel('Cell count/ml')    


plt.tight_layout(rect = [0.02,0,0.97,0.95])   
    
fig.text(0.97, 0.8, '1:1', ha='right', va='top', rotation=-90)
fig.text(0.97, 0.4875, '1:10', ha='right', va='center', rotation=-90)
fig.text(0.97, 0.165, '1:100', ha='right', va='bottom', rotation=-90)

fig.text(0.99, 0.4875, 'Dilutions', ha='right', va='center', rotation=-90, size = 13)

fig.text(0.02, 0.83, 'Cell count/ml', ha='right', va='top', rotation=90, size = 11)
fig.text(0.02, 0.4875, 'Cell count/ml', ha='right', va='center', rotation=90, size = 11)
fig.text(0.02, 0.135, 'Cell count/ml', ha='right', va='bottom', rotation=90, size = 11)


plt.show()
    