# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 16:34:49 2017

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
EnrichDate = 'Dec4'    #input('Enter string of enrichment batch desired with no space ei. "May16":')

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
        
data = data[data['Name'].str.contains('sip')]

#        
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
       
        for dilu in [1, 10]:
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
#                m = m[m['PI-A']>np.exp(5.8)]

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

           
plt.figure(0)
mng = plt.get_current_fig_manager()
mng.window.showMaximized()

#plot time series  
matplotlib.style.use('ggplot')                 
fig = plt.figure(0)
plt.clf()
plt.suptitle(r'Time Series of Flow Cytometry Data Filtered by Sybrsafe', fontsize = 20)

#countData = countData[countData['Name'] != 'Neg0sip']

countData = countData[~countData['Name'].str.contains('live')]
countData = countData[~countData['Name'].str.contains('test')]
countData = countData[~countData['Name'].str.contains('Test')]

nameList = list(set(countData['Name']))
nameList.sort()
#nameList = [nameList[0]] + nameList[2:7]
#countData1 = countData[countData['Name'] != '0mM20bsip']
#countData1 = countData1[countData1['Name'] != '0mM20bsip2x']

i = 1

MaxData1 = countData[countData['Dilution'] == 1]
Max1 = max(MaxData1['Count'])/50*1000
MaxDay1 = max(MaxData1['Day'])
#MaxData10 = countData[countData['Dilution'] == 10]
#Max10 = max(MaxData10['Count'])/50*1000
#MaxDay10 = max(MaxData10['Day'])


#######get rid of numbers in middle plots (axis), change backgroud to white
#######make it pretty
#nameList = ['Neg0sip', '0mM20bbsip', '0mM20basip']

for samp in nameList:
    
    plotData = countData[countData['Name'] == samp]
    plotData1 = plotData[plotData['Dilution'] == 1]
    ax = plt.subplot(1,len(nameList), i)
    plt.plot(plotData1['Day'], plotData1['Count']/50*1000, color = 'k')
    plt.plot(plotData1['Day'], plotData1['Count']/50*1000, 'ok')
    plt.axis([0, MaxDay1, 0, Max1+(Max1/10)])
    plt.title(r'%s'% samp, fontsize = 20)
    plt.rc('font', size = 12)
    plt.xticks(fontsize = 10)
    plt.yticks(fontsize = 10)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useMathText=True)
    ax.yaxis.major.formatter._useMathText = True    

#    plotData10 = plotData[plotData['Dilution'] == 10]
#    plt.subplot(2,len(nameList), i+len(nameList))
#    plt.plot(plotData10['Day'], plotData10['Count']/50*1000, color = 'k')
#    plt.axis([0, MaxDay10, 0, Max10+(Max10/10)])
#    
    plt.xlabel('Day', fontsize = 18)
    
    i += 1
   
    
plt.subplot(1,len(nameList), 1) 
plt.ylabel('Cell count/ml', fontsize = 18) 
#plt.subplot(2,len(nameList), 1+len(nameList))
#plt.ylabel('Cell count') 

plt.tight_layout(rect = [0.0,0,0.99,0.95])   

plt.show()
    
           

#
#matplotlib.style.use('ggplot')                 
#fig = plt.figure(1)
#plt.suptitle(r'Time series of Flow cytometry Data filtered by PI and sybr', fontsize = 18)
#
##nameList = list(set(countData['Name']))
##nameList = [nameList[1], nameList[7]]
#i = 1
#MaxData1 = countData[countData['Dilution'] == 1]
#Max1 = max(MaxData1['Count'])/50*1000
#MaxDay1 = max(MaxData1['Day'])
##MaxData10 = countData[countData['Dilution'] == 10]
##Max10 = max(MaxData10['Count'])/50*1000
##MaxDay10 = max(MaxData10['Day'])
#
#
########get rid of numbers in middle plots (axis), change backgroud to white
########make it pretty
#
#
#for samp in nameList:
#    
#    plotData = countData[countData['Name'] == samp]
#    plotData1 = plotData[plotData['Dilution'] == 1]
#    plt.subplot(1,len(nameList), i)
#    plt.plot(plotData1['Day'], plotData1['Count']/50*1000, color = 'k')
#    plt.axis([0, MaxDay1, 500000, Max1+(Max1/10)])
#    plt.title(r'%s'% samp)
#
##    plotData10 = plotData[plotData['Dilution'] == 10]
##    plt.subplot(2,len(nameList), i+len(nameList))
##    plt.plot(plotData10['Day'], plotData10['Count']/50*1000, color = 'k')
##    plt.axis([0, MaxDay10, 0, Max10+(Max10/10)])
##    
#    plt.xlabel('Day')
#    
#    i += 1
#   
#    
#plt.subplot(1,len(nameList), 1) 
#plt.ylabel('Cell count/ml') 
##plt.subplot(2,len(nameList), 1+len(nameList))
##plt.ylabel('Cell count') 
#
#
#fig.text(0.93, 0.5, '1:1', ha='right', va='center', rotation=-90)