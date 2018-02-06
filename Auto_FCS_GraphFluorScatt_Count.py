# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 15:22:45 2017

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
EnrichDate = 'May29'    #input('Enter string of enrichment batch desired with no space ei. "May16":')


#####MAKE FILE LOCATION MORE GENERAL?############
Files = g.glob('C:\Users\Lana\Desktop\FACS\*\*'+ EnrichDate + '*.csv')

data = pd.DataFrame(columns = ['File','Name','Dilution','Day', 'Stain'])

# Collect file info into data frame
i = 0 
for csv in Files:
    data.loc[i,'File'] = csv
    name = csv.split("\\")[-1]    
    data.loc[i,'Name'] =  capitalize(name.split("_")[0])
    
    csv = csv.lower()
    date= csv.split("\\")[-1]
    date= date.split('day',1)[1]
    date = date.split('_',1)[0]
    data.loc[i,'Day'] = int(date.split('.',1)[0])

#####ALLOW CHOOSING FOR DYE?#####
    if 'pi' in csv:
        data.loc[i,'Stain'] = True
    else:
        data.loc[i,'Stain'] = False
    
    if '100' in csv:
        data.loc[i,'Dilution'] = 100
    elif '10' in csv:
        data.loc[i,'Dilution'] = 10
    else:
        data.loc[i,'Dilution'] = 1
        
    i += 1
        
##########GIVE OPTION TO DO IT FOR SPECIFIC GROUPS##################
##########(especially in graphing)    
     
for day in list(set(data['Day'])):   
    #Filter for dataframe for this day 
    dataDay = data[data['Day'] == day]
    
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
            
            s_dye = pd.read_csv(samp_stain_file)
            m_dye = pd.read_csv(medium_stain_file)
            s_nodye = pd.read_csv(samp_nostain_file)
            m_nodye = pd.read_csv(medium_nostain_file)
            
            s_dye = s_dye.drop(['VioBlue-A', 'PE-A','PE-Cy7-A','APC-A', 'APC-Cy7-A'], axis=1)
            m_dye = m_dye.drop(['VioBlue-A', 'PE-A','PE-Cy7-A','APC-A', 'APC-Cy7-A'], axis=1)
            m_nodye = m_nodye.drop(['VioBlue-A', 'PE-A','PE-Cy7-A','APC-A', 'APC-Cy7-A'], axis=1)
            s_nodye = s_nodye.drop(['VioBlue-A', 'PE-A','PE-Cy7-A','APC-A', 'APC-Cy7-A'], axis=1)            
            
            log_s_dye = np.log(s_dye)
            log_m_dye = np.log(m_dye)
            log_s_nodye = np.log(s_nodye)
            log_m_nodye = np.log(m_nodye)
            
            matplotlib.style.use('ggplot')
            
            plt.figure(0)
            plt.subplot(2,1,1)
            a = log_s_dye['FITC-A']
            b = log_m_dye['FITC-A']
            d = log_s_nodye['FITC-A']
            e = log_m_nodye['FITC-A']
            a.plot.kde(label='s dye', legend = True, lw=2)
            b.plot.kde(label='m dye', legend = True, lw=2)
            d.plot.kde(label='s no dye', legend = True, lw=2)
            e.plot.kde(label='m no dye', legend = True, lw=2)
            plt.axvline(6.2, color='k', linestyle='--', lw=2)
            plt.axis([-1, 12, 0, 1.5])
            plt.xlabel(r'FITC-A')
            plt.title(r'Sybrsafe Flourescence vs Particle Density Plot')
            
            plt.subplot(2,1,2)
            a = log_s_dye['PI/PE-Cy5.5-A']
            b = log_m_dye['PI/PE-Cy5.5-A']
            d = log_s_nodye['PI/PE-Cy5.5-A']
            e = log_m_nodye['PI/PE-Cy5.5-A']
            a.plot.kde(label='s dye', legend = True, lw=2)
            b.plot.kde(label='m dye', legend = True, lw=2)
            d.plot.kde(label='s no dye', legend = True, lw=2)
            e.plot.kde(label='m no dye', legend = True, lw=2)
            plt.axvline(5.8, color='k', linestyle='--', lw=2)
            plt.axis([-1, 12, 0, 1.5])
            plt.xlabel(r'FITC-A')
            plt.title(r'Propidium Iodide Flourescence vs Particle Density Plot')
            
            
            
            #FILTER##########################################
            s = log_s_dye[log_s_dye['FITC-A']>6.2]
            m = log_m_dye[log_m_dye['FITC-A']>6.2]
            s_bg = log_s_dye[log_s_dye['FITC-A']<6.2]
            m_bg = log_m_dye[log_m_dye['FITC-A']<6.2]
            s = s[s['PI/PE-Cy5.5-A']>5.8]
            m = m[m['PI/PE-Cy5.5-A']>5.8]
            s_bg_bg = s_bg[s_bg['PI/PE-Cy5.5-A']<5.8]
            m_bg = m_bg[m_bg['PI/PE-Cy5.5-A']<5.8]
            #################################################
            
            
            plt.figure(1)
            plt.suptitle(r'Scatter plot of Flow cytometry Data filtered by PI and sybr', fontsize = 18)
            plt.subplot(2,2,1)
            plt.scatter(s['FSC-A'], s['SSC-A'], color = 'r', alpha = 0.1, label = 'filtered sample')
            plt.scatter(m['FSC-A'], m['SSC-A'], color = 'b', alpha = 0.1, label = 'filtered medium')
            plt.axis([5, 11, 3, 11])
            plt.xlabel('FSC-A')
            plt.ylabel('SSC-A')
            plt.legend(loc ='lower left')
            plt.subplot(2,2,2)
            plt.scatter(m_bg['FSC-A'], m_bg['SSC-A'], color = 'y', alpha = 0.1, label = 'Background medium')
            plt.scatter(s_bg['FSC-A'], s_bg['SSC-A'], color = 'g', alpha = 0.1, label = 'Background sample')
            plt.axis([5, 11, 3, 11])
            plt.xlabel('FSC-A')
            plt.ylabel('SSC-A')
            plt.legend(loc ='lower left')
            plt.subplot(2,2,3)
            plt.scatter(m['FSC-A'], m['SSC-A'], color = 'b', alpha = 0.1, label = 'filtered medium')
            plt.scatter(m_bg['FSC-A'], m_bg['SSC-A'], color = 'c', alpha = 0.1, label = 'Background medium')
            plt.axis([5, 11, 3, 11])
            plt.xlabel('FSC-A')
            plt.ylabel('SSC-A')
            plt.legend(loc ='lower left')
            plt.subplot(2,2,4)
            plt.scatter(s['FSC-A'], s['SSC-A'], color = 'r', alpha = 0.1, label = 'filtered sample')
            plt.scatter(s_bg['FSC-A'], s_bg['SSC-A'], color = 'g', alpha = 0.1, label = 'Background sample')
            plt.axis([5, 11, 3, 11])
            plt.xlabel('FSC-A')
            plt.ylabel('SSC-A')
            plt.legend(loc ='lower left')
            
            print 'Current particle count in experimental sample %s: ' % dataDayDilu.loc[i,'Name'],len(s.index)
            print 'Current particle count in Medium: ',len(m.index)
                       
    
