# -*- coding: utf-8 -*-
"""
Created on Wed May 24 21:34:31 2017

@author: Lana
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd


def neg(df):
    i = False
    for d in list(df):
        i |= any(n < 0 for n in df[d])
    return i

samp_stain_file = 'C:/Users/Lana/Desktop/FACS/May29th_Day3/0mM_1-100_PI_sybr_day3_may29.0030.fcs.csv'
negctl_stain_file = 'C:/Users/Lana/Desktop/FACS/May29th_Day3/negg_1-10_PI_sybr_day3_may29.0010.fcs.csv'
medium_stain_file = 'C:/Users/Lana/Desktop/FACS/May29th_Day3/Medium_PI_sybr_day3_may29.0008.fcs.csv'
samp_nostain_file = 'C:/Users/Lana/Desktop/FACS/May29th_Day3/ExpP_1-10_nodye_day3_may29.0006.fcs.csv'
medium_nostain_file = 'C:/Users/Lana/Desktop/FACS/May29th_Day3/Medium_nodye_day3_may29.0001.fcs.csv'

s_dye = pd.read_csv(samp_stain_file)
m_dye = pd.read_csv(medium_stain_file)

#be able to choose which?
s_dye = s_dye.drop(['VioBlue-A', 'PE-A','PE-Cy7-A','APC-A', 'APC-Cy7-A'], axis=1)
m_dye = m_dye.drop(['VioBlue-A', 'PE-A','PE-Cy7-A','APC-A', 'APC-Cy7-A'], axis=1)

#TRY TO MAKE LOGGING FASTER
#while neg(s_dye):
#    s_dye += 10
#while neg(ng_dye):
#    ng_dye += 1
#while neg(m_dye):
#    m_dye += 10
#while neg(s_nodye):
#    s_nodye += 1
#while neg(m_nodye):
#    m_nodye += 1


log_s_dye = np.log(s_dye)
log_m_dye = np.log(m_dye)

matplotlib.style.use('ggplot')

plt.figure(0)
plt.subplot(2,1,1)
a = log_s_dye['FITC-A']
b = log_m_dye['FITC-A']

a.plot.kde(label='s dye', legend = True, lw=2)
b.plot.kde(label='m dye', legend = True, lw=2)
ax = plt.gca()

#####CAN WE GET THE SAME X VALUES?
line = ax.lines[0]
line2 = ax.lines[1]
dfS = pd.DataFrame({'x' : line.get_xdata(), 'y': line.get_ydata()})
dfM = pd.DataFrame({'x' : line2.get_xdata(), 'y': line2.get_ydata()})

sMaxX = dfS.loc[dfS['y'].idxmax(),'x']
mMaxX = dfM.loc[dfM['y'].idxmax(),'x']

#plt.axvline(sMaxX, color='k', linestyle='--', lw=2)
#plt.axvline(mMaxX, color='k', linestyle='--', lw=2)

dfS = dfS[dfS['x']<sMaxX]
dfS = dfS[dfS['x']>mMaxX]
dfS.reset_index(inplace = True)
dfM = dfM[dfM['x']<sMaxX]
dfM = dfM[dfM['x']>mMaxX]
dfM.reset_index(inplace = True)

dfDiff = pd.DataFrame({'x' : dfM['x'], 'y': np.abs(dfM['y']-dfS['y'])})

thresholdS = dfDiff.loc[dfDiff['y'].idxmin(),'x']

plt.axvline(thresholdS, color='k', linestyle='--', lw=2)
plt.axis([-1, 12, 0, 1.5])
plt.xlabel(r'FITC-A')
plt.title(r'Sybrsafe Flourescence vs Particle Density Plot')

plt.subplot(2,1,2)
a = log_s_dye['PI/PE-Cy5.5-A']
b = log_m_dye['PI/PE-Cy5.5-A']
a.plot.kde(label='s dye', legend = True, lw=2)
b.plot.kde(label='m dye', legend = True, lw=2)

ax1 = plt.gca()
#####CAN WE GET THE SAME X VALUES?
line = ax1.lines[0]
line2 = ax1.lines[1]
dfS = pd.DataFrame({'x' : line.get_xdata(), 'y': line.get_ydata()})
dfM = pd.DataFrame({'x' : line2.get_xdata(), 'y': line2.get_ydata()})


sMaxX = dfS.loc[dfS['y'].idxmax(),'x']
mMaxX = dfM.loc[dfM['y'].idxmax(),'x']

#plt.axvline(sMaxX, color='k', linestyle='--', lw=2)
#plt.axvline(mMaxX, color='k', linestyle='--', lw=2)

dfS = dfS[dfS['x']<sMaxX]
dfS = dfS[dfS['x']>mMaxX]
dfS.reset_index(inplace = True)
dfM = dfM[dfM['x']<sMaxX]
dfM = dfM[dfM['x']>mMaxX]
dfM.reset_index(inplace = True)

dfDiff = pd.DataFrame({'x' : dfM['x'], 'y': np.abs(dfM['y']-dfS['y'])})

thresholdP = dfDiff.loc[dfDiff['y'].idxmin(),'x']

plt.axvline(thresholdP, color='k', linestyle='--', lw=2)
plt.axis([-1, 12, 0, 1.5])
plt.xlabel(r'FITC-A')
plt.title(r'Propidium Iodide Flourescence vs Particle Density Plot')


s = log_s_dye[log_s_dye['FITC-A']>thresholdS]
m = log_m_dye[log_m_dye['FITC-A']>thresholdS]
s_bg = log_s_dye[log_s_dye['FITC-A']<thresholdS]
m_bg = log_m_dye[log_m_dye['FITC-A']<thresholdS]
s = s[s['PI/PE-Cy5.5-A']>thresholdP]
m = m[m['PI/PE-Cy5.5-A']>thresholdP]
s_bg_bg = s_bg[s_bg['PI/PE-Cy5.5-A']<thresholdP]
m_bg = m_bg[m_bg['PI/PE-Cy5.5-A']<thresholdP]


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

print 'Particle count in experimental sample before scatter filter: ',len(s.index)
print 'Particle count in Medium before scatter filter: ',len(m.index)


#go through each point left in med, look to see if there is anything within a 
#small range from it, take down index # in list. 
# delete collected index numbers, 
overlapArea = 0.1

for particle in m.index:
    done = False
    for i in s.index:
        if done==False and (m.loc[particle,'FSC-A']-s.loc[i,'FSC-A']) < overlapArea:
            if (m.loc[particle,'SSC-A']-s.loc[i,'SSC-A']) < overlapArea:
                s = s.drop(i)
                m = m.drop(particle)   #append to m_bg instead?
                done = True
        

#visually this can be shown but for count calc
#this is the same as subtraction the 2?
#actually it will not delete relevant data.

plt.figure(2)
plt.suptitle(r'Scatter plot of Flow cytometry Data filtered by PI and sybr and overlap', fontsize = 18)
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


#another option to consider savefig('foo.png')

print 'Current particle count in experimental sample: ',len(s.index)
print 'Current particle count in Medium: ',len(m.index)