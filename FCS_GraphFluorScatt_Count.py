# -*- coding: utf-8 -*-
"""
Created on Wed May 24 21:34:31 2017

@author: Lana
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd


samp_stain_file = 'C:/Users/Lana/Desktop/FACS/July4_day3/6mM_10_pi_sybr_july4_day3.0026.fcs.csv'
negctl_stain_file = 'C:/Users/Lana/Desktop/FACS/July4_day3/negg8mM_10_pi_sybr_july4_day3.0020.fcs.csv'
medium_stain_file = 'C:/Users/Lana/Desktop/FACS/July4_day3/medium8mM_pi_sybr_july4_day3.0012.fcs.csv'
samp_nostain_file = 'C:/Users/Lana/Desktop/FACS/July4_day3/exp_10_no_dye_july4_day3.0008.fcs.csv'
medium_nostain_file = 'C:/Users/Lana/Desktop/FACS/July4_day3/medium_no_dye_july4_day3.0002.fcs.csv'
mediumfil_nostain_file = 'C:/Users/Lana/Desktop/FACS/July4_day3/medium8mM_no_dye_july4_day3.0003.fcs.csv'

s_dye = pd.read_csv(samp_stain_file)
ng_dye = pd.read_csv(negctl_stain_file)
m_dye = pd.read_csv(medium_stain_file)
s_nodye = pd.read_csv(samp_nostain_file)
m_nodye = pd.read_csv(medium_nostain_file)
munfil_dye = pd.read_csv(medium_nostain_file)

s_dye = s_dye.drop(['VioBlue-A', 'PE-A','PE-Cy7-A','APC-A', 'APC-Cy7-A'], axis=1)
m_dye = m_dye.drop(['VioBlue-A', 'PE-A','PE-Cy7-A','APC-A', 'APC-Cy7-A'], axis=1)
ng_dye = ng_dye.drop(['VioBlue-A', 'PE-A','PE-Cy7-A','APC-A', 'APC-Cy7-A'], axis=1)
m_nodye = m_nodye.drop(['VioBlue-A', 'PE-A','PE-Cy7-A','APC-A', 'APC-Cy7-A'], axis=1)
s_nodye = s_nodye.drop(['VioBlue-A', 'PE-A','PE-Cy7-A','APC-A', 'APC-Cy7-A'], axis=1)
munfil_dye = munfil_dye.drop(['VioBlue-A', 'PE-A','PE-Cy7-A','APC-A', 'APC-Cy7-A'], axis=1)



log_s_dye = np.log(s_dye)
log_ng_dye = np.log(ng_dye)
log_m_dye = np.log(m_dye)
log_s_nodye = np.log(s_nodye)
log_m_nodye = np.log(m_nodye)
log_munfil_dye = np.log(munfil_dye)

matplotlib.style.use('ggplot')

plt.figure(0)
plt.subplot(2,1,1)
a = log_s_dye['FITC-A']
b = log_m_dye['FITC-A']
c = log_ng_dye['FITC-A']
d = log_s_nodye['FITC-A']
e = log_m_nodye['FITC-A']
f = log_munfil_dye['FITC-A']
a.plot.kde(label='s dye', legend = True, lw=2)
b.plot.kde(label='m dye', legend = True, lw=2)
c.plot.kde(label='ng dye', legend = True, lw=2)
d.plot.kde(label='s no dye', legend = True, lw=2)
e.plot.kde(label='m no dye', legend = True, lw=2)
f.plot.kde(label='m fil dye', legend = True, lw=2)
plt.axvline(6.2, color='k', linestyle='--', lw=2)
plt.axis([-1, 12, 0, 1.5])
plt.xlabel(r'FITC-A')
plt.title(r'Sybrsafe Flourescence vs Particle Density Plot')

plt.subplot(2,1,2)
a = log_s_dye['PI/PE-Cy5.5-A']
b = log_m_dye['PI/PE-Cy5.5-A']
c = log_ng_dye['PI/PE-Cy5.5-A']
d = log_s_nodye['PI/PE-Cy5.5-A']
e = log_m_nodye['PI/PE-Cy5.5-A']
f = log_munfil_dye['PI/PE-Cy5.5-A']
a.plot.kde(label='s dye', legend = True, lw=2)
b.plot.kde(label='m dye', legend = True, lw=2)
c.plot.kde(label='ng dye', legend = True, lw=2)
d.plot.kde(label='s no dye', legend = True, lw=2)
e.plot.kde(label='m no dye', legend = True, lw=2)
f.plot.kde(label='m fil dye', legend = True, lw=2)
plt.axvline(5.8, color='k', linestyle='--', lw=2)
plt.axis([-1, 12, 0, 1.5])
plt.xlabel(r'FITC-A')
plt.title(r'Propidium Iodide Flourescence vs Particle Density Plot')

s = log_s_dye[log_s_dye['FITC-A']>6.2]
m = log_m_dye[log_m_dye['FITC-A']>6.2]
s_bg = log_s_dye[log_s_dye['FITC-A']<6.2]
m_bg = log_m_dye[log_m_dye['FITC-A']<6.2]
s = s[s['PI/PE-Cy5.5-A']>5.8]
m = m[m['PI/PE-Cy5.5-A']>5.8]
s_bg_bg = s_bg[s_bg['PI/PE-Cy5.5-A']<5.8]
m_bg = m_bg[m_bg['PI/PE-Cy5.5-A']<5.8]

plt.figure(1)
plt.suptitle(r'Scatter plot of Flow cytometry Data filtered by PI and sybr', fontsize = 18)
plt.subplot(2,2,1)
plt.scatter(s['FSC-A'], s['SSC-A'], color = 'r', alpha = 0.1, label = 'filtered sample')
plt.scatter(m['FSC-A'], m['SSC-A'], color = 'b', alpha = 0.1, label = 'filtered medium')
plt.axis([5, 13, 3, 13])
plt.xlabel('FSC-A')
plt.ylabel('SSC-A')
plt.legend(loc ='lower left')
plt.subplot(2,2,2)
plt.scatter(m_bg['FSC-A'], m_bg['SSC-A'], color = 'y', alpha = 0.1, label = 'Background medium')
plt.scatter(s_bg['FSC-A'], s_bg['SSC-A'], color = 'g', alpha = 0.1, label = 'Background sample')
plt.axis([5, 13, 3, 13])
plt.xlabel('FSC-A')
plt.ylabel('SSC-A')
plt.legend(loc ='lower left')
plt.subplot(2,2,3)
plt.scatter(m['FSC-A'], m['SSC-A'], color = 'b', alpha = 0.1, label = 'filtered medium')
plt.scatter(m_bg['FSC-A'], m_bg['SSC-A'], color = 'c', alpha = 0.1, label = 'Background medium')
plt.axis([5, 13, 3, 13])
plt.xlabel('FSC-A')
plt.ylabel('SSC-A')
plt.legend(loc ='lower left')
plt.subplot(2,2,4)
plt.scatter(s['FSC-A'], s['SSC-A'], color = 'r', alpha = 0.1, label = 'filtered sample')
plt.scatter(s_bg['FSC-A'], s_bg['SSC-A'], color = 'g', alpha = 0.1, label = 'Background sample')
plt.axis([5, 13, 3, 13])
plt.xlabel('FSC-A')
plt.ylabel('SSC-A')
plt.legend(loc ='lower left')

print 'Current particle count in experimental sample: ',len(s.index)
print 'Current particle count in Medium: ',len(m.index)