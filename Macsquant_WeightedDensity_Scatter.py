# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 23:32:33 2017

@author: Lana
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import scipy.stats



sybr = 'FITC-A'
pi = 'PI/PE-Cy5.5-A'





def NoNanInf(pd):
    pd = pd.dropna()
    
samp_stain_file = 'C:/Users/Lana/Desktop/FACS/June12_day2_and_day3/0mMC_1_pi_sybr_june12_day3.0023.fcs.csv'
negctl_stain_file = 'C:/Users/Lana/Desktop/FACS/June12_day2_and_day3/negg_100_pi_sybr_june12_day3.0018.fcs.csv'
medium_stain_file = 'C:/Users/Lana/Desktop/FACS/June12_day2_and_day3/medim_pi_sybr_june12_day3.0005.fcs.csv'
samp_nostain_file = 'C:/Users/Lana/Desktop/FACS/June12_day2_and_day3/Exp_100_dye_june12_day3.0002.fcs.csv'
medium_nostain_file = 'C:/Users/Lana/Desktop/FACS/June12_day2_and_day3/medium_no_dye_june12_day3.0001.fcs.csv'
medium_PIstain_file = 'C:/Users/Lana/Desktop/FACS/July4_day3/medium_pi_sybr_july4_day3.0011.fcs.csv'
samppi_file = 'C:/Users/Lana/Desktop/FACS/July4_day3/medium_pi_sybr_july4_day3.0011.fcs.csv'

s_dye = pd.read_csv(samp_stain_file)
ng_dye = pd.read_csv(negctl_stain_file)
m_dye = pd.read_csv(medium_stain_file)
s_nodye = pd.read_csv(samp_nostain_file)
m_nodye = pd.read_csv(medium_nostain_file)
mpi_dye = pd.read_csv(medium_PIstain_file)
samppi_dye = pd.read_csv(samppi_file)

#s_dye = s_dye.drop(['BL2-A', 'FSC-H', 'SSC-H', 'SYBR Green-H', 'BL2-H', 'PI-H', 'PI_2-H', 'FSC-W', 'SSC-W', 'SYBR Green-W', 'BL2-W', 'PI-W', 'PI_2-W'], axis=1)
#m_dye = m_dye.drop(['BL2-A', 'FSC-H', 'SSC-H', 'SYBR Green-H', 'BL2-H', 'PI-H', 'PI_2-H', 'FSC-W', 'SSC-W', 'SYBR Green-W', 'BL2-W', 'PI-W', 'PI_2-W'], axis=1)
#ng_dye = ng_dye.drop(['BL2-A', 'FSC-H', 'SSC-H', 'SYBR Green-H', 'BL2-H', 'PI-H', 'PI_2-H', 'FSC-W', 'SSC-W', 'SYBR Green-W', 'BL2-W', 'PI-W', 'PI_2-W'], axis=1)
#m_nodye = m_nodye.drop(['BL2-A', 'FSC-H', 'SSC-H', 'SYBR Green-H', 'BL2-H', 'PI-H', 'PI_2-H', 'FSC-W', 'SSC-W', 'SYBR Green-W', 'BL2-W', 'PI-W', 'PI_2-W'], axis=1)
#s_nodye = s_nodye.drop(['BL2-A', 'FSC-H', 'SSC-H', 'SYBR Green-H', 'BL2-H', 'PI-H', 'PI_2-H', 'FSC-W', 'SSC-W', 'SYBR Green-W', 'BL2-W', 'PI-W', 'PI_2-W'], axis=1)
#samppi_dye = samppi_dye.drop(['BL2-A', 'FSC-H', 'SSC-H', 'SYBR Green-H', 'BL2-H', 'PI-H', 'PI_2-H', 'FSC-W', 'SSC-W', 'SYBR Green-W', 'BL2-W', 'PI-W', 'PI_2-W'], axis=1)
#mpi_dye = mpi_dye.drop(['BL2-A', 'FSC-H', 'SSC-H', 'SYBR Green-H', 'BL2-H', 'PI-H', 'PI_2-H', 'FSC-W', 'SSC-W', 'SYBR Green-W', 'BL2-W', 'PI-W', 'PI_2-W'], axis=1)

#s_dye = pd.DataFrame(s_dye[sybr])
#m_dye = m_dye[sybr]
#ng_dye = ng_dye[sybr]
#m_nodye = m_nodye[sybr]
#s_nodye = s_nodye[sybr]
#munfil_dye = munfil_dye[sybr]


log_s_dye = np.log(s_dye)
log_ng_dye = np.log(ng_dye)
log_m_dye = np.log(m_dye)
log_s_nodye = np.log(s_nodye)
log_m_nodye = np.log(m_nodye)
log_mpi_dye = np.log(mpi_dye)
log_samppi_dye = np.log(samppi_dye)


#SYBR
a = log_s_dye[sybr]
b = log_m_dye[sybr]
c = log_ng_dye[sybr]
d = log_s_nodye[sybr]
e = log_m_nodye[sybr]
#f = log_munfil_dye[sybr]
g = log_samppi_dye[sybr]
a = a.replace([np.inf, -np.inf], np.nan).dropna()
b = b.replace([np.inf, -np.inf], np.nan).dropna()
c = c.replace([np.inf, -np.inf], np.nan).dropna()
d = d.replace([np.inf, -np.inf], np.nan).dropna()
e = e.replace([np.inf, -np.inf], np.nan).dropna()
#f = f.replace([np.inf, -np.inf], np.nan).dropna()
g = g.replace([np.inf, -np.inf], np.nan).dropna()


kdea = scipy.stats.gaussian_kde(a)
kdeb = scipy.stats.gaussian_kde(b)
kdec = scipy.stats.gaussian_kde(c)
kded = scipy.stats.gaussian_kde(d)
kdee = scipy.stats.gaussian_kde(e)
#kdef = scipy.stats.gaussian_kde(f)
kdeg = scipy.stats.gaussian_kde(g)



total = 0
for i in [a,b,c,d,e,g]:
    total += len(i)

grid = np.arange(0,12,0.1)



#weighted kde curves
wa = kdea(grid)*(len(a)/float(total))
wb = kdeb(grid)*(len(b)/float(total))
wc = kdec(grid)*(len(c)/float(total))
wd = kded(grid)*(len(d)/float(total))
we = kdee(grid)*(len(e)/float(total))
#wf = kdef(grid)*(len(f)/float(total))
wg = kdeg(grid)*(len(g)/float(total))


#wa = kdea(grid)
#wb = kdeb(grid)
#wc = kdec(grid)
#wd = kded(grid)
#we = kdee(grid)
#wf = kdef(grid)
#wg = kdeg(grid)

fig, ax = plt.subplots(1)
ax.plot(grid, wa, lw=1, label = 's dye (Sybr)')
ax.plot(grid, wb, lw=1, label = 'm dye (Sybr)')
ax.plot(grid, wc, lw=1, label = 'ng dye (Sybr)')
ax.plot(grid, wd, lw=1, label = 's no dye')
ax.plot(grid, we, lw=1, label = 'm no dye')
#ax.plot(grid, wg, lw=1, label = 's dye (PI)')
ax.axvline(6.6, color='k', linestyle='--', lw=2)
ax.set_title('Sybrsafe Flourescence vs Particle Density Plot')
ax.set_xlabel('Sybr')
ax.set_ylabel('Density')
ax.legend()


## PI
#a = log_s_dye[pi]
#b = log_mpi_dye[pi]
#c = log_ng_dye[pi]
#d = log_s_nodye[pi]
#e = log_m_nodye[pi]
##f = log_munfil_dye[pi]
#g = log_samppi_dye[pi]
#a = a.replace([np.inf, -np.inf], np.nan).dropna()
#b = b.replace([np.inf, -np.inf], np.nan).dropna()
#c = c.replace([np.inf, -np.inf], np.nan).dropna()
#d = d.replace([np.inf, -np.inf], np.nan).dropna()
#e = e.replace([np.inf, -np.inf], np.nan).dropna()
##f = f.replace([np.inf, -np.inf], np.nan).dropna()
#g = g.replace([np.inf, -np.inf], np.nan).dropna()
#kdea = scipy.stats.gaussian_kde(a)
#kdeb = scipy.stats.gaussian_kde(b)
#kdec = scipy.stats.gaussian_kde(c)
#kded = scipy.stats.gaussian_kde(d)
#kdee = scipy.stats.gaussian_kde(e)
##kdef = scipy.stats.gaussian_kde(f)
#kdeg = scipy.stats.gaussian_kde(g)
#
#total = 0
#for i in [a,b,c,d,e,g]:
#    total += len(i)
#
#grid = np.arange(0,15,0.1)
#
##weighted kde curves
#wa = kdea(grid)*(len(a)/float(total))
#wb = kdeb(grid)*(len(b)/float(total))
#wc = kdec(grid)*(len(c)/float(total))
#wd = kded(grid)*(len(d)/float(total))
#we = kdee(grid)*(len(e)/float(total))
##wf = kdef(grid)*(len(f)/float(total))
#wg = kdeg(grid)*(len(g)/float(total))
#
##wa = kdea(grid)
##wb = kdeb(grid)
##wc = kdec(grid)
##wd = kded(grid)
##we = kdee(grid)
##wf = kdef(grid)
##wg = kdeg(grid)
#
#ax2.plot(grid, wa, lw=1, label = 's dye (Sybr)')
#ax2.plot(grid, wb, lw=1, label = 'm dye (PI)')
#ax2.plot(grid, wc, lw=1, label = 'ng dye (Sybr)')
#ax2.plot(grid, wd, lw=1, label = 's no dye')
#ax2.plot(grid, we, lw=1, label = 'm no dye')
##ax2.plot(grid, wg, lw=1, label = 's dye (PI)')
#ax2.axvline(7, color='k', linestyle='--', lw=2)
#ax2.set_title(r'Propidium Iodide Flourescence vs Particle Density Plot')
#ax2.set_xlabel('PI')
#ax2.legend()
plt.show()




ssyb = log_s_dye[log_s_dye[sybr]>6.6]
msyb = log_m_dye[log_m_dye[sybr]>6.6]
#spi = log_samppi_dye[log_samppi_dye[pi]>7]
#mpi = log_mpi_dye[log_mpi_dye[pi]>7]

#
#plt.figure(3)
#plt.suptitle(r'Scatter plot of Flow cytometry Data filtered by PI and sybr', fontsize = 18)
#plt.subplot(2,2,1)
#plt.scatter(ssyb['FSC-A'], ssyb['SSC-A'], color = 'r', alpha = 0.1, label = 'syb filtered sample')
#plt.scatter(msyb['FSC-A'], msyb['SSC-A'], color = 'b', alpha = 0.1, label = 'syb filtered medium')
#plt.axis([7, 14, 7, 14])
#plt.xlabel('FSC-A')
#plt.ylabel('SSC-A')
#plt.legend(loc ='lower right')
#plt.subplot(2,2,2)
#plt.scatter(log_s_dye['FSC-A'], log_s_dye['SSC-A'], color = 'r', alpha = 0.1, label = 'syb sample')
#plt.scatter(log_m_dye['FSC-A'], log_m_dye['SSC-A'], color = 'b', alpha = 0.1, label = 'syb medium')
#plt.axis([7, 14, 7, 14])
#plt.xlabel('FSC-A')
#plt.ylabel('SSC-A')
#plt.legend(loc ='lower right')
#plt.subplot(2,2,3)
#plt.scatter(spi['FSC-A'], spi['SSC-A'], color = 'y', alpha = 0.1, label = 'PI filtered sample')
#plt.scatter(mpi['FSC-A'], mpi['SSC-A'], color = 'g', alpha = 0.1, label = 'PI filtered medium')
#plt.axis([7, 14, 7, 14])
#plt.xlabel('FSC-A')
#plt.ylabel('SSC-A')
#plt.legend(loc ='lower right')
#plt.subplot(2,2,4)
#plt.scatter(log_samppi_dye['FSC-A'], log_samppi_dye['SSC-A'], color = 'y', alpha = 0.1, label = 'PI sample')
#plt.scatter(log_mpi_dye['FSC-A'], log_mpi_dye['SSC-A'], color = 'g', alpha = 0.1, label = 'PI medium')
#plt.axis([7, 14, 7, 14])
#plt.xlabel('FSC-A')
#plt.ylabel('SSC-A')
#plt.legend(loc ='lower right')


plt.figure(3)
plt.suptitle(r'Scatter plot of Flow cytometry Data filtered by PI and sybr', fontsize = 18)
plt.subplot(2,1,1)
plt.scatter(ssyb['FSC-A'], ssyb['SSC-A'], color = 'r', alpha = 0.1, label = 'syb filtered sample')
plt.scatter(msyb['FSC-A'], msyb['SSC-A'], color = 'b', alpha = 0.1, label = 'syb filtered medium')
plt.axis([5, 11, 0, 13])
plt.xlabel('FSC-A')
plt.ylabel('SSC-A')
plt.legend(loc ='lower right')
plt.subplot(2,1,2)
plt.scatter(log_s_dye['FSC-A'], log_s_dye['SSC-A'], color = 'r', alpha = 0.1, label = 'syb sample')
plt.scatter(log_m_dye['FSC-A'], log_m_dye['SSC-A'], color = 'b', alpha = 0.1, label = 'syb medium')
plt.axis([5, 11, 0, 13])
plt.xlabel('FSC-A')
plt.ylabel('SSC-A')
plt.legend(loc ='lower right')


print 'Current particle count in pi: ',len(spi.index)
print 'Current particle count in syb: ',len(ssyb.index)