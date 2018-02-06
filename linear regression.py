# -*- coding: utf-8 -*-
"""
Created on Sat Dec 16 15:33:54 2017

@author: Lana
"""


import glob as g
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import lmfit
import random
import scipy.stats


#DATA = pd.DataFrame(columns = ['Name','Day', 'Cell', 'Nitrate'])
#
#DATA['Name'] = Cell['Name']
#DATA['Day'] = Cell['Day']
#DATA['Cell'] = Cell['Count']
#
#Nitrate['Name'].replace(to_replace='Neg0', value='Negg0mM', inplace = True)
#Nitrate['Name'].replace(to_replace='Neg8', value='Negg8mM', inplace = True)


#for i in Nitrate.index:
#    for n in DATA.index:
#        if DATA.loc[n,'Name'] == Nitrate.loc[i,'Name'] and DATA.loc[n,'Day'] == Nitrate.loc[i,'Day']:
#            DATA.loc[n,'Nitrate'] = Nitrate.loc[i,'Nitrate']
#           
            

nList = list(DATAa['Nitrate'])
cList = list(DATAa['Cell'])
#core = np.corrcoef(cList,nList)
#print core[0,1]**2
print scipy.stats.pearsonr(cList,nList)
    